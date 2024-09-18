## these code are used to generate data for plotting figure 2 and 3 the only two new functions comapred to simulation.R are true_distance and true_CMDS

library(segmented)
library(igraph)
library(irlba)
library(locfit)
library(dplyr)
library(ggplot2)
library(plyr)
library(doParallel)
library(foreach)
library(broom)
library(vegan)
library(ggrepel)
library(tidyr)
library(tidyverse)
library(lattice)
library(Matrix)
library("corrplot")

rdpg.sample <- function(X, rdpg_scale=FALSE) {
  P <- X %*% t(X)
  if (rdpg_scale) {
    P <- scales::rescale(P)
  }
  n <-  nrow(P)
  U <- matrix(0, nrow = n, ncol = n)
  U[col(U) > row(U)] <- runif(n*(n-1)/2)
  U <- (U + t(U))
  diag(U) <- runif(n)
  A <- (U < P) + 0 ;
  diag(A) <- 0
  return(graph.adjacency(A,"undirected"))
}
full.ase <- function(A, d, diagaug=TRUE, doptr=FALSE) {
  #    require(irlba)
  
  # doptr
  if (doptr) {
    g <- ptr(A)
    A <- g[]
  } else {
    A <- A[]
  }
  
  # diagaug
  if (diagaug) {
    diag(A) <- rowSums(A) / (nrow(A)-1)
  }
  
  A.svd <- irlba(A,d)
  Xhat <- A.svd$u %*% diag(sqrt(A.svd$d))
  Xhat.R <- NULL
  
  if (!isSymmetric(A)) {
    Xhat.R <- A.svd$v %*% diag(sqrt(A.svd$d))
  }
  
  return(list(eval=A.svd$d, Xhat=Matrix(Xhat), Xhat.R=Xhat.R))
}
doSim <- function(n=300, tmax=40, delta=0.1, p=0.4, q=0.9, tstar=20)
{    
  glist <- NULL
  Xlist <- NULL
  Xt <- matrix(0,n,(tmax+1))
  
  for (t in 2:(tstar+1)) {
    tmp <- runif(n) < p
    Xt[,t] <- Xt[,t] + Xt[,t-1]
    Xt[tmp,t] <- Xt[tmp,t] + delta
  }
  
  for (t in (tstar+2):(tmax+1)) {
    tmp <- runif(n) < q
    Xt[,t] <- Xt[,t] + Xt[,t-1]
    Xt[tmp,t] <- Xt[tmp,t] + delta
  }
  Xt <- Xt[,-1]
  Xt <- Xt+0.1
  
  df <- tibble(time=1:tmax) %>%
    mutate(Xt = map(time, function(x) matrix(Xt[,x],n,1)  )) %>%
    mutate(g = map(Xt, ~rdpg.sample(.)))
  
  df
}

procrustes2 <- function(X, Y) {
  tmp <- t(X) %*% Y
  tmp.svd <- svd(tmp)
  W <- tmp.svd$u %*% t(tmp.svd$v)
  newX <- X %*% W
  return(list(newX = newX, error = norm(newX-Y, type="F"), W = W))
}
getD <- function(Xlist, k=0, etype="proc") {
  if (k==0) {
    ind <- 1:n
  } else {
    ind <- which(Yhat==k)
  }
  comb <- combn(m,2)
  Dout <- foreach (k = 1:ncol(comb), .combine='rbind') %dopar% {
    i <- comb[1,k]
    j <- comb[2,k]
    #cat("i = ", i, ", j = ", j, "\n")
    
    if (etype == "proc") {
      Xhati <- Xlist[[i]][ind,] # 32277 x Khat
      Xhatj <- Xlist[[j]][ind,]
      proc <- procrustes2(as.matrix(Xhati), as.matrix(Xhatj))
      Xhati <- Xhati %*% proc$W
    } else {
      Xhati <- Xlist[[i]][ind,] # 32277 x Khat
      Xhatj <- Xlist[[j]][ind,]
    }
    
    D <- norm(Xhati - Xhatj, type="2")^2/n
    tibble(i=i, j=j, D=D)
  }
  D2 <- matrix(0,m,m)
  D2[t(comb)] <- Dout$D
  D2 <- (D2 + t(D2)) / 1
  #as.dist(D2)
  D2 <- sqrt(D2)
  D2
}
doMDS <- function(D, doplot=TRUE)
{
  mds <- cmdscale(D, m-1)
  df.mds <- tibble(ind=1:tmax, time=sprintf("%2d",1:tmax), x=mds[,1], y=mds[,2], z=mds[,3], w=mds[,4])
  
  if (doplot) {
    plot(apply(mds,2,sd), type="b", main="", xlab="dimension", ylab="column stdev")
    
    p1 <- df.mds %>% ggplot(aes(x=ind, y=x, color=time, group=1)) +
      geom_point(size=3) + geom_line() +
      geom_vline(xintercept = tstar, linetype="dashed") +
      theme(legend.position = "none") + labs(x="time",y="mds1") #+
    print(p1)
    
    p1 <- df.mds %>% ggplot(aes(x=ind, y=y, color=time, group=1)) +
      geom_point(size=3) + geom_line() +
      geom_vline(xintercept = tstar, linetype="dashed") +
      theme(legend.position = "none") + labs(x="time",y="mds2") #+
    print(p1)
    
    p1 <- df.mds %>% ggplot(aes(x=ind, y=z, color=time, group=1)) +
      geom_point(size=3) + geom_line() +
      geom_vline(xintercept = tstar, linetype="dashed") +
      theme(legend.position = "none") + labs(x="time",y="mds3") #+
    print(p1)
    
    p2 <- df.mds %>% ggplot(aes(x=x, y=y, color=time)) +
      geom_point(size=3) +
      geom_label_repel(aes(label=time), size=2) +
      theme(legend.position = "none") + labs(x="mds1",y="mds2") #+
    print(p2)
  }
  
  return(list(mds=mds, df.mds=df.mds))
}    
doIso <- function(mds, mdsd=2, isod=1, doplot=F)
{    
  df.iso <- NULL
  dis <- vegdist(mds[,1:mdsd,drop=F], "euclidean")
  knn <- 1
  success <- FALSE
  while(!success) {
    tryCatch({
      iso = isomap(dis, k=knn, ndim=isod, path="shortest")$points
      success <- TRUE
    },
    error = function(e) {
      knn <<- knn + 1
    })
  }
  iso2 <- tibble(iso=iso[,1]) %>% mutate(i=1:nrow(mds), ind=df.mds$ind, time=df.mds$time, knn=knn)
  df.iso <- rbind(df.iso, cbind(iso2, mdsd=mdsd))
  df.iso <- df.iso %>% group_by(mdsd) %>% mutate(iso = if(iso[1] > 0) {-iso} else {iso}) %>% ungroup()
  
  if (doplot) {
    p <- df.iso %>% filter(mdsd==mdsd) %>%
      ggplot(aes(x=ind, y=iso, color=time, group=1)) +
      geom_point(size=3) + geom_line() +
      theme(legend.position = "none") +
      labs(x="time", y="isomap embedding") +
      # scale_x_date(breaks = scales::breaks_pretty(8), labels=label_date_short()) +
      theme(axis.text.x=element_text(hjust=0.7))
    # theme(axis.text.x = element_text(size = 12, angle = 90, vjust = 0.3),
    #       axis.text.y = element_text(size = 12),
    #       axis.title = element_text(size = 14, face="bold"))
    #    p <- p + scale_x_date(breaks = scales::breaks_pretty(8), labels=label_date_short())
    print(p)
    
    df.isok <- df.iso %>% filter(mdsd==mdsd) #%>% mutate(date2 = format(ymd(paste0(date,"-01")),"%m/%y"))
    row.names(df.isok) <- df.isok$time
    fit <- lm(iso ~ i, data=df.isok)
    # print(tidy(fit))
    # print(glance(fit))
    myfor <- augment(fit)
    myfor2 <- myfor %>% mutate(date=.rownames, 
                               ranks = rank(.sigma),                                       
                               mycol=sprintf("%2d",rank(.fitted)))
    p <- myfor2 %>%
      ggplot(aes(.fitted, .resid)) +
      geom_point(aes(color=mycol)) + 
      geom_hline(yintercept = 0, linetype="dashed", color="grey") +
      geom_smooth(method="loess", se=FALSE) +
      labs(x="Fitted Values", y="Residuals") +
      theme(legend.position = "none",
            axis.title = element_text(size=14, face="bold"))
    p <- p + geom_label_repel(aes(label=date), data=myfor2 %>% filter(ranks %in% 1:3))
    print(p)
  }
  
  return(df.iso)
}

## true dmv distance
true_distance=function(m,p){
  return( (m*(m-1)*p^2+m*p)*delta^2 )
}

##CMDS on the true distance matrix for model with change point                  
true_cmds=function(tt,t0,p,q,cdim){
  
  D=matrix(0,tt,tt)
  for (i in 1:t0) {
    for (j in (i):t0) {
      D[i,j]=sqrt(true_distance(abs(i-j),p))
    }
  }
  
  for (i in t0:tt) {
    for (j in (i):tt) {
      D[i,j]=sqrt(true_distance(abs(i-j),q))
    }
  }
  
  for (i in 1:(t0-1)) {
    for (j in (t0+1):tt) {
      m=t0-i
      n=j-t0
      D[i,j]=sqrt(true_distance(m,p)+true_distance(n,q)+2*m*n*p*q*delta^2)
    }
  }
  
  D=D+t(D)
  D2=D^2
  
  P=diag(rep(1,tt))-rep(1,tt)%*%t(rep(1,tt))/tt
  B=(-1/(2))*P%*%D2%*%P
  
  svdb= irlba(B,cdim)
  df_true=as.data.frame(svdb$u[,1:cdim]%*%sqrt(diag(svdb$d[1:cdim])))
  
  if(df_true$V1[1]>0){
    df_true$V1=-df_true$V1
  }
  
  df_true
}


pacman::p_load("doParallel")
registerDoParallel(detectCores()-1)

## mds_nochange and df_true_nochange generation
set.seed(2)
n=1500
tmax <- m <- 40
p <- q <- 0.4
delta <- (1-0.1)/tmax
tstar <- 20

df <- doSim(n,tmax,delta,p,q,tstar)
Dhat <- getD(df$Xt) 
df <- df %>% mutate(Xhat = map(g, function(x) full.ase(x,2)$Xhat[,1,drop=F]))
D2 <- getD(df$Xhat) 
df.mds <- doMDS(D2,doplot = F)
mds_nochange <- df.mds$mds
df_true_nochange=true_cmds(tmax,tstar,p,q,3)


#mds_change and df_true_change generation
set.seed(2)
n=1500
tmax <- m <- 40
p <- 0.4
q <- 0.2
delta <- (1-0.1)/tmax
tstar <- 20

df <- doSim(n,tmax,delta,p,q,tstar)
Dhat <- getD(df$Xt) 
df <- df %>% mutate(Xhat = map(g, function(x) full.ase(x,2)$Xhat[,1,drop=F]))
D2 <- getD(df$Xhat) 
df.mds <- doMDS(D2,doplot = F)
mds_change <- df.mds$mds
df_true_change=true_cmds(tmax,tstar,p,q,3)
## code above are used to generate mds_change/nochange and df_true_change/nochange data, if you just want to replicates the data, you can only use code below 



                               
##Please load the figure_data, although the above code can generate mds_change/nochange and df_true_change/nochange, the real data and the numerical experiment result is not included.  
load("~/Figure_data.RData")

pacman::p_load(tidyverse)

theme_update(axis.text = element_text(size=12),
             strip.text = element_text(size = 12, face="bold"),
             axis.title = element_text(size=14, face="bold"))

#### Figure 2
par(mfrow=c(1,3))
tt=m
x=1:m
plot(x,df_true_nochange$V1,xlab = '',main = 'dim 1' ,ylab = 'multidimensional scaling',cex.lab=1.5,cex.axis=1.5)
points(x,mds_nochange[,1],col='blue')
plot(x,mds_nochange[,2],col='blue',xlab = '',main = 'dim 2',ylab='',ylim=c(-0.02,0.02),cex.lab=1.5,cex.axis=1.5)
points(x,df_true_nochange$V2)
## figure size is 9**3.53
plot(x,mds_nochange[,3],col='blue',xlab = '',main = 'dim 3',ylab='',ylim=c(-0.02,0.02),cex.lab=1.5,cex.axis=1.5)
points(x,df_true_nochange$V3)
mtext('time',side = 1, line = -2, outer = TRUE,cex=1)

## ggplot version
df.fig2 = rbind(tibble(x=x, y=df_true_nochange$V1, dim="dim 1", type="numerical"),
                tibble(x=x, y=mds_nochange[,1], dim="dim 1", type="mirror"),
                tibble(x=x, y=df_true_nochange$V2, dim="dim 2", type="numerical"),
                tibble(x=x, y=mds_nochange[,2], dim="dim 2", type="mirror"),
                tibble(x=x, y=df_true_nochange$V3, dim="dim 3", type="numerical"),
                tibble(x=x, y=mds_nochange[,3], dim="dim 3", type="mirror"))
df.fig2 |> ggplot(aes(x=x, y=y, color=type)) +
  geom_point(alpha=0.5) +
  facet_wrap(~dim, scale="free") +
  labs(x="time", y="multidimensional scaling") +
  scale_color_manual(values=c("numerical"="black", "mirror"="blue")) +
  theme_bw() +
  theme(legend.position = "none",
        axis.text = element_text(size=12),
        strip.text = element_text(size = 12, face="bold"),
        axis.title = element_text(size=14, face="bold"))
ggsave("Fig2-ggplot.pdf", width=9, height=3, units="in")



####New fig 2 without out the estimated mirror 
df.fig2 = rbind(tibble(x=x, y=df_true_nochange$V1, dim="dim 1", type="numerical"),
                #tibble(x=x, y=mds_nochange[,1], dim="dim 1", type="mirror"),
                tibble(x=x, y=df_true_nochange$V2, dim="dim 2", type="numerical"),
                #tibble(x=x, y=mds_nochange[,2], dim="dim 2", type="mirror"),
                tibble(x=x, y=df_true_nochange$V3, dim="dim 3", type="numerical"))
                #tibble(x=x, y=mds_nochange[,3], dim="dim 3", type="mirror")
                
df.fig2 |> ggplot(aes(x=x, y=y, color=type)) +
  geom_point(alpha=0.5) +
  facet_wrap(~dim, scale="free") +
  labs(x="time", y="multidimensional scaling") +
  scale_color_manual(values=c("numerical"="black")) +
  theme_bw() +
  theme(legend.position = "none",
        axis.text = element_text(size=12),
        strip.text = element_text(size = 12, face="bold"),
        axis.title = element_text(size=14, face="bold"))
ggsave("Fig2-ggplot-no-mirror.pdf", width=9, height=3, units="in")


###Figure 3 this is for change point model
par(mfrow=c(1,3))
x=1:tmax
plot(x,df_true_change$V1,xlab = '',main = 'dim 1' ,ylab = 'multidimensional scaling',cex.lab=1.5,cex.axis=1.5
     ,ylim=c(-0.15,0.15))
points(x,mds_change[,1],col='blue')
abline(v=x[tstar],lty=2)
plot(x,mds_change[,2],col='blue',xlab = '',main = 'dim 2',ylab='',ylim=c(-0.02,0.02),cex.lab=1.5,cex.axis=1.5)
points(x,df_true_change$V2)
abline(v=x[tstar],lty=2)
plot(x,mds[,3],col='blue',xlab = '',main = 'dim 3',ylab='',ylim=c(-0.02,0.02),cex.lab=1.5,cex.axis=1.5)
points(x,df_true$V3)
abline(v=x[tstar],lty=2)
mtext('time',side = 1, line = -2, outer = TRUE,cex=1)

df.fig3 = rbind(tibble(x=x, y=df_true_change$V1, dim="dim 1", type="numerical"),
                tibble(x=x, y=mds_change[,1], dim="dim 1", type="mirror"),
                tibble(x=x, y=df_true_change$V2, dim="dim 2", type="numerical"),
                tibble(x=x, y=mds_change[,2], dim="dim 2", type="mirror"),
                tibble(x=x, y=df_true_change$V3, dim="dim 3", type="numerical"),
                tibble(x=x, y=mds_change[,3], dim="dim 3", type="mirror"))
df.fig3 |> ggplot(aes(x=x, y=y, color=type)) +
  geom_point(alpha=0.5) +
  geom_vline(xintercept = x[tstar], linetype="dashed") +
  facet_wrap(~dim, scale="free") +
  labs(x="time", y="multidimensional scaling") +
  scale_color_manual(values=c("numerical"="black", "mirror"="blue")) +
  theme_bw() +
  theme(legend.position = "none",
        axis.text = element_text(size=12),
        strip.text = element_text(size = 12, face="bold"),
        axis.title = element_text(size=14, face="bold"))
ggsave("Fig3-ggplot.pdf", width=9, height=3, units="in")


### fig 3 without mirror 
df.fig3 = rbind(tibble(x=x, y=df_true_change$V1, dim="dim 1", type="numerical"),
#tibble(x=x, y=mds_change[,1], dim="dim 1", type="mirror"),
tibble(x=x, y=df_true_change$V2, dim="dim 2", type="numerical"),
#tibble(x=x, y=mds_change[,2], dim="dim 2", type="mirror"),
tibble(x=x, y=df_true_change$V3, dim="dim 3", type="numerical"))
#tibble(x=x, y=mds_change[,3], dim="dim 3", type="mirror"))
df.fig3 |> ggplot(aes(x=x, y=y, color=type)) +
  geom_point(alpha=0.5) +
  geom_vline(xintercept = x[tstar], linetype="dashed") +
  facet_wrap(~dim, scale="free") +
  labs(x="time", y="multidimensional scaling") +
  scale_color_manual(values=c("numerical"="black")) +
  theme_bw() +
  theme(legend.position = "none",
        axis.text = element_text(size=12),
        strip.text = element_text(size = 12, face="bold"),
        axis.title = element_text(size=14, face="bold"))
ggsave("Fig3-ggplot-no-mirror.pdf", width=9, height=3, units="in")






### Figure 4
#m=12 fixed
plottt1<- ggplot(sm_mse_change_n, aes(x=V4, y=V2)) +
  geom_line() +geom_errorbar(aes(ymin=V1, ymax=V3))+
  ylim(0,0.0115)+
  scale_x_continuous(breaks = nn)+
  labs(y='relative MSE',x='n')+
  theme(axis.text=element_text(size=25),
        axis.title=element_text(size=25,face="bold"))
print(plottt1)

#n=800 fixed
plottt2<- ggplot(sm_mse_change_m, aes(x=V4, y=V2)) +
  geom_line() +geom_errorbar(aes(ymin=V1, ymax=V3))+
  ylim(0,0.0115)+
  scale_x_continuous(breaks = mm)+
  labs(y='',x='m')+
  theme(axis.text=element_text(size=25),
        axis.title=element_text(size=25,face="bold"))
print(plottt2)

# # using facets
df.fig4 = rbind(sm_mse_change_n |> mutate(fix="fixed m"),
                sm_mse_change_m |> mutate(fix="fixed n"))
df.fig4 |> filter(fix=="fixed m") |>
  ggplot(aes(x=V4, y=V2)) +
  geom_line() +
  geom_errorbar(aes(ymin=V1, ymax=V3)) +
  facet_wrap(fix ~ ., scale="free") + #, labeller = label_both) +
  scale_x_continuous(breaks = nn) +
  labs(y="relative MSE", x="n") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text = element_text(size=12),
        strip.text = element_text(size = 12, face="bold"),
        axis.title = element_text(size=14, face="bold"))
ggsave("Fig4a-ggplot.pdf", width=5, height=4, units="in")
df.fig4 |> filter(fix=="fixed n") |>
  ggplot(aes(x=V4, y=V2)) +
  geom_line() +
  geom_errorbar(aes(ymin=V1, ymax=V3)) +
  facet_wrap(fix ~ ., scale="free") + #, labeller = label_both) +
  scale_x_continuous(breaks = mm) +
  labs(y="relative MSE", x="m") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text = element_text(size=12),
        strip.text = element_text(size = 12, face="bold"),
        axis.title = element_text(size=14, face="bold"))
ggsave("Fig4b-ggplot.pdf", width=5, height=4, units="in")


### Figure 5
ggplot(all_mse, aes(x=V4, y=V2)) +
  geom_line() +geom_errorbar(aes(ymin=V1, ymax=V3),width=0.5)+
  facet_grid(m_group ~ n_group)+
  scale_x_continuous(breaks = all_mse$V4)+
  labs(y='',x='',size=10)+
  xlab(expression("MDS embedding dimension " * d * " for the " * (d %->% 1) * "-iso-mirror"))+
  ylab('relative MSE')+
  theme_bw() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12,face="bold"),
        strip.text.x = element_text(size=12, face="bold"),
        strip.text.y = element_text(size=12, face="bold"))
ggsave("Fig5-ggplot.pdf", width=9, height=9, units="in")

##Fig 6
# (a)
days=as.numeric(well_34_result[1,])
isomirror=as.numeric(-well_34_result[2,])

## function for l_inf change point localization
linf_cp=function(t,y,cp){
  n=length(t)
  nl=sum(t<cp)+1
  XL=matrix(1,nrow = nl,ncol=4)
  XL[,4]=0
  XL[,3]=t[1:nl]-cp

  #XL; y[1:nl]

  XL2=XL
  XL2[,2:3]=-XL[,2:3]

  #rbind(XL,XL2); c(y[1:nl],-y[1:nl])

  XR=matrix(1,nrow = n-nl,ncol = 4)
  XR[,3]=0
  XR[,4]=t[(nl+1):n]-cp

  XR2=XR
  XR2[,c(2,4)]=-XR[,c(2,4)]


  X=rbind(XL,XR,XL2,XR2)
  Y=c(y,-y)

  library(lpSolveAPI)
  lprec <- make.lp(0,4)
  set.objfn(lprec,c(1,0,0,0))
  for (i in 1:(nrow(X)) ) {
    add.constraint(lprec, X[i,], ">=", Y[i])
  }

  set.bounds(lprec, lower = c(0,-Inf,-Inf,-Inf), columns = c(1,2,3,4))
  ColNames <- c('Z', "alpha", "bl","br")
  dimnames(lprec)[[2]] <- ColNames
  solve(lprec)
  return(get.variables(lprec))
}

## find the point which minimize the obj func, that is the change point
obf=NULL
for (nk in 2:(length(days)-1)) {
  obf[nk]=linf_cp(days,isomirror,days[nk])[1]
}

## output the piecewise linear fitting corresponding to the change point found from the above step
aa=linf_cp(days,isomirror,days[which(obf==min(obf[-1]))] )

pl=function(x,alpha,bl,br,tstar){
  if(x<tstar){
    y=bl*(x-tstar)+alpha
  } else {
    y=br*(x-tstar)+alpha
  }
}
vpl=Vectorize(pl)

plot(days,isomirror,xlab='time',ylab='iso-mirror')
lines(days,vpl(days,aa[2],aa[3],aa[4],days[23]),lty=2)
abline(v=days[23],lty=2)

## ggplot
df.fig6a = rbind( tibble(x=days, y=isomirror, type="numerical"),
                  tibble(x=days, y=vpl(days,aa[2],aa[3],aa[4],days[23]), type="mirror"))
df.fig6a |> ggplot(aes(x=x, y=y, color=type)) +
  geom_point(data=df.fig6a |> filter(type=="numerical"), alpha=0.7) +
  geom_line(data=df.fig6a |> filter(type=="mirror"), linetype="dashed") +
  geom_vline(xintercept = days[23], linetype="dashed") +
  labs(x="day", y="iso-mirror") +
  scale_color_manual(values=c("numerical"="black", "mirror"="blue")) +
  theme_bw() +
  theme(legend.position = "none",
        axis.text = element_text(size=12),
        strip.text = element_text(size = 12, face="bold"),
        axis.title = element_text(size=14, face="bold"))
ggsave("Fig6a-ggplot.pdf", width=5, height=4, units="in")

# (b)
days=c(30,as.numeric(well_34_controll[1,]))
frobdiff=as.numeric(well_34_controll[2,])

# Generate labels automatically
labels <- paste(days[-length(days)], days[-1], sep=":")
labels
# Create the plot with custom x-axis
plot(days[-1], # Exclude the last point for plotting to match the labels
     frobdiff, # Exclude the last point for y values as well
     type = "b",
     xaxt = 'n', # No default x-axis
     ylab = expression("||" * A[i] - A[i-1] * "||"["F"]),
     xlab = ""
)

# Now, add the main x-axis label
mtext("day", side = 1, line = 3.4)

# Add the custom axis labels
axis(1,
     at = days[-1], # Position labels between ticks
     labels = labels,
     padj = 1,
     cex.axis = 0.8,
     tck = 0.02,
     las=2
)

# ggplot
labels2 <- paste(sprintf("%3d",days[-length(days)]), sprintf("%3d",days[-1]), sep=":")
df.fig6b = tibble(x=days[-1], y=frobdiff, x2=labels2)
df.fig6b |> ggplot(aes(x=x2, y=y, group=1)) +
  geom_point() + geom_line() +
  labs(x="day", y=expression("||" * A[i] - A[i-1] * "||"["F"])) +
  # scale_x_discrete(breaks = seq(1,nrow(df.fig6b),2)) +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle=90, vjust=.5, hjust=1, size=8),
        axis.text = element_text(size=12),
        strip.text = element_text(size = 12, face="bold"),
        axis.title = element_text(size=14, face="bold"))
ggsave("Fig6b-ggplot.pdf", width=5, height=4, units="in")
