pacman::p_load(segmented, igraph, irlba, locfit, tidyverse, doParallel, broom, vegan, Matrix)

registerDoParallel(detectCores()-1)

##generate RDPG from given latent position X
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

##ASE for a network A with embedding dimension d
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

## Generate latent positions process from Model3 Constant drift with jump probability first-order changepoint 
doSim <- function(n=300, tmax=40, delta=0.1, p=0.4, q=0.9, tstar=20,startpoint=0.1)
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
  Xt <- Xt+startpoint
  
  df <- tibble(time=1:tmax) %>%
    mutate(Xt = map(time, function(x) matrix(Xt[,x],n,1)  )) %>%
    mutate(g = map(Xt, ~rdpg.sample(.)))
  
  df
}

## Procruste/i.e gain W in the estimated d_MV distance. Note in our case latent positions are 1d so this is trivial with w=1or -1
procrustes2 <- function(X, Y) {
  tmp <- t(X) %*% Y
  tmp.svd <- svd(tmp)
  W <- tmp.svd$u %*% t(tmp.svd$v)
  newX <- X %*% W
  return(list(newX = newX, error = norm(newX-Y, type="F"), W = W))
}

## Get distance matrix
getD <- function(Xlist, k=0, etype="proc") {
  m <- length(Xlist)
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

## Apply CMDS on distance matrix 
doMDS <- function(D, doplot=TRUE)
{
  tmax <- m <- nrow(D)
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

#apply ISOMAP on the CMDS result with chosen dimension mdsd from CMDS step defaultly it always embeds to 1 
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

## This is another slope change point algorithm called segmented that is not used in the paper 
break_point_dection=function(D,k){
  tmax <- nrow(D)
  df.mds <- doMDS(D,doplot = F)
  mds <- df.mds$mds
  df.iso <- doIso(mds, mdsd=k)
  x=1:tmax/tmax
  y1=df.iso$iso
  os1<-segmented(lm(y1~x),psi=c(0.2))
  result=as.data.frame(matrix(0,1,3))
  result[1,1:4]=c((os1$psi[1,2]-1.95*os1$psi[1,3])*tmax,(os1$psi[1,2])*tmax,(os1$psi[1,2]+1.95*os1$psi[1,3])*tmax,os1$psi[1,3]*tmax)
  return(result)
  ## this result returns you the confidence interval of point estimation and the standard deviation
}


## Implementation of the 3rd step in Algorithm 2 by recasting it as a linear programming problem.
## This function will return the objective function value Sk in the paper for a given change point t 
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



#### fix n change m

set.seed(2)
df.mds=NULL
n=800
p <- 0.4
q <- 0.3
nmc=2000
mm=c(16,24,32,40)


pacman::p_load("doParallel")
registerDoParallel(detectCores()-4)

out <- foreach (mc = 1:nmc) %do% {
  tmp1 <- tmp2 <- rep(0,length(mm))
  for(i in 1:length(mm)){
    m <- tmax <- mm[i]
    delta <- (1-0.1)/tmax
    tstar <- tmax/2
    df <- doSim(n,tmax,delta,p,q,tstar,startpoint = 0.1)
    df <- df %>% mutate(Xhat = map(g, function(x) full.ase(x,2)$Xhat[,1,drop=F]))
    D2 <- getD(df$Xhat) 
    df.mds <- doMDS(D2,doplot = F)
    mds <- df.mds$mds
    df.iso <- doIso(mds, mdsd=1)
    obf=NULL
    for (nk in 2:(tmax-1)) { ## find the point which minimize the obj func Sk, that is the change point 
      obf[nk]=linf_cp(1:tmax,df.iso$iso,nk)[1]
    }
    ecp=min(which(obf==min(obf[-1])))
    tmp1[i]=(ecp-tstar)/tmax
    bpd=break_point_dection(D2,1)
    round_bpd=which(abs(bpd[1,2]-1:tmax)==min(abs(bpd[1,2]-1:tmax)))
    tmp2[i]=(round_bpd-tstar)/tmax
  }
  list(tmp1, tmp2)
}

maehat1=Reduce('cbind', lapply(out, "[[", 1))
#maehat2=Reduce('cbind', lapply(out, "[[", 2))

maehat1[1,]

sm_mae1=as.data.frame(matrix(0,length(mm),4))
sm_mae1[,2]=apply(abs(maehat1), 1, mean)
sm_mae1[,1]=apply(abs(maehat1), 1, mean)-apply(abs(maehat1), 1, sd)*1.96/sqrt(nmc)  
sm_mae1[,3]=apply(abs(maehat1), 1, mean)+apply(abs(maehat1), 1, sd)*1.96/sqrt(nmc) 
sm_mae1

msehat1=maehat1^2
sm_mse1=as.data.frame(matrix(0,length(mm),4))
sm_mse1[,2]=apply(msehat1, 1, mean)
sm_mse1[,1]=apply(msehat1, 1, mean)-apply(msehat1, 1, sd)*1.96/sqrt(nmc)  
sm_mse1[,3]=apply(msehat1, 1, mean)+apply(msehat1, 1, sd)*1.96/sqrt(nmc) 

sm_mse1[,4]=mm
plottt2<- ggplot(sm_mse_change_m, aes(x=V4, y=V2)) + 
  geom_line() +geom_errorbar(aes(ymin=V1, ymax=V3))+
  ylim(0,0.0115)+
  scale_x_continuous(breaks = mm)+labs(y='',x='m')+
  theme(axis.text=element_text(size=25),axis.title=element_text(size=25,face="bold"))
print(plottt2)


###############   fix m change n 
set.seed(3)
df.mds=NULL
tmax <- m <- 12
p <- 0.4
q <- 0.3
delta <- (1-0.1)/tmax
tstar <- tmax/2
nmc=2000
nn=c(200,800,1600)
msehat1=msehat2=matrix(0,length(nn),nmc)

pacman::p_load("doParallel")
registerDoParallel(detectCores()-4)


out=foreach (mc=1:nmc) %do% {
  tmp1 <- tmp2 <- rep(0,length(nn))
  for(i in 1:length(nn)){
    n=nn[i]
    df <- doSim(nn[i],tmax,delta,p,q,tstar,startpoint = 0.1)
    df <- df %>% mutate(Xhat = map(g, function(x) full.ase(x,2)$Xhat[,1,drop=F]))
    D2 <- getD(df$Xhat) 
    df.mds <- doMDS(D2,doplot = F)
    mds <- df.mds$mds
    df.iso <- doIso(mds, mdsd=1)
    ## l-inf regression change point localization
    obf=NULL
    for (nk in 2:(tmax-1)) {
      obf[nk]=linf_cp(1:tmax,df.iso$iso,nk)[1]
    }
    ecp=min(which(obf==min(obf[-1])))
    tmp1[i]=(ecp-tstar)/tmax
    ##the break point algorithm but i change it so that it only detect the change point that is one of the time point 
    bpd=break_point_dection(D2,1)
    round_bpd=which(abs(bpd[1,2]-1:tmax)==min(abs(bpd[1,2]-1:tmax)))
    tmp2[i]=(round_bpd-tstar)/tmax
  }
  list(tmp1,tmp2)
}

msehat1=Reduce('cbind', lapply(out, "[[", 1))

## summary matrix
sm_mse1=as.data.frame(matrix(0,length(nn),4))
sm_mse1[,2]=apply(abs(msehat1)^2, 1, mean)
sm_mse1[,1]=apply(abs(msehat1)^2, 1, mean)-apply(abs(msehat1)^2, 1, sd)*1.96/sqrt(nmc)  
sm_mse1[,3]=apply(abs(msehat1)^2, 1, mean)+apply(abs(msehat1)^2, 1, sd)*1.96/sqrt(nmc) 
sm_mse1[,4]=nn

## make the figure in paper 
plottt1<- ggplot(sm_mse_change_n, aes(x=V4, y=V2)) + 
  geom_line() +geom_errorbar(aes(ymin=V1, ymax=V3))+
  ylim(0,0.0115)+
  scale_x_continuous(breaks = nn)+
  labs(y='MSE',x='n')+
  theme(axis.text=element_text(size=25),axis.title=element_text(size=25,face="bold"))
print(plottt1)





############# fix m n change CMDS embedding d 
set.seed(3)
df.mds=NULL
#n=800 best
tmax <- m <- 12
p <- 0.4
q <- 0.3
delta <- (1-0.1)/tmax
tstar <- tmax/2
#nmc = 150 best 
nmc=2000
n=800
msehat1=msehat2=matrix(0,10,nmc)

pacman::p_load("doParallel")
registerDoParallel(detectCores()-4)

tmp1=NULL
msehat <- foreach (mc = 1:nmc, .combine='cbind') %do% {
  df <- doSim(n,tmax,delta,p,q,tstar,startpoint = 0.1)
  df <- df %>% mutate(Xhat = map(g, function(x) full.ase(x,2)$Xhat[,1,drop=F]))
  D2 <- getD(df$Xhat) 
  df.mds <- doMDS(D2,doplot = F)
  mds <- df.mds$mds
  for(dd in 1:10){
    df.iso <- doIso(mds, mdsd=dd)
    ## l-inf regression change point localization
    obf=NULL
    for (nk in 2:(tmax-1)) {
      obf[nk]=linf_cp(1:tmax,df.iso$iso,nk)[1]
    }
    ecp=min(which(obf==min(obf[-1])))
    tmp1[dd]=(ecp-tstar)/tmax
  }
  tmp1
}

msehat

sm_mse=as.data.frame(matrix(0,10,4))
sm_mse[,2]=apply(abs(msehat)^2, 1, mean)
sm_mse[,1]=apply(abs(msehat)^2, 1, mean)-apply(abs(msehat)^2, 1, sd)*1.96/sqrt(nmc)  
sm_mse[,3]=apply(abs(msehat)^2, 1, mean)+apply(abs(msehat)^2, 1, sd)*1.96/sqrt(nmc) 
sm_mse[,4]=1:10

sm_mse

plottt1<- ggplot(sm_mse[1:5,], aes(x=V4, y=V2)) + 
  geom_line() +geom_errorbar(aes(ymin=V1, ymax=V3))+
  ylim(-0.0001,0.01)+
  scale_x_continuous(breaks = sm_mse$V4)+labs(y=expression(paste('MSE(',frac(hat(t),m),')')),x='MDS embedding dimension d for the $(d \to 1)$-iso-mirror'
                                              ,title=paste('p=',p,'q=',q, 't*=',tstar, 'm=',tmax,'n=',n),size=15  )+
  theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold"))
print(plottt1)









