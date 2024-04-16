#setwd("~/Dropbox/Mirror/Paper/Tianyi")

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
