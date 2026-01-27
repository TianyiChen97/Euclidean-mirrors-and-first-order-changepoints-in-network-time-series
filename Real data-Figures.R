# --- Libraries ---
pacman::p_load(segmented, igraph, RSpectra, locfit, tidyverse, doParallel, broom, vegan, Matrix)
library(iGraphMatch)
library(transport)
library(gridExtra) 
library(grid) 
library(lpSolveAPI)
registerDoParallel(detectCores()-1)

# --- Embedding and Geometry Functions ---

full.ase <- function(A, d, diagaug=TRUE, doptr=FALSE) {
  # Performs Adjacency Spectral Embedding on a network matrix A with embedding dimension d.
  require(irlba)
  if (doptr) {
    g <- ptr(A)
    A <- g[]
  } else {
    A <- A[]
  }
  if (diagaug) {
    diag(A) <- rowSums(A) / (nrow(A)-1)
  }
  A.svd <- svds(A,k=d)
  Xhat <- A.svd$u %*% diag(sqrt(A.svd$d))
  Xhat.R <- NULL
  if (!isSymmetric(A)) {
    Xhat.R <- A.svd$v %*% diag(sqrt(A.svd$d))
  }
  return(list(eval=A.svd$d, Xhat=Matrix(Xhat), Xhat.R=Xhat.R))
}

procrustes2 <- function(X, Y) {
  # Computes the Procrustes transformation to optimally align matrix X to matrix Y.
  tmp <- t(X) %*% Y
  tmp.svd <- svd(tmp)
  W <- tmp.svd$u %*% t(tmp.svd$v)
  newX <- X %*% W
  return(list(newX = newX, error = norm(newX-Y, type="F"), W = W))
}

getD <- function(Xlist, k=0, etype="proc") {
  # Computes a pairwise distance matrix between embeddings in a list, optionally using Procrustes alignment.
  m <- length(Xlist)
  n <- nrow(Xlist[[1]])
  
  if (k==0) {
    ind <- 1:n
  } else {
    ind <- which(Yhat==k)
  }
  comb <- combn(m,2)
  Dout <- foreach (k = 1:ncol(comb), .combine='rbind') %dopar% {
    i <- comb[1,k]
    j <- comb[2,k]
    
    if (etype == "proc") {
      Xhati <- Xlist[[i]][ind,] 
      Xhatj <- Xlist[[j]][ind,]
      proc <- procrustes2(as.matrix(Xhati), as.matrix(Xhatj))
      Xhati <- Xhati %*% proc$W
    } else {
      Xhati <- Xlist[[i]][ind,] 
      Xhatj <- Xlist[[j]][ind,]
    }
    
    D <- norm(Xhati - Xhatj, type="2")^2/n
    tibble(i=i, j=j, D=D)
  }
  D2 <- matrix(0,m,m)
  D2[t(comb)] <- Dout$D
  D2 <- (D2 + t(D2)) / 1
  D2 <- sqrt(D2)
  D2
}

doMDS <- function(D, doplot=TRUE) {
  # Applies Classical Multidimensional Scaling (CMDS) to the distance matrix and optionally plots results.
  tmax <- m <- nrow(D)
  mds <- cmdscale(D, m-1)
  
  df.mds <- tibble(
    ind = 1:tmax,
    time = sprintf("%2d", 1:tmax),
    x = mds[, 1],  
    y = if (ncol(mds) >= 2) mds[, 2] else NA, 
    z = if (ncol(mds) >= 3) mds[, 3] else NA,
    w = if (ncol(mds) >= 4) mds[, 4] else NA
  )
  
  if (doplot) {
    matrix_name <- deparse(substitute(D))
    big_title <- textGrob(paste("CMDS on", matrix_name), gp=gpar(fontsize=18, fontface="bold"))
    
    sd_data <- data.frame(
      dimension = 1:ncol(mds), 
      column_sd = apply(mds, 2, sd)
    )
    
    p0 <- ggplot(sd_data, aes(x = dimension, y = column_sd)) +
      geom_point() + geom_line() +
      labs(x = "dimension", y = "column stdev") + theme_minimal()
    
    p1 <- df.mds %>% ggplot(aes(x=ind, y=x, color=time, group=1)) +
      geom_point(size=3) + geom_line() +
      geom_vline(xintercept = tstar, linetype="dashed") +
      theme(legend.position = "none") + labs(x="time", y="mds1")
    
    p2 <- df.mds %>% ggplot(aes(x=ind, y=y, color=time, group=1)) +
      geom_point(size=3) + geom_line() +
      geom_vline(xintercept = tstar, linetype="dashed") +
      theme(legend.position = "none") + labs(x="time", y="mds2")
    
    p3 <- df.mds %>% ggplot(aes(x=ind, y=z, color=time, group=1)) +
      geom_point(size=3) + geom_line() +
      geom_vline(xintercept = tstar, linetype="dashed") +
      theme(legend.position = "none") + labs(x="time", y="mds3")
    
    p4 <- df.mds %>% ggplot(aes(x=x, y=y, color=time)) +
      geom_point(size=3) + geom_label_repel(aes(label=time), size=2) +
      theme(legend.position = "none") + labs(x="mds1", y="mds2")
    
    grid.arrange(big_title, arrangeGrob(p0, p1, p2, p3, ncol=4), ncol=1, heights=c(0.3, 4))    
  }
  return(list(mds=mds, df.mds=df.mds))
}

doIso <- function(mds, mdsd=2, isod=1, doplot=F) {
  # Applies ISOMAP to the MDS results to estimate intrinsic geometry and optionally plots residuals.
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
  iso2 <- tibble(iso=iso[,1]) %>% mutate(i=1:nrow(mds), ind=1:nrow(mds), time=1:nrow(mds), knn=knn)
  df.iso <- rbind(df.iso, cbind(iso2, mdsd=mdsd))
  df.iso <- df.iso %>% group_by(mdsd) %>% mutate(iso = if(iso[1] > 0) {-iso} else {iso}) %>% ungroup()
  
  if (doplot) {
    p <- df.iso %>% filter(mdsd==mdsd) %>%
      ggplot(aes(x = 1: length(iso) , y=iso, color=time, group=1)) +
      geom_point(size=3) + geom_line() +
      theme(legend.position = "none") +
      labs(x="time", y="isomap embedding") +
      theme(axis.text.x=element_text(hjust=0.7))
    print(p)
    
    df.isok <- df.iso %>% filter(mdsd==mdsd) 
    row.names(df.isok) <- df.isok$time
    fit <- lm(iso ~ i, data=df.isok)
    myfor <- augment(fit)
    myfor2 <- myfor %>% mutate(date=.rownames, ranks = rank(.sigma), mycol=sprintf("%2d",rank(.fitted)))
    p <- myfor2 %>% ggplot(aes(.fitted, .resid)) +
      geom_point(aes(color=mycol)) +
      geom_hline(yintercept = 0, linetype="dashed", color="grey") +
      geom_smooth(method="loess", se=FALSE) +
      labs(x="Fitted Values", y="Residuals") +
      theme(legend.position = "none", axis.title = element_text(size=14, face="bold"))
    p <- p + geom_label_repel(aes(label=date), data=myfor2 %>% filter(ranks %in% 1:3))
    print(p)
  }
  return(df.iso)
}

# --- Network Metric Functions ---

calculate_edge_density <- function(mat) {
  # Calculates the density of edges in a directed graph, excluding self-loops.
  adj <- (mat > 0) 
  n <- nrow(adj)
  actual_edges <- sum(adj) - sum(diag(adj))
  possible_edges <- n * (n - 1)
  return(actual_edges / possible_edges)
}

calculate_reciprocity <- function(mat) {
  # Calculates the reciprocity (ratio of mutual edges to total edges) of the graph.
  adj <- (mat > 0) + 0
  diag(adj) <- 0
  total_edges <- sum(adj)
  if (total_edges == 0) return(0)
  mutual_edges <- sum(adj * t(adj))
  return(mutual_edges / total_edges)
}

calculate_path_length <- function(mat) {
  # Calculates the average geodesic path length for reachable pairs in the graph.
  adj <- (mat > 0) + 0 
  g <- graph_from_adjacency_matrix(adj, mode = "directed", diag = FALSE)
  avg_path <- mean_distance(g, directed = TRUE, unconnected = TRUE)
  return(avg_path)
}

check_weights_is_nine <- function(mat) {
  # Verifies if all non-zero weights in a matrix are exactly 9.
  non_zero_values <- mat[mat != 0]
  if (length(non_zero_values) == 0) return(TRUE)
  return(all(non_zero_values == 9))
}

# --- Change Point Detection Functions ---

linf_cp=function(t,y,cp){
  # Solves a linear programming problem to localize a change point based on L-infinity norm.
  n=length(t)
  nl=sum(t<cp)+1
  XL=matrix(1,nrow = nl,ncol=4)
  XL[,4]=0
  XL[,3]=t[1:nl]-cp
  
  XL2=XL
  XL2[,2:3]=-XL[,2:3]
  
  XR=matrix(1,nrow = n-nl,ncol = 4)
  XR[,3]=0
  XR[,4]=t[(nl+1):n]-cp
  
  XR2=XR
  XR2[,c(2,4)]=-XR[,c(2,4)]
  
  X=rbind(XL,XR,XL2,XR2)
  Y=c(y,-y)
  
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

pl=function(x,alpha,bl,br,tstar){
  # Piecewise linear function that calculates a value based on a change point tstar.
  if(x<tstar){
    y=bl*(x-tstar)+alpha
  } else {
    y=br*(x-tstar)+alpha
  }
}
vpl=Vectorize(pl)










## Data loading 
days = c(30, 41, 44, 48, 51, 59, 65, 66, 72, 76,
         79, 86, 93, 99, 107, 113, 125, 128, 135, 
         136, 142, 149, 156, 163, 170, 177, 191,
         228, 232, 246)

# Load data (ensure this path is correct)
load("reproduce/glist_matrices.RData")


### This is to check that all weights in all graphs are only nine, 
### thus the graphs are effectively unweighted as claimed in the paper.
weight_check_results <- sapply(glist_matrices, check_weights_is_nine)




## 3. Summary Statistics

# Calculate metrics across all matrices
edge_densities <- sapply(glist_matrices, calculate_edge_density)
reciprocities <- sapply(glist_matrices, calculate_reciprocity)
path_lengths <- sapply(glist_matrices, calculate_path_length)

edge_densities + path_lengths ## confirm the special structure of these graphs

# Combine into DataFrame
df_all <- tibble(
  day = days,
  `average path length` = path_lengths, 
  `edge density` = edge_densities,      
  reciprocity = reciprocities
) |> 
  pivot_longer(cols = -day, names_to = "metric", values_to = "value")

# Plot combined statistics
p_combined <- df_all |> 
  ggplot(aes(x = day, y = value)) +
  geom_point(alpha = 0.7, color = "black") +
  geom_vline(xintercept = days[23], linetype = "dashed") +
  facet_wrap(~metric, ncol = 1, scales = "free_y", strip.position = "left") +
  labs(x = "Day", y = NULL) + 
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 12),
    strip.text = element_text(size = 12, face = "bold"),
    strip.background = element_blank(), 
    strip.placement = "outside",        
    axis.title = element_text(size = 14, face = "bold")
  )

print(p_combined)



## 4. Frobenius Control Chart

# Validate weights
weight_check_results <- sapply(glist_matrices, check_weights_is_nine)

# Calculate Frobenius differences
frobdiff <- numeric(length(glist_matrices))
for (i in 2:length(glist_matrices)) {
  diff_matrix <- glist_matrices[[i]] - glist_matrices[[i-1]]
  frobdiff[i] <- norm(diff_matrix, type = "F")
}

# Plot Control Chart
labels2 <- paste(sprintf("%3d",days[-length(days)]), sprintf("%3d",days[-1]), sep=":")
df.fig6b = tibble(x=days[-1], y=frobdiff[-1], x2=labels2)

df.fig6b |> ggplot(aes(x=x2, y=y, group=1)) +
  geom_point() + geom_line() +
  labs(x="day", y=expression("||" * A[i] - A[i-1] * "||"["F"])) +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle=90, vjust=.5, hjust=1, size=8),
        axis.text = element_text(size=12),
        strip.text = element_text(size = 12, face="bold"),
        axis.title = element_text(size=14, face="bold"))


## 5. ISOmirror
# Perform ASE on list
Xhat_list <- lapply(glist_matrices, function(g) {
  g = (g > 0) + 0 
  ase <- full.ase(g, 2)
  Xhat <- cbind(as.matrix(ase$Xhat), ase$Xhat.R)
  return(Xhat)
})

# Calculate Distance and MDS
D = getD(Xhat_list)
MDS = doMDS(D, doplot = F)
plot(-MDS$df.mds$x, -MDS$df.mds$y)


## Plot the Scree Plot

scree_data <- tibble(
  dimension = 1:ncol(MDS$mds),
  std_dev = apply(MDS$mds, 2, sd)
)
scree_data |> 
  ggplot(aes(x = dimension, y = std_dev)) +
  geom_point(size = 3) +            
  scale_x_continuous(breaks = scree_data$dimension) + 
  labs(
    title = "Scree Plot fo CMDS",
    x = "Dimension",
    y = "Standard Deviation"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title = element_text(face = "bold"),
    panel.grid.minor.x = element_blank() # Clean up vertical grid lines
  )

## Choose the elbow to be 2

# Perform ISOMAP
iso = doIso(MDS$mds, mdsd = 2)
isomirror = iso$iso
plot(days, iso$iso)

# Find Change Point minimizing objective function
obf=NULL
for (nk in 2:(length(days)-1)) {
  obf[nk]=linf_cp(days,isomirror,days[nk])[1]
}

# Generate Piecewise Linear Fit
aa=linf_cp(days,isomirror,days[which(obf==min(obf[-1]))] )

# Plot ISOMirror with Fit
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


