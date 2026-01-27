This repository contains code supporting the text "Euclidean mirrors and first-order changepoints in network time series" by _Tianyi Chen, Zachary Lubberts, Avanti Athreya, Youngser Park, Carey E. Priebe_.  
https://arxiv.org/abs/2405.11111

## Abstract: 
We describe a model for a class of network time series whose evolution is governed by an underlying stochastic process, known as the latent position process, in which network evolution can be represented in Euclidean space by a curve, called the Euclidean mirror. We define the notion of a first-order changepoint for a time series of networks, and construct a family of latent position process networks with underlying first-order changepoints. We prove that a spectral estimate of the associated Euclidean mirror localizes these changepoints, even when the graph distribution evolves continuously, but at a rate that changes. Simulated and real data examples on organoid networks show that this localization identifies empirically significant shifts in network evolution.  
## Data

Data including the simulation results and the real data. We make our numerical experiment results (Section 3.3) available and contains all of them in `Figure_data.RData`. The numerical experiment results are summarized in the matrices named with form: "sm_mse_...". As for the real data the time series of graphs for brain organoid is available, it is stored in "glist_matrices.RData". 

## Code
`Simulation.R` includes all the code needed to reproduce the numerical experiment for localization performance result "sm_mse_...".
Then `Simulation-Figures.R` includes the code to plot the simulation results figures using "sm_mse_..."s. Also it includes the code we need to reproduce the figs 2 3 4 7 (simulation for Model with changepoint or without changepoint and the real data example). 
Then `Real data-Figures.R` includes the codes that read in the TSG with 44 nodes and 30 time points directed graphs and then perform summary statistics analysis, Frobenius norm control chart analysis and isomirror analysis.  