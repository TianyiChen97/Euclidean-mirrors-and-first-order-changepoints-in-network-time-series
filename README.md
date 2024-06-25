This repository contains code supporting the text "Euclidean mirrors and first-order changepoints in network time series" by _Tianyi Chen, Zachary Lubberts, Avanti Athreya, Youngser Park, Carey E. Priebe_.  
https://arxiv.org/abs/2405.11111

## Abstract: 
We describe a model for a class of network time series whose evolution is governed by an underlying stochastic process, known as the latent position process, in which network evolution can be represented in Euclidean space by a curve, called the Euclidean mirror. We define the notion of a first-order changepoint for a time series of networks, and construct a family of latent position process networks with underlying first-order changepoints. We prove that a spectral estimate of the associated Euclidean mirror localizes these changepoints, even when the graph distribution evolves continuously, but at a rate that changes. Simulated and real data examples on organoid networks show that this localization identifies empirically significant shifts in network evolution.  
## Data

Data including the simulation results and the real data. We make our numerical experiement results (Section 3.3) avaliable and contians all of them in Figure_data.RData. The numerical experiement results are summarized in the matrices named with form: "sm_mse_...". As for the real data only the iso-mirror results of the time series of brain organoid connectome is avaliable, it is stored in "well_34_controll". 

## Code
`Simulation.R` includes all the code needed to reproduce the numerical experiement result "sm_mse_..." for Figs 5 and 6. Then `All-figures+real_data.R` includes the code to plot the figures 5 6 using "sm_mse_..."s. Also `All-figures+real_data.R` incldues all the code we need to reproduce the figs 2 3 4 7 (simulation for Model 2.25 and the real data example). 
