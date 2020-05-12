To get SS values for each tumor-
Call_NRF2p53_SSlastval.m, which uses TCGA_FCs_medians.csv and calls ‘NRF2p53noH2O2_DDE.m’ —> runs simulations for 5000 minutes to reach steady state, and then saves the last value for every species per tumor in ‘SS_lastvalues.mat’

For 100 simulations with increased ROS generation rate per tumor, using SS values found above as initial conditions, each simulation is 1000 minutes (increased ROS generation rate for first 2hrs)-
Call_NRF2p53_TCGA100sims.m, which uses SS_lastvalues.mat and calls ‘NRF2p53plusH2O2_DDE.m’ —>  saves all species at all times as a mat file per tumor