This repository contains the code for replication of the results in the paper "Efficient Nonparametric Inference for Mediation Analysis with Nonignorable Missing Confounders".

- The sub-directory `main_simulation`  contains the codes that produce the simulation results reported in the main text. Files beginning with `1.x` correspond to DGP1, and those beginning with `2.x` correspond to DGP2. These codes generate Tables 2 and 3 and Figure 2 in the main text, and Figure G.1 in the the Supplementary Material.
- The subdirectory `supp_dml` contains the codes that produce the simulation implementing debiased machine learning with higher-dimensional covariates in Supplementary Appendix G.2. These codes generate Tables G.1 in the the Supplementary Material.
- The subdirectory `supp_MAR` contains the codes that produce the simulation with missing-at-random mechanisms in Supplementary Appendix G.3. These codes generate Tables G.2 and G.3 and Figure G.2 in the the Supplementary Material.
- The subdirectory `supp_trim` contains the codes that produce the simulation implementing the trimming approach in Supplementary Appendix G.4. These codes generate Figure G.3 in the the Supplementary Material.
- The subdirectory `realdata` contains the code that implements our empirical application in the main text. Due to copyright restrictions, we do not include the original data in this folder. Readers can download the original data by applying at https://www.isss.pku.edu.cn/cfps/en/.  All other codes are included.

