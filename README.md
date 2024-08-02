# longitudinal snSMART

## Simulation studies `R` code for paper "Integrating inverse probability weighted external control data into a Bayesian longitudinal small sample, sequential, multiple assignment randomized trial (snSMART) design"

Our proposed approach represents an enhancement of the robust MAC-snSMART method (Wang et al. (2023)), offering major improvements to the existing snSMART methods: 1) enabling the analysis of longitudinal data, 2) allowing for the inclusion of patient baseline characteristics, 3) handling missing data through multiple imputation, 4) limiting the heterogeneity between data sources through propensity score (PS), and 5) addressing the possible stage-wise treatment effect non-exchangeability caused by unexpected disease progression or ineffective washout period in-between stages. These enhancements significantly broaden the applicability of the innovative snSMART design and increase efficiency in rare disease drug development. 

Six files are provided:
- `external_data_simulation.R` - simulate sample external control data;
- `current_data_simulation.R` - simulate sample current trial data;
- `BLPM_test.bug` - `JAGS` code of our prposed BLPM method;
- `BLPM_robust_test.bug` - `JAGS` code of our proposed robust BLPM method;
- `JointStageBayes_mixture.bug` - `JAGS` code of the Bayesian joint stage model proposed by [Fang Fang](https://www.tandfonline.com/doi/abs/10.1080/19466315.2022.2118162);
- `test.R` - the original `R` code used to generate the simulation studies results presented in the `Result` section of our paper.

To conduct simultion studies yourself, please put all files under one folder, set working directory to that folder, and run `current_data_simulation.R` and `external_data_simulation.R` to generate simulated datasets. Then, run `test.R` to generate rMSE, bias, credible interval width, and coverage rate of different simulation scenarios.

Contact: sidiwang@umich.edu
Reference: 
Sidi Wang, Kelley M. Kidwell, Satrajit Roychoudhury, Dynamic Enrichment of Bayesian Small-Sample, Sequential, Multiple Assignment Randomized Trial Design Using Natural History Data: A Case Study from Duchenne Muscular Dystrophy, Biometrics, Volume 79, Issue 4, December 2023, Pages 3612–3623, https://doi-org.proxy.lib.umich.edu/10.1111/biom.13887 
