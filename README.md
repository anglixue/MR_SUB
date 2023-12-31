# Causal association between substance use behaviors and common diseases
This is the analysis code repository for Mendelian Randomisation (MR) study between substance use behaviors (SUB) and common diseases, as part of the manuscript "Dissecting the complexity of causal effect estimation for substance use behaviours on common diseases"

In this study, we systematically evaluated the causal relationship between 7 SUB traits and a wide range of health outcomes. We compared 11 MR methods both in simulation and real data. We also dissect the putative dosage-dependent effects of coffee and tea intake by integrating stratified regression, genetic correlation, and MR analyses.

The improved GSMR method (GSMR2) is formally released along with the paper. We have implemented the C++ version [GCTA-GSMR](https://yanglab.westlake.edu.cn/software/gcta/index.html#GSMR) within the GCTA software and also updated the [gsmr](https://yanglab.westlake.edu.cn/software/gsmr/) R package. The GitHub repository of GSMR can be found [here](https://github.com/jianyanglab/gsmr2).

Scripts are listed by the order in the methods section of the manuscript:

1. GWAS of 7 SUB traits
2. Genetic correlation between SUB and common diseases
3. Mendelian Randomisation (MR) between SUB and common diseases
4. Validation using data from published studies
5. Dosage-dependent analysis
6. mtCOJO adjustment for covariates
7. Bi-directional analysis
8. Simulation


# Citation
Angli Xue, Zhihong Zhu, Huanwei Wang, Longda Jiang, Peter M. Visscher, Jian Zeng, Jian Yang. Unravelling the complex causal effects of substance use behaviours on common diseases. ***Communications Medicine***. 2023. [preprint](https://www.researchsquare.com/article/rs-3465061/v1)

For questions, please email us at Angli Xue (a.xue@garvan.org.au) or Jian Yang (jian.yang@westlake.edu.cn)
