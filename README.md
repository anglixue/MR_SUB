# Causal association between substance use behaviors and common diseases
This is the analysis code repository for Mendelian Randomisation (MR) study between substance use behaviors (SUB) and common diseases, as part of the manuscript "Unravelling the complex causal effects of substance use behaviours on common diseases"

In this study, we systematically evaluated the causal relationship between 7 SUB traits and a wide range of health outcomes. We compared 11 MR methods both in simulation and real data. We also dissect the putative dosage-dependent effects of coffee and tea intake by integrating stratified regression, genetic correlation, and MR analyses.

Scripts are listed by the order in the methods section of the manuscript:

1. GWAS of 7 SUB traits
2. Genetic correlation between SUB and common diseases
3. Mendelian Randomisation (MR) between SUB and common diseases
4. Validation using data from published studies
5. Dosage-dependent analysis
6. mtCOJO adjustment for covariates
7. Bi-directional analysis
8. Simulation

# GSMR2 software
The improved GSMR method (GSMR2) is formally released along with the paper. 

C++ version: [GCTA-GSMR](https://yanglab.westlake.edu.cn/software/gcta/index.html#GSMR) within the GCTA software 

R package: [GSMR2](https://github.com/jianyanglab/gsmr2) GitHub repository.

# Download GWAS summary statistics
We have provided GWAS summary statistics for seven SUB traits (tobacco smoking, alcohol consumption, coffee and tea intake). Please find the download link [here](https://yanglab.westlake.edu.cn/pub_data.html).

# Citation
Angli Xue, Zhihong Zhu, Huanwei Wang, Longda Jiang, Peter M. Visscher, Jian Zeng, Jian Yang. Unravelling the complex causal effects of substance use behaviours on common diseases. ***Communications Medicine***. 2024. [[Full text](https://doi.org/10.1038/s43856-024-00473-3)] [[Preprint](https://www.researchsquare.com/article/rs-3465061/v1)]

For questions, please email us at Angli Xue (a.xue@garvan.org.au) or Jian Yang (jian.yang@westlake.edu.cn)
