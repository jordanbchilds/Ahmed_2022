Data analysis to support the manuscript "Estimating the allele frequency threshold of the pathogenic mtDNA variant m.3243A>G tolerated by human myofibres" by Ahmed et al. (2022)

* Open R  
* Install the required packages: ```install.packages(c("colortools","data.table","fitdistrplus","gmp","mclust","Hotelling","rootSolve","beanplot"))```
* Set this directory as working directory
* Open and run the code in analyse.R (e.g. type ```source("analyse.R")``` at the prompt and hit enter)

The script will generate a variety of reports, including Report3.pdf which represents the clustering of myofibres into OXPHOS normal and OXPHOS deficient sub-populations and the estimation of the threshold mutation level above which OXPHOS defect occurs in single cells.  The script also outputs TestStatistic.pdf, which shows the result of permutation tests to look for pairs of patients with significantly different thresholds.
