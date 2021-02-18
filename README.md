# PORCUPINE

PORCUPINE (PCA to Obtain Regulatory Contributions Using Pathway-based Interpretation of Network Estimates). The purpose of PORCUPINE is to detect statistically significant key regulatory pathways which classify the population of patient gene regulatory networks into biologically meaningful subtypes.This approach determines if a subset of variables, in this case a set of genes in a specific pathway, can separate samples into biologically meaningful groups and therefore allows to identify differentially regulated pathways among the individuals. PORCUPINE performs a PCA analysis on the scaled TF-gene edge weight matrix extracted for genes in each pathway.The amount of variance explained by PC1 is then compared to the amount of variance explained in random data. For randomization, a set of 1000 gene subtsets of the equal size as the pathway of interest is created. If the amount of variance explained by a principal component of a pathway is significantly higher than expected by generated 1000 random gene sets of the same size, a pathway is considered as significant contributor to heterogeneity. To identify significant pathways, one sided t-test p-value and the effect size are calculated.

Here we provide example how we analysed heterogeneity among single gene regulatory sample networks in Leiomyosarcomas (LMS). Networks were obtained with PANDA and LIONESS algorithms. Our example dataset contains data for 20 LMS patient specific gene regulatory networks. To run PORCOPINE analysis one needs dataset with networks and gmt pathway file.

First, we load the network data. Network file is quite big and it takes time to read it in. Our example contains 20 netowks represented by 11,151,077 edges between 623 Tfs and 17,899 target genes. it is improtant to note, that the first three columns in network file are “reg”, “tar”, “prior”.


install.packages("devtools")
library(devtools)
devtools::install_github("aiwa1986/PORCUPINE")

require(porcupine)
lms <- readRDS(file.path(data_dir, "20_LMS.Rdata"))
lms[1:5, 1:5]

##    reg  tar prior 0E244FE2-7C17-4642-A51F-2CCA796D9C70
## 1 AIRE A1BG     0                                 0.76
## 2 ALX1 A1BG     0                                 0.94
## 3 ALX3 A1BG     0                                 1.09
## 4 ALX4 A1BG     0                                 1.13
## 5   AR A1BG     0                                -0.71
##   75435ED8-93E8-45FB-8480-98D8EB2EF8CB
## 1                                 0.10
## 2                                 1.43
## 3                                 2.78
## 4                                 2.60
## 5                                -1.42
