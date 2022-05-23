# PORCUPINE

PORCUPINE 
**P**rincipal Components Analysis to **O**btain **R**egulatory **C**ontributions **U**sing **P**athway-based Interpretation of **N**etwork **E**stimates 
An R package to identify biological pathway which drive inter-tumour heterogeneity in a population of gene regulatory networks. It is a Principal Components Analysis (PCA)-based approach that can be used to identify key pathways that drive heterogeneity among individuals in a dataset to determine whether a specific set of variables—for example a set of genes in a specific pathway—have coordinated variability in their regulation.

## Method
PORCUPINE uses as input individual patient networks, for example networks modeled using PANDA and LIONESS, as well as a .gmt file that includes biological pathways and the genes belonging to them. For each pathway, it extracts all edges connected to the genes belonging to that pathway and scales each edge across individuals. It then performs a PCA analysis on these edge weights, as well as on a null background that is based on random pathways. To identify significant pathways, PORCUPINE applies a one-tailed t-test and calculates the effect size (ES). 
Here we provide example how we analysed heterogeneity among single gene regulatory sample networks in Leiomyosarcomas (LMS). Networks were obtained with PANDA and LIONESS algorithms. Our example dataset contains data for 20 LMS patient specific gene regulatory networks. 

## Usage
```{r}
library(remotes)
install_github('kuijjerlab/PORCUPINE')
library("PORCUPINE")
```{r}
net <- fread(file.path(data_dir, "20LMS_nets.txt"))
net[1:5, 1:5]
#    0E244FE2-7C17-4642-A51F-2CCA796D9C70 75435ED8-93E8-45FB-8480-98D8EB2EF8CB
# 1:                                 0.76                                 0.10
# 2:                                 0.94                                 1.43
# 3:                                 1.09                                 2.78
# 4:                                 1.13                                 2.60
# 5:                                -0.71                                -1.42
#    B6D11678-15A9-4F43-A0A2-225067DCAF1C B7F5A41E-9559-4329-81F5-1B88A74730B7
# 1:                                -1.27                                 0.01
# 2:                                 0.30                                 0.91
# 3:                                 1.01                                 2.13
# 4:                                 1.66                                 1.71
# 5:                                 0.02                                 0.27
#    04823F53-A12D-4852-8F34-77B9DCBB7DF0
# 1:                                -7.18
# 2:                                -5.69
# 3:                                -6.09
# 4:                                -6.04
# 5:                                 5.31
```
edges <- fread(file.path(data_dir, "edges.txt"))
head(edges)