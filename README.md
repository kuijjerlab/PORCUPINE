# PORCUPINE
**P**rincipal Components Analysis to **O**btain **R**egulatory **C**ontributions **U**sing **P**athway-based Interpretation of **N**etwork **E**stimates is a an R package to identify biological pathway which drive inter-tumour heterogeneity in a population of gene regulatory networks. It is a Principal Components Analysis (PCA)-based approach that can be used to determine whether a specific set of variables—for example a set of genes in a specific pathway—have coordinated variability in their regulation.

## Method
PORCUPINE uses as input individual patient networks, for example networks modeled using PANDA and LIONESS, as well as a .gmt file that includes biological pathways and the genes belonging to them. For each pathway, it extracts all edges connected to the genes belonging to that pathway and scales each edge across individuals. It then performs a PCA analysis on these edge weights, as well as on a null background that is based on random pathways. To identify significant pathways, PORCUPINE applies a one-tailed t-test and calculates the effect size (ES). 
Here we provide example how we analysed heterogeneity among single gene regulatory sample networks in Leiomyosarcomas (LMS). Networks were obtained with PANDA and LIONESS algorithms. Our example dataset contains data for 20 LMS patient specific gene regulatory networks. 

## Usage
```{r}
library(remotes)
install_github('kuijjerlab/PORCUPINE')
library("PORCUPINE")
```
For this analysis, we need to load the following packages:
```{r}
require("data.table")
require("fgsea")
require("dplyr")
require("plyr")
require("purrr")
require("stats"
require("parallel")
require("lsr")
```
First, we load the network data. Network file is quite big and it takes time to read it in. Our example contains 20 networks, represented by 11,151,077 edges.

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

dim(net)
# [1] 11151077       20
```

Next, we need to load in edges information, corresponding to our network data. Edges table should contain columns "tar" and "reg". In our case, edge table contains 623 Tfs and 17,899 target genes.

```{r}
edges <- fread(file.path(data_dir, "edges.txt"))
head(edges)
#      reg  tar prior
# 1:  AIRE A1BG     0
# 2:  ALX1 A1BG     0
# 3:  ALX3 A1BG     0
# 4:  ALX4 A1BG     0
# 5:    AR A1BG     0
# 6: ARID2 A1BG     0

length(unique(edges$reg))
# [1] 623
length(unique(edges$tar))
# [1] 17899
```
Then, we need to load in pathway file (.gmt file). Gmt files can be downloaded from http://www.gsea-msigdb.org/gsea/msigdb/collections.jsp
```{r}
pathways <- load_gmt(file.path(data_dir, "c2.cp.reactome.v7.1.symbols.gmt"))
length(pathways)
# [1] 1532
```
We need to filter pathways in order to include only pathways with genes that are present in our network file. Additionally, we filter pathways based on their size, and all pathways with less than 10 and more than 150 genes are filtered out. This filtering leaves 1128 pathways.
```{r}
pathways <- filter_pathways(pathways, edges,  minSize = 10, maxSize = 200)
length(pathways)
# [1] 1128
```
We select the top 10 pathways for the analysis (just to reduce computational time).
```{r}
pathways_to_use <- pathways[1:10]
```
Next, we perform PCA analysis for the 10 selected pathawys and extact the information on the variance explained by the first principal component in each pathway.
```{r}
pca_res_pathways <- pca_pathway(pathways_to_use, net, edges, ncores = 10)
head(pca_res_pathways)
#                                                                   pathway
# 1 REACTOME_A_TETRASACCHARIDE_LINKER_SEQUENCE_IS_REQUIRED_FOR_GAG_SYNTHESIS
# 2                               REACTOME_ABACAVIR_TRANSPORT_AND_METABOLISM
# 3                          REACTOME_ABC_FAMILY_PROTEINS_MEDIATED_TRANSPORT
# 4                                       REACTOME_ABC_TRANSPORTER_DISORDERS
# 5                           REACTOME_ABC_TRANSPORTERS_IN_LIPID_HOMEOSTASIS
# 6   REACTOME_ABORTIVE_ELONGATION_OF_HIV_1_TRANSCRIPT_IN_THE_ABSENCE_OF_TAT
#      pc1 n_edges pathway_size
# 1 34.897   16198           26
# 2 36.113    6230           10
# 3 29.476   60431           97
# 4 31.932   44233           71
# 5 32.061   11214           18
# 6 31.695   13083           21
```
Then we perform a PCA analysis based on random gene sets. In this case we create 500 random gene sets for each pathway and run PCA for each gene set.

```{r}
pca_res_random <- pca_random(net, edges, pca_res_pathways, pathways, n_perm = 500, ncores = 20)
head(pca_res_random)

#   pathway    pc1 n_edges pathway_size
# 1       1 28.071   16198           26
# 2       2 32.989   16198           26
# 3       3 25.577   16198           26
# 4       4 28.452   16198           26
# 5       5 22.744   16198           26
# 6       6 37.569   16198           26

```
Then to identify significant pathways we run PORCUPINE, which compares the observed PCA score for a pathway to a set of PCA scores of random gene sets of the same size as pathway. Calculates p-value and effect size.
```{r}
res_porcupine <- PORCUPINE(pca_res_pathways, pca_res_random)
res_porcupine$p.adjust <- p.adjust(res_porcupine$pval, method = "fdr")
res_porcupine

#                                                                    pathway
# 1                    REACTOME_ACTIVATED_NTRK2_SIGNALS_THROUGH_FRS2_AND_FRS3
# 2                     REACTOME_ACETYLCHOLINE_NEUROTRANSMITTER_RELEASE_CYCLE
# 3                      REACTOME_ACETYLCHOLINE_BINDING_AND_DOWNSTREAM_EVENTS
# 4    REACTOME_ABORTIVE_ELONGATION_OF_HIV_1_TRANSCRIPT_IN_THE_ABSENCE_OF_TAT
# 5                            REACTOME_ABC_TRANSPORTERS_IN_LIPID_HOMEOSTASIS
# 6                                        REACTOME_ABC_TRANSPORTER_DISORDERS
# 7                           REACTOME_ABC_FAMILY_PROTEINS_MEDIATED_TRANSPORT
# 8                                REACTOME_ABACAVIR_TRANSPORT_AND_METABOLISM
# 9  REACTOME_A_TETRASACCHARIDE_LINKER_SEQUENCE_IS_REQUIRED_FOR_GAG_SYNTHESIS
# 10                REACTOME_ACTIVATED_NOTCH1_TRANSMITS_SIGNAL_TO_THE_NUCLEUS
   pathway_size         pval         es     p.adjust
# 1            11 1.000000e+00 1.22831105 1.000000e+00
# 2            15 4.860446e-20 0.42427094 1.215112e-19
# 3            12 1.000000e+00 1.13698534 1.000000e+00
# 4            21 8.951552e-02 0.06017867 1.491925e-01
# 5            18 1.623396e-01 0.04408903 2.319137e-01
# 6            71 2.429933e-39 0.63868153 1.214967e-38
# 7            97 1.000000e+00 0.48791384 1.000000e+00
# 8            10 8.087782e-23 0.45858385 2.695927e-22
# 9            26 4.592328e-70 0.93116198 4.592328e-69
# 10           26 1.478399e-02 0.09758639 2.956799e-02
```