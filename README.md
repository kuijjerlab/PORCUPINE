# PORCUPINE

PORCUPINE (Principal Components Analysis to Obtain Regulatory Contributions Using Pathway-based Interpretation of Network Estimates) is a Principal Components Analysis (PCA)-based approach that can be used to identify key pathways that drive heterogeneity among individuals in a dataset. It determines whether a specific set of variables—for example a set of genes in a specific pathway—have coordinated variability in their regulation.
PORCUPINE uses as input individual patient networks, for example networks modeled using PANDA and LIONESS, as well as a .gmt file that includes biological pathways and the genes belonging to them. For each pathway, it extracts all edges connected to the genes belonging to that pathway and scales each edge across individuals. It then performs a PCA analysis on these edge weights, as well as on a null background that is based on random pathways. To identify significant pathways, PORCUPINE applies a one-tailed t-test and calculates the effect size (ES). 
Here we provide example how we analysed heterogeneity among single gene regulatory sample networks in Leiomyosarcomas (LMS). Networks were obtained with PANDA and LIONESS algorithms. Our example dataset contains data for 20 LMS patient specific gene regulatory networks. Please check vignettes for example of how to run PORCUPINE.


