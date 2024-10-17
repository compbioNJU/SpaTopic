# SpaTopic

## A Statistical Learning Framework for Exploring Tumor Spatial Architecture from Spatially Resolved Transcriptomic Data.

![](https://github.com/compbioNJU/SpaTopic/blob/main/Fig/Fig1.png)

Here, we introduce SpaTopic, a versatile statistical learning framework that harmonizes spatial clustering and cell-type deconvolution by integrating single-cell transcriptomics and SRT data. The objective of SpaTopic is to identify spatial clusters of spots within the SRT data, characterized by homogenous gene expression and cell-type organization. Each unique cluster of spots is viewed as a distinct spatial domain, presenting discernible patterns that set it apart from other clusters. SpaTopic significantly aids in the characterization of spatial domains, capitalising on the contributions from their corresponding cell-type topics. This enables the quantitative comparison of spatial domains and the identification of spatial regions prevalent across various SRT datasets. SpaTopic is implemented as an open-source R package, freely available at <https://github.com/compbioNJU/SpaTopic>.

## Citation

Yuelei Zhang et al. ,SpaTopic: A statistical learning framework for exploring tumor spatial architecture from spatially resolved transcriptomic data. Sci. Adv.10,eadp4942(2024). DOI:[10.1126/sciadv.adp4942](https://www.science.org/doi/10.1126/sciadv.adp4942)

## Installation

SpaTopic is implemented as an R package, which can be installed from GitHub.

### Dependencies

-   R (≥ 3.5.0)
-   testthat (\>= 3.0.0)
-   R packages: modeltools, slam, stats, topicmodels

\*\* Install devtools if necessary \*\*

``` r
install.packages('devtools')
```

Install SpaTopic

``` r
devtools::install_github('compbioNJU/SpaTopic')
```

Load package

``` r
 library(SpaTopic)
```

## Issues

All feedback, bug reports and suggestions are warmly welcomed! Please make sure to raise issues with a detailed and reproducible exmple and also please provide the output of your sessionInfo() in R!

## How to use SpaTopic

This tutorial is the example analysis with SpaTopic on the human pancreatic ductal adenocarcinomas data from [Moncada et al, 2020](https://www.nature.com/articles/s41587-019-0392-8?proof=t). Before runing the tutorial, make sure that the SpaTopic package is installed. Installation instructions see the Installation. See more at <https://compbioNJU.github.io/SpaTopic/>.

``` r
library(SpaTopic)
```

### load the data

``` r
load("data/spot_clusters.rda")
load("data/spot_celltype.rda")
```

``` r
spot_clusters[1:5,1:5]
      row col sizeFactor cluster.init spatial.cluster
10x10  10  10  4.7761108            1               2
10x13  10  13  1.0052199            2               2
10x14  10  14  0.8106812            2               2
10x15  10  15  0.4987377            2               2
10x16  10  16  0.4346143            2               2

spot_celltype[1:5,1:5]
      Acinar_cells Ductal_cells Cancer_clone_A Cancer_clone_B         DCs
10x10 5.838572e-02    0.2349066   1.365076e-03   3.892868e-04 0.165860789
10x13 4.807943e-05    0.9984677   1.654640e-06   9.032885e-06 0.001244634
10x14 4.701190e-02    0.8373601   4.846860e-03   9.009235e-04 0.003541947
10x15 5.047613e-02    0.8020465   1.911570e-04   3.325224e-02 0.084113110
10x16 4.694120e-03    0.9718078   1.719378e-06   6.266388e-04 0.007665514
```

### use SpaTopic to the spot_celltype and spot_clusters

``` r
#result_list: A list with three data frame and one vector. 
#MetaTopic is a data frame which can be add to a Seurat object. 
#The domain_topic is a data frame, row is CellTopic. and col is domain.
#The celltype_topic is a data frame, row is celltype and col is CellTopic. 
#Cell_topic is a vector of which topic be chosen in each CellTopic. 
#If meta.cell = TRUE, one more result will be given in result list, MetaTopic is a data frame of the cluster result of CellTopic.

result_list <- CellTopic(spot_celltype,spot_clusters,cluster = "spatial.cluster", num_topics = 13,percent = 0.7,
            Binarization = FALSE, meta.cell = FALSE, k = NULL)
```

``` r
#show the result
head(result_list[["CellTopic"]])
       CellTopic        CellTopic1        CellTopic2         CellTopic3         CellTopic4
10x10 CellTopic2 0.577382618544802 0.787303032098654 0.0080170243865711 0.0853445821596965
10x13 CellTopic2 0.577382618544802 0.787303032098654 0.0080170243865711 0.0853445821596965
10x14 CellTopic2 0.577382618544802 0.787303032098654 0.0080170243865711 0.0853445821596965
10x15 CellTopic2 0.577382618544802 0.787303032098654 0.0080170243865711 0.0853445821596965
10x16 CellTopic2 0.577382618544802 0.787303032098654 0.0080170243865711 0.0853445821596965
10x17 CellTopic2 0.577382618544802 0.787303032098654 0.0080170243865711 0.0853445821596965

head(result_list[["domain_topic"]])
           spot_domain_1 spot_domain_2 spot_domain_3 spot_domain_4
CellTopic1    0.78207686   0.577382619   0.174953872    0.10799194
CellTopic2    0.44211741   0.787303032   0.007438506    0.06603564
CellTopic3    0.12712585   0.008017024   0.787181571    0.03422577
CellTopic4    0.05105707   0.085344582   0.018005238    0.78840065

head(result_list[["celltype_topic"]])
               CellTopic1 CellTopic2   CellTopic3 CellTopic4
Acinar_cells   0.04503436 0.03515437 5.404895e-02 0.17033144
Ductal_cells   0.11062714 0.14213899 3.273553e-06 0.02809054
Cancer_clone_A 0.03090381 0.02157599 1.770954e-01 0.01540925
Cancer_clone_B 0.02943409 0.01671419 1.662732e-01 0.01086386
DCs            0.07268275 0.06557971 2.767503e-02 0.14179855
Tuft_cells     0.06113322 0.04374492 5.155769e-02 0.14100323

head(result_list[["Cell_topic"]])
    CellTopic1     CellTopic2     CellTopic3     CellTopic4 
"3_11_4_5_7_2"   "2_8_1_11_3"         "9_12"        "13_10" 
```
