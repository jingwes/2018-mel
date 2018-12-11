# Simultaneous measurement of RNA and Protein in single-cells

We want to develop a workflow for the processing of high-dimensional data obtained from single-cell measurements.
For this example, we used 3647 Peripheral Mononuclear Blood Cells and the data can be downloaded [here](https://www.dropbox.com/sh/3sfe8a0nyv1xn0t/AAAe79KVWUmgwJ76Ka0vTAxCa?dl=0).

## Setting up R
```r
library(RColorBrewer)
library(scales)
library(readr)
library(Rtsne)
library(ggalt)
library(ggfortify)
library(plyr)
library(factoextra)
library(cluster)
library(ggplot2)
library(umap)
library(tidyverse)
library(gplots)
```
We also define a theme for plotting later, for aesthetic reason.
```r
  our_theme <- theme(legend.position = "none",
                    panel.grid = element_blank(),
                    axis.title = element_blank(),
                    axis.text = element_blank(),
                    axis.ticks = element_blank(),
                    panel.background = element_blank()
```

The dimensionality reduction tools t-SNE and UMAP already have a package imported, but not PCA. Hence, we need to define PCA with the following line of code:
```r
calcpc <- function(variables, loadings)
  {
  as.data.frame(variables)
  numsamples <- nrow(variables)
  pc <- numeric(numsamples)
  numvariables <- length(variables)
  for (i in 1:numsamples)
  {
    valuei <- 0
    for (j in 1:numvariables)
    {
      valueij <- variables[i, j]
      loadingj <- loadings[j]
      valuei <- valuei + (valueij * loadingj)
    }
    pc[i] <- valuei
  }
  return(pc)
```
Although we run PCA more than once, we only need to run the chunk of code above once before running PCA.

## Importing data
From the .csv file, we can see that column 1 is used to store the unique cell barcode for each of the 3647 cells, and row 1 is used to list the variables measured. Hence, we would want to indicate `header = T` to remove the first row from the list of observations.
```r
adt <- read.csv("common_adt.csv", header = T)
row.names(adt) <- adt[, 1]
adt <- adt[, 2:53]                                             # formatting cell barcodes as observation names
adt_matrix <- as.matrix(adt)                                   # converting to matrix class
adt_matrix_2 <- adt_matrix                                     # making duplicate matrix to compare clusterings later

rna <- read.csv("common_rna.csv", header = T)
row.names(rna) <- rna[,1]
rna <- rna[,2:1087]
rna_matrix <- as.matrix(rna)
rna_matrix_2 <- rna_matrix
```
## Dimensionality reduction
We try the various techniques and compare the results. It is observed that PCA is the fastest algorithm, followed by UMAP and lastly t-SNE.

### Principal Component Analysis (PCA)
PCA is known to be the 'standard' method for such analysis protocols, as results are consistent every time with no random variables involved. 

#### PCA on Protein data

```r
adt_matrix.pca <- prcomp(adt_matrix, retx = TRUE)              
screeplot(adt_matrix.pca, type = "lines")                      # calculate variance contribution across increasing number of PCs
```
We observe a decreasing trend and 2 Principal Components (PCs) is generally sufficient to account for most of the variance. Hence, we plot the data onto the new eigenvectors, PC1 and PC2 and scale the data to the regular x-y axis.
```r
adt_pca <- data.frame(prcomp(adt_matrix)$x[, 1:2])
```

#### PCA on RNA data

The above is repeated for RNA, which also uses 2 PCs.
```r
rna_matrix.pca <- prcomp(rna_matrix, retx = TRUE)              
  screeplot(rna_matrix.pca, type = "lines")
  rna_pca <- data.frame(prcomp(rna_matrix)$x[, 1:2])
```

PCA however favours the global variance and may lead to loss of information within smaller clusters. Hence, we try other methods to compare our findings.

### t-distributed Stochastic Neighbour Embedding (t-SNE)
t-SNE was [proposed in 2008](http://www.jmlr.org/papers/volume9/vandermaaten08a/vandermaaten08a.pdf) and has been favoured by researchers for its ability to preserve relationships between neighbouring points within clusters while reducing dimensionality.

#### t-SNE on Protein data
We run t-SNE on protein data:
```
set.seed(1)                                                    # for reproducibilty as t-SNE results are not perfectly identical         
adt_matrix_tsne <- read.csv("common_adt.csv", header = T)
adt_tsne <- Rtsne(rna_matrix_tsne[, -1],
                  dims = 2,
                  initial_dims = 52,
                  perplexity = 35,
                  verbose = TRUE,
                  max_iter = 1000,
                  pca = FALSE)
}
```
#### t-SNE on RNA data
We repeat the above for RNA.

```
set.seed(1)                                             
rna_matrix_tsne <- read.csv("common_rna.csv", header = T)
rna_tsne <- Rtsne(rna_matrix_tsne[, -1],
                  dims = 2,
                  initial_dims = 1086,
                  perplexity = 35,
                  verbose = TRUE,
                  max_iter = 1000,
                  pca = FALSE)
```

### Uniform Manifold Approximation and Projection
UMAP is a recent, alternative to t-SNE proposed with (notable decrease in runtime and consistency) [https://www.biorxiv.org/content/early/2018/04/10/298430].

#### UMAP on Protein data
```r
adt_umap <- umap(adt)
```

#### UMAP on RNA data
```r
rna_umap <- umap(rna)
```

## Cell clustering algorithms
After reducing dimensionality, cell clustering techniques allow us to identify similar populations of cells within the new 2-dimensional space and deduce possible cell types.

### _k_-means clustering
k-means clustering is used to group points into similar clusters based on similar features. Points are defined on n-axes based on the number of measurements. Based on hypothesis input, random centroids are generated, and points are considered to ‘belong’ to the centroid-defined cluster. The centroids are then repositioned to the centre position of all points in that cluster, and re-clustered until there is no change to attribution of points to clusters. Thus, points belonging to the same group are more similar to each other than to points in other groups.
#### _k_-means onto PCA-reduced datasets

We estimate less than 15 clusters, so we generate clustering for between 2 to 15 clusters and calculate Within-Group Sum of Squares (**WSS**), the average distance within points of a cluster to its centroid. We take the number of clusters, _k_ to be when **WSS** shows asymptotic behaviour and is minimised.

##### Protein
```r
wss <- (nrow(adt_pca) - 1) * sum(apply(adt_pca, 2, var))
for (i in 2:15)
  wss[i] <- sum(kmeans(adt_pca, centers = i)$withinss)
plot(1:15, wss, type = "b",
     xlab = "Number of Clusters",
     ylab = "Within groups sum of squares")
```
##### RNA


### Hierarchical clustering
