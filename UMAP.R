library(umap)
library(tidyverse)

adt_umap <- umap (adt)
rna_umap <- umap (rna)

  
  wss <- (nrow(adt_umap$layout) - 1) * sum(apply(adt_umap$layout, 2, var))
  for (i in 2:15)
    wss[i] <- sum(kmeans(adt_umap$layout, centers = i)$withinss)
  plot(1:15, wss, type = "b",
       xlab = "Number of Clusters",
       ylab = "Within groups sum of squares")
  
  kmeans_adt_umap <- kmeans(adt_umap$layout, 3, nstart = 25, iter.max = 1000)     # cluster based on selected number of clusters
  palette(alpha(brewer.pal(9, 'Set1'), 0.5))
  
  
  wss <- (nrow(rna_umap$layout) - 1) * sum(apply(rna_umap$layout, 2, var))
  for (i in 2:15)
    wss[i] <- sum(kmeans(rna_umap$layout, centers = i)$withinss)
  plot(1:15, wss, type = "b",
       xlab = "Number of Clusters",
       ylab = "Within groups sum of squares")
  
  kmeans_rna_umap <- kmeans(rna_umap$layout, 2, nstart = 25, iter.max = 1000)     # cluster based on selected number of clusters
  palette(alpha(brewer.pal(9, 'Set1'), 0.5))




ggplot(data.frame(x = adt_umap$layout[,1],
                  y = adt_umap$layout[,2]),
       aes(x, y)) + geom_point(aes(colour = factor(kmeans_adt_umap$clust))) + our_theme

ggplot(data.frame(x = rna_umap$layout[,1],
                  y = rna_umap$layout[,2]),
       aes(x, y)) + geom_point(aes(colour = factor(kmeans_rna_umap$clust))) + our_theme

