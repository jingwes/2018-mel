hclust_adt <- hclust(dist(adt))
hclust_rna <- hclust(dist(rna))

plot(hclust_adt, labels=FALSE)
plot(hclust_rna, labels=FALSE)


adt_clusterCut <- cutree(hclust_adt, 6)
rna_clusterCut <- cutree(hclust_rna, 9)


hclust_adt <- hclust(dist(adt_umap$layout))
hclust_rna <- hclust(dist(rna_umap$layout))

hclust_adt <- hclust(dist(adt_pca))
hclust_rna <- hclust(dist(rna_pca))

plot(hclust_adt, labels=FALSE)
plot(hclust_rna, labels=FALSE)


adt_clusterCut <- cutree(hclust_adt, 4)
rna_clusterCut <- cutree(hclust_rna, 3)

{
  ggplot(data.frame(x = (adt_umap$layout[,1]),
                    y = (adt_umap$layout[,2])),
         aes(x, y)) + geom_point(aes(colour = factor(adt_clusterCut))) + our_theme
}

{
  ggplot(data.frame(x = (rna_umap$layout[,1]),
                    y = (rna_umap$layout[,2])),
         aes(x, y)) + geom_point(aes(colour = factor(rna_clusterCut))) + our_theme
}
