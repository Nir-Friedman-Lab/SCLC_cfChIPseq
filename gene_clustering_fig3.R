# prepare data ------------------------------------------------------------
library(dendextend)
s = s.samples[estimated.tumor[s.samples, "SCLC.n"] > .3]
# s = high.ctDNA.samples
selected_genes_file = paste0(baseDir, "high_cor_genes2.csv")
# genes = read.table(selected_genes_file)$V1
genes = high.sig.genes; outname = "high_cor_genes_high_samples"# genes that had high correlation to RNA in high samples
# genes = diff.genes.common
clus_matrix = log(1+chip_data_all[genes, s])
clus_matrix = sweep(clus_matrix, 1, rowMeans(clus_matrix), "-")
hist(rowSds(clus_matrix), breaks = "fd")
clus_matrix = clus_matrix[(rowSds(clus_matrix) > 1),]
hr.clus = PClusterMatrixRows(clus_matrix) # cluster rows
hc.clus = PClusterMatrixRows(t(clus_matrix)) # cluster columns
# writeCluster(clus_matrix, row.clust = hr.clus, col.clust = hc.clus,
#              ClusterFilename = "~/treeview_files/SCLC/clusters") # for viewing in tree view
r.ord = hr.clus$labels[hr.clus$order]
c.ord = hc.clus$labels[hc.clus$order]

pheatmap(clus_matrix[r.ord, c.ord], cluster_rows = F, cluster_cols = F, show_rownames = F, 
         show_colnames = F, breaks = cfm.breaks, color = cfm.color,
         annotation_col = estimated.tumor[,"SCLC.n", drop = F])

clusters = data.frame(cluster = as.factor(cutree(hr.clus, k = 2)))
# write.csv(clusters, paste0(figDirPaper, "gene_clusters.csv"))
samp.clusters = data.frame(cluster = as.factor(cutree(hc.clus, k = 4)))
ca = data.frame(SCLC_score = estimated.tumor[rownames(samp.clusters),"SCLC.n"], samp.clusters)
c1.genes = rownames(clusters)[clusters$cluster == 1]
c2.genes = rownames(clusters)[clusters$cluster == 2]
c3.genes = rownames(clusters)[clusters$cluster == 3]
c4.genes = rownames(clusters)[clusters$cluster == 4]
# c5.genes = rownames(clusters)[clusters$cluster == 5]
pheatmap(clus_matrix[r.ord, c.ord], cluster_rows = F, cluster_cols = F, fontsize_col = 4, 
         show_rownames = F, breaks = cfm.breaks, color = cfm.color, 
         annotation_row = clusters, annotation_col = estimated.tumor[,"SCLC.n", drop = F], 
         filename = paste0(figDirPaper, "figure4/", outname, ".pdf"))


c1 = colSums(chip_data_all[c1.genes,s])
c2 = colSums(chip_data_all[c2.genes,s])
c3 = colSums(chip_data_all[c3.genes,s])
c4 = colSums(chip_data_all[c4.genes,s])
c5 = colSums(chip_data_all[c5.genes,s])

cont.var = list(ne_score = "NE_SCORE", sclc = "SCLC_score", ascl1 = "ASCL1", 
                neurod1 = "NEUROD1", yap1 = "YAP1", pou2f3 = "POU2F3")
catg.var = list(timepoint = "Timepoint", protocol = "protocol", ne = "NE_status")


data.cluster = data.frame(metadata[s, c(unlist(cont.var), unlist(catg.var))], c1, c2)#, c3, c4, c5) #, c4, c5, c6, c7)
data.cluster[,unlist(catg.var)] = lapply(data.cluster[,unlist(catg.var)], factor)
data.cluster[,unlist(cont.var)] = lapply(data.cluster[,unlist(cont.var)], as.numeric)
cluster_names = c("cluster_1", "cluster_2")#, "cluster_3", "cluster_4", "cluster_5")
names(data.cluster) = c(names(cont.var), names(catg.var), cluster_names)

boxFn = function(data, mapping, ...) {
  p = ggplot(data = data, mapping = mapping) +  geom_point(size = .2, alpha = .7)
  p = p + theme(axis.text = element_blank(), strip.text.x = element_text(size = 8, face = "bold"))
  p
}
p = ggpairs(data.cluster[,c(names(catg.var), cluster_names)], 
        lower = list(continuous = "density", combo = "box_no_facet"))
ggsave(paste0(figDirPaper, "figure4/clusters_catg.png"), p, width = 10, height = 7)

p = ggpairs(data.cluster[,c(names(cont.var), cluster_names)])
ggsave(paste0(figDirPaper, "figure4/clusters_cont.png"), p, width = 10, height = 7)

data.protocol = reshape2::melt(data.cluster, id.vars = "protocol", measure.var = cluster_names, 
               variable.name = "cluster")
p = ggplot(data.protocol, aes(cluster, value, fill = protocol, color = protocol)) 
p = p + geom_boxplot(alpha = .8) + geom_point(size = 1, position = position_jitterdodge())
p = p + labs(ylab = "log2(reads/cluster)")

p = ggplot(data.cluster, aes(cluster_1, ne_score)) + geom_point() + geom_smooth(method = "lm", se = F) + stat_cor()
ggplot(data.cluster, aes(cluster_2, ne_score)) + geom_point() + geom_smooth(method = "lm", se = F) + stat_cor()
ggplot(data.cluster, aes(cluster_3, ne_score)) + geom_point() + geom_smooth(method = "lm", se = F) + stat_cor()
