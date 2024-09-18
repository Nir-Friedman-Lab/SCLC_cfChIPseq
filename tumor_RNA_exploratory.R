# define liver signature. strongest signal in many of the sample --------
liver.ind = which(names(rna_roadmap) %in% c("Adult_Liver", "HEPG2"))
other.ind = which(!names(rna_roadmap) %in% c("Adult_Liver", "HEPG2"))
liver.genes = which(rna_roadmap[,"Adult_Liver"] / rowMaxs(as.matrix(rna_roadmap[,other.ind])) > 2 & 
                      rna_roadmap[,"Adult_Liver"] > 8)
liver.genes = sort(rownames(rna_roadmap)[liver.genes])
liver.genes.more = c("AMBP", "APOB", "APOC3", "FGB", "GC", "RBP4")
liver.genes = c(liver.genes, liver.genes.more)
liver.genes = liver.genes[liver.genes %in% rownames(rna_data)]

metadata.matching %>% 
  filter(Sample_id %in% colnames(rna_data)) %>% 
  left_join(metadata.biopsies %>% select(-c(PatientID, `biopsy site`)), 
                                         by = "RNA_id") %>% 
  mutate(liver.sum = colSums(rna_data[liver.genes,Sample_id]), 
         biopsy_site = if_else(Site.type == "Liver", "liver", "other")) %>% 
  column_to_rownames("Sample_id") -> ca

liver.data = log2(1+rna_data[liver.genes,])
col.ord = rankSubgroups(ca[colnames(liver.data),], "biopsy_site", "liver.sum", 
                        decreasing = T)
pheatmap::pheatmap(liver.data[,col.ord], fontsize = 4,
         annotation_col = ca %>% select(biopsy_site),
         show_colnames = F, treeheight_row = 0, treeheight_col = 0,
         # filename = paste0(figDirPaper, "figureS3/rna_liver.pdf"),
         cluster_cols = F, width = 3, height = 3)

HeatmapAnnotation(df = ca %>% select(biopsy_site), 
                  show_annotation_name = T, 
                  simple_anno_size = unit(1, "mm"), 
                  annotation_label = "biopsy site", 
                  annotation_name_gp = gp,
                  show_legend = F) -> top.a

data.frame(ind = 1:nrow(liver.data), 
           gene = rownames(liver.data)) %>% 
  slice_sample(n = 30) -> data.gene 
HeatmapAnnotation(foo = anno_mark(labels_gp = gp, 
                                  at = data.gene$ind, 
                                  labels = data.gene$gene,
                                  lines_gp = gpar(lwd = .2), 
                                  link_width = unit(2, "mm")
                                  ), which = "row", 
                  annotation_name_gp = gp, 
                  gp = gp, 
                  show_legend = T, 
                  annotation_legend_param = list(legend_gp = gp, 
                                                 labels_gp = gp,
                                                 title_gp = gp,
                                                 direction = "horizontal",
                                                 # legend_label_gp= gp,
                                                 legend_height = unit(4, "mm"))) -> gene.a 
pdf(paste0(figDirPaper, "figureS3/rna_liver_v1.pdf"), width = 4, height = 4)
Heatmap(liver.data[,rownames(ca)], 
        column_order = col.ord, 
        heatmap_width = unit(80, "mm"),
        heatmap_height = unit(90, "mm"),
        name = " ",
        column_split = ca %>% select(biopsy_site),
        show_row_names = F,
        show_row_dend = F,
        show_column_names = F,
        top_annotation = top.a, 
        right_annotation = gene.a, 
        column_title_gp = gp, 
        column_dend_gp = gp, 
        column_gap = unit(1, "pt"),
        use_raster = T,
        heatmap_legend_param = list(title_gp = gp, 
                                    legend_gp = gp,
                                    labels_gp = gp,
                                    legend_gp = gp,
                                    legend_label_gp = gp, 
                                    legend_height = unit(4, "mm")))
dev.off()

high.val = round(max(liver.data))
data.frame(liver.data[,col.ord], check.names = F) %>% 
  mutate(gene = rownames(.)) %>%
  melt(id.vars = "gene", variable.name = "samp") %>%
  heatmap(x = "samp", y = "gene", fill = "value", clow = "white", chigh = "red",
          divergent = F, breaks = c(0,high.val), lim = c(0,high.val), 
          leglab = "log2\n(1+RPKM)", leg = c("0", high.val)) + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
        axis.text.y = element_text(size = 5)) +
  scale_y_discrete(guide = guide_axis(n.dodge = 2)) -> p.hm

ord.samps = colnames(liver.data[,col.ord])
ca[ord.samps,] %>% 
  mutate(samp = factor(ord.samps, levels = ord.samps)) %>%
  ggplot() + 
  geom_tile(aes(samp, "biopsy site", fill = biopsy_site)) +
  labs(x = "", y = "", fill = "") + 
  scale_fill_aaas() + 
  theme(axis.text.x = element_blank(), axis.ticks = element_blank(),
        axis.line = element_blank(), legend.key.size = unit(5, "mm")) -> p.bs


p = plot_grid(p.bs, p.hm, ncol = 1, rel_heights = c(1/10,9/10), align = "hv")
p = rasterise(p)
ggsave(paste0(figDirPaper, "figureS3/rna_liver.pdf"), p, width = 90, 
       height = 100, units = "mm")

pheatmap::pheatmap(rna_roadmap[liver.genes,], 
                   labels_col = substr(colnames(rna_roadmap), 1, 12), 
         # filename = paste0(figDirPaper, "figureS3/rna_liver_roadmap.pdf"), 
         fontsize_row = 5, fontsize_col = 7, width = 3, height = 3)



pheatmap::pheatmap(log2(1+chip_data_all[liver.genes[liver.genes %in% rownames(chip_data_all)],c(s.samples, h.samples)]), 
         # filename = paste0(figDirPaper, "temp/chip_liver_roadmap.png"), width = 12, height = 7, treeheight_row = 0, 
         labels_col = substr(c(s.samples, h.samples), 1, 12), annotation_col = sample.annotation, fontsize_row = 5, fontsize_col = 5)

hist(colSums(rna_data[liver.genes,]), nbins = 20)
high.liver.samples = colnames(rna_data)[colSums(rna_data[liver.genes,]) >200]
low.liver.samples = colnames(rna_data)[colSums(rna_data[liver.genes,]) < 50]

high.liver.samples[high.liver.samples %in% high.ctDNA.samples.wRNA]
# consider removing these samples from RNA-ChIP comparison

matched.rna.roadmap.genes = intersect(rownames(rna_data), rownames(rna_roadmap))
pca.rna = rna_data[matched.rna.roadmap.genes,rna.samples]
pca.rna = cbind(pca.rna, roadmap_liver = rna_roadmap[matched.rna.roadmap.genes,"Adult_Liver"])
res.pca = prcomp(t(pca.rna))
data.pca.rna = data.frame(samp = substr(names(pca.rna),1,11), pc1 = res.pca$x[,1], pc2 = res.pca$x[,2], 
                          liver = colSums(pca.rna[liver.genes,]))
p = ggplot(data.pca.rna, aes(pc1, pc2, label = samp)) 
p = p + geom_point(shape = 16, aes(alpha = log10(liver))) + geom_text_repel(size = 1) 
p = p + labs(x = "PC1", y = "PC2", alpha = "log10\n(counts/liver-genes)") 
ggsave(paste0(figDirPaper, "figureS3/pca_rna_liver.pdf"), p, width = 90, height = 60, units = "mm")

p = scatter.plot(data.pca.rna, "pc1", "liver", xlab = "PC1", ylab = "counts/liver-genes")
ggsave(paste0(figDirPaper, "figureS3/pc1_vs_liver_genes.pdf"), p, width = 60, height = 60, units = "mm")
# TF signatures -----------------------------------------------------------
# effect of liver fraction on TF 
data.rna = data.frame(samp = rna.samples, 
                      row.names = rna.samples, 
                      liver = colSums(rna_data[liver.genes,]), 
                      t(rna_data[sclc.subtypes,]))
data.rna.m = melt(data.rna, id.vars = c("samp", "liver"), variable.name = "tf", value.name = "val")
head(data.rna.m)
ggplot(data.rna.m, aes(liver, val, fill = tf)) + geom_point() + stat_cor() + facet_wrap(~tf, scales = "free")
ggplot(data.rna.m, aes(x = val)) + geom_histogram(color = "white", fill = "black", binwidth = .5) +
  facet_wrap(~tf, scale = "free")

ca = data.rna[,c("liver", "ASCL1", "NEUROD1", "YAP1", "POU2F3")]
# ascl1 targets
low.ascl1 = data.rna$samp[data.rna$ASCL1 < 1]
high.ascl1 = data.rna$samp[data.rna$ASCL1 > 4]
ascl1.targets = rownames(rna_data)[which(rowMedians(as.matrix(rna_data[,high.ascl1])) / 
                                           rowMedians(as.matrix(rna_data[,low.ascl1])) > 2 &
                                           rowMedians(as.matrix(rna_data[,high.ascl1])) > 5)]
pheatmap(rna_data[ascl1.targets,], annotation_col = ca, show_colnames = F)

low.pou2f3 = data.rna$samp[data.rna$POU2F3 < 1]
high.pou2f3 = data.rna$samp[data.rna$POU2F3 > 2]
pou2f3.targets = rownames(rna_data)[which(rowMins(as.matrix(rna_data[,high.pou2f3])) > 
                                            rowMaxs(as.matrix(rna_data[,low.pou2f3])) &
                                            rowMedians(as.matrix(rna_data[,high.pou2f3])) > 3)]
pheatmap(rna_data[pou2f3.targets,], annotation_col = ca, show_colnames = F)