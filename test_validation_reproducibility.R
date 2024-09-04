estimated.tumor %>% 
  filter(rownames(.) %in% train.samples, SCLC.n > .7) %>% 
  rownames(.) -> high.train

estimated.tumor %>% 
  filter(rownames(.) %in% validation.samples, SCLC.n > .7) %>% 
  rownames(.) -> high.validation


data.frame(gene = rownames(chip_data_all), 
           train = rowMeans(chip_data_all[,high.train]), 
           validation = rowMeans(chip_data_all[,high.validation])) %>%
  mutate(diff.genes = gene %in% diff.genes.common) -> data.test.validation 
  
data.test.validation %>% 
  ggplot(aes(train, validation, color = diff.genes)) +
  geom_point(data = data.test.validation %>% filter(!diff.genes), size = .1, 
             color = "gray") +
  geom_point(data = data.test.validation %>% filter(diff.genes), size = .1, 
             color = "darkred") +
  scale_x_continuous(trans = "log1p", breaks=c(0,1,3,10,30,100,300)) + 
  scale_y_continuous(trans = "log1p", breaks=c(0,1,3,10,30,100,300)) + 
  labs(x = "train cohort", y = "validation cohort") + 
  geom_abline(slope = 1, intercept = 0, color="black", linewidth = .2) + 
  geom_density2d(colour="black", show.legend = FALSE, bins=20, size = .2) +
  stat_cor(color = "black", size = base_size/.pt) -> p
p = rasterize(p, dpi = 500)
ggsave(paste0(figDirPaper, "validation/train_vs_validation_scatter.pdf"), p, 
       width = 50, height = 50, units = "mm")  
ggsave(paste0(figDirPaper, "validation/train_vs_validation_scatter.png"), p, 
       width = 50, height = 50, units = "mm")  


chip_data_all[diff.genes.common, c(train.samples, validation.samples, h.samples, 
                                   c.samples)] -> d
d = log2(1+d)
d = sweep(d, 1, rowMeans(d), "-")
d[d > 3] = 3
d[d < -3] = -3
col.annotation %>% 
  mutate(group = if_else(sample %in% validation.samples, "validation", 
                         group), 
         sclc_score = estimated.tumor[sample, "SCLC.n"]) %>% 
  select(group, sclc_score) -> col.a
col.ord = rankSubgroups(col.a[colnames(d),], "group", "sclc_score",decreasing = T, 
                        group_order = c("SCLC", "validation", "NSCLC", "CRC", 
                                        "Healthy"))

pheatmap::pheatmap(d[,col.ord], color = cfm.color,
                   annotation_col = col.a, 
                   cluster_cols = F, 
                   show_rownames = F, show_colnames = F, 
                   treeheight_row = 0, fontsize = 6,
                   filename = paste0(figDirPaper, "validation/SCLC_signature.pdf"),
                   width = 5, height = 4)

pdf(paste0(figDirPaper, "validation/SCLC_signature.pdf"), width = 4, height = 3)
# png(paste0(figDirPaper, "validation/SCLC_signature.png"), width = 100, 
    # height = 60, unit = "mm", res = 100)
Heatmap(d, col = cfm.color, 
        column_title_gp = gpar(fontsize = base_size),
        show_row_names = F, 
        show_column_names = F, show_column_dend = F, 
        name = "gene counts", use_raster = T, 
        column_split = col.a[colnames(d),"group"], show_row_dend = F,
        )
dev.off()
