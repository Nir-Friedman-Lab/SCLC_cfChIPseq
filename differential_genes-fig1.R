# SCLC diffGenes --------------------------------------------------
# diff.genes are defined in linear-regression-SCLC-score.R
diff.genes.common = c()
for (s in diff.genes.high.samples) {
  d = read.csv(paste0(baseDir, "Output/H3K4me3/DiffGenes/", s, ".csv"))
  diff.genes.common = c(diff.genes.common, d$X[d$FoldChange..log2. > 0])
}
diff.genes.common = table(diff.genes.common)
diff.genes.common = names(which(diff.genes.common > 5))
write.table(diff.genes.common, paste0(tableDir, "common.diff.genes.txt"), 
            quote = F, row.names = F, col.names = F)

# histogram of number of differential genes -------------------------------
lab = paste0("high SCLC\n(n=", length(diff.genes.high.samples),")"); 

diff.genes %>% 
  mutate(pretreatment = case_match(
    rownames(.),
    pre.samples ~ "Pre",
    post.samples ~ "Post", 
    .default = NA
  )) %>% 
  filter(!is.na(pretreatment)) %>% 
  ggplot(aes(x = n.genes, fill = pretreatment)) +
  geom_histogram(binwidth = 100, linewidth = .2, color = "white") + 
  scale_fill_aaas() +
  geom_vline(xintercept = diff.genes.cutoff, linetype="dashed", color = "red", size = .2) + 
  labs(x = "# of genes", y = "# of SCLC samples", fill = "treatment") +
  annotate("text", x = 3000, y = 30, label = lab, size = base_size/.pt) +
  theme(legend.position = c(.8, .8), legend.key.size = unit(5, 'pt')) -> p
ggsave(paste0(figDirPaper, "figureS1/diffGenes_counts.pdf"), p, 
       width = 50, height = 50, units = "mm")


# scatter of mean high SCLC vs healthy ------------------------------------
lab.genes = c("NFIB", "CHGA", "INSM1", "DLL3", "FOXA2", "EGFR")
data.frame(sclc = rowMeans(chip_data_all[,diff.genes.high.samples]), 
           healthy = rowMeans(chip_data_all[,h.samples]), 
           diffgene = NA) %>% 
  mutate(diffgene = if_else(diffgene %in% diff.genes.common, "SCLC", diffgene)) %>% 
  rownames_to_column("gene") -> 
           data.high_sclc_healthy

mVal = ceiling(max(quantile(data.high_sclc_healthy$sclc,0.995), quantile(data.high_sclc_healthy$healthy, 0.995)))
scatter.plot(data.high_sclc_healthy, x = "healthy", y = "sclc", 
             xlab = "Healthy baseline", ylab = "High SCLC samples", 
             title = paste(length(diff.genes.common), "significant"), 
             color = "diffgene", size = 0.15, cor = F) +
  geom_point(data = data.high_sclc_healthy %>% 
               filter(gene %in% diff.genes.common), shape = 16, 
             color = group.colors$SCLC,size = .15) + 
  coord_fixed(ratio=1,xlim=c(0,mVal), ylim=c(0,mVal)) + 
  scale_x_continuous(trans = "log1p", breaks=c(0,1,3,10,30,100,300)) + 
  scale_y_continuous(trans = "log1p", breaks=c(0,1,3,10,30,100,300)) + 
  geom_abline(slope = 1, intercept = 0, color="black", linewidth = base_line_size) + 
  geom_density2d(colour="black", show.legend = FALSE, bins=20, size = .2) + 
  geom_label_repel(data = data.high_sclc_healthy %>% 
                     filter(gene %in% lab.genes), color = "black", 
                         mapping = aes(label = gene), size = base_size/.pt, 
                         label.size = NA, label.padding = unit(1, "pt")) -> p
p = rasterize(p, dpi = 500)
ggsave(paste0(figDirPaper, "figure1/sclc_healthy_scatter.pdf"), p, width = 50,
       height = 50, units = "mm")

# boxplot SCLC score in plasma groups  ----------------------------
estimated.tumor %>% 
  filter(sample %in% c(s.samples, h.samples, l.samples, c.samples)) %>% 
  mutate(group = case_match(
    sample,
    pre.samples ~ "SCLC\npre",
    post.samples ~ "SCLC\npost",
    .default = group), 
    group = factor(group, levels = c("SCLC\npre", "SCLC\npost", "Healthy", 
                                     "NSCLC", "CRC"))) -> data.score.groups

comparisons = list(c("SCLC\npre", "SCLC\npost"), c("SCLC\npre", "Healthy"), 
                   c("SCLC\npre", "NSCLC"), c("SCLC\npre", "CRC"))

data.score.groups %>% 
  filter(!is.na(group)) %>% 
  boxplotWpoints(x = "group", y = "SCLC.n", fill = "group", 
                 ylab = "SCLC score", comparisons = comparisons) +
  scale_y_continuous(breaks = sclc.breaks, labels = sclc.lab) -> p
ggsave(paste0(figDirPaper, "figure1/tumor_fraction_boxplot.pdf"), p, width = 50,
       height = 60, units = "mm") 

# median of groups (add to text in manuscript)
data.score.groups %>% 
  group_by(group) %>% 
  summarise(median = median(SCLC.n))



