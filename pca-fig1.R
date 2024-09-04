# PCA healthy SCLC - all genes  ------------------------------------------------------------
pca.samples = c(h.samples, s.samples)
res.pca = prcomp(t(chip_data_all[,pca.samples]))
fviz_eig(res.pca, geom = "bar", barfill = "black", barcolor = "white", 
         ggtheme = base_theme) + 
  theme(plot.title = element_blank()) + 
  labs(x = "Principal components", y = "% explained variance") -> p
ggsave(paste0(figDirPaper, "figureS1/pca_scree_plot.pdf"), p, width = 50, height = 50, units = "mm")

data.frame(pc1 = res.pca$x[,1], 
           pc2 = res.pca$x[,2], 
           estimated.tumor[pca.samples,]) -> data.pca
data.pca %>% 
  ggplot(aes(pc1, pc2, color = group, text = rownames(data.pca))) +
  geom_point(aes(color = group, alpha = SCLC.n), shape = 16, size = 1) +
  geom_point(shape = 1, aes(color = group), stroke = .3, size = 1) + 
  scale_color_manual(values = group.colors, limits = force) + 
  labs(x = "PC1", y = "PC2", alpha = "SCLC score", color = "") +
  theme(legend.key.size = unit(2, "mm")) -> p
ggsave(paste0(figDirPaper, "figure1/pca_sclc_healthy.pdf"), p, width = 70, 
       height = 50, units = "mm")
ggplotly(p)

data.pca %>% 
  scatter.plot("pc1", "SCLC.n", color = "group", xlab = "PC1", 
             ylab = "SCLC score", alpha = .8) +
  scale_y_continuous(breaks = sclc.breaks, labels = sclc.lab) +
  theme(panel.background = element_blank(), panel.border = element_blank(), 
        panel.ontop = T) -> p
ggsave(paste0(figDirPaper, "figureS1/pca_tumor_frac.pdf"), p, width = 53, 
       height = 53, units = "mm")
