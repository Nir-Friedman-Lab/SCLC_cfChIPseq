## compare ChIP and RNA in two patients with high estimated tumor and different ChIP signal
# prepare data ------------------------------------------------------------
samp1 = "SCLC0147-1529"
samp2 = "SCLC0030-435"
data.comp = data.frame(chip1 = chip_data_all[matched.genes, samp1],
                       chip2 = chip_data_all[matched.genes, samp2], 
                       chipH = healthy.ref[matched.genes],
                       rna1 = rna_data[matched.genes, samp1],
                       rna2 = rna_data[matched.genes, samp2])
# use gene body coordinates of pou2f3 and atoh1
data.comp["POU2F3", "chip1"] = colSums(as.matrix(win_data_all[pou2f3.wins,samp1]))
data.comp["POU2F3", "chip2"] = colSums(as.matrix(win_data_all[pou2f3.wins,samp2]))
data.comp["ATOH1", "chip1"] = win_data_all[atoh1.wins,samp1]
data.comp["ATOH1", "chip2"] = win_data_all[atoh1.wins,samp2]
high.chip1 = which(data.comp$chip1/data.comp$chip2 > 2 & data.comp$chip1 > 5 & 
                     data.comp$chipH < .5)
high.chip2 = which(data.comp$chip2/data.comp$chip1 > 2 & data.comp$chip2 > 5 & 
                     data.comp$chipH < .5)
notsig.chip = setdiff(1:nrow(data.comp), c(high.chip1, high.chip2))
data.comp$gene = rownames(data.comp)

# statistics - scatter plot  ---------------------------------
# count gene in the 4 quarters 
n.s1.p = sum(log2((1+data.comp[high.chip1,"rna1"]) / 
                    (1+data.comp[high.chip1,"rna2"])) > 1)
n.s1.n = sum(log2((1+data.comp[high.chip1,"rna1"]) / 
                    (1+data.comp[high.chip1,"rna2"])) < -1)
n.s2.p = sum(log2((1+data.comp[high.chip2,"rna1"]) / 
                    (1+data.comp[high.chip2,"rna2"])) > 1)
n.s2.n = sum(log2((1+data.comp[high.chip2,"rna1"]) / 
                    (1+data.comp[high.chip2,"rna2"])) < -1)

ggplot(data.comp, aes(log2((1+chip1)/(1+chip2)), log2((1+rna1)/(1+rna2)), 
                      alpha = ifelse(rna1+rna2 > 3,1,.7), label = gene)) + 
  geom_point(shape = 16, data = data.comp[high.chip1, ], color = "red", size = .25) +
  geom_point(shape = 16, data = data.comp[high.chip2, ], color = "blue", size =.25) +
  geom_vline(size = base_line_size, xintercept = 0) + 
  geom_hline(size = base_line_size, yintercept = 0) +
  geom_hline(size = base_line_size, yintercept = -1, linetype = "dashed") + 
  geom_hline(size = base_line_size, yintercept = 1, linetype = "dashed") +
  labs(x = "ChIP - log2 ratio", y = "RNA - log2 ratio") +
  ylim(c(-6,6)) + 
  guides(alpha = "none") + 
  annotate(geom = "text", x=c(-4,-4,4,4), y=c(6,-6,-6,6), 
           label = paste(c(n.s2.p,n.s2.n,n.s1.n,n.s1.p), "genes"), 
           vjust="inward",hjust="inward", color="black", size = base_size/.pt) -> p
ggsave(paste0(figDirPaper, "figureS3/rna_chip_pair_ratio_v1.pdf"), p, 
       width = 60, height = 60, units = "mm")

# ChIP and RNA comparison - scatter plot  ---------------------------------
# add names to SCLC cannonical genes
data.comp[sclc.genes, ] # consider removing genes that are high in healthy
selected.genes = c("POU2F3", "ASCL1", "CHGA", "DLL3", "SOX2", "FOXA2", "BCL2")

mValrna = ceiling(max(quantile(data.comp$rna1,0.98), quantile(data.comp$rna2, 0.98)))
data.comp %>%
  filter(rna1 > 0 | rna2 > 0) %>%
  ggplot(aes(rna1, rna2, label = gene)) + 
  geom_point(data = data.comp[notsig.chip,], shape = 16, size=0.2, color="gray80") +
  geom_point(data = data.comp[high.chip1,], shape = 16, size=0.3, color="red") + 
  geom_point(data = data.comp[high.chip2,], shape = 16, size=0.3, color="blue") +
  coord_fixed(ratio=1,xlim=c(0,mValrna), ylim=c(0,mValrna)) + 
  scale_x_continuous(trans = "log1p", breaks=c(0,1,3,10,30,100,300,1000)) +
  scale_y_continuous(trans = "log1p", breaks=c(0,1,3,10,30,100,300,1000)) +
  geom_abline(slope = 1, intercept = 0, color="black") +
  geom_density2d(colour="black", show.legend = FALSE, bins=10, size = .1) + 
  geom_label_repel(data = data.comp[selected.genes,], size = base_size/.pt, 
                   label.padding = unit(1, "pt"), label.size = NA) +
  labs(x = "", y = samp2, title = "tumor RNA-seq") + 
  theme(plot.title = element_text(size = base_size)) -> p.rna

mValchip = ceiling(max(quantile(data.comp$chip1,0.995), quantile(data.comp$chip2, 0.995)))
data.comp %>%
  filter(chip1 > 0 | chip2 > 0) %>%
  ggplot(aes(chip1, chip2, label = gene)) +
  geom_point(data = data.comp[notsig.chip,], shape = 16, size=0.2, color="gray80") +
  geom_point(data = data.comp[high.chip1,], shape = 16, size=0.3, color="red") +
  geom_point(data = data.comp[high.chip2,], shape = 16, size=0.3, color="blue") +
  coord_fixed(ratio=1,xlim=c(0,mValchip), ylim=c(0,mValchip)) +
  scale_x_continuous(trans = "log1p", breaks=c(0,1,3,10,30,100,300)) +
  scale_y_continuous(trans = "log1p", breaks=c(0,1,3,10,30,100,300)) +
  geom_abline(slope = 1, intercept = 0, color="black") + 
  geom_density2d(colour="black", show.legend = FALSE, bins=10, size = .1) +
  labs(x = samp1, y = samp2, title = "plasma ChIP-seq") + 
  geom_label_repel(data = data.comp[selected.genes,], size = base_size/.pt, 
                   label.padding = unit(1, "pt"), label.size = NA) + 
  theme(plot.title = element_text(size = base_size)) -> p.chip

p = plot_grid(p.rna, p.chip, ncol = 1, align = "v")
p = rasterise(p)
ggsave(paste0(figDirPaper, "figure3/rna_chip_pair_comparison.pdf"), p, 
       width = 45, height = 90, units = "mm")

