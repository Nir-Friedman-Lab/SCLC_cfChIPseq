## RNA-ChIP correlation  
# compute correlation  ----------------------------------------------------
recalculate = F
if (recalculate) {
  chip_data_matched = chip_data_all[,rna.samples.chip.passQC]
  rna_data_matched = as.matrix(rna_data[,rna.samples.chip.passQC])
  # chip_data_matched = chip_data_all[,rna.samples.chip.passQC]
  # rna_data_matched = as.matrix(2^rna_data[,rna.samples.chip.passQC]-1)
  gene.cor = sapply(matched.genes, function(g) cor.test(t(rna_data_matched[g, ]), 
                                                        t(chip_data_matched[g, ])))
  gene.cor.rand = sapply(matched.genes, function(g) 
    cor.test(t(rna_data_matched[g,]), 
             t(chip_data_matched[g, sample(colnames(chip_data_matched))])))
  gene.cor.high = sapply(matched.genes, function(g) 
    cor.test(t(rna_data_matched[g, high.ctDNA.samples.wRNA]),
             t(chip_data_matched[g, high.ctDNA.samples.wRNA])))
  gene.cor.rand.high = sapply(matched.genes, function(g) 
    cor.test(t(rna_data_matched[g, high.ctDNA.samples.wRNA]),
             t(chip_data_matched[g, sample(high.ctDNA.samples.wRNA)])))
  
  data.cor = data.frame(random = unlist(gene.cor.rand["estimate",]), all = unlist(gene.cor["estimate",]),
                        high = unlist(gene.cor.high["estimate",]), random.high = unlist(gene.cor.rand.high["estimate",]),
                        random.p = unlist(gene.cor.rand["p.value",]), all.p = unlist(gene.cor["p.value",]),
                        high.p = unlist(gene.cor.high["p.value",]), random.high.p = unlist(gene.cor.rand.high["p.value",]))
  rownames(data.cor) = sub(".cor", "", rownames(data.cor))
  data.cor = data.cor[rowSums(is.na(data.cor)) == 0,]
  data.cor$gene = rownames(data.cor)
  write.csv(data.cor, paste0(paperDir, "Figures/RNA_ChIP_correlation_v1.csv")) # changed to _v1
} else {
  data.cor = read.csv(paste0(paperDir, "Figures/RNA_ChIP_correlation_v1.csv"), row.names = 1)
}

# define high dynamic range genes and/or genes low in healthy   ------------------------------
s = rna.samples.chip.passQC
cor.g = data.cor$gene
# genes that are low in healthy and high dynamic range
low.healthy.genes = names(which(rowMeans(chip_data_all[cor.g,h.samples]) <= 5))
high.healthy.genes = names(which(rowMeans(chip_data_all[cor.g,h.samples]) > 5))
chip.high = log2(1+rowQuantiles(as.matrix(chip_data_all[cor.g,s]), probs = .95))
chip.low = log2(1+rowQuantiles(as.matrix(chip_data_all[cor.g,s]), probs = .05))
chip.range = chip.high - chip.low
# hist(chip.range, breaks = "fd", xlab = "ChIP - dynamic_range")
chip.cutoff = 2
rna.high = log2(1 + rowQuantiles(as.matrix(rna_data[cor.g,s]), probs = .95))
rna.low = log2(1 + rowQuantiles(as.matrix(rna_data[cor.g,s]), probs = .05))
rna.range = rna.high - rna.low
hist(rna.range, breaks = "fd", xlab = "dynamic_range_rna")
# rna.cutoff = log2(3) # gene where high rna is 3-fold more than low rna
rna.cutoff = 3 # changed May 2024
high.range.genes = names(which(chip.range > chip.cutoff & rna.range > rna.cutoff))
low.range.chip = names(which(chip.range < chip.cutoff))

# collect statistics  ------------------------------------------
data = data.cor; outname = "all_genes"; outw = 90; outh = 50
data = data.cor[high.range.genes,]; outname = "high_range"; outw = 140; outh = 70
# data = data.cor[intersect(high.range.genes, low.healthy.genes),]; outname = "high_range_low_h"
# data = data.cor[low.healthy.genes,]; outname = "low_h"

random.q = p.adjust(data[,"random.p"], method = "fdr")
random.sig.genes = data$gene[random.q < .05 & data$random > 0]
random.high.q = p.adjust(data[,"random.high.p"], method = "fdr")
random.high.sig = data$gene[random.high.q < .05 & data$random.high > 0]

all.q = p.adjust(data[,"all.p"], method = "fdr")
all.q.min = min(data$all[all.q < .05 & data$all > 0])
all.q.max = max(data$all[all.q < .05 & data$all > 0])
all.sig.genes = data$gene[all.q < .05 & data$all > 0]
high.q = p.adjust(data[,"high.p"], method = "fdr")
high.q.min = min(data$high[high.q < .05 & data$high > 0])
high.q.max = max(data$high[high.q < .05 & data$high > 0])
high.sig.genes = data$gene[high.q < .05 & data$high > 0]
c(length(high.sig.genes), length(all.sig.genes), length(intersect(high.sig.genes, all.sig.genes)))
# write.table(unique(c(high.sig.genes, all.sig.genes)),
# "~/BloodChIP/Analysis/Analysis-Projects/NIH_SCLC/high_cor_genes2.csv", quote = F,
# row.names = F, col.names = F)
# plot histogram of correlation -------------------------------------------
df.m = melt(data[grep("\\.p", names(data), invert = T)], id.vars = "gene", 
            variable.name = "group", value.name = "cor")
levels(df.m$group) = c("random\npermutation", "SCLC", "high SCLC samples", "random\nhigh")

# lab = data.frame(x = c(.9,-.5), y = c(Inf, Inf), label = c(paste0("FDR < 0.05\n", length(all.sig.genes), " genes"), 
#                                                      paste0(nrow(data), " genes")))
lab = data.frame(x = .8, y = Inf, 
                 label = c(paste0(length(all.sig.genes), "\ngenes")))

df.m %>% 
  ggplot(aes(x = cor)) +
  geom_histogram(data = df.m[df.m$group == "random\npermutation",], 
                 position="identity", binwidth = .03, 
                 size = base_line_size, fill = "gray", color = "black") + 
  geom_histogram(data = df.m[df.m$group == "SCLC",], alpha=0.5, 
                 position="identity", binwidth = .03, 
                 size = base_line_size, color = "black", fill = "blue") + 
  labs(x = "", y = "", fill = "",
       title = paste0("matched samples (n = ", length(rna.samples.chip.passQC), ")")) + 
  xlim(c(-1,1)) + 
  geom_vline(xintercept = all.q.min, linetype="dotted", color = "red") + 
  theme(plot.title = element_text(size = base_size)) +
  geom_text(data = lab, aes(x,y,label=label), 
            size = base_size/.pt, color="black",vjust="inward") + 
  geom_vline(xintercept = 0, size = .3, color = "black") -> p

df.m %>% 
  ggplot(aes(x = cor)) + 
  xlim(c(-1,1)) +
  geom_histogram(data = df.m[df.m$group == "random\nhigh",], position="identity", 
                 binwidth = .03, size = base_line_size, fill = "gray",
                 color = "black") + 
  geom_histogram(data = df.m[df.m$group == "high SCLC samples",], alpha=0.5, 
                 position="identity", binwidth = .03, 
                 size = base_line_size, fill = "blue", color = "black") + 
  geom_vline(xintercept = high.q.min, linetype="dotted", color = "red") +#+ geom_vline(xintercept = 0)
# p.high = p.high + labs(x = "correlation RNA (tumor) H3K4me3 (plasma)", y = "# genes", title = paste0("High SCLC samples (n = ", length(high.ctDNA.samples.wRNA), ")"), fill = "")
  labs(x = "", y = "", fill = "", 
       title = paste0("High SCLC samples (n = ", length(high.ctDNA.samples.wRNA), ")")) + 
    geom_text(data = data.frame(x =.8, y = Inf, 
                                label = paste0(length(high.sig.genes), "\ngenes")), 
                            aes(x,y,label = label), size = base_size/.pt, 
              color="black", vjust="inward") + 
    geom_vline(xintercept = 0, size = .3, color = "black") + 
    theme(plot.title = element_text(size = base_size)) -> p.high
# p0 = ggplot(df.m, aes(x = cor))
# p0 = p0 + geom_histogram(data = df.m[df.m$group == "random\nhigh",], position="identity", binwidth = .03, fill = "transparent")
# p0 = p0 + geom_histogram(data = df.m[df.m$group == "high SCLC samples",], position="identity", binwidth = .03, fill = "transparent")
# p0 = p0 + geom_vline(xintercept = 0, size = .3, color = "black") + xlim(c(-1,1))
# p0 = p0 + theme(axis.text=element_text(colour = "transparent"), axis.line=element_blank(), 
#                 axis.ticks = element_blank(), axis.title = element_text(colour = "transparent"), 
#                 panel.background = element_rect(fill='transparent'),
#                 plot.background = element_rect(fill='transparent', color=NA),
#                 panel.grid.major = element_blank(),
#                 panel.grid.minor = element_blank(),
#                 legend.background = element_rect(fill='transparent'),
#                 legend.box.background = element_rect(fill='transparent'))
# p.both = plot_grid(p, p.high, nrow = 1, labels = "temp", label_x = 1, label_y = 0.1)
p.both = grid.arrange(p, p.high, nrow = 1, 
                      left = textGrob("# genes", gp=gpar(fontsize = base_size), rot = 90), 
                      bottom = textGrob("correlation", gp=gpar(fontsize=base_size)), 
                      top = textGrob(paste0("RNA-seq (tumor) vs. cfChIP-seq (plasma). ", nrow(data), " genes"), gp=gpar(fontsize=base_size)))
# p.all = ggdraw(p.both) + draw_plot(p0)
# ggsave(paste0(figDirPaper, "figure3/RNA_ChIP_both-", outname, ".pdf"), p, width = 75, height = 120, units = "mm")
ggsave(paste0(figDirPaper, "figure3/RNA_ChIP_both-", outname, ".pdf"), 
       p.both, width = outw, height = outh, units = "mm")



# low correlation genes exploration  --------------------------------------
# attempt to understand the low correlation genes
data.low.cor = data.frame(gene = data.cor$gene, 
                          chip.dr = chip.range, 
                          rna.dr = rna.range, 
                          cor.observed = data.cor$all, 
                          cor.random = data.cor$random, 
                          healthy_mean = rowMeans(chip_data_all[data.cor$gene,h.samples]),
                          sclc_mean = rowMeans(chip_data_all[data.cor$gene, high.ctDNA.samples.wRNA]), 
                          rna_liver = rna_roadmap[data.cor$gene, "Adult_Liver"])

# data.low.cor = data.low.cor[high.range.genes,]
lab = c(cor.observed = "observed correlation", cor.random = "random permutation")
# data.low.cor$sh_ratio = log2(data.low.cor$sclc_mean/data.low.cor$healthy_mean)
# summary(lm(cor.observed ~ chip.dr + rna.dr + sh_ratio, data = data.low.cor))
# chip dynamic range is significantly possitivly correlated and healthy/sclc ratio is negativly correlated

data.low.cor %>% 
  filter(gene %in% high.range.genes) %>%
  mutate(sh_ratio = log2(sclc_mean/healthy_mean)) %>%
  melt(id.vars = c("gene", "chip.dr", "rna.dr"), 
       measure.vars = c("cor.observed", "cor.random"), variable.name = "type", 
       value.name = "cor") -> data.range.m
outname = "high_range_genes"

data.range.m %>% 
  ggplot(aes(chip.dr, cor)) + 
  geom_point(shape = 16, size = 1, alpha = .5) +
  geom_density_2d(size = 3*base_line_size) + 
  stat_cor(size = base_size/.pt, digits = .001, label.y.npc = 1) + 
  labs(x = "ChIP - dynamic range (log2 reads)", y = "correlation") + 
  geom_smooth(method = "lm", size = 3*base_line_size) + 
  facet_wrap(~type, labeller = labeller(type = lab)) + 
  theme(aspect.ratio = 1) -> p
ggsave(paste0(figDirPaper, "figureS3/dynamic_range_vs_cor_chip_", outname, ".pdf"), 
       p, width = 80, height = 45, units = "mm")

data.range.m %>%
  ggplot(aes(rna.dr, cor)) + 
  geom_point(shape = 16, size = 1, alpha = .5) + 
  geom_density_2d(bins = 10, size = 3*base_line_size) + 
  stat_cor(size = base_size/.pt, digits = .001) + 
  labs(x = "RNA - dynamic range (log2 CPM)", y = "correlation") + 
  geom_smooth(method = "lm", size = 3*base_line_size) + 
  facet_wrap(~type, labeller = labeller(type = lab)) + 
  theme(aspect.ratio = 1) -> p
ggsave(paste0(figDirPaper, "figureS3/dynamic_range_vs_cor_RNA_", outname, ".pdf"), 
       p, width = 80, height = 45, units = "mm")
# low dynamic range genes are less correlated

# data.ratio.m = melt(data.low.cor[high.range.genes,], id.vars = c("gene", "sh_ratio"), 
#                     measure.vars = c("cor.observed", "cor.random"), variable.name = "type", 
#                     value.name = "cor"); 
outname = "high_range_genes"

data.low.cor %>%
  mutate(sh_ratio = log2(sclc_mean/healthy_mean)) %>% 
  filter(gene %in% high.range.genes) %>%
  melt(id.vars = c("gene", "sh_ratio"), 
       measure.vars = c("cor.observed", "cor.random"), variable.name = "type", 
       value.name = "cor") %>%
  ggplot(aes(sh_ratio,cor)) + 
  geom_point(size = .2) + 
  geom_density_2d(bins = 10, size = 3*base_line_size) +
  stat_cor(size = base_size/.pt, digits = .001) +
  geom_smooth(method = "lm", size = 3*base_line_size) +
  facet_wrap(~type, labeller = labeller(type = lab)) +
# p = p + labs(x = expression(paste("log2 ", frac(healthy, SCLC))), y = "correlation")
  labs(x = "log2 (SCLC / healthy)", y = "correlation") + 
  theme(aspect.ratio = 1) -> p
ggsave(paste0(figDirPaper, "figureS3/sclc_healthy_ratio_vs_cor.pdf"), p, 
       width = 80, height = 45, units = "mm")

# write.table(all.genes, paste0(figDirPaper, "figureS3/high.cor.genes.csv"), quote = F, row.names = F, col.names = F)
# write.table(low.cor.genes, paste0(figDirPaper, "figureS3/low.cor.genes.csv"), quote = F, row.names = F, col.names = F)
# write.table(high.genes, paste0(figDirPaper, "figureS3/high.cor.genes.high_samps.csv"), quote = F, row.names = F, col.names = F)


# data.liver.m = reshape2::melt(data.low.cor[high.range.genes,], id.vars = c("gene", "rna_liver"), 
#                     measure.vars = c("cor.observed", "cor.random"), variable.name = "type", 
#                     value.name = "cor"); 
outname = "liver_genes"

data.low.cor %>%
  filter(gene %in% high.range.genes) %>%
  melt(id.vars = c("gene", "rna_liver"), 
       measure.vars = c("cor.observed", "cor.random"), variable.name = "type",
       value.name = "cor") %>%
  ggplot(aes(rna_liver,cor)) + 
  geom_point(size = .2) + 
  geom_density_2d(bins = 10, size = 3*base_line_size) +
  stat_cor(size = base_size/.pt, digits = .001) +
  geom_smooth(method = "lm", size = 3*base_line_size) + 
  facet_wrap(~type, labeller = labeller(type = lab)) + 
# p = p + labs(x = expression(paste("log2 ", frac(healthy, SCLC))), y = "correlation")
  labs(x = "expression in Roadmap liver (log2 RPKM)", y = "correlation") + 
  theme(aspect.ratio = 1) -> p
ggsave(paste0(figDirPaper, "figureS3/liver_vs_cor.pdf"), p, 
       width = 80, height = 45, units = "mm")

if (0) {
  g = "G6PC"
  df = data.frame(t(rna_data[g,high.ctDNA.samples.wRNA]), 
                  log2(1+t(chip_data_all[g, high.ctDNA.samples.wRNA])), liver = colSums(rna_data[liver.genes,high.ctDNA.samples.wRNA]))
  names(df) = c("rna", "chip", "liver")
  p = ggplot(df, aes(rna, chip, alpha = liver)) + geom_point() + labs(title = g, x = "RNA", y = "ChIP")
  p = p + stat_cor()
  ggsave(paste0(figDirPaper, "temp/", g, "_chip_rna_liver.png"),p)
} # examples for genes where liver contamination disterbs correlation 

# RNA - ChIP scatter of canonical genes  --------------------------------
genes.selected = c("BCL2", "NFIB", "SOX2"); outname = "_cannonical_genes.pdf"
# genes.selected = data.cor$gene[data.cor$high.p < .0001 & data.cor$high > .85 &
#     rowQuantiles(as.matrix(chip_data_all[data.cor$gene, high.ctDNA.samples.wRNA]), probs = .5) > 5]
ggplot(data.frame(gene = genes.selected,
                  h = rowMeans(chip_data_all[genes.selected,h.samples]),
                  s = rowMeans(chip_data_all[genes.selected,high.ctDNA.samples.wRNA])),
       aes(s,h, label = gene)) + geom_point() + geom_abline(slope = 1, intercept = 0) + geom_text_repel()
# genes.selected = c("SFTA3", "MYO10", "TSPYL5", "TIMM10"); outname = "high_cor_genes.pdf"

samp = high.ctDNA.samples.wRNA

# samp = samp[!samp %in% high.liver.samples] # TODO based on rna data exploratory - 
df.rna = data.frame(samp = samp, t(rna_data[genes.selected,samp]))
df.chip = data.frame(samp = samp, t(chip_data_all[genes.selected,samp]))
df.chip.m = melt(df.chip, id.vars = "samp", variable.name = "gene", value.name = "chip")
df.rna.m = melt(df.rna, id.vars = "samp", variable.name = "gene", value.name = "rna")
df.selected = cbind(df.chip.m, df.rna.m)
df.selected$chip = df.selected$chip
# df.selected$rna = 2^df.selected$rna
df.selected = df.selected[!duplicated(colnames(df.selected))]
df.selected %>% 
  ggplot(aes(rna, chip, label = substr(samp,3,11))) +
  geom_point(size = .5, shape = 16) + 
  labs(x = "tumor RNA (CPM)", y = "plasma cfChIP (promoter reads)") +
  # scale_y_continuous(trans = "log1p", breaks = c(5,10,50,100,300,500)) + 
  # scale_x_continuous(trans = "log1p", breaks = c(5,10,50,100,300,500)) +
  facet_wrap(~gene, scales = "free", ncol = 1) + 
  theme(aspect.ratio = 1) +
  stat_cor(color = "black", method = "spearman", size = base_size/.pt, 
           p.accuracy = 0.001, position = position_nudge_repel()) -> p
ggsave(paste0(figDirPaper, "figure3/rna_chip", outname), p, width = 45, 
       height = 120, units = "mm")



# compare between RNA-cfChIP correlation and roadmap RNA-ChIP correlation ----------------
chip_roadmap = readRDS("~/BloodChIP/Data/ExternalChIP/Analysis-Roadmap/Atlas.rds")
chip_roadmap = chip_roadmap$GeneCounts.QQnorm
colnames(chip_roadmap) = roadmap.names[colnames(chip_roadmap), "Universal_Human_Reference"]
chip_roadmap = chip_roadmap[,colnames(rna_roadmap)]
roadmap.cor = sapply(matched.genes, function(g) cor(t(rna_roadmap[g, ]), chip_roadmap[g, ]))

roadmap.cor = roadmap.cor[names(roadmap.cor) %in% rownames(data.cor)]
data.roadmap.cor = data.frame(roadmap = roadmap.cor, cfchip = data.cor$high); outname = "roadmap_vs_cfchip.png"
high.range.roadmap = matched.genes[which(rowMaxs(as.matrix(log2(1+chip_roadmap[matched.genes,]))) > chip.cutoff & 
                     rowMaxs(as.matrix(rna_roadmap[matched.genes,])) > rna.cutoff)]
high.range.roadmap = high.range.roadmap[high.range.roadmap %in% high.range.genes]
data.roadmap.cor = data.roadmap.cor[high.range.roadmap,]; outname = "roadmap_vs_cfchip_high_range.png"
p = ggplot(data.roadmap.cor, aes(roadmap, cfchip)) + geom_point(size = .3) + geom_density2d() + geom_abline()
p = p + labs(x = "Roadmap", y = "SCLC cfChIP", title = "RNA ChIP Correlation")
p = p + theme(aspect.ratio = 1)
ggsave(paste0(figDirPaper, "figureS3/SCLC_vs_Roadmap_cor/", outname), p, width = 5, height = 5)


# genes that have higher correlation in cfChIP 
G = rownames(data.roadmap.cor)[which(data.roadmap.cor$roadmap < 0.3 & data.roadmap.cor$cfchip > .8)]
G = G[which(rowMaxs(chip_roadmap[G,]) > 20 & rowMaxs(as.matrix(rna_roadmap[G,])) > 2)]

for (g in G) {
  print(g)
  sclc.g = data.frame(t(rna_data[g,high.ctDNA.samples.wRNA]), chip_data_all[g, high.ctDNA.samples.wRNA])
  names(sclc.g) = c("rna", "chip")
  tiss.g = data.frame(t(rna_roadmap[g, ]), chip_roadmap[g,])
  names(tiss.g) = c("rna", "chip")
  sclc.p = ggplot(sclc.g, aes(rna, chip)) + geom_point() + stat_cor(p.digits = 0) +
    labs(x = "RNA", y = "ChIP", title = paste(g, "- SCLC"))
  tiss.p = ggplot(tiss.g, aes(rna, chip)) + geom_point() + stat_cor() +
     labs(x = "RNA", y = "ChIP", title = paste(g, "- Roadmap"))
  p = plot_grid(sclc.p, tiss.p)
  ggsave(paste0(figDirPaper, "figureS3/SCLC_vs_Roadmap_cor/", g, ".png"), width = 6, height = 3)

}


