# load data ---------------------------------------------------------------
atacDir = "~/Documents/SCLC_data_from_NIH_box/ATAC-seq/"
histDir = paste0(atacDir, "Length_hist/")
atac.samples = read.table(paste0(atacDir, "atac_samples.txt"))$V1
atac.samples.all = read.table(paste0(atacDir, "atac_samples_all.txt"))$V1
atacFigDir = paste0(figDirPaper, "ATAC_seq/")
atac.qc = read.csv(paste0(atacDir, "Output/H3K4me3/atac_qc.csv"))
sort(atac.qc$X.Signal.at.TSS)
atac.samples.passQC = atac.qc$X[atac.qc$X.Signal.at.TSS > 30]
# atac.pairs = data.frame(cfchip = c("SCLC0043-518_NA_H3K4me3-38_26032021-13", 
#                                    "SCLC0046-572_NA_H3K4me3-38_26032021-13",
#                                    "SCLC0019-235_NA_H3K4me3-38_26032021-13", 
#                                    "SCLC0059-764_NA_H3K4me3-39_26032021-13"), 
#                         atac = c("2781080_S30_L003", 
#                                  "2792010_S0_L001", 
#                                  "2792020_S32_L003", 
#                                  "3049430_S42_L003"), 
#                         comments = c("A_N", "A", "Y", "NA"))
atac.pairs = read.csv(paste0(atacDir, "atac_chip_mapping.csv"))
atac.pairs %>%
  mutate(cfchip = sub("_.*", "", cfchip)) %>%
  filter(cfchip %in% s.samples, atac %in% atac.samples.all, 
         atac %in% atac.samples.passQC) %>%
  mutate(sclc.score = estimated.tumor[cfchip, "SCLC.n"], 
         subtype = metadata[cfchip, "NAPY_allcluster_expression"]) -> atac.pairs


# atac.pairs = atac.pairs[which(atac.pairs$cfchip %in% rownames(metadata) &
#                                 atac.pairs$atac %in% atac.samples.all &
#                                 atac.pairs$atac %in% atac.samples.passQC), ]
# atac.pairs$sclc.score = metadata[atac.pairs$cfchip, "SCLC_score"]
# atac.pairs$subtype = metadata[atac.pairs$cfchip, "NAPY_allcluster_expression"]
# atac.pairs$short.name = sub("_.*", "", atac.pairs$cfchip)


atac.data.genes = sapply(atac.samples.all, function(s) readRDS(
  paste0(atacDir, "Samples/H3K4me3/", s, ".rdata"))$GeneCounts.norm) # GeneCounts.QQnorm
# atac.data.wins = sapply(atac.samples.all, function(s) readRDS(
#   paste0(atacDir, "Samples/H3K4me3/", s, ".rdata"))$Counts.QQnorm)
# chip.data.normgenes = sapply(unique(c(h.samples, atac.pairs$cfchip)), function(s) readRDS(
#   paste0(rawdataDirNew, s, ".rdata"))$GeneCounts.norm)
atac.data.genes = atac.data.genes[Genes.notexcluded,]
chip.data.normgenes = chip_data_all[,unique(c(h.samples, atac.pairs$cfchip))]

healthy.normgenes = rowMeans(chip.data.normgenes[,h.samples])
point.size = .03

# compare ATACseq and cfChIPseq - global correlation ---------------------------
# g = sclc.genes
# g = diff.genes.common
# atac.pairs$cor.sclc = sapply(1:nrow(atac.pairs), function(i) 
#   cor(chip.data.normgenes[g,atac.pairs$cfchip[i]], atac.data.genes[g,atac.pairs$atac[i]]))
atac.pairs$cor = sapply(1:nrow(atac.pairs), function(i) 
  cor(chip.data.normgenes[,atac.pairs$cfchip[i]], 
      atac.data.genes[,atac.pairs$atac[i]]))
atac.pairs$cor.healthy = sapply(1:nrow(atac.pairs), function(i) 
  cor(healthy.normgenes, 
      atac.data.genes[,atac.pairs$atac[i]]))

atac.pairs %>% 
  ggplot(aes(sclc.score, cor)) + 
  geom_point(size = 1) +
  geom_point(mapping = aes(sclc.score, cor.healthy), color = "#bdbdbd", 
             size = 1) + 
  geom_smooth(mapping = aes(sclc.score, cor.healthy), color = "#bdbdbd", 
              method = "rlm", se = F, linetype = "dashed", linewidth = .5) +
  labs(x = "SCLC score", y = "correlation (cfChIP-seq vs. ATAC-seq)") -> p
ggsave(paste0(atacFigDir, "chip_atac_correlation.pdf"), p, 
       width = 50, height = 50, units = "mm")


# pairwise comparisson ----------------------------------------------------
atac.pairs[which.min(atac.pairs$cor.healthy),]
qplot(atac.data.genes[, "2792040_S0_L001"], healthy.ref)
which(atac.data.genes[, "2792040_S0_L001"] > 3000 & healthy.ref < 200)
df = data.frame(cbind(chip.data.normgenes[,unique(atac.pairs$cfchip)], 
                      atac.data.genes[,unique(atac.pairs$atac)], 
                      healthy.normgenes),
                check.names = F)
df$sclc_gene = factor(1*rownames(chip_data_all) %in% diff.genes.common)
rownames(df)[df$`SCLC0043-518_NA_H3K4me3-38_26032021-13` > 50 & df$`2781080_S30_L003` < 1]
sort(rownames(df)[df$`2781080_S30_L003` > 50 & df$`SCLC0043-518_NA_H3K4me3-38_26032021-13` < 1])


for (s in 1:nrow(atac.pairs)) {
  print(s)
  name = atac.pairs$cfchip[s]
  ctDNA = round(atac.pairs$sclc.score[s], digits = 2)
  scatter.plot(df,x = atac.pairs$cfchip[s], y = atac.pairs$atac[s], 
               color = "sclc_gene", trans.x = "log1p", trans.y = "log1p", 
               shape = 19, xlab = "cfChIP-seq (sample)", 
               ylab = "ATAC-seq (sample)", size = point.size, 
               title = paste0(name, ": ", "SCLC score = ", ctDNA)) + 
    geom_point(data = df[diff.genes.common, ], size = point.size) + 
    scale_color_manual(values = list("1" = group.colors$SCLC, "0" = "gray")) +
    scale_x_continuous(trans = "log1p", breaks = c(0,1,3,10,30,100,300)) + 
    scale_y_continuous(trans = "log1p", breaks = c(0,1,3,10,30,100,300)) + 
    geom_density2d(color = "blue",  linewidth = .2, bins = 20) + 
    geom_abline(color = "blue", linewidth = .2) -> p
  p = rasterize(p, dpi = 500)
  ggsave(paste0(atacFigDir, name, "_cfchip_vs_atac.pdf"), p, width = 50, height = 50, units = "mm")
}

scatter.plot(df, x = "healthy.normgenes", y = "2781080_S30_L003", 
             color = "sclc_gene", trans.x = "log1p", trans.y = "log1p", 
             shape = 19, xlab = "cfChIP-seq (healthy baseline)", 
             ylab = "ATAC-seq (sample)", size = point.size, 
             title = paste0("SCLC0043-518: SCLC score = 1")) + 
  geom_point(data = df[diff.genes.common, ], size = point.size) + 
  scale_color_manual(values = list("1" = group.colors$SCLC, "0" = "gray")) +
  scale_x_continuous(trans = "log1p", breaks = c(0,1,3,10,30,100,300)) +
  scale_y_continuous(trans = "log1p", breaks = c(0,1,3,10,30,100,300)) +
  geom_density2d(color = "blue",  linewidth = .2, bins = 20) +
  geom_abline(color = "blue", linewidth = .2) -> p
p = rasterize(p, dpi = 500)
ggsave(paste0(atacFigDir, "SCLC0043-518_atac_vs_healthy.pdf"), p, width = 50, 
       height = 50, units = "mm")

# diff chip vs. diff atac -------------------------------------------------
# high.1 = rownames(chip_data_all)[which(df[,1] > 30 & df[,2] < 1)]
# high.2 = rownames(chip_data_all)[which(df[,2] > 30 & df[,1] < 1)]
# 

samp.1.ind = grep("46-572", atac.pairs$cfchip)
samp.2.ind = grep("43-518", atac.pairs$cfchip)
samp.1.chip = atac.pairs$cfchip[samp.1.ind]
samp.2.chip = atac.pairs$cfchip[samp.2.ind]
samp.1.atac = atac.pairs$atac[samp.1.ind]
samp.2.atac = atac.pairs$atac[samp.2.ind]
name.1 = atac.pairs$short.name[samp.1.ind]
name.2 = atac.pairs$short.name[samp.2.ind]

high.1 = rownames(chip_data_all)[which(df[,samp.1.chip] > 30 & 
                                         df[,samp.1.chip] / df[,samp.2.chip] > 3)]
high.2 = rownames(chip_data_all)[which(df[,samp.2.chip] > 30 & 
                                         df[,samp.2.chip] / df[,samp.1.chip] > 3)]

high.1.selected = rownames(chip_data_all)[which(df[,samp.1.chip] > 20 & 
                                                  df[,samp.2.chip] < 2 &
                                                  df[,samp.1.atac] > 20 &
                                                  df[,samp.2.atac] < 2)]
high.2.selected = rownames(chip_data_all)[which(df[,samp.2.chip] > 20 & 
                                                  df[,samp.1.chip] < 2 &
                                                  df[,samp.2.atac] > 20 &
                                                  df[,samp.1.atac] < 2)]
sink(paste0(atacFigDir, "genes.high.in.samples.txt"))
print(list(samples = c(samp.1.chip, samp.2.chip), samp.1.chip = high.1.selected, 
           samp.2.chip = high.2.selected))
sink()
gr = read.csv(paste0(atacFigDir, "ATAC_diffpeaks.csv"), comment.char = "#")
atac.genes = setdiff(gr$name, "ACTB")
# atac.genes = c(high.1.selected, high.2.selected)

# p = scatter.plot(df, samp.1.chip, samp.2.chip, xlab = name.1, ylab = name.2, size = point.size, 
#                  title = "cfChIP-seq", cor = F)

df %>%
  ggplot(aes(.data[[samp.1.chip]], .data[[samp.2.chip]])) + 
  geom_point(color = "gray", size = point.size) +
  scale_x_continuous(trans = "log1p", breaks = c(0,1,3,10,30,100,300)) +
  scale_y_continuous(trans = "log1p", breaks = c(0,1,3,10,30,100,300)) +
  geom_density2d(color = "blue", linewidth = .2, bins = 30) + 
  geom_abline(color = "blue", linewidth = .2) +
  geom_point(data = df[high.1,], size = point.size, mapping = 
               aes(.data[[samp.1.chip]], .data[[samp.2.chip]]), color = "red") +
  geom_point(data = df[high.2,], size = point.size, mapping = 
               aes(.data[[samp.1.chip]], .data[[samp.2.chip]]), color = "blue") +
  geom_label_repel(data = df[atac.genes,], size = base_size/.pt,
                  mapping = aes(.data[[samp.1.chip]], .data[[samp.2.chip]],
                                label = atac.genes), 
                  label.padding = unit(.5, "pt"), 
                  label.size = NA) +
  labs(x = samp.1.chip, y = samp.2.chip, title = "cfChIP-seq") -> p
p = rasterise(p, dpi = 500)
ggsave(paste0(atacFigDir, "cfchip.diff.pdf"), p, width = 50, height = 50, 
       units = "mm")

df %>%
  ggplot(aes(.data[[samp.1.atac]], .data[[samp.2.atac]])) + 
  geom_point(color = "gray", size = point.size) + 
  scale_x_continuous(trans = "log1p", breaks = c(0,1,3,10,30,100)) +
  scale_y_continuous(trans = "log1p", breaks = c(0,1,3,10,30,100,300)) +
  geom_density2d(color = "blue", linewidth = .2, bins = 30) +
  geom_abline(color = "blue", linewidth = .2) +
  geom_point(data = df[high.1,], size = point.size,
             mapping = aes(.data[[samp.1.atac]], .data[[samp.2.atac]]), 
             color = "red") + 
  geom_point(data = df[high.2,], size = point.size, 
             mapping = aes(.data[[samp.1.atac]], .data[[samp.2.atac]]), 
             color = "blue") + 
  labs(x = samp.1.chip, y = samp.2.chip, title = "ATAC-seq") +
  geom_label_repel(data = df[atac.genes,], size = base_size/.pt,
                   mapping = aes(.data[[samp.1.chip]], .data[[samp.2.chip]],
                                 label = atac.genes), 
                   label.padding = unit(.5, "pt"), 
                   label.size = NA) -> p
p = rasterise(p, dpi = 500)
ggsave(paste0(atacFigDir, "atac.diff.pdf"), p, width = 50, 
       height = 50, units = "mm")


