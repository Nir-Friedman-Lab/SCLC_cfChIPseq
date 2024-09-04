# NE score and signature --------------------------------------------------

# compute NE score --------------------------------------------------------
# function obtained from the ssGSEA package 
# https://github.com/alwaysblack777/Therapeutic-targeting-of-ATR-yields-durable-regressions-in-high-replication-stress-tumors/tree/main/ssGSEA
NE = read_xlsx(paste0(baseDir, "data_from_NIH/NCI_SCLC_transcriptome_subtype_20200910.xlsx"), 
               sheet = "NE signature genes", skip = 2)


gsets = list(SCLC_Neuroendocrine = NE$`High NE genes`, 
             SCLC_Non_Neuroendocrine = NE$`Low NE genes`)
# remove genes with high signal in healthy plasma 
gsets$SCLC_Neuroendocrine = 
  gsets$SCLC_Neuroendocrine[healthy.ref[gsets$SCLC_Neuroendocrine] < 10]
gsets$SCLC_Non_Neuroendocrine = 
  gsets$SCLC_Non_Neuroendocrine[healthy.ref[gsets$SCLC_Non_Neuroendocrine] < 30]
  
data.frame(t(GSVA::gsva(expr = SCLC_ChIP_fixed, gset.idx.list = gsets,
                        method = 'ssgsea', kcdf = "Gaussian", verbose=FALSE))) %>%
  mutate(NE_score = (SCLC_Neuroendocrine - SCLC_Non_Neuroendocrine)) / 2 ->
  data.ne.score.chip
data.frame(t(GSVA::gsva(expr = SCLC_RNA_fixed, gset.idx.list = gsets, # changed to rna.cpm from rna_data
                        method = 'ssgsea', kcdf = "Gaussian", verbose=FALSE))) %>%
  mutate(NE_score = (SCLC_Neuroendocrine - SCLC_Non_Neuroendocrine)) / 2 ->
  data.ne.score.rna

metadata %>%
  mutate(NE_SCORE = as.numeric(NE_SCORE), 
         NE_SCORE_chip = data.ne.score.chip[Sample_id,"NE_score"], 
         NE_SCORE_RNA = data.ne.score.rna[Sample_id, "NE_score"]) -> metadata

metadata %>% 
  ggplot(aes(NE_SCORE, NE_SCORE_RNA, text = Sample_id)) + 
  geom_point() + 
  geom_abline() -> p
ggplotly(p)
# metadata$NE_SCORE_chip = NA
# metadata[s.samples, "NE_SCORE_chip"] = data.ne.score.chip[s.samples, "NE_score"]

# NE score ChIP vs RNA - scatter  -----------------------------------------
# samps = rna.samples.chip.passQC; outname = "all_samples"
# samps = high.ctDNA.samples.wRNA; outname = "high_samples"
# samps = high.ctDNA.samples.wRNA.matched.time; outname = "high_samples_matched_time"
samps = matched.sclc_l.high.samples; outname = "fixed"
# samps = not_matched.sclc_l.all.samples

metadata %>% 
  filter(Sample_id %in% samps) %>% 
  scatter.plot("NE_SCORE_RNA", "NE_SCORE_chip", xlab = "tumor RNA", 
               ylab = "plasma cfChIP", cor = T) -> p
p
# validation.rna %>%
#   inner_join(patient.sample.table, 
#              join_by("patient_id" == "Collaborators.Subject.ID")) %>%
metadata %>%
  # filter(SCLC_score > low.score.cutoff) %>% 
  filter(Sample_id %in% samps) %>%
  # filter(!grepl("33-489", Sample_id)) %>%
  # filter(`Tumor_ctDNA_timing_match (determined based on treatment timing standpoint)` == "Yes") %>%
  ggplot(aes(NE_SCORE_RNA, NE_SCORE_chip)) +
  # geom_point(aes(alpha = SCLC_score, color = cohort), shape = 16, size = 1) +
  geom_point(shape = 16, size = 1) +
  labs(x = "tumor RNA", y = "plasma cfChIP", title = "NE score") + 
  # geom_text_repel(size = 2) + 
  scale_color_aaas() + 
  stat_cor(label.y.npc = 1, size = base_size/.pt) +
  geom_smooth(formula = y~x-1, method = "rlm", se = F, linewidth = .5) + 
  # geom_vline(xintercept = 0, linetype = "dashed", linewidth = .2) + 
  # geom_hline(yintercept = 0, linetype = "dashed", linewidth = .2) + 
  theme(aspect.ratio = 1, legend.key.size = unit(2, "mm"), 
        legend.position = c(.9, .2)) -> p
ggsave(paste0(figDirPaper, "figure4/ne.score_RNA_ChIP_", outname, "v1.pdf"), 
       width = 50, height = 50, units = "mm")
ggplotly(p)

 
## old version NE score in RNA and ChIP - based on Zhang 2018 https://www.mendeley.com/catalogue/29d708a9-9dac-3797-a991-e3307a62f78d/?utm_source=desktop&utm_medium=1.19.8&utm_campaign=open_catalog&userDocumentId=%7B728683f7-9890-3c91-a97f-7e3610740f4c%7D
# ne.score = read.csv(paste0(baseDir, "cfChIP-paper/Figures/NE_score_Zhang.csv"))
# rna.ne.score = sapply(samps, function(i) (cor(rna_data[ne.score$gene,i], ne.score$mean_exp_NE) -
#                                             cor(rna_data[ne.score$gene,i], ne.score$mean_exp_nonNE))/2)
# chip.ne.score = sapply(samps, function(i) (cor(log2(1+chip_data_all[ne.score$gene,i]), ne.score$mean_exp_NE) -
#                                              cor(log2(1+chip_data_all[ne.score$gene,i]), ne.score$mean_exp_nonNE))/2)
# data.ne = data.frame(samp = samps, rna = rna.ne.score, chip = chip.ne.score)
# p = scatter.plot(data.ne, x = "rna", y = "chip", xlab = "tumor RNA", ylab = "plasma cfChIP", cor = T)
# ggplot(data.ne, aes(rna, chip, label = substr(samp, 5, 13))) + geom_point() + geom_text_repel(size = 2)
ggsave(paste0(figDirPaper, "figure4/ne.score_RNA_ChIP_", outname, ".pdf"), width = 55, height = 55, units = "mm")


# NE signature heatmap ----------------------------------------------------
## NE signature in RNA and ChIP
NE = read_xlsx(paste0(baseDir, "data_from_NIH/NCI_SCLC_transcriptome_subtype_20200910.xlsx"), 
               sheet = "NE signature genes", skip = 2)
NE.genes = c(NE$`High NE genes`, NE$`Low NE genes`)
data.NE.genes = data.frame(row.names = NE.genes, NE = rep(c("high", "low"), 1, each = 25))

data.rna.ne = sweep(rna_data[NE.genes,], 1, rowMeans(rna_data[NE.genes,]), "-")
# row.ord = hclust(dist(data.rna.ne))$order
# col.ord = hclust(dist(data.rna.ne))$order
row.ord = ClusterMatrixRows(as.matrix(data.rna.ne))$order
# col.ord = ClusterMatrixRows(as.matrix(t(data.rna.ne)))$order
col.ord = order(as.numeric(metadata[colnames(rna_data), "NE_SCORE"]))
g.max = 2
g.min = -3
data.rna.ne[data.rna.ne > g.max] = g.max
data.rna.ne[data.rna.ne < g.min] = g.min
data.ne.m = reshape2::melt(t(data.rna.ne[row.ord, col.ord]))
names(data.ne.m) = c("name", "gene", "rna")


p = heatmap(data.ne.m, "name", "gene", "rna", lim = c(g.min,g.max), breaks = c(g.min,0,g.max), 
            leglab = "log2\n(1+RPKM)", leg = c("X1/8", "mean", "X4"))
p = p + labs(x = "tumor RNA samples", y = "NE signature genes")
ne.break = data.ne.m$gene[min(which(data.ne.m$gene %in% NE$`Low NE genes`))]
# p = p + geom_hline(yintercept = ne.break, color = "white")
p = p + theme(axis.text.x =  element_blank(), axis.ticks.x = element_blank(), 
              axis.text.y = element_text(size = 5), 
              legend.key.height = unit(3, 'mm'),legend.key.width = unit(2,"mm"),
              legend.title = element_text(size=base_size), legend.text = element_text(size=base_size))
ggsave(paste0(figDirPaper, "figureS4/ne_signature_rna.pdf"), p, width = 80, height = 80, units = "mm")

pheatmap(data.rna.ne[row.ord, col.ord], cluster_rows = F, cluster_cols = F, show_colnames = F, 
         annotation_row = data.NE.genes, annotation_col = metadata[,"NE_SCORE",drop = F], 
         fontsize = 4, color = cfm.color, border_color = NA, width = 4, height = 3,
         filename = paste0(figDirPaper, "figureS4/ne_signature_rna_v2.pdf"))


data.chip.ne = sweep(log2(1+chip_data_all[NE.genes,high.ctDNA.samples]), 1, 
                     rowMeans(log2(1+chip_data_all[NE.genes,high.ctDNA.samples])), "-")
row.ord = ClusterMatrixRows(as.matrix(data.chip.ne))$order
col.ord = order(as.numeric(metadata[high.ctDNA.samples, "NE_SCORE_chip"]))
data.chip.ne[data.chip.ne > 4] = 4
data.chip.ne[data.chip.ne < -4] = -4
pheatmap(data.chip.ne[row.ord, col.ord], cluster_rows = F, cluster_cols = F, fontsize_col = 2, 
         annotation_row = data.NE.genes, annotation_col = metadata[,"NE_SCORE_chip",drop = F], 
         fontsize= 4, color = cfm.color, border_color = NA, width = 4, height = 3,
         filename = paste0(figDirPaper, "figureS4/ne_signature_chip_v2.pdf"))



hist(healthy.ref$Gene.avg[NE.genes], breaks = "fd")
NE.genes.lowH = NE.genes[healthy.ref$Gene.avg[NE.genes] < 20]
d = log2(1+chip_data_all[NE.genes.lowH,high.ctDNA.samples])
data.chip.ne = sweep(d, 1, rowMeans(d), "-")
pheatmap(data.chip.ne, show_colnames = F, annotation_row = data.NE.genes, 
         breaks = seq(-2,2, length.out = 50))
row.ord = ClusterMatrixRows(as.matrix(data.chip.ne))$order
col.ord = ClusterMatrixRows(as.matrix(t(data.chip.ne)))$order
data.ne.m = reshape2::melt(t(data.chip.ne[row.ord, col.ord]))
names(data.ne.m) = c("name", "gene", "chip")
g.max = 2
g.min = -3
data.ne.m$chip[data.ne.m$chip > g.max] = g.max
data.ne.m$chip[data.ne.m$chip < g.min] = g.min
p = heatmap(data.ne.m, "name", "gene", "chip", lim = c(-3,2), breaks = c(-3,0,2), 
            leglab = "log2\n(1+counts)", leg = c("X1/8", "mean", "X4"))
p = p + labs(x = "plasma cfChIP samples", y = "NE signature genes")
p = p + theme(axis.text.x = element_text(angle = 90))
p = p + theme(axis.text.x =  element_blank(), axis.ticks.x = element_blank(), 
              axis.text.y = element_text(size = 5), 
              legend.key.height = unit(3, 'mm'),legend.key.width = unit(2,"mm"),
              legend.title = element_text(size=base_size), legend.text = element_text(size=base_size))
ggsave(paste0(figDirPaper, "figureS4/ne_signature_chip.pdf"), p, width = 80, height = 80, units = "mm")


