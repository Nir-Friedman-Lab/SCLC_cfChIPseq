# subtype labeling --------------------------------------------------------
# subtype.cutoff are set in SCLC_transcription_factors-fig4.R
metadata.matching %>% 
  mutate(is.ascl1 = rna_data["ASCL1", Sample_id] > subtype.cutoff$ascl1, 
         sclc_high = SCLC > .5) %>% 
  count(sclc_high, cohort, is.ascl1)

metadata.matching %>% 
  filter(SCLC > .05, cohort == "train") %>% 
  pull(Sample_id) -> samples.train

metadata.matching %>% 
  filter(SCLC > .05, cohort == "validation") %>% 
  pull(Sample_id) -> samples.validation

metadata.matching %>% 
  filter(rna_data["ASCL1", Sample_id] > subtype.cutoff$ascl1, 
         SCLC > .05, cohort == "train") %>% 
  pull(Sample_id) -> ascl1.samples.train

metadata.matching %>% 
  filter(rna_data["ASCL1", Sample_id] > subtype.cutoff$ascl1, 
         SCLC > .05, cohort == "validation") %>% 
  pull(Sample_id) -> ascl1.samples.validation

metadata.matching %>% 
  filter(rna_data["ASCL1", Sample_id] <= subtype.cutoff$ascl1, 
         SCLC > .05, cohort == "train") %>% 
  pull(Sample_id) -> nascl1.samples.train

metadata.matching %>% 
  filter(rna_data["ASCL1", Sample_id] > subtype.cutoff$ascl1, 
         SCLC > .05, cohort == "train") %>% 
  pull(Sample_id) -> ascl1.samples.train

metadata.matching %>% 
  filter(rna_data["NEUROD1", Sample_id] > subtype.cutoff$neurod1, 
         SCLC > .05, cohort == "train") %>% 
  pull(Sample_id) -> neurod1.samples.train

metadata.matching %>% 
  filter(rna_data["NEUROD1", Sample_id] > subtype.cutoff$neurod1, 
         SCLC > .05, cohort == "validation") %>% 
  pull(Sample_id) -> neurod1.samples.validation

metadata.matching %>% 
  filter(rna_data["NEUROD1", Sample_id] <= subtype.cutoff$neurod1, 
         SCLC > .05, cohort == "train") %>% 
  pull(Sample_id) -> nneurod1.samples.train

colnames(chip_data_all)
rowMeans(win_data_all[,sub("_.*", "", h.samples)]) %>% quantile()


# single gene classifier --------------------------------------------------
g = "ASCL1"; g.tit =  bquote(ASCL1^"+"); cu = subtype.cutoff$ascl1
g = "NEUROD1"; g.tit =  bquote(NEUROD1^"+"); cu =  subtype.cutoff$neurod1 # cu = 60;
g = "POU2F3"; g.tit =  bquote(POU2F3^"+"); cu = subtype.cutoff$pou2f3
g = "ATOH1"; g.tit =  bquote(ATOH1^"+"); cu = subtype.cutoff$atoh1 # cu = 4 ; subtype.cutoff$atoh1 = 4

metadata.matching %>% 
  filter(SCLC.state != "other") %>%
  mutate(rna = rna_data[g, Sample_id], 
         chip = chip.matched[g, Sample_id], 
         positive = factor(rna > cu, levels = c(T,F))) -> df.gene

df.gene %>% 
  ggplot(aes(SCLC, chip, color = positive)) +
  geom_point(size = .5) + 
  geom_smooth(formula = y~x-1, method = "rlm", linewidth = .2, se = F) +
  theme(legend.position = c(.2,.7), legend.key.size = unit(3, "mm")) + 
  labs(x = "SCLC-score", y = g, color = g.tit) + 
  scale_x_continuous(breaks = sclc.breaks, labels = sclc.lab) + 
  scale_color_aaas() -> p
ggsave(paste0(figDirPaper, "validation/", g, "_scatter.pdf"), p, width = 50, 
       height = 50, units = "mm")


metadata.matching %>% 
  mutate(ascl1.rna = rna_data["ASCL1", Sample_id], 
         neurod1.rna = rna_data["NEUROD1", Sample_id], 
         pou2f3.rna = rna_data["POU2F3", Sample_id], 
         atoh1.rna = rna_data["ATOH1", Sample_id], 
         ascl1.chip = chip.matched["ASCL1", Sample_id], 
         ascl1.score = colSums(win_data_all[ascl1.ttest.wins, Sample_id]), # ascl1.ttest.wins defined below
         neurod1.chip = chip.matched["NEUROD1", Sample_id], 
         pou2f3.chip = chip.matched["POU2F3", Sample_id], 
         atoh1.chip = chip.matched["ATOH1", Sample_id], 
         ascl1.pos = factor(ascl1.rna > subtype.cutoff$ascl1, levels = c(T,F)), 
         neurod1.pos = factor(neurod1.rna > subtype.cutoff$neurod1, levels = c(T,F)),
         pou2f3.pos = factor(pou2f3.rna > subtype.cutoff$pou2f3, levels = c(T,F)),
         atoh1.pos = factor(atoh1.rna > subtype.cutoff$atoh1, levels = c(T,F))) -> fig4

fig4 %>% 
  select(Sample_id, SCLC, Matched, SCLC.state, 
         ascl1.pos, ascl1.chip, neurod1.pos, neurod1.chip, pou2f3.pos, pou2f3.chip, atoh1.pos, atoh1.chip) %>% 
  rename(SCLC_score = SCLC) %>% 
  filter(SCLC.state != "other") -> fig4ab

fig4 %>% 
  filter(SCLC.state != "other", SCLC > 0.1, Matched != "No") -> data.roc


ascl1.roc = roc(data.roc$ascl1.pos, data.roc$ascl1.chip)
ascl1.auc = round(ascl1.roc$auc, digits = 2)
ascl1.score.roc = roc(data.roc$ascl1.pos, data.roc$ascl1.score)
ascl1.score.auc = round(ascl1.score.roc$auc, digits = 2)
neurod1.roc = roc(data.roc$neurod1.pos, data.roc$neurod1.chip)
neurod1.auc = round(neurod1.roc$auc, digits = 2)
pou2f3.roc = roc(data.roc$pou2f3.pos, data.roc$pou2f3.chip)
pou2f3.auc = round(pou2f3.roc$auc, digits = 2)
atoh1.roc = roc(data.roc$atoh1.pos, data.roc$atoh1.chip)
atoh1.auc = round(atoh1.roc$auc, digits = 2)

roc.labels = list(ascl1 = paste0("ASCL1 (AUC = ", round(ascl1.auc, digits = 2), ")"), 
                  neurod1 = paste0("NEUROD1 (AUC = ", round(neurod1.auc, digits = 2), ")"), 
                  pou2f3 = paste0("POU2F3 (AUC = ", round(pou2f3.auc, digits = 2), ")"), 
                  atoh1 = paste0("ATOH1 (AUC = ", round(atoh1.auc, digits = 2), ")"))
ggroc(list(ascl1 = ascl1.roc, neurod1 = neurod1.roc, 
           pou2f3 = pou2f3.roc, atoh1 = atoh1.roc)) +
  geom_abline(intercept = 1, slope = 1, lwd = base_line_size, 
              linetype = "dashed") +
  scale_color_manual(values = c(group.colors$ASCL1, group.colors$NEUROD1,
                                group.colors$POU2F3, "purple"), 
                     labels = roc.labels) +
  labs(color = "")  +
  theme(legend.position = c(.7, .2), aspect.ratio = 1, 
        legend.key.size = unit(3, "mm"), legend.background = element_blank()) -> p
ggsave(paste0(figDirPaper, "validation/ROC_all.pdf"), p, width = 60, height = 60, 
       units = "mm")

roc_df <- bind_rows(
  ggroc(ascl1.roc)$data %>% mutate(model = paste0("ASCL1 (AUC=", ascl1.auc, ")")),
  ggroc(neurod1.roc)$data %>% mutate(model = paste0("NEUROD1 (AUC=", neurod1.auc, ")")),
  ggroc(pou2f3.roc)$data %>% mutate(model = paste0("POU2F3 (AUC=", pou2f3.auc, ")")),
  ggroc(atoh1.roc)$data %>% mutate(model = paste0("ATOH1 (AUC=", atoh1.auc, ")")))

ggplot(roc_df, aes(x = specificity, y = sensitivity)) +
  geom_line(lwd = .2) +
  coord_equal() +
  scale_x_reverse(limits = c(1, 0)) +
  facet_wrap(~model, ncol = 2) +
  labs(x = "Specificity", y = "Sensitivity") +
  theme(strip.text = element_text(face = "bold"), 
        strip.background = element_blank()) -> p
ggsave(paste0(figDirPaper, "validation/ROC_facet.pdf"), p, width = 80, height = 80, 
       units = "mm")

data.roc %>% 
  filter(cohort == "train") -> data.roc.t
ascl1.roc.t = roc(data.roc.t$ascl1.pos, data.roc.t$ascl1.score)
ascl1.auc.t = round(ascl1.roc.t$auc, digits = 2)
neurod1.roc.t = roc(data.roc.t$neurod1.pos, data.roc.t$neurod1.chip)
neurod1.auc.t = round(neurod1.roc.t$auc, digits = 2)
pou2f3.roc.t = roc(data.roc.t$pou2f3.pos, data.roc.t$pou2f3.chip)
pou2f3.auc.t = round(pou2f3.roc.t$auc, digits = 2)
atoh1.roc.t = roc(data.roc.t$atoh1.pos, data.roc.t$atoh1.chip)
atoh1.auc.t = round(atoh1.roc.t$auc, digits = 2)

data.roc %>% 
  filter(cohort == "validation") -> data.roc.v
ascl1.roc.v = roc(data.roc.v$ascl1.pos, data.roc.v$ascl1.score)
ascl1.auc.v = round(ascl1.roc.v$auc, digits = 2)
neurod1.roc.v = roc(data.roc.v$neurod1.pos, data.roc.v$neurod1.chip)
neurod1.auc.v = round(neurod1.roc.v$auc, digits = 2)
pou2f3.roc.v = roc(data.roc.v$pou2f3.pos, data.roc.v$pou2f3.chip)
pou2f3.auc.v = round(pou2f3.roc.v$auc, digits = 2)
atoh1.roc.v = roc(data.roc.v$atoh1.pos, data.roc.v$atoh1.chip)
atoh1.auc.v = round(atoh1.roc.v$auc, digits = 2)

g = "ASCL1"; gene.t = ascl1.roc.t; gene.v = ascl1.roc.v; gene.t.auc = ascl1.auc.t; gene.v.auc = ascl1.auc.v
g = "NEUROD1"; gene.t = neurod1.roc.t; gene.v = neurod1.roc.v; gene.t.auc = neurod1.auc.t; gene.v.auc = neurod1.auc.v
g = "POU2F3"; gene.t = pou2f3.roc.t; gene.v = pou2f3.roc.v; gene.t.auc = pou2f3.auc.t; gene.v.auc = pou2f3.auc.v
g = "ATOH1"; gene.t = atoh1.roc.t; gene.v = atoh1.roc.v; gene.t.auc = atoh1.auc.t; gene.v.auc = atoh1.auc.v
roc.labels = list(gene.t = paste0("test (AUC = ", gene.t.auc, ")"), 
                  gene.v = paste0("validation (AUC = ", gene.v.auc, ")"))
ggroc(list(gene.t = gene.t, gene.v = gene.v)) +
  geom_abline(intercept = 1, slope = 1, size = base_line_size, 
              linetype = "dashed") +
  scale_color_manual(values = c("darkblue", "darkred"), 
                     labels = roc.labels) +
  labs(color = "", title = g)  +
  theme(legend.position = c(.7,.3), aspect.ratio = 1, 
        legend.key.size = unit(3, "mm"), legend.background = element_blank()) -> p
ggsave(paste0(figDirPaper, "validation/", g, "_ROC_test_validation.png"), p, width = 60, 
       height = 60, units = "mm", dpi = 500)



# leave one out cross validation classifier signature --------------------------
SCLCWins = readRDS("~/gavriel.fialkoff@mail.huji.ac.il - Google Drive/Shared drives/Friedman Lab Shared Drive/BloodChIP/Analysis/Projects/NIH_SCLC/new_data/SCLC_ChIPSamples.rds")

neurod1.ttest = readRDS(paste0(baseDir, "Classify/class-NEUROD1-50.rdata"))
neurod1.ttest.wins = as.integer(neurod1.ttest$Sig)
ascl1.ttest = readRDS(paste0(baseDir, "Classify/class-ASCL1-100.rdata"))
ascl1.ttest.wins = as.integer(ascl1.ttest$Sig)

g = "ASCL1"; w = ascl1.ttest.wins; s = c(ascl1.samples.train, ascl1.samples.validation)
g = "NEUROD1"; w = neurod1.ttest.wins; s = c(neurod1.samples.train, neurod1.samples.validation)

# gene.mat = sapply(names(SCLCWins), function(s) SCLCWins[[s]]$Counts.QQnorm)[w,]
gene.mat = win_data_all[w,s]
rownames(gene.mat) = TSS.windows[w]$name
rownames(gene.mat)[rownames(gene.mat) == "."] = ""
rownames(gene.mat)[is.na(rownames(gene.mat))] = ""
d = log2(1+gene.mat)
data.frame(Sample_id = colnames(d)) %>% 
  left_join(metadata.matching %>% select(SCLC, cohort, Sample_id), 
            by = "Sample_id") %>% 
  mutate(is.gene = factor(1* Sample_id %in% s)) %>% 
  column_to_rownames("Sample_id") -> col.a
col.ord = clusterSubgroups(d, col.a, "is.gene")
col.ord = rankSubgroups(col.a, "is.gene", "SCLC", decreasing = T, group_order = c("1","0"))
pheatmap::pheatmap(d[,col.ord], fontsize = 5, annotation_col = col.a, 
                   cluster_cols = F)

data.frame(Sample_id = colnames(d)) %>% 
  left_join(metadata.matching %>% select(SCLC, cohort, Sample_id), 
            by = "Sample_id") %>% 
  rename("SCLC-score" = SCLC) %>% 
  mutate(pos = factor(Sample_id %in% s,
                      # labels = c("ASCL1" %p% supsc('+'), "ASCL1" %p% supsc('-')),
                      labels = c(paste0(g, "high"), paste0(g,"low")),
                      levels = c(T, F)), 
         "ASCL1 (RNA)" = log2(1+rna_data["ASCL1", Sample_id]), 
         "ASCL1 (ChIP)" = chip.matched["ASCL1", Sample_id],
         "ASCL1 score (ChIP)" = colSums(win_data_all[ascl1.ttest.wins,
                                                     Sample_id])) %>% 
  column_to_rownames("Sample_id") -> col.a
top.a = HeatmapAnnotation(df = col.a %>% select(-c(pos, cohort)) %>% 
                            select(c("ASCL1 (RNA)", "ASCL1 (ChIP)", 
                                     "ASCL1 score (ChIP)", "SCLC-score")), 
                          annotation_name_gp = gp, 
                          gp = gp, 
                          show_legend = T, 
                          annotation_legend_param = list(legend_gp = gp, 
                                                         labels_gp = gp,
                                                         title_gp = gp,
                                                         direction = "horizontal",
                                                         # legend_label_gp= gp,
                                                         legend_height = unit(4, "mm"))) 

data.frame(ind = 1:nrow(d), 
           gene = rownames(d)) %>% 
  filter(!is.na(gene), !gene == ".", !gene == "") %>% 
  group_by(gene) %>% 
  slice_sample(n = 1) %>% 
  ungroup() %>% 
  slice_sample(n = 25) -> data.gene 
# row.lab = rownames(d)
# row.ind = which(!is.na(row.lab) & row.lab != "." & row.lab != "")
# row.ind = sample(row.ind, 20)
# row.lab[row.ind]
gene.a = HeatmapAnnotation(foo = anno_mark(labels_gp = gp, 
                                           at = data.gene$ind, 
                                           labels = data.gene$gene,
                                           lines_gp = gpar(lwd = .2), 
                                           link_width = unit(2, "mm")
), 
which = "row")

pdf(file = paste0(figDirPaper, "validation/", g, "_sig_heatmap.pdf"), width = 6, 
    height = 5)
Heatmap(d,
        heatmap_width = unit(60,"mm"), 
        heatmap_height = unit(100, "mm"),
        name = " ",
        column_title_gp = gp, 
        column_dend_gp = gp, 
        column_gap = unit(1, "pt"),
        use_raster = T,
        show_row_dend = F,
        show_row_names = F,
        right_annotation = gene.a,
        heatmap_legend_param = list(legend_gp = gp,
                                    labels_gp = gp,
                                    legend_gp = gp,
                                    legend_label_gp = gp, 
                                    legend_height = unit(4, "mm")),
        show_column_names = F,
        column_split = col.a[,c("cohort", "pos")],
        row_names_gp = gpar(fontsize = 5),
        column_names_gp = gp,
        column_order = col.ord,
        cluster_columns = F,
        top_annotation = top.a)
dev.off()

g = "ASCL1"
train = ascl1.samples.train
validation = ascl1.samples.validation
sig = ascl1.ttest
g.tit = bquote(ASCL1^"+")

g = "NEUROD1"
train = neurod1.samples.train
validation = neurod1.samples.validation
sig = neurod1.ttest
g.tit = bquote(NEUROD1^"+")

metadata.matching %>% 
  left_join(data.frame(score = c(sig$XX, sig$YY)) %>% 
              mutate(Sample_id = rownames(.)), by = "Sample_id") %>% 
  filter(SCLC > 0.05) %>% 
  mutate(gene.pos = factor(Sample_id %in% c(train, validation),
                           levels = c(T, F))) ->
  gene.data  

gene.data %>%
  filter(cohort == "train") %>% 
  ggplot(aes(SCLC, score, color = gene.pos)) + 
  geom_point(size = .5) +
  scale_color_aaas() +
  scale_x_continuous(breaks = sclc.breaks, labels = sclc.lab) + 
  geom_smooth(method = "rlm", linewidth = .2, se = F) +
  theme(legend.position = c(.2,.7), legend.key.size = unit(3, "mm"), 
        legend.background = element_blank()) + 
  labs(x = "SCLC-score", y = paste(g, "score"), color = g.tit) -> p
ggplotly(p)
ggsave(paste0(figDirPaper, "validation/", g, ".scatter.pdf"), width = 50, 
       height = 50, units = "mm")

gene.data %>% 
  filter(Sample_id %in% samples.train) -> roc.train
roc.t = roc(roc.train$gene.pos, roc.train$score/roc.train$SCLC)

gene.data %>% 
  filter(Sample_id %in% samples.validation) -> roc.validation
roc.v = roc(roc.validation$gene.pos, roc.validation$score/roc.validation$SCLC)

roc.labels = list(train = paste0("train (AUC = ", round(roc.t$auc, digits = 2), ")"), 
                  validation = paste0("validation (AUC = ", round(roc.v$auc, digits = 2), ")"))

ggroc(list(train = roc.t, validation = roc.v)) +
  geom_abline(intercept = 1, slope = 1, size = base_line_size, 
              linetype = "dashed") +
  scale_color_manual(values = c("darkblue", "darkred"), 
                     labels = roc.labels) +
  labs(color = "", title = paste(g, "score"))  +
  theme(legend.position = c(.7, .2), aspect.ratio = 1, 
        legend.key.size = unit(3, "mm")) -> p
ggsave(paste0(figDirPaper, "validation/ROC_", g, ".pdf"), p, width = 50, height = 50, 
       units = "mm")  


# enrichment of signature genes -------------------------------------------
ascl1.sig.genes = unique(TSS.windows[ascl1.ttest.wins]$name)
ascl1.sig.genes = ascl1.sig.genes[ascl1.sig.genes != "." & !is.na(ascl1.sig.genes)]
write.table(ascl1.sig.genes, paste0(figDirPaper, "validation/ascl1.sig.genes.txt"), 
            quote = F, row.names = F, col.names = F)
enricher(genes = ascl1.sig.genes, output.path = paste0(figDirPaper, "validation/"), 
         outputfilename = "ascl1.enrichment.csv")


# heatmap binary RNA and ChIP score --------------------------------------------
subtypes.four = c("ASCL1", "NEUROD1", "POU2F3", "ATOH1")
metadata.matching %>% 
  filter(SCLC > 0.05) %>% 
  pull(Sample_id) -> samples.above.5pct

data.frame(t(rna_data[subtypes.four,])) %>% 
  filter(rownames(.) %in% samples.above.5pct) %>% 
  mutate(ASCL1 = factor(ASCL1 > subtype.cutoff$ascl1, levels = c(T,F), labels = c(1,0)),
         NEUROD1 = factor(NEUROD1 > subtype.cutoff$neurod1, levels = c(T,F), labels = c(1,0)),
         POU2F3 = factor(POU2F3 > subtype.cutoff$pou2f3, levels = c(T,F), labels = c(1,0)),
         ATOH1 = factor(ATOH1 > subtype.cutoff$atoh1, levels = c(T,F), labels = c(1,0))) ->
  # ASCL1 = 1*(ASCL1 > subtype.cutoff$ascl1),
  # NEUROD1 1*(NEUROD1 > subtype.cutoff$neurod1), 
  # POU2F3 = 1*(POU2F3 > subtype.cutoff$pou2f3), 
  # ATOH1 = 1*(ATOH1 > subtype.cutoff$atoh1) -> 
  rna.binary
row.ord = hclust(dist(rna.binary))$order
col.ord = hclust(dist(t(rna.binary)))$order


subtypes.four = c("ASCL1", "NEUROD1", "POU2F3", "ATOH1")
data.frame(t(log2(1+rna_data[subtypes.four,]))) %>%
  filter(rownames(.) %in% samples.above.5pct) %>% 
  sweep(2, log2(1+c(subtype.cutoff$ascl1, subtype.cutoff$neurod1, 
                    subtype.cutoff$pou2f3, subtype.cutoff$atoh1)), "-") -> 
  rna.binary.n

metadata.matching %>% 
  column_to_rownames("Sample_id") %>% 
  select(SCLC) -> row.a
right.a = HeatmapAnnotation(df = row.a[rownames(rna.binary[,col.ord]),], 
                            which = "row", annotation_label = "SCLC score", 
                            show_legend = F)

data.frame(t(chip.matched[subtypes.four,])) %>% 
  filter(rownames(.) %in% samples.above.5pct) %>% 
  mutate(ASCL1 = colMeans(sapply(rownames(.), function(s) 
    SCLCWins[[s]]$Counts.QQnorm[ascl1.ttest.wins]))) -> chip.tf 

chip.tf %>% 
  rownames_to_column("Sample_id") %>% 
  left_join(metadata.matching %>% select(Sample_id, SCLC), 
            by = "Sample_id") %>% 
  filter(SCLC > .05) %>% 
  select(-Sample_id) %>% 
  mutate(ASCL1 = ASCL1/SCLC, 
         NEUROD1 = NEUROD1/SCLC, 
         POU2F3 = POU2F3/SCLC, 
         ATOH1 = ATOH1/SCLC) %>% 
  select(-SCLC) -> chip.tf.n

h1 = Heatmap(rna.binary[row.ord, col.ord],
             right_annotation = right.a, 
             cluster_rows = F, 
             cluster_columns = F, 
             show_row_names = F, 
             name = "RNA", 
             col = colorRampPalette(c("#00A7FE", "black", "#E54C00"))(2))
h2 = Heatmap(rna.binary.n[row.ord, col.ord],
             right_annotation = right.a, 
             cluster_rows = F, 
             cluster_columns = F, 
             show_row_names = F, 
             name = "RNA norm", 
             col = colorRampPalette(c("#00A7FE", "black", "#E54C00"))(2))
h3 = Heatmap(log2(1+chip.tf[row.ord, col.ord]),
             cluster_rows = F, 
             cluster_columns = F, 
             show_row_names = F, 
             name = "ChIP", 
             col = colorRampPalette(c("#00A7FE", "#E54C00"))(50))
h4 = Heatmap(log2(1+chip.tf.n[row.ord, col.ord]),
             cluster_rows = F, 
             cluster_columns = F, 
             show_row_names = F, 
             name = "ChIP norm", 
             col = colorRampPalette(c("#00A7FE", "#E54C00"))(50))
pdf(file = paste0(figDirPaper, "validation/rna_vs_chip.pdf"), width = 12, 
    height = 6)
# pdf(file = paste0("~/Downloads/rna_vs_chip.pdf"), width = 8, height = 6)
h1 + h2 + h3 + h4
dev.off()



