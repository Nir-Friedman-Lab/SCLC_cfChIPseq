# SCLC TFs (ASCL1, NEUROD1, YAP1, POU2F3, (ATOH1)--------------------------------------------------------------
sclc.subtypes = c("ASCL1", "NEUROD1", "YAP1", "POU2F3", "ATOH1")

# TF rna cutoffs  ---------------------------------------------------------
rna_data[sclc.subtypes,] %>% 
  as.data.frame() -> df
df[df > 150] = 150

df %>% 
  data.frame(check.names = F) %>% 
  rownames_to_column("gene") %>% 
  melt(id.vars = "gene") %>% 
  ggplot(aes(x = value)) +
  geom_histogram(bins = 50, color = "white", linewidth = base_line_size) +
  facet_wrap(~gene, scales = "free") -> p
ggsave("~/Downloads/TF_RNA_hist.png", p, width = 100, 
       height = 60, units = "mm", dpi = 500)
subtype.cutoff = list(ascl1 = 10, neurod1 = 60, yap1 = 50, pou2f3 = 40, atoh1 = 4)
# heatmap of TF in RNA (show that they are more of less mutually exclusive) ------------
metadata.matching %>% 
  filter(SCLC.state != "other") %>% 
  pull(Sample_id) -> s
rna.subtypes = log2(1+t(rna_data[sclc.subtypes, s]))
# rna.subtypes = sweep(rna.subtypes, 1, rowMedians(rna.subtypes), "-")
rna.subtypes = sweep(rna.subtypes, 2, log2(1+unlist(subtype.cutoff)), "-")

rna.subtypes %>% 
  as_tibble() %>% 
  mutate(Sample_id = rownames(rna.subtypes)) %>% 
  select(Sample_id, everything()) -> fig3a

metadata.matching %>% 
  filter(Sample_id %in% s) %>% 
  select(SCLC.state, cohort) -> rna.groups

hlp = heatmap_legend_param = list(
  at = c(-3, 0, 3),
  labels = c("x1/8", "median", "x8"),
  legend_gp = gp, 
  labels_gp = gp,
  title = "relative \n CPM",
  title_gp = gp,
  legend_height = unit(5, "mm"))

HeatmapAnnotation(df = rna.groups[,"SCLC.state", drop = F], which = "row",
                  show_annotation_name = F,
                  annotation_legend_param = list(legend_gp = gp, 
                                                 labels_gp = gp,
                                                 legend_width = unit(1, "mm"),
                                                 title = "tumor type",
                                                 title_gp = gp,
                                                 col_title = NULL)) -> ra

pdf(paste0(figDirPaper, "figure4/TF_RNA_wATOH1_v0.pdf"), width = 2.5, height = 3.1)
Heatmap(rna.subtypes, show_row_names = F, 
        # split = rna.groups[,"cohort"], 
        right_annotation = ra,
        cluster_columns = F,
        show_column_dend = F, show_row_dend = F, 
        heatmap_legend_param = hlp, cluster_row_slices = F,
        col=colorRamp2(c(-3, 0, 3),c(group.colors$hm.low, group.colors$hm.mid, 
                                     group.colors$hm.high)),
        name = "relative \n CPM", gap = unit(2,"pt"),
        column_names_gp = gp, row_title_gp = gp)
dev.off()

# boxplot of 5 TFs in healthy and SCLC ChIP (to demonstrate that  --------
bind_cols(sample.annotation,
          t(chip_data_all[sclc.subtypes, sample.annotation$Sample_id])) %>% 
  filter(Sample_id %in% sh.samples) %>% 
  pivot_longer(cols = all_of(sclc.subtypes), 
               values_to = "value", 
               names_to = "gene") %>% 
  mutate(gene = factor(gene, levels = sclc.subtypes)) -> fig3b

fig3b %>% 
  boxplotWOpoints("gene", "value", "group", ylab = "cfChIP reads", guides = T, 
                  stat_test = "wilcox.test") +
  labs(fill = "") + 
  theme(legend.position = c(.85, .7)) -> p
ggsave(paste0(figDirPaper, "figure4/TF_ChIP_boxplot.pdf"), p, width = 60, 
       height = 70, units = "mm")

# RNA ChIP scatter - SCLC transcription factors ---------------------------
chip.matched[sclc.subtypes,] %>% 
  data.frame(check.names = F) %>% 
  rownames_to_column("gene") %>% 
  melt(id.vars = "gene", value.name = "chip", 
       variable.name = "Sample_id") -> data.tf.chip 

rna_data[sclc.subtypes,] %>% 
  data.frame(check.names = F) %>% 
  rownames_to_column("gene") %>% 
  melt(id.vars = "gene", value.name = "rna", 
       variable.name = "Sample_id") -> data.tf.rna 

data.tf.chip %>% 
  left_join(data.tf.rna, join_by(gene, Sample_id)) %>% 
  left_join(metadata.matching, join_by(Sample_id)) %>% 
  filter(SCLC.state != "other", SCLC > 0.5, Matched != "No") %>%
  mutate(gene = factor(gene, levels = sclc.subtypes)) -> data.tf.both

data.tf.both %>% 
  ggplot(aes(rna, chip)) +
  geom_point(aes(color = SCLC.state)) +
  labs(x = "tumor RNA", y = "plasma cfChIP", color = "", shape = "") +
  facet_wrap(~gene, scales = "free")  +
  stat_cor(color = "black", size = base_size/.pt, label.y.npc = 1) + 
  scale_x_continuous(trans = "log1p", breaks = c(0,1,3,10,30,100),
                     expand = expansion(mult = .1)) +
  scale_y_continuous(trans = "log1p", breaks = c(0,1,3,10,30,100),
                     expand = expansion(mult = c(0.1, 0.1))) +
  theme(strip.background = element_blank(), 
        legend.position = c(.8, .3), aspect.ratio = 1,
        strip.text = element_text(face = "bold")) +
  scale_color_aaas() -> p
ggplotly(p)
ggsave(paste0(figDirPaper, "figure4/rna_chip_TF_scatter_high.pdf"), p, 
       width = 160, height = 90, units = "mm")

data.tf.both %>% 
  select(Sample_id, gene, chip, rna, SCLC.state) -> fig3c
  
data.tf.chip %>% 
  left_join(data.tf.rna, join_by(gene, Sample_id)) %>% 
  left_join(metadata.matching, join_by(Sample_id)) %>% 
  filter(SCLC.state != "other", SCLC > 0.1, Matched != "No") %>%
  mutate(gene = factor(gene, levels = c("ASCL1", "NEUROD1", "POU2F3", 
                                        "ATOH1", "YAP1"))) %>% 
  ggplot(aes(rna * SCLC, chip, color = cohort)) +
  geom_point(shape = 16) +
  labs(x = "normalized tumor RNA", y = "plasma cfChIP", color = "") +
  facet_wrap(~gene, scales = "free") + theme(aspect.ratio = 1) +
  stat_cor(color = "black", size = base_size/.pt, label.y.npc = 1) + 
  # scale_x_continuous(trans = "log1p", breaks = c(0,1,3,10,30,100),
  #                    expand = expansion(mult = .1)) +
  # scale_y_continuous(trans = "log1p", breaks = c(0,1,3,10,30,100),
  #                    expand = expansion(mult = .1)) +
  theme(strip.background = element_blank(), 
        legend.position = c(.8, .3),
        strip.text = element_text(face = "bold")) +
  scale_color_aaas() -> p
ggsave(paste0(figDirPaper, "figure4/rna_chip_TF_scatter_all", ".pdf"), p, 
       width = 160, height = 90, units = "mm")
ggplotly(p)
# pou2f3, atoh1 ChIP-RNA cor genebody vs. promoter ------------------------

lab.gene = c(pou2f3 = "POU2F3", atoh1 = "ATOH1")
lab.var = c(promoter = "promoter", genebody = "genebody")

data.frame(promoter = c(chip_data_all["POU2F3",s], 
                        atoh1.promoter = chip_data_all["ATOH1",s]), 
           genebody = c(chip.matched["POU2F3",s], 
                        chip.matched["ATOH1",s]),
           rna = c(t(rna_data["POU2F3",s]), t(rna_data["ATOH1",s])), 
           samp = s, 
           gene = c(rep("pou2f3", length(s)), rep("atoh1", length(s)))) %>%
  # mutate(rna = (2^rna) -1) %>%
  melt(id.vars = c("samp", "gene", "rna")) -> figS4b

figS4b %>% 
  ggplot(aes(rna, value)) +
  stat_cor(color = "black", size = base_size/.pt, label.y.npc = 1, 
           label.x.npc = .2) +
  geom_point(shape = 16) +
  labs(x = "tumor RNA (CPM)", y = "plasma cfChIP (reads)") +
  # p = p + facet_wrap(~gene, scales = "free") + theme(aspect.ratio = 1)
  scale_x_continuous(trans = "log1p", breaks = c(0,1,3,10,30,100)) + 
  scale_y_continuous(trans = "log1p", breaks = c(0,1,3,10,30,100)) +
  facet_grid(gene ~ variable, scale = "free_y", 
             labeller = labeller(variable = lab.var,
                                 gene = lab.gene)) +
  theme(aspect.ratio = 1, strip.text = element_text(face = "bold"),
        panel.border = element_rect(color = "black", fill = NA),
        panel.background = element_blank(), 
        strip.background = element_blank()) -> p
ggsave(paste0(figDirPaper, "figureS4/pou2f3_atoh1_promoter_vs_genebody.pdf"),p, 
       width = 70, height = 70, units = "mm")

# multiple subtypes -------------------------------------------------------
metadata.matching %>% 
  filter(SCLC.state != "other") %>% 
  pull(Sample_id) -> s

chip.matched["ATOH1", s] %>% quantile(seq(0,1,.1))
subtype.cutoff.chip = list(ASCL1 = 2, NEUROD1 = 4, POU2F3 = 3, ATOH1 = 3)
subtype.cutoff.rna = list(ASCL1 = 10, NEUROD1 = 10, POU2F3 = 8, ATOH1 = 1)
mult_subtypes = c("ASCL1", "NEUROD1", "POU2F3", "ATOH1")  
table(sweep(chip.matched[mult_subtypes, s], 1, unlist(subtype.cutoff.chip), ">") %>% 
  as.table() %>% 
  colSums())
table(sweep(rna_data[mult_subtypes,s], 1, unlist(subtype.cutoff.rna), ">") %>% 
  as.table() %>% 
  colSums())


pheatmap::pheatmap(log2(1+rna_data[mult_subtypes, s]), 
                   width = 8, height = 4, fontsize_col = 6,
                   treeheight_row = 0,
                   filename = paste0(figDirPaper, "JCI_resubmission/multiple_subtypes_rna.png"))
pheatmap::pheatmap(log2(1+chip.matched[mult_subtypes, s]), 
                   annotation_col = metadata.matching %>% 
                     column_to_rownames("Sample_id") %>% 
                     select(SCLC), 
                   treeheight_row = 0,
                   width = 10, height = 4, fontsize_col = 6,
                   filename = paste0(figDirPaper, "JCI_resubmission/multiple_subtypes_chip.png"))


