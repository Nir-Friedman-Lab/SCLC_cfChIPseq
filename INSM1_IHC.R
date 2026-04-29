# metadata %>%
#   filter(!is.na(INSM1_Hscore)) %>%
#   # filter(`cfDNA Patient identical/duplicate` == "Identical") %>%
#   filter(Sample_id %in% rna.samples) %>%
#   select(Sample_id, INSM1_Hscore, SCLC_score) %>%
#   mutate(insm1.chip = chip_data_all["INSM1", Sample_id], 
#          INSM1 = t(rna_data["INSM1", Sample_id])) %>%
#   # filter(SCLC_score > low.score.cutoff) %>%
#   ggplot(aes(INSM1_Hscore, insm1.chip, label = Sample_id)) +
#   geom_point(position = position_jitter(width = 10)) +
#   stat_cor(size = base_size/.pt) + 
#   # geom_text_repel(size = 2) +
#   labs(x = "INSM1 H-score (tumor)", 
#        y = expression(~italic(INSM1)~' reads (plasma)')) + 
#   geom_smooth(method = "rlm", se = F, linewidth = .5) -> p
# ggsave(paste0(figDirPaper, "figureS4/INSM1.pdf"), width = 45, height = 45, 
#        units = "mm")

# the correlation of INSM1_Hscore and INSM1 tumor RNA is not great

metadata.matching %>% 
  inner_join(metadata.biopsies) %>% 
  filter(!is.na(INSM1_Hscore)) %>% 
  mutate(insm1.chip = chip.matched["INSM1", Sample_id], 
         insm1.rna = rna_data["INSM1", Sample_id]) -> data_insm1

data_insm1 %>% 
  ggplot(aes(INSM1_Hscore, insm1.chip)) +
  geom_point(position = position_jitter(width = 10)) +
  stat_cor(size = base_size/.pt) + 
  # geom_text_repel(size = 2) +
  geom_smooth(method = "rlm", se = F, linewidth = .5) +
  labs(x = "INSM1 H-score (tumor)",
       y = expression(~italic(INSM1)~' cfChIP (plasma)')) +
  theme() -> p
ggsaveG(p, "~/Downloads/insm1_IHC_vs_cfchip", 60, 60, format = ".pdf")

data_insm1 %>% 
  ggplot(aes(INSM1_Hscore, insm1.rna)) +
  geom_point(position = position_jitter(width = 10)) +
  stat_cor(size = base_size/.pt) + 
  # geom_text_repel(size = 2) +
  geom_smooth(method = "rlm", se = F, linewidth = .5) +
  labs(x = "INSM1 H-score (tumor)",
       y = expression(~italic(INSM1)~'- RNA')) +
  theme() -> p
ggsaveG(p, "~/Downloads/insm1_IHC_vs_rna", 60, 60, format = ".pdf")

data_insm1 %>% 
  ggplot(aes(insm1.chip, insm1.rna)) +
  geom_point(position = position_jitter(width = 10)) +
  stat_cor(size = base_size/.pt) + 
  # geom_text_repel(size = 2) +
  geom_smooth(method = "rlm", se = F, linewidth = .5) +
  labs(x = expression(~italic(INSM1)~'- cfChIP'),
       y = expression(~italic(INSM1)~'- RNA')) +
  theme() -> p
ggsaveG(p, "~/Downloads/insm1_rna_vs_cfchip", 60, 60, format = ".pdf")


insm1_wins = which(TSS.windows$name == "INSM1")
d = log2(1+win_data_all[insm1_wins, ])
rownames(d) = insm1_wins
pheatmap::pheatmap(d, show_colnames = F, cluster_rows = F,
                   annotation_col = metadata %>% 
                     column_to_rownames("Sample_id") %>% 
                     select(SCLC.n), 
                   filename = "~/Downloads/insm1_wins.pdf", 
                   width = 10)

data_insm1 %>% 
  mutate(win1 = win_data_all[insm1_wins[1], Sample_id], 
         win2 = win_data_all[insm1_wins[2], Sample_id],
         win3 = win_data_all[insm1_wins[3], Sample_id],
         win4 = win_data_all[insm1_wins[4], Sample_id],
         win5 = win_data_all[insm1_wins[5], Sample_id]) -> data_insm1

cor(data_insm1 %>% 
      select(INSM1_Hscore, insm1.chip, win1, win2, win3, win4, win5)) -> insm1_wins_cor

insm1_wins_cor[insm1_wins_cor == 1] = NA
pheatmap::pheatmap(insm1_wins_cor, cluster_cols = F, cluster_rows = F, 
                   filename = "~/Downloads/insm1_wins_cor.pdf")
