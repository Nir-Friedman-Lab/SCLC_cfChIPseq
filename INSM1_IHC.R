metadata %>%
  filter(!is.na(INSM1_Hscore)) %>%
  # filter(`cfDNA Patient identical/duplicate` == "Identical") %>%
  filter(Sample_id %in% rna.samples) %>%
  select(Sample_id, INSM1_Hscore, SCLC_score) %>%
  mutate(insm1.chip = chip_data_all["INSM1", Sample_id], 
         INSM1 = t(rna_data["INSM1", Sample_id])) %>%
  # filter(SCLC_score > low.score.cutoff) %>%
  ggplot(aes(INSM1_Hscore, insm1.chip, label = Sample_id)) +
  geom_point(position = position_jitter(width = 10)) +
  stat_cor(size = base_size/.pt) + 
  # geom_text_repel(size = 2) +
  labs(x = "INSM1 H-score (tumor)", 
       y = expression(~italic(INSM1)~' reads (plasma)')) + 
  geom_smooth(method = "rlm", se = F, linewidth = .5) -> p
ggsave(paste0(figDirPaper, "figureS4/INSM1.pdf"), width = 45, height = 45, 
       units = "mm")

# the correlation of INSM1_Hscore and INSM1 tumor RNA is not great
