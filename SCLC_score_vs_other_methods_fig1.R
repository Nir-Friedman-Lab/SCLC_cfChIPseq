# SCLC score vs. ctDNA, CT volume and RECIST ------------------------------
metadata %>% 
  left_join(read_xlsx(paste0(paperDir, "data/SCLC_cfChIP_DELFI_TF.xlsx")) %>% 
              filter(Timepointmatch == "Yes") %>%
              select(cfChiPseq_ID, DELFI_TF, Timepointmatch) %>% 
              rename(Sample_id = cfChiPseq_ID, 
                     ctDNA_delfi = DELFI_TF)) %>% 
  # filter(Sample_id %in% s.samples) %>% 
  select(Sample_id, SCLC.n, ctDNA, cfDNA, CT_volume_ALL, CTC, RECIST, 
         ctDNA_delfi, Timepointmatch) -> 
  data.score

data.score %>% 
  scatter.plot("ctDNA_delfi", "ctDNA", 
               xlab = "ctDNA fraction (delfi)",
               ylab = "ctDNA fraction (ichorCNA)", 
               title = paste0("n=", sum(!is.na(data.score$ctDNA) &
                                                 !is.na(data.score$ctDNA_delfi)))) +
  # geom_text_repel(aes(label = Sample_id), size = 4/.pt) +
  geom_abline(lwd = .2) + 
  theme(aspect.ratio = 1) -> p
ggsave(paste0(figDirPaper, "figure1/delfi_vs_ctDNA.pdf"),  p, 
       width = 40, height = 40, units = "mm")
ggsave(paste0(figDirPaper, "figure1/delfi_vs_ctDNA.png"),  p, 
       width = 40, height = 40, units = "mm", dpi = 500)

data.score %>% 
  scatter.plot("ctDNA_delfi", "SCLC.n",
               ylab = "SCLC score",
               xlab = "ctDNA fraction (delfi)", 
               title = paste0("n=", sum(!is.na(data.score$SCLC.n) &
                                                 !is.na(data.score$ctDNA_delfi)))) + 
  scale_y_continuous(breaks = sclc.breaks, labels = sclc.lab) -> p
ggsave(paste0(figDirPaper, "figure1/delfi_vs_est_tumor.pdf"),  p, 
       width = 40, height = 40, units = "mm")
ggsave(paste0(figDirPaper, "figure1/delfi_vs_est_tumor.png"),  p, 
       width = 40, height = 40, units = "mm", dpi = 500)
  
data.score %>% 
  scatter.plot("SCLC.n", "ctDNA_delfi", xlab = "", 
             ylab = "ctDNA fraction (DELFI)", 
             title = paste0("n=", sum(!is.na(data.score$ctDNA)))) + 
  theme(axis.text.x = element_blank(), 
        plot.margin = margin(b = -2, unit = "cm")) -> p1

data.score %>% 
  scatter.plot("SCLC.n", "CT_volume_ALL", xlab = "", 
             ylab = expression(CT~volume~(cm^{3})), trans.y = "log10", 
             title = paste0("n=", sum(!is.na(data.score$CT_volume_ALL)))) + 
  theme(axis.text.x = element_blank(), 
        plot.margin = margin(b = -2, t = -2, unit = "cm")) -> p2

data.score %>% 
  scatter.plot("SCLC.n", "RECIST", xlab = "SCLC score", 
             ylab = "RECIST (cm)", trans.y = "log10", 
             title = paste0("n=", sum(!is.na(data.score$RECIST)))) + 
  scale_x_continuous(breaks = sclc.breaks, labels = sclc.lab) + 
  theme(plot.margin = margin(t =-2, unit = "cm")) -> p3
p = plot_grid(p1, p2, p3, ncol = 1, align = "hv")
ggsave(paste0(figDirPaper, "figure1/est_tumor_vs_ctDNA_CT_RECIST.pdf"),  p, 
       width = 40, height = 120, units = "mm")

# SCLC score vs. cfDNA and CTC --------------------------------------------
data.score %>% 
  scatter.plot("SCLC.n", "cfDNA", xlab = "", 
             ylab = "cfDNA (ng/ml)", 
             title = paste0("n=", sum(!is.na(data.score$cfDNA)))) + 
  theme(axis.text.x = element_blank()) -> p1

data.score %>% 
  scatter.plot("SCLC.n", "CTC", xlab = "SCLC score", 
             ylab = "EpCAM+ CTC", trans.y = "log10", 
             title = paste0("n=", sum(!is.na(data.score$CTC)))) + 
  scale_x_continuous(breaks = sclc.breaks, labels = sclc.lab) -> p2
p = plot_grid(p1, p2, ncol = 1, align = "hv")
ggsave(paste0(figDirPaper, "figureS1/est_tumor_vs_cfDNA_CTC.pdf"),  p, 
       width = 45, height = 90, units = "mm")

# ctDNA and cfDNA  --------------------------------------------------------
data.score %>% 
  scatter.plot("cfDNA", "ctDNA", xlab = "cfDNA (ng/ml)", 
                 ylab = "ctDNA fraction", 
               title = paste0("n=", sum(!is.na(data.score$cfDNA) &
                                          !is.na(data.score$ctDNA)))) -> p
ggsave(paste0(figDirPaper, "figureS1/ctDNA_vs_cfDNA.pdf"), 
       p, width = 45, height = 45, units = "mm")

data.score %>% 
  rename(SCLC_score = SCLC.n) -> fig1eS1de
# SCLC score vs. response -------------------------------------------------
metadata %>%
  filter(Sample_id %in% s.samples & Responder_NonResponder != "NA") %>%
  mutate(response = factor(Responder_NonResponder, 
                           levels = c("NonResponder", "Responder"), 
                           labels = c("non-responder", "responder"))) -> data.resp

data.resp %>% 
  select(Sample_id, PatientID, Responder_NonResponder, SCLC.n) %>% 
  rename(SCLC_score = SCLC.n) -> figS1g

data.resp %>% 
  boxplotWpoints("response", "SCLC.n", fill = "response", xlab = "", 
                   ylab = "SCLC score", plot_stat = F) + 
  scale_y_continuous(breaks = sclc.breaks, labels = sclc.lab) + 
  scale_color_aaas() + 
  scale_fill_aaas() -> p
ggsave(paste0(figDirPaper, "figureS1/tumor_load_vs_response.pdf"), p, width = 50, 
       height = 50, units = "mm")


# SCLC score vs. timepoint ------------------------------------------------
selected.timepoints = c("Pre-treatment", "Post-treatment", 
                        "Ongoing-response", "Disease-progression")
metadata %>%
  filter(Sample_id %in% s.samples) %>%
  filter(Timepoint %in% selected.timepoints) %>%
  select(c(SCLC.n, Timepoint, PatientID)) %>%
  mutate(Timepoint = factor(Timepoint, levels = selected.timepoints, 
                            labels = c("Pre", "Post", "Ongoing", "Progression"))) -> figS1f

figS1f %>% 
  mutate(SCLC.n = jitter(SCLC.n, 100)) %>%
  mutate(SCLC.n = if_else(SCLC.n > 1, 1, SCLC.n), 
         SCLC.n = if_else(SCLC.n < 0, 0, SCLC.n)) %>% 
  group_by(PatientID, Timepoint) %>%
  slice_max(n = 1, order_by = SCLC.n) %>%
  group_by(PatientID) %>% 
  filter(n() > 1) %>% 
  ggplot(aes(Timepoint, SCLC.n, group = PatientID, text = PatientID)) + 
  geom_line(color = "gray", linewidth = .2) +
  geom_point(shape = 16, alpha = .7, size = .8) +
  scale_y_continuous(breaks = sclc.breaks, labels = sclc.lab) + 
  labs(x = "treatment stage", y = "SCLC score") -> p
ggsave(paste0(figDirPaper, "figureS1/tumor_load_vs_progression_v1.pdf"), p, 
       width = 50, height = 50, units = "mm")
ggplotly(p)
