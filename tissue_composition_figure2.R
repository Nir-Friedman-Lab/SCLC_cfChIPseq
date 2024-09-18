# tissue compositoin based on deconvolution with roadmap atlas 
# prepare data ------------------------------------------------------------
signatures = read.csv(paste0(baseDir, "Output/H3K4me3/signatures.csv"), 
                      row.names = 1)
############### temporary patch (till I rerun the signatures on the full data set) ##############
signatures.healthy = read.csv(paste0("~/gavriel.fialkoff@mail.huji.ac.il - Google Drive/Shared drives/Friedman Lab Shared Drive/BloodChIP/Analysis/Projects/NIH_SCLC/Output/H3K4me3/SCLC_Healthy_March_2022_sig-Counts.csv"), row.names = 1)
signatures.healthy = signatures.healthy[-grep("PS0011-1", rownames(signatures.healthy))[2],]
# rownames(signatures.healthy) = sub("_.*", "", rownames(signatures.healthy))
signatures.healthy = signatures.healthy[h.samples,]
qc.healthy = read.csv(paste0("~/gavriel.fialkoff@mail.huji.ac.il - Google Drive/Shared drives/Friedman Lab Shared Drive/BloodChIP/Analysis/Projects/NIH_SCLC/Output/H3K4me3/all_Mar-2022_qc.csv"), row.names = 1)
qc.healthy = qc.healthy[-grep("PS0011-1", rownames(qc.healthy))[2],]

############################ till here ############################

signatures = rbind(signatures, signatures.healthy)
# remove redundant and problematic tissues and rename some others
signatures %>%
  mutate(concentration = rbind(qc[,"Total.uniq.est", drop = F], 
                               qc.healthy[,"Total.uniq.est", drop = F])[rownames(.),] / 1e7) %>%
  mutate(sample = sub("_.*", "", rownames(.))) %>%
  filter(sample %in% c(s.samples, h.samples)) %>%
  left_join(estimated.tumor, join_by(sample)) %>% 
  # remove redundant and problematic tissues and rename some others
  select(!c(Leukocytes, Lymphocytes, Placenta)) %>%
  rename("B cell" = "B.Cells", "Megakaryocytes" = "Megakaryocyte") %>%
  # mutate(group = sample.annotation[sample, "group"],
  #        SCLC = estimated.tumor[sample, "SCLC.n"]) %>%
  melt(id.vars = c("sample", "group", "concentration", "SCLC.n"), 
       variable.name = "tissue") %>%
  mutate(abs.value = value * concentration) -> data.tiss.abs

# boxplot all tissues absolute counts (normalized to estimated concentrati --------
boxplotWOpoints(data.tiss.abs, "tissue", "abs.value", "group") + 
  labs(x = "", y = "reads/kb") + 
  coord_flip()  + 
  ylim(c(0,80)) -> p
ggsave(paste0(figDirPaper, "figureS2/tissue_signatures_abs_boxplot.pdf"), p, 
       width = 80, height = 85, units = "mm")

# boxplot selected high tissues -----------------------------------------
high.tissues = c("Neutrophils", "Monocytes", "Megakaryocytes", "Lung", "Brain", 
                 "B cell")
boxplotWOpoints(data.tiss.abs[data.tiss.abs$tissue %in% high.tissues,], 
                "tissue", "abs.value", "group") + 
  labs(x = "", y = "reads/kb") + 
  coord_flip() + 
  ylim(c(0,80)) -> p
ggsave(paste0(figDirPaper, "figure2/high_tissue_signatures_abs_boxplot.pdf"), p, 
       width = 55, height = 55, units = "mm")

# tissue signature vs. SCLC score -----------------------------------------
data.tiss.abs %>% 
  filter(tissue %in% high.tissues) %>%
  mutate(tissue = factor(tissue, levels = high.tissues)) %>% 
  ggplot(aes(x = SCLC.n, y = value, colour = group)) + 
  geom_point(shape = 16, size = .25) + 
  scale_color_manual(values = group.colors, limits = force) + 
  scale_x_continuous(breaks = sclc.breaks, labels = sclc.lab) + 
  # scale_y_continuous(trans = "log1p", breaks = c(0,2,4,6,8,10,12,14,16,18)) + 
  guides(color = "none") +
  facet_wrap(~tissue, scales='free_y', nrow = 3, dir = "v", ) + 
  labs(x = "SCLC score", y = "reads/signature", color = "") + 
  theme(aspect.ratio = 1, strip.background = element_blank(), 
        plot.background = element_blank(), panel.spacing=unit(1,"mm")) +
  stat_cor(color = "black", size = base_size/.pt, label.y.npc = 1) -> p
ggsave(paste0(figDirPaper, "figureS2/tumor_vs_tissue_signature.pdf"), p, 
       width = 70, height = 110, units = "mm")


# p.margin = margin(t=1, b=-3,r=-2,l=-3)
# p1 = scatter.plot(data.tiss[data.tiss$tissue == "Neutrophils",], "SCLC", "value", 
#                   color = "group", ylab = "reads/signature", size = .25, title = "Neutrophils") +
#   theme(plot.margin = margin(t=1, b=-3,r=-2,l=0), axis.text.x = element_blank())
# p2 = scatter.plot(data.tiss[data.tiss$tissue == "Monocytes",], "SCLC", "value", 
#                   color = "group", size = .25, title = "Monocytes") +
#   theme(plot.margin = p.margin, axis.text.x = element_blank())
# p3 = scatter.plot(data.tiss[data.tiss$tissue == "Megakaryocytes",], "SCLC", "value", 
#                   color = "group", size = .25, title = "Megakaryocytes") +
#   theme(plot.margin = margin(t=1, b=-3,r=-1,l=0), axis.text.x = element_blank())
# p4 = scatter.plot(data.tiss[data.tiss$tissue == "Lung",], "SCLC", "value", 
#                   ylab = "reads/signature", color = "group", size = .25, title = "Lung") +
#   theme(plot.margin = margin(t=1, b=-3,r=-2,l=0))
# # p4 = p4 + scale_x_continuous(breaks = sclc.breaks, labels = sclc.lab)
# p5 = scatter.plot(data.tiss[data.tiss$tissue == "Brain",], "SCLC", "value", 
#                   color = "group", size = .25, title = "Brain") +
#   theme(plot.margin = p.margin)
# # p5 = p5 + scale_x_continuous(breaks = sclc.breaks, labels = sclc.lab)
# p6 = scatter.plot(data.tiss[data.tiss$tissue == "B cell",], "SCLC", "value", 
#                   color = "group", size = .25, title = "B cell") +
#   theme(plot.margin = margin(t=1, b=-3,r=-2,l=0))
# # p6 = p6 + scale_x_continuous(breaks = sclc.breaks, labels = sclc.lab)
# p = plot_grid(p1, p2, p3, p4, p5, p6, nrow = 2, align = "hv")
# ggsave(paste0(figDirPaper, "figureS2/tumor_vs_tissue_signature_v2.pdf"), p, width = 80, height = 55, units = "mm")

# effect of brain metastisis on signature (not is use) -------------------------
data.tiss.abs %>% 
  left_join(metadata %>% 
              rename(Brain_met = names(.)[grepl("BrainMet", names(.))]) %>% 
              select(Brain_met, Sample_id), 
            join_by("sample" == "Sample_id")) %>% 
  filter(tissue == "Brain") %>% 
  mutate(Brain_met = factor(Brain_met, labels = c("yes", "no"))) %>% 
  filter(!is.na(Brain_met)) %>% 
  ggplot(aes(Brain_met, value)) +
  geom_boxplot(outliers = F, lwd = base_line_size) +
  geom_point(size = .5, position = position_jitter(width = .2)) + 
  stat_compare_means(label = "p.signif", comparisons = list(c("yes", "no")), 
                     size = base_size/.pt, lwd = base_line_size) + 
  labs(x = "Brain Metastasis", y = "Brain signature (reads/kb)") -> p
ggsave(paste0(figDirPaper, "figureS2/brain_metastasis_sig.png"), p, width = 50, 
       height = 70, unit = "mm")

data.tiss.abs %>% 
  left_join(metadata %>% 
              rename(liver_met = names(.)[grepl("LiverMet", names(.))]) %>% 
              select(liver_met, Sample_id), 
            join_by("sample" == "Sample_id")) %>% 
  filter(tissue == "Liver") %>% 
  mutate(liver_met = factor(liver_met, labels = c("yes", "no"))) %>% 
  filter(!is.na(liver_met)) %>% 
  ggplot(aes(liver_met, value)) +
  geom_boxplot(outliers = F, lwd = base_line_size) +
  geom_point(size = .5, position = position_jitter(width = .2)) + 
  stat_compare_means(label = "p.signif", comparisons = list(c("yes", "no")), 
                     size = base_size/.pt, lwd = base_line_size) + 
  labs(x = "Liver Metastasis", y = "Liver signature (reads/kb)") -> p
ggsave(paste0(figDirPaper, "figureS2/liver_metastasis_sig.png"), p, width = 50, 
       height = 70, unit = "mm")
