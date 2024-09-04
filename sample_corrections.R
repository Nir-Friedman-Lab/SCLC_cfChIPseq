chip_data_all = readRDS(paste0(baseDir, "chip_data_all.rds"))
win_data_all = readRDS(paste0(baseDir, "win_data_all.rds"))

rawdataDirNew = "~/Documents/SCLC_data/RDS/H3K4me3/"
qc.new = read.csv("~/Documents/SCLC_data/Output/H3K4me3/qc_all.csv", row.names = 1)
low.qc.samples = rownames(qc.new)[qc.new$TSS < cutoff.yeild | 
                                qc.new$X.Signal.at.TSS < cutoff.signal]
low.qc.samples.short = sub("_.*", "", low.qc.samples)
passed.qc.samples = setdiff(rownames(qc.new), low.qc.samples)
passed.qc.samples.short = sub("_.*", "", passed.qc.samples)


# samples in metadata but not in atlas - all good
metadata %>% 
  filter(!Sample_id %in% sub("_.*", "", c(low.qc.samples, passed.qc.samples))) %>% 
  # filter(!Sample_id %in% fixed_names$Old_name) %>% 
  select(Sample_id)

# sample in atlas but not in metadata
colnames(chip_data_all)[!colnames(chip_data_all) %in% rownames(metadata) &
                          grepl("SCLC", colnames(chip_data_all))] 
# 200-100 -> 19-100; 200-133 -> 19-133; 243-283 -> 50-283; 243-443 -> 50-443;
# SCLC0253-736 and SCLC0341-265 ??


# colnames(rna_data) = sub("SCLC0019-372", "SCLC0019-100", colnames(rna_data))

chip_data_updated = sapply(passed.qc.samples, function(s)
  readRDS(paste0(rawdataDirNew, s, ".rdata"))$GeneCounts.QQnorm)

chip_data_all = chip_data_all[, grep("SCLC", colnames(chip_data_all), invert = T)]
chip_data_all = cbind(chip_data_all, chip_data_updated)

chip_data_all = chip_data_all[,grep("PS0011-1_NA_H3K4me3-39_26032021-13", colnames(chip_data_all), invert = T)]
any(sort(table(sub("_.*", "", colnames(chip_data_all)))) > 1)
colnames(chip_data_all) = sub("_.*", "", colnames(chip_data_all))

chip_data_all = chip_data_all[Genes.notexcluded,]

win_data_updated = sapply(passed.qc.samples, function(s)
  readRDS(paste0(rawdataDirNew, s, ".rdata"))$Counts.QQnorm)
win_data_all = win_data_all[, grep("SCLC", colnames(win_data_all), invert = T)]
win_data_all = cbind(win_data_all, win_data_updated)
win_data_all = win_data_all[,grep("PS0011-1_NA_H3K4me3-39_26032021-13", colnames(win_data_all), invert = T)]
colnames(win_data_all) = sub("_.*", "", colnames(win_data_all))

############################# backup for the drive #############################
saveRDS(win_data_all, "~/Documents/SCLC_data_from_NIH_box/win_data_all.rds")
saveRDS(chip_data_all, "~/Documents/SCLC_data_from_NIH_box/chip_data_all.rds")
write.table(colnames(win_data_all), "~/Documents/SCLC_data_from_NIH_box/all_samples.txt", 
            quote = F, row.names = F, col.names = F)
win_data_all = readRDS("~/Documents/SCLC_data_from_NIH_box/win_data_all.rds")
chip_data_all = readRDS("~/Documents/SCLC_data_from_NIH_box/chip_data_all.rds")
################################################################################


metadata %>% mutate(short = sub("_.*", "", Sample_id)) %>%
  filter(! short  %in% colnames(chip_data_all)) %>%
  filter(! short %in% sub("_.*", "", low.qc.samples)) %>% 
  select(Sample_id)

metadata %>% 
  mutate(short = sub("_.*", "", Sample_id)) %>%
  filter(short %in% sub("_.*", "", passed.qc.samples)) %>%
  filter(Histology %in% sclc.histology, cohort == "train") %>%
  select(Sample_id) %>% 
  .$Sample_id -> s.samples

metadata %>% 
  mutate(short = sub("_.*", "", Sample_id)) %>%
  filter(short %in% sub("_.*", "", passed.qc.samples)) %>%
  filter(Histology %in% sclc.histology, cohort == "validation") %>%
  select(Sample_id) %>% 
  .$Sample_id -> v.samples

metadata %>% 
  mutate(short = sub("_.*", "", Sample_id)) %>%
  filter(short %in% sub("_.*", "", passed.qc.samples)) %>%
  filter(Histology %in% c("NSCLC")) %>%
  select(Sample_id) %>% 
  .$Sample_id -> l.samples

metadata %>% 
  mutate(short = sub("_.*", "", Sample_id)) %>%
  filter(short %in% sub("_.*", "", passed.qc.samples)) %>%
  filter(!Histology %in% c(sclc.histology, "NSCLC")) %>%
  select(Sample_id) %>% 
  .$Sample_id -> o.samples

metadata %>% filter(Sample_id %in% o.samples) %>% count(cohort) 
# all other samples are from the validation cohort 

sv.samples = c(s.samples, v.samples)
h.samples = colnames(chip_data_all)[grep("PS", colnames(chip_data_all))]
sh.samples = c(s.samples, h.samples)
a.samples = colnames(chip_data_all)
c.samples = colnames(chip_data_all)[grep("CRC", colnames(chip_data_all))]
l.samples = c(l.samples, colnames(chip_data_all)[grep("^LC", colnames(chip_data_all))])
test.samples = c(s.samples, h.samples, l.samples, c.samples)
test.groups = c("SCLC", "Healthy", "NSCLC", "CRC")

col.annotation = data.frame(row.names = a.samples, 
                            sample = a.samples)
col.annotation$group = NA
col.annotation[s.samples,"group"] = "SCLC"
col.annotation[v.samples,"group"] = "validation"
col.annotation[h.samples,"group"] = "Healthy"
col.annotation[c.samples,"group"] = "CRC"
col.annotation[l.samples,"group"] = "NSCLC"
col.annotation[o.samples,"group"] = "other"
write.csv(col.annotation, quote = F, row.names = F,
          file = paste0(paperDir, "sample_annotations_2024.csv"))

# quality control validation metadata  ------------------------------------
metadata %>% 
  ggplot(aes(x = NE_SCORE, fill = cohort, color = cohort)) +
  # geom_histogram(color = "white") +
  geom_density(alpha = .7) +
  # stat_ecdf() 
  scale_color_aaas() +
  scale_fill_aaas() +
  theme(legend.position.inside = c(.8,.8), legend.key.size = unit(3, "mm")) -> p
ggsave(paste0(figDirPaper, "validation/NE_score_density.pdf"), p, width = 60, 
       height = 60, units = "mm")

metadata %>% 
  mutate(POU2F3 = as.numeric(POU2F3)) %>%
  ggplot(aes(x = POU2F3, fill = cohort, color = cohort)) +
  # geom_histogram(color = "white", position = "stacked") +
  # geom_density(alpha = .7) +
  # geom_density(alpha = .7) +
  stat_ecdf() + 
  scale_color_aaas() +
  scale_fill_aaas() +
  theme(legend.position = c(.8,.3), legend.key.size = unit(3, "mm")) -> p
ggsave(paste0(figDirPaper, "validation/POU2F3_ecdf.pdf"), p, width = 60, 
       height = 60, units = "mm")

metadata %>% 
  mutate(ASCL1 = as.numeric(ASCL1)) %>%
  ggplot(aes(x = ASCL1, fill = cohort, color = cohort)) +
  # geom_histogram(color = "white", position = "stacked") +
  # geom_density(alpha = .7) +
  # geom_density(alpha = .7) +
  stat_ecdf() + 
  scale_color_aaas() +
  scale_fill_aaas() +
  theme(legend.position = c(.8,.3), legend.key.size = unit(3, "mm")) -> p
ggsave(paste0(figDirPaper, "validation/ASCL1_ecdf.pdf"), p, width = 60, 
       height = 60, units = "mm")


metadata %>% 
  ggplot(aes(x = SCLC_score, fill = cohort, color = cohort)) +
  # geom_histogram(color = "white", position = "stacked") +
  # geom_density(alpha = .7) +
  # geom_density(alpha = .7) +
  stat_ecdf() + 
  scale_color_aaas() +
  scale_fill_aaas() +
  theme(legend.position = c(.8,.3), legend.key.size = unit(3, "mm")) -> p
ggsave(paste0(figDirPaper, "validation/SCORE_ecdf.pdf"), p, width = 60, 
       height = 60, units = "mm")

metadata %>% 
  # melt(id.vars = "NE_status", measure.vars = "cohort") %>%
  filter(NE_status != ".") %>%
  count(NE_status, cohort) %>%
  ggplot(aes(NE_status, n, fill = cohort)) +
  geom_bar(stat="identity", position=position_dodge()) +
  scale_fill_aaas() +
  labs(x = "", y = "# samples") +
  theme(legend.position = c(.8,.7), legend.key.size = unit(3, "mm")) -> p
ggsave(paste0(figDirPaper, "validation/NE_bar.pdf"), p, width = 60, 
       height = 60, units = "mm")

metadata %>% 
  # melt(id.vars = "NE_status", measure.vars = "cohort") %>%
  filter(NAPY_allcluster_expression != ".") %>%
  count(NAPY_allcluster_expression, cohort) %>%
  ggplot(aes(NAPY_allcluster_expression, n, fill = cohort)) +
  geom_bar(stat="identity", position=position_dodge()) +
  scale_fill_aaas() +
  labs(x = "", y = "# samples") +
  theme(legend.position = c(.8,.7), legend.key.size = unit(3, "mm")) -> p
ggsave(paste0(figDirPaper, "validation/NAPY_bar.pdf"), p, width = 80, 
       height = 60, units = "mm")
# RNA of validation set ---------------------------------------------------

v.samples

validation.file = paste0(baseDir, "data_from_NIH/nature_geentics_validation_cohort.xlsx")
validation.rna = data.frame(read_xlsx(validation.file)) 
validation.rna %>%
  mutate(plasma_sample = sub("_.*", "", plasma_sample), 
         INSM1_Hscore = as.numeric(INSM1_Hscore),
         biopsy2blood = round(difftime(biopsy.date, blood_draw_date, units = "d")), 
         matched_time = abs(biopsy2blood) < 10, 
         train_validation = if_else(plasma_sample %in% s.samples, 
                                    "test", "validation")) -> validation.rna
rownames(validation.rna) = validation.rna$sample_id

patient.sample.table = read.csv("~/BloodChIP/Analysis/Projects/NIH_SCLC/SCLC-NIH-All_samples_Israa.csv")

validation.rna %>%
  filter(!patient_id %in% patient.sample.table$Collaborators.Subject.ID) %>%
  filter(!sample_id %in% metadata$Primary_tumor_sample) %>% 
  filter(!sample_id %in% metadata$Tumor_sample_1) %>% 
  filter(!sample_id %in% metadata$Tumor_sample_2) %>%
  filter(!sample_id %in% metadata$Tumor_sample_3) %>%
  filter(!sample_id %in% metadata$Tumor_sample_4) %>%
  filter(!sample_id %in% metadata$Tumor_sample_5) %>%
  select(patient_id, sample_id) -> unkown.samples
  

validation.rna %>%
inner_join(patient.sample.table, 
             join_by("patient_id" == "Collaborators.Subject.ID")) %>%
  # select(Sample.ID, patient_id, sample_id) %>%
  filter(!Sample.ID %in% metadata$Sample_id) %>%
  # filter(!Sample.ID %in% s.samples) %>%
  select(Sample.ID, patient_id) %>%
  .$Sample.ID

  # left_join(patient.sample.table, c("patient_id", "Collaborators.Subject.ID"))
  # filter(!patient_id %in% metadata$Patient_id) %>%


write.table(setdiff(validation.rna$sample_id, colnames(x)), 
            "~/Downloads/validation_only.txt", quote = F, row.names = F, col.names = F)

validation.rna %>%
  filter(plasma_sample %in% s.samples) %>%
  mutate(INSM1.chip = chip_data_all["INSM1", plasma_sample],
        sclc_score = metadata[plasma_sample, "SCLC_score"]) %>%
  # filter(sclc_score > low.score.cutoff) %>%
  ggplot(aes(INSM1_Hscore, (2^INSM1)-1)) +
  geom_point(aes(alpha = sclc_score)) +
  stat_cor() 
  


validation.rna %>% 
  select(ASCL1, NEUROD1, POU2F3, YAP1, ATOH1) %>% 
  mutate(ATOH1 = log2(1+ATOH1)) %>%
  round(digits = 2) %>% 
  unique() -> validation.tf

metadata %>%
  filter(Sample_id %in% colnames(rna_data)) %>%
  mutate(ATOH1 = t(rna_data["ATOH1", Sample_id])) %>%
  select(ASCL1, NEUROD1, POU2F3, YAP1, ATOH1) %>% 
  mutate_if(is.character, as.numeric) %>% 
  round(digits = 2) %>%
  unique() -> test.tf

dim(validation.tf)
dim(test.tf)
all.tf = rbind(validation.tf, test.tf)
dim(all.tf)
dim(unique(all.tf))
all.tf = unique(all.tf)
pheatmap::pheatmap(all.tf, fontsize_row = 4)


# old stuff ---------------------------------------------------------------



# old
fixed_names = read.csv("~/BloodChIP/Analysis/Projects/NIH_SCLC/Sample Correction/fixing_runs_names.csv",
                       sep = "\t")
matched.ind = match(metadata$Sample_id, fixed_names$Old_name)
# sanity check
all(metadata$Sample_id[which(!is.na(matched.ind))] == 
  fixed_names$Old_name[matched.ind[which(!is.na(matched.ind))]])

metadata$Sample_id[which(!is.na(matched.ind))] = 
  fixed_names$New_name[matched.ind[which(!is.na(matched.ind))]]

passed.qc.samples[which(!sub("_.*", "", passed.qc.samples) %in% sub("_.*", "", metadata$Sample_id))]

sort(table(sub("_.*", "", metadata$Sample_id)))

duplicated.samples = names(which(table(sub("_.*", "", metadata$Sample_id)) > 1))
metadata %>%
  filter(grepl("Ali", Sample_id)) %>%
  dplyr::select(c(Sample_id))

  