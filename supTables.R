# supplementary tables  ---------------------------------------------------
# table 1 
metadata %>% 
  filter(SCLC.state != "other") %>% 
  select(-c(Collaborators.Sample.ID)) %>% 
  rename(SCLC_score = SCLC.n) %>% 
  write.xlsx(file = paste0(tableDir, "Table-S1_sample_metadata.xlsx"),
             sheetName = "plasma samples", append = FALSE)

metadata.patients %>% 
  filter(SCLC.state != "other") %>% 
  select(-c(NIH_PatientID)) %>% 
  write.xlsx(file = paste0(tableDir, "Table-S1_sample_metadata.xlsx"), 
             sheetName="patients", append=TRUE)
  
metadata.biopsies %>% 
  filter(!PatientID %in% sub("-.*", "", o.samples)) %>% 
  select(-c(RNA_id, `Surgical Case`)) %>% 
  write.xlsx(file = paste0(tableDir, "Table-S1_sample_metadata.xlsx"), 
           sheetName="biopsy samples", append=TRUE)

metadata.matching %>%
  filter(SCLC.state != "other") %>% 
  select(c(Sample_id, Matched)) %>% 
  write.xlsx(file = paste0(tableDir, "Table-S1_sample_metadata.xlsx"), 
             sheetName="plasma-biopsy time match", append=TRUE)

metadata.timeline %>% 
  filter(!PatientID %in% sub("-.*", "", o.samples)) %>% 
  write.xlsx(file = paste0(tableDir, "Table-S1_sample_metadata.xlsx"), 
             sheetName="patient timelines", append=TRUE)

# table 2
sample.annotation %>% 
  filter(group != "other") %>% 
  select(-patient) %>% 
  arrange(Sample_id) %>% 
  write.csv(paste0(tableDir, "Table-S2_sample_annotations.csv"), 
          quote = F, row.names = F)

# table 3
# data_qc is created at SequencingQC.R
data_qc %>% 
  rename(Sample_id = sample_name, 
         "unique reads" = total_unique, 
         "% signal" = tss_pct_signal) %>% 
  mutate(Sample_id = sub("_.*", "", Sample_id)) %>% 
  filter(Sample_id %in% (sample.annotation %>% 
                           filter(group != "other") %>% 
                           pull(Sample_id))) %>% 
  select(Sample_id, "unique reads", "% signal") %>% 
  write.csv(paste0(tableDir, "Table-S3_sequencing_statistics.csv"), row.names = F)

# table 4
# defined in differential_genes-fig1.R
write.table(diff.genes.common, paste0(tableDir, "Table-S4_SCLC_differential_genes.txt"), 
            quote = F, row.names = F, col.names = F)

# table 5
# defined in lung_celltypes-fig2.R
write.csv(data.cluster.genes, paste0(tableDir, 'Table-S5_lung_cell-type_marker_genes.csv'), 
          row.names = F)

# table 6
# liver.genes are defined in tumor_RNA_exploratory.R
write.table(liver.genes, paste0(tableDir, "Table-S6_liver_specific_genes.txt"), 
            quote = F, row.names = F, col.names = F)

# table 7
data.frame(TSS.windows[ascl1.ttest.wins]) %>% 
  select(c(seqnames, start, end, name)) %>% 
  mutate(name = case_when(
    is.na(name) ~ "",
    name == "." ~ "",
    .default  = name
  )) %>% 
  write.csv(paste0(tableDir, "Table-S7_ascl1_siganture.csv"), row.names = F)


# group numbers for text  -------------------------------------------------
# total number of samples (excluding other)
sample.annotation %>% 
  group_by(group) %>%
  filter(group != "other") %>%
  summarise(
    n_samples = n_distinct(Sample_id),
    n_patients = n_distinct(patient)
  )

metadata.matching %>% 
  filter(SCLC.state != "other") %>%
  count(Matched)

# number of samples by timepoint (figure 1A)
metadata %>% 
  filter(SCLC.state == "SCLC") %>% 
  count(Timepoint) 

# number of train patients and samples by groups 
sample.annotation %>% 
  left_join(metadata %>% 
              select(Sample_id, cohort), 
            join_by("sample" == "Sample_id")) %>% 
  # filter(cohort != "validation" | is.na(cohort)) %>%
  filter(sample %in% ssl.samples) %>%
  mutate(patient = sub("-.*", "", sample)) %>% 
  group_by(group, cohort) %>% 
  summarise(n_samples = length(patient))

# number of validation patients by groups 
sample.annotation %>% 
  left_join(metadata %>% 
              select(Sample_id, cohort), 
            join_by("sample" == "Sample_id")) %>% 
  filter(cohort == "validation", sample %in% ssl.samples) %>%
  mutate(patient = sub("-.*", "", sample)) %>% 
  group_by(group) %>% 
  summarise(n_patients = length(unique(patient)), 
            n_samples = length(patient))

metadata %>% 
  inner_join(metadata.matching %>% 
               mutate(PatientID = sub("-.*", "", Sample_id)) %>% 
               select(PatientID, RNA_id), 
             by = "PatientID") %>% 
  filter(SCLC.state %in% c("SCLC", "SCLC-like")) %>% #select(Sample_id, cohort, RNA_id)
  count(cohort)


  


# distribution of TF in cfChIP samples
rna_data[c("ASCL1", "NEUROD1", "POU2F3", "ATOH1"),] %>% 
  t() %>% 
  sweep(2, unlist(chip.tf.cutoff), ">") %>% 
  as.data.frame() %>% 
  summarise(a.n = sum(ASCL1)/length(ASCL1), 
            n.n = sum(NEUROD1)/length(ASCL1), 
            p.n = sum(POU2F3)/length(ASCL1), 
            at.n = sum(ATOH1)/length(ASCL1))


mutate(rna = SCLC_RNA_fixed[g, SampleShort_ID], 
       chip = SCLC_ChIP_fixed[g, SampleShort_ID], 
       positive = factor(rna > cu, levels = c(T,F))) -> df.gene

# data for upload to Zenodo -----------------------------------------------
zenodo.samples = c(s.samples, nec.samples, l.samples, h.samples, c.samples)
write.csv(chip_data_all[,zenodo.samples], 
          gzfile(paste0(paperDir, "Zenodo/GeneCounts.gz")), quote = F, 
          row.names = T)

zenodo.rna.samples = colnames(rna_data)[!colnames(rna_data) %in% o.samples]
  
write.csv(rna_data[,zenodo.rna.samples], 
          gzfile(paste0(paperDir, "Zenodo/RNAcpm.gz")), quote = F, 
          row.names = T)


# supporting data tables (data of all figure panels) ---------------------------
library("openxlsx")

sup_data_values = list(
  fig1b = fig1b,
  fig1c = fig1c,
  fig1dS1c = fig1dS1c,
  fig1eS1de = fig1eS1de,
  fig1g = fig1g,
  figS1a = figS1a,
  figS1b = figS1b,
  figS1f = figS1f,
  figS1g = figS1g, 
  fig2a = fig2a, 
  fig2b = fig2b,
  fig2d = fig2d,
  fig2eS3ab = fig2eS3ab,
  fig2fS3c1 = fig2fS3c1,
  fig2g = fig2g,
  figS2 = figS2,
  fig3a = fig3a, 
  fig3b = fig3b, 
  fig3c = fig3c,
  figS3c2 = figS3c2,
  figS3d = figS3d,
  fig4ab = fig4ab,
  figS4b = figS4b,
  figS5ab = figS5ab)
write.xlsx(sup_data_values, file = paste0(tableDir, "Supporting data values.xlsx"))



