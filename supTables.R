# supplementary tables  ---------------------------------------------------
# table 1 
metadata %>% 
  filter(SCLC.state != "other") %>% 
  select(-c(Collaborators.Sample.ID)) %>% 
  write.xlsx(file = paste0(tableDir, "Table-S1_sample_metadata.xlsx"),
             sheetName = "plasma samples", append = FALSE)

metadata.patients %>% 
  filter(SCLC.state != "other") %>% 
  select(-c(NIH_PatientID)) %>% 
  write.xlsx(file = paste0(tableDir, "Table-S1_sample_metadata.xlsx"), 
             sheetName="patients", append=TRUE, row.names = F)
  
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

# table 2
sample.annotation %>% 
  filter(group != "other") %>% 
  write.csv(paste0(tableDir, "Table-S2_sample_annotations.csv"), 
          quote = F, row.names = F)

# table 3
qc %>% 
  mutate(Sample_ID = sub("_.*", "", rownames(qc))) %>% 
  filter(Sample_ID %in% s.samples) %>% # TODO - use also other samples?
  # pull(Total.uniq) %>%
  # mean() / 1e6
  `rownames<-`( NULL ) %>% 
  column_to_rownames("Sample_ID") %>% 
  rename("unique reads" = Total.uniq) %>% 
  select("unique reads") %>% # TODO - add other fields?
  write.csv(paste0(tableDir, "Table-S3_sequencing_statistics.csv"))

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
  left_join(metadata %>% 
              select(Sample_id, cohort), 
            join_by("sample" == "Sample_id")) %>% 
  filter(group != "other") %>%
  nrow()

# total number of training set
sample.annotation %>% 
  left_join(metadata %>% 
              select(Sample_id, cohort), 
            join_by("sample" == "Sample_id")) %>% 
  filter(cohort != "validation" | is.na(cohort)) %>% 
  nrow()

# number of train patients and samples by groups 
sample.annotation %>% 
  left_join(metadata %>% 
              select(Sample_id, cohort), 
            join_by("sample" == "Sample_id")) %>% 
  filter(cohort != "validation" | is.na(cohort)) %>%
  mutate(patient = sub("-.*", "", sample)) %>% 
  group_by(group) %>% 
  summarise(n_patients = length(unique(patient)), 
            n_samples = length(patient))

# number of test patients by groups 
sample.annotation %>% 
  left_join(metadata %>% 
              select(Sample_id, cohort), 
            join_by("sample" == "Sample_id")) %>% 
  filter(cohort == "validation") %>%
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

# number of samples by timepoint (figure 1A)
metadata %>% 
  filter(cohort != "validation" | is.na(cohort), SCLC.state != "other") %>% 
  # dim()
  count(Timepoint) 
  


# data for upload to Zenodo -----------------------------------------------
zenodo.samples = c(s.samples, nec.samples, l.samples, h.samples, c.samples)
write.csv(chip_data_all[,zenodo.samples], 
          gzfile(paste0(paperDir, "Zenodo/GeneCounts.gz")), quote = F, 
          row.names = T)

zenodo.rna.samples = colnames(rna_data)[!colnames(rna_data) %in% o.samples]
  
write.csv(rna_data[,zenodo.rna.samples], 
          gzfile(paste0(paperDir, "Zenodo/RNAcpm.gz")), quote = F, 
          row.names = T)

