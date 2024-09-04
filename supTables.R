# RNA samples high/low matched/non-matched  -------------------------------
qc %>% 
  mutate(Sample_ID = sub("_.*", "", rownames(qc))) %>% 
  filter(Sample_ID %in% s.samples) %>% # TODO - use also other samples?
  # pull(Total.uniq) %>% 
  # mean()
  `rownames<-`( NULL ) %>% 
  column_to_rownames("Sample_ID") %>% 
  rename("unique reads" = Total.uniq) %>% 
  select("unique reads") %>% # TODO - add other fields?
  write.csv(paste0(tableDir, "TableS2_sequencing_statistics.csv"))

write.table(diff.genes.common, paste0(tableDir, "TableS3_SCLC_differential_genes.txt"), 
            quote = F, row.names = F, col.names = F)

write.csv(data.cluster.genes, paste0(tableDir, 'TableS4-lung_cell-type_marker_genes.csv'), 
                                     row.names = F)

# liver.genes are defined in tumor_RNA_exploratory.R
write.table(liver.genes, paste0(tableDir, "TableS5-liver_specific_genes.txt"), 
            quote = F, row.names = F, col.names = F)

data.frame(TSS.windows[ascl1.ttest.wins]) %>% 
  select(c(seqnames, start, end)) %>% 
  write.csv(paste0(tableDir, "TableS6_ascl1_siganture.csv"), row.names = F)


# OLD ---------------------------------------------------------------------
# metadata %>% 
filter(`Tumor_ctDNA_timing_match (determined based on treatment timing standpoint)` == "Yes") %>% 
  filter(Sample_id %in% sv.samples) %>%
  .$Sample_id %in% rna.samples.matched.time


# total number of samples by groups (samples and patients)
sample.annotation %>% 
  left_join(metadata %>% select(Sample_id, cohort, Primary_tumor_sample), 
            join_by(sample == Sample_id)) %>% 
  filter(cohort == "train" | 
           !Primary_tumor_sample %in% c("No tumor RNA", "No tumor RNAseq") | 
           sample %in% c(c.samples, l.samples)) %>% 
  count(group)

sample.annotation %>% 
  left_join(metadata %>% select(Sample_id, cohort, Primary_tumor_sample), 
            join_by(sample == Sample_id)) %>% 
  filter(cohort == "train" | sample %in% rna.samples | 
           sample %in% c(c.samples, l.samples)) %>% 
  mutate(patient = sub("-.*", "", sample)) %>% 
  group_by(patient) %>% 
  slice_sample(n = 1) %>% 
  group_by(group) %>% 
  count(group)

# number of tumors by group
sample.annotation %>% 
  filter(sample %in% colnames(rna_data)) %>% 
  group_by(group) %>% 
  count()
