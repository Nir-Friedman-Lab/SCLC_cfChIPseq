# table1 ------------------------------------------------------------------
library(DBI)
dbcon <- dbConnect(RPostgres::Postgres(),
                   user = "gavriel.fialkoff@senseerahealth.com",
                   password = system(paste0("gcloud sql generate-login-token"), 
                                     intern = TRUE),
                   host = "34.165.66.84",
                   dbname = "senseera-prod",
                   sslmode='require')
dbGetQuery(dbcon, "SELECT *
                    FROM sample_qc_view;") %>% 
  mutate(sample_name_short = sub("_.*", "", sample_name)) -> all.qc
dbGetQuery(dbcon, "SELECT analysis_id, sex, confidence
                    FROM analysis_sex_pred") %>% 
  inner_join(all.qc %>% 
               select(analysis_id, sample_name, sample_name_short)) -> all_sex

read.csv("~/gavriel.fialkoff@mail.huji.ac.il - Google Drive/My Drive/CollatedFolders/Projects/NIH_SCLC/cfChIP-paper/Tables/Table1_information.csv") %>% 
  left_join(all_sex %>% 
               group_by(sample_name_short) %>% 
               slice_max(confidence, with_ties = F) %>% 
               ungroup() %>% 
              mutate(batch = sub(".*H3K4me3-([^_]+)_.*", "\\1", sample_name)) %>% 
               select(sample_name_short, sex, batch), 
            by = "sample_name_short") %>% 
  mutate(PatientID = if_else(PatientID != "",
                             PatientID,
                             sub("-.*", "", sample_name_short)),
         sex = if_else(is.na(sex.y), sex.x, sex.y)) %>% 
  select(-sex.x, -sex.y) -> 
  table1_info

source("~/gavriel.fialkoff@senseerahealth.com - Google Drive/My Drive/Code/Projects/GlobalAnalysis/agePrediction.R")
runAgePrediction = F
if (runAgePrediction) {
  agePredication(all.qc %>% 
                   filter(grepl("_H3K4me3", sample_name), total_unique > 5e5) %>% 
                   filter(sample_name_short %in% table1_info$sample_name_short) %>% 
                   select(sample_name, sample_name_short, tss_win_count) %>% 
                   group_by(sample_name_short) %>% 
                   slice_max(tss_win_count) %>% 
                   pull(sample_name)
  ) -> 
    sample_age
}


sample_age %>% 
  mutate(predicted_age = case_when(
    predicted_age < 20 ~ 20, 
    predicted_age > 90 ~ 90, 
    .default = predicted_age
  )) %>% 
  mutate(PatientID = sub("-.*", "", sample_name)) %>% 
  inner_join(metadata.patients %>% 
               select(PatientID, `age at diagnosis`, `Date of Birth`)) %>% 
  mutate(age = 2015 - year(`Date of Birth`)) %>% 
  # filter(n_overexp < 100) %>% 
  ggplot(aes(predicted_age, age)) +
  geom_point() +
  geom_abline() +
  lims(x = c(0,NA), y = c(0, NA))

table1_info %>%
  group_by(group) %>%
  summarise(
    n_patients = n_distinct(PatientID),
    n_samples = n()
  )
  
table1_info %>%
  left_join(metadata.patients %>% 
               select(PatientID, `age at diagnosis`, `Date of Birth`, `Race(s)`)) %>% 
  mutate(age = 2015 - year(`Date of Birth`)) %>% 
  left_join(sample_age) %>% 
  # ggplot(aes(age, predicted_age)) +
  # geom_point() +
  # geom_abline() +
  # lims(x = c(0,NA), y = c(0, NA))
  mutate(Age = if_else(!is.na(age), age, predicted_age), 
         Race = factor(if_else(!is.na(`Race(s)`), `Race(s)`, Race))) %>% 
  select(sample_name_short, group, Race, sex, Age, batch, PatientID) ->
table1_info_age_sex_race

library(tableone)
library(gtsummary)
library(gt)
library(flextable)
# 1. Prepare your unique-patient dataset
table1_info_age_sex_race %>%
  group_by(PatientID) %>%
  slice_head(n = 1) %>%
  ungroup() -> 
  metadata_unique

myVars <- c("sex", "Age", "Race", "batch") 
catVars <- c("sex", "Race") 

tab1 <- metadata_unique %>%
  select(group, all_of(myVars)) %>%
  tbl_summary(
    by = group, 
    missing = "no",
    type = list(sex ~ "dichotomous", 
                Race ~ "categorical"),
    value = list(sex ~ "female"), # Displays only the Female row
    label = list(
      Age ~ "Age (years)",
      sex ~ "Sex (Female)"
      ),
    statistic = list(
      all_continuous() ~ "{median} ({p25}, {p75})",
      # Changed {p_formatted} to {p}
      all_categorical() ~ "{n} / {N} ({p}%)" 
    )
  ) %>%
  add_p() %>%
  bold_labels()
tab1 

tab1 %>%
  as_flex_table() %>%
  save_as_docx(path = paste0(figDirPaper, "JCI_resubmission/Table1.docx"))

# batch effect ------------------------------------------------------------
library(ggcorrplot)
conflicts_prefer(base::as.factor)

# pca.samples = c(h.samples, s.samples)
pca.samples = colnames(chip_data_all)
res.pca = prcomp(t(chip_data_all[,pca.samples]))
fviz_eig(res.pca, geom = "bar", barfill = "black", barcolor = "white", 
         ggtheme = base_theme) + 
  theme(plot.title = element_blank()) + 
  labs(x = "Principal components", y = "% explained variance")
data.frame(pc1 = res.pca$x[,1], 
           pc2 = res.pca$x[,2], 
           pc3 = res.pca$x[,3],
           pc4 = res.pca$x[,4],
           pc5 = res.pca$x[,5],
           PatientID = sub("-.*", "", pca.samples), 
           sample_name_short = pca.samples) %>% 
  left_join(estimated.tumor %>% 
              rownames_to_column("sample_name_short") %>% 
              select(sample_name_short, SCLC.n)) %>% 
  inner_join(table1_info_age_sex_race, by = "sample_name_short") -> data.pca



data.pca %>%
  mutate(
    Sex = as.numeric(as.factor(sex)),
    Age = Age,
    Race = as.numeric(as.factor(Race)),
    TumorFraction = SCLC.n
  ) %>%
  select(pc1, pc2, pc3, pc4, pc5, Age, Sex, Race, TumorFraction, sample_name_short) -> 
  data.pca_numeric

# 2. Calculate the correlation matrix
cor_mat <- cor(data.pca_numeric %>% 
                 select(-sample_name_short), use = "complete.obs")

# 3. Filter the matrix to show only PCs vs Metadata
# (Rows = PCs, Columns = Clinical Variables)
plot_mat <- cor_mat[1:5, 6:(ncol(data.pca_numeric)-1)]

plot_mat %>% 
  as.data.frame() %>%
  rownames_to_column("Var1") %>%
  pivot_longer(-Var1, names_to = "Var2", values_to = "correlation") %>% 
  ggplot(aes(Var1, Var2, fill = correlation)) +
  geom_tile(color = "black") +
  geom_text(aes(label = round(correlation, 2)), size = 4/.pt) + 
  scale_fill_gradient2(low = "#6D9EC1", high = "#E46726", mid = "white", 
                       midpoint = 0, limit = c(-1,1)) +
  labs(x = "", y = "") + 
  coord_fixed() +
  theme(legend.key.size = unit(2,"mm"),
        # legend.position = "bottom",
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        legend.text = element_text(size = 4),
        legend.title = element_text(size = 4),
        # axis.text.x = element_text(size = base_size/.pt, angle = 45, vjust = 1, hjust = 1),
        # axis.text.y = element_text(size = base_size/.pt)
  ) ->
  p
ggsaveG(p, paste0(figDirPaper, "JCI_resubmission/cor_pc_vars"), 50, 50)


# List of clinical columns you want to test
clinical_vars <- c("Age", "Sex", "Race", "TumorFraction")

results <- lapply(clinical_vars, function(var) {
  # Convert column to numeric (handles factors/characters automatically)
  y <- as.numeric(as.factor(data.pca_numeric[[var]]))
  x <- data.pca_numeric$pc3
  
  test <- cor.test(x, y, use = "complete.obs")
  
  data.frame(
    Variable = var,
    Correlation = test$estimate,
    P_Value = test$p.value
  )
})

# Combine and print
bind_rows(results)

# Linear model to see which variable actually "wins"
summary(lm(pc1 ~ TumorFraction + Age + Race + Sex, data = data.pca_numeric))
summary(lm(pc2 ~ TumorFraction + Age + Race + Sex, data = data.pca_numeric))
summary(lm(pc3 ~ TumorFraction + Age + Race + Sex, data = data.pca_numeric))
summary(lm(pc4 ~ TumorFraction + Age + Race, data = data.pca_numeric))

# replicates --------------------------------------------------------------
library(patchwork)
qc.signal.yeild %>% 
  filter(grepl("_H3K4me3", sample_name)) %>% 
  filter(grepl("SCLC", sample_name)) %>%
  filter(grepl("Beads", sample_name)) %>% 
  filter(grepl("-11", sample_name)) %>% 
  select(sample_name_short, sample_name) -> replicates

qc.signal.yeild %>% 
  filter(grepl("_H3K4me3", sample_name)) %>% 
  filter(grepl("SCLC", sample_name)) %>%
  filter(tss_pct_signal > 30, total_unique > 5e5) %>% 
  filter(grepl("Aliquot", sample_name)) %>%
  group_by(sample_name_short) %>% 
  pull(sample_name_short) -> tech_replicates

qc.signal.yeild %>% 
  filter(grepl("_H3K4me3", sample_name)) %>% 
  filter(sample_name_short %in% tech_replicates) %>% 
  filter(!grepl("M-5ulBeads", sample_name)) %>% 
  filter(tss_pct_signal > 30, total_unique > 5e5) %>% 
  group_by(sample_name_short) %>% 
  filter(tss_pct_signal > 55) %>% 
  # filter((max(tss_pct_signal) - min(tss_pct_signal)) < 10) %>% 
  filter(n() > 1) %>% 
  select(sample_name, sample_name_short) -> tech_replicates

SourceDIR = "~/Documents/cfChIP_Pipeline/cfChIP-development/"
source("~/Documents/cfChIP_Pipeline/cfChIP-development/Core/cfChIP-Cmd-Common.R")

tiss_sigs = c("Erythroblast", "Neutrophils", "Megakaryocyte", "Liver",
              "Lung", "Brain", "GI", "SCLC", "Loycocyte"
              # "Monocytes","Macrophage", "Colon",
              )
plot_list <- list() 
plot_tiss_list <- list() 
data_gene_list <- list()
data_tiss_list <- list()
for (samp in unique(tech_replicates$sample_name_short)) {
  idx = which(unique(tech_replicates$sample_name_short) == samp)
  message("Processing: ", samp, " (", idx, " out of ", length(unique(tech_replicates$sample_name_short)), ")")
  # Identify the two replicate names
  indices = which(qc.signal.yeild$sample_name_short == samp)
  samp1 = qc.signal.yeild$sample_name[indices[1]]
  samp2 = qc.signal.yeild$sample_name[indices[2]]
  
  # Load data
  dat1 = cfChIP.open.file(paste0("gs://senseera-main/data/rds/H3K4me3/", samp1, ".rdata"))
  dat2 = cfChIP.open.file(paste0("gs://senseera-main/data/rds/H3K4me3/", samp2, ".rdata"))
  
  # --- Gene Plot Section ---
  g_ind = names(dat1$GeneCounts) %in% Genes.notexcluded
  
  # --- CHANGE: Save the dataframe to a variable first ---
  df_gene <- data.frame(samp1 = dat1$GeneCounts.QQnorm[g_ind], 
                        samp2 = dat2$GeneCounts.QQnorm[g_ind])
  data_gene_list[[samp]] <- df_gene 
  # ------------------------------------------------------
  
  p_gene <- df_gene %>% 
    ggplot(aes(samp1, samp2)) + 
    geom_point(size = .1, alpha = 0.2) + 
    scale_x_continuous(transform = "log1p", breaks = c(0,1,3,10,30,100,300)) + 
    scale_y_continuous(transform = "log1p", breaks = c(0,1,3,10,30,100,300)) + 
    geom_density2d(lwd = .2, color = "blue") + 
    geom_abline(lwd = .3, linetype = "dashed", color = "red") + 
    labs(x = "Repeat 1", y = "Repeat 2", title = samp) + 
    stat_cor(label.x = 0.5) + 
    theme(base_size = base_size, aspect.ratio = 1) 
  
  plot_list[[samp]] <- p_gene
  
  # --- Tissue Plot Section ---
  # --- CHANGE: Optimized data extraction for tissue signatures ---
  tiss_dat1 <- dat1$Signatures$vs.Background.Corrected %>% filter(Signature %in% tiss_sigs)
  tiss_dat2 <- dat2$Signatures$vs.Background.Corrected %>% filter(Signature %in% tiss_sigs)
  
  df_tiss <- data.frame(
    tiss  = tiss_dat1$Signature,
    samp1 = as.numeric(tiss_dat1$NormalizedCounts),
    samp2 = as.numeric(tiss_dat2$NormalizedCounts)
  )
  
  data_tiss_list[[samp]] <- df_tiss # Save values
  # ---------------------------------------------------------------
  
  p_tiss <- df_tiss %>% 
    ggplot(aes(samp1, samp2)) + 
    geom_point(size = 2) + 
    scale_x_continuous(transform = "log1p", breaks = c(0,1,3,10,30,100,300)) +
    scale_y_continuous(transform = "log1p", breaks = c(0,1,3,10,30,100,300)) +
    geom_abline(lwd = .3, linetype = "dashed", color = "red") + 
    labs(x = "Repeat 1", y = "Repeat 2", title = samp) + 
    stat_cor(label.x = 0.5) + 
    theme(base_size = base_size, aspect.ratio = 1) 
  
  plot_tiss_list[[samp]] <- p_tiss
}

# Use patchwork to wrap all 6 plots into one
final_plot <- wrap_plots(plot_list, ncol = 3)
final_plot_tiss <- wrap_plots(plot_tiss_list, ncol = 3)
print(final_plot_tiss)
ggsaveG(final_plot, paste0(figDirPaper, "JCI_resubmission/tech_replicates_v1"), 130, 130)
ggsaveG(final_tiss_plot, paste0(figDirPaper, "JCI_resubmission/tech_replicates_tiss_v1"), 220, 220)
bind_rows(data_tiss_list, .id = "sample_id") -> 
  figS2


qc.signal.yeild %>%
  filter(grepl("_H3K4me3", sample_name)) %>% 
  filter(sample_name_short %in% tech_replicates$sample_name_short) %>% 
  pull(sample_name) -> tech_samples
LL_rep = lapply(tech_samples, function(s) {
  print(grep(s, tech_samples))
  cfChIP.open.file(paste0("gs://senseera-main/data/rds/H3K4me3/", s, ".rdata"))})
names(LL_rep) = sapply(LL_rep, function(s) s$Name)
rep_sigs = sapply(LL_rep, function(s) {
  s$Signatures$vs.Background.Corrected %>% 
    filter(Signature %in% tiss_sigs) %>% 
    pull(NormalizedCounts) %>% 
    as.numeric()})
colnames(rep_sigs) = names(LL_rep)
rownames(rep_sigs) = dat$Signatures$vs.Background.Corrected %>% 
  filter(Signature %in% tiss_sigs) %>% 
  pull(Signature)

rep_sigs %>% 
  as.data.frame() %>% 
  rownames_to_column("tissue") %>% 
  pivot_longer(
    cols = -tissue, 
    names_to = "sample_name", 
    values_to = "value"
  ) %>% 
  mutate(sample_name_short = sub("_.*", "", sample_name)) %>% 
  mutate(rep_id = ifelse(grepl("Aliquot2", sample_name), 2, 1)) %>%
  select(tissue, sample_name_short, rep_id, value) %>%
  pivot_wider(names_from = rep_id, values_from = value, values_fn = mean) %>% 
  ggplot(aes(`1`, `2`)) +
  geom_point() +
  geom_abline() + 
  facet_wrap(~tissue, scale = "free")
  
# depth -------------------------------------------------------------------
metadata %>% 
  select(Sample_id, SCLC.n) %>% 
  inner_join(data_qc, join_by("Sample_id" == "sample_name_short")) %>% 
  ggplot(aes(SCLC.n, tss_pct_signal)) + 
  geom_point() + 
  scale_x_continuous(breaks = sclc.breaks, labels = sclc.lab) +
  labs(x = "SCLC score", y = "% Signal reads") + 
  stat_cor() -> p
ggsaveG(p, paste0(figDirPaper, "JCI_resubmission/signal_vs_sclc_score"), 50, 50)

data_qc %>% 
  select(sample_name, total_unique, est_library_size) %>% 
  dplyr::rename(Sample_id = sample_name, 
                unique_reads = total_unique, 
                estimated_library = est_library_size) -> 
  figS1a

data_qc %>% 
  ggplot(aes(total_unique, est_library_size)) +
  geom_point() +
  scale_x_continuous(transform = "log10") +
  scale_y_continuous(transform = "log10") +
  labs(x = "Unique Reads", y = "Estimated Library") +
  geom_abline() -> 
  p
ggsaveG(p, paste0(figDirPaper, "JCI_resubmission/unique_vs_estimated"), 50, 50)



# SCLC score reproducibility ----------------------------------------------
# reg.atlas are from the liner-regression.R script
oldDataDir = "gs://senseera-archive/historical_backups/historical_09072025_1505/senseera-main/data/rds/H3K4me3/"
rep_samples = tech_replicates$sample_name
rep_score = matrix(data = 0, nrow = length(rep_samples), ncol = ncol(reg.atlas))
rownames(rep_score) = rep_samples; colnames(rep_score) = colnames(reg.atlas)

for (s in rep_samples) {
  print(grep(s, rep_samples))
  dat = cfChIP.open.file(paste0(oldDataDir, s, ".rdata"))
  names(dat$GeneCounts.QQnorm) = names(dat$GeneCounts)
  b = dat$GeneCounts.QQnorm[rownames(reg.atlas)]
  dim(reg.atlas)
  res = round(nnls(A = reg.atlas, b = b)$x, digits = 3)
  names(res) = colnames(reg.atlas)[grep(s, colnames(reg.atlas), invert = T)]
  rep_score[s,names(res)] = res
}

data.frame(healthy = rowSums(rep_score[,h.samples]), 
           SCLC = rep_score[,"sclc.arctype"]) %>% 
  mutate(SCLC.n = SCLC/(SCLC + healthy), 
         healthy.n = healthy/(SCLC + healthy)) -> estimated_tumor_rep

estimated_tumor_rep %>% 
  rownames_to_column("sample_name") %>% 
  mutate(sample_name_short = sub("_.*", "", sample_name)) %>% 
  mutate(rep_id = ifelse(grepl("Aliquot2", sample_name), 2, 1)) %>%
  select(sample_name_short, rep_id, SCLC.n) %>%
  pivot_wider(names_from = rep_id, values_from = SCLC.n) %>%
  ggplot(aes(100*`1`, 100*`2`)) +
  geom_point(size = 2) +
  lims(x = c(0,100), y = c(0,100)) + 
  labs(x = "Repeat 1", y = "Repeat 2", title = "SCLC score", 
       subtitle = "Biological Repeats") + 
  stat_cor() + 
  geom_abline() -> p
ggsaveG(p, paste0(figDirPaper, "JCI_resubmission/replicates_sclc"), 50, 50)

# subtype reproducibility -------------------------------------------------
# Map over the list and bind the results into a single table
all_reps = unique(c(tech_samples, replicates$sample_name))
sort(all_reps)
LL_rep_old = lapply(all_reps, function(s) {
  print(grep(s, all_reps))
  cfChIP.open.file(paste0(oldDataDir, s, ".rdata"))})
names(LL_rep_old) = sapply(LL_rep_old, function(s) s$Name)

tf_reps <- bind_rows(lapply(names(LL_rep_old), function(sn) {
  s <- LL_rep_old[[sn]]
  data.frame(
    sample_name = sn,
    ascl1   = sum(s$Counts.QQnorm[ascl1.wins], na.rm = TRUE), 
    neurod1 = sum(s$Counts.QQnorm[neurod1.wins], na.rm = TRUE), 
    pou2f3  = sum(s$Counts.QQnorm[pou2f3.wins], na.rm = TRUE), 
    atoh1   = sum(s$Counts.QQnorm[atoh1.wins], na.rm = TRUE)
  )
}))

tf_reps %>% 
  inner_join(qc.signal.yeild) %>% 
  # filter(tss_pct_signal > 40) %>% 
  mutate(sample_name_short = sub("_.*", "", sample_name)) %>% 
  group_by(sample_name_short) %>% 
  filter(max(tss_pct_signal) - min(tss_pct_signal) < 20) %>%
  select(-tss_pct_signal, -total_unique, -est_library_size) %>%
  # 1. Pivot to long first
  pivot_longer(cols = -c(sample_name, sample_name_short), 
               names_to = "gene", 
               values_to = "value") %>% 
  # 2. Group by the sample and gene to count the occurrences
  group_by(sample_name_short, gene) %>% 
  mutate(rep_id = row_number()) %>% 
  ungroup() %>% 
  # 3. Remove the unique full name
  select(-sample_name) %>% 
  # 4. Now pivot wider will be perfectly unique
  pivot_wider(
    names_from = rep_id, 
    names_prefix = "rep_", 
    values_from = value
  ) %>% 
  ggplot(aes(rep_1, rep_2, text = sample_name_short)) +
  geom_point(size = 2) +
  stat_cor(size = 5/.pt) + 
  labs(x = "Repeat 1", y = "Repeat 2") + 
  geom_abline()  +
  # scale_x_continuous(transform = "log1p", breaks = c(0,3,10,30)) +
  # scale_y_continuous(transform = "log1p", breaks = c(0,3,10,30)) + 
  facet_wrap(~gene) -> p
ggplotly(p)
ggsaveG(p, paste0(figDirPaper, "JCI_resubmission/replicates_subtype"), 70, 70)

# high pou2f3 example
a = "SCLC0030-435_Aliquot2_H3K4me3-728_06072023-96"
b = "SCLC0030-435_NA_H3K4me3-39_26032021-13"
outname = "POU2F3"

a = "SCLC0034-462_Aliquot2_H3K4me3-728_06072023-96"
b = "SCLC0034-462_NA_H3K4me3-39_26032021-13"
outname = "ASCL1"

a = "SCLC0051-612_Aliquot2_H3K4me3-728_06072023-96"
b = "SCLC0051-612_NA_H3K4me3-39_26032021-13"
outname = "NEUROD1"

a = "SCLC0057-676_NA_H3K4me3-38_26032021-13"
b = "SCLC0057-676_Aliquot2_H3K4me3-728_06072023-96"
outname = "ASCL1_2"


data_frame(gene = names(LL_rep_old[[a]]$GeneCounts), 
           rep1 = LL_rep_old[[a]]$GeneCounts.QQnorm, 
           rep2 = LL_rep_old[[b]]$GeneCounts.QQnorm) %>% 
  mutate(rep1 = case_when(
    gene == "ASCL1" ~ sum(LL_rep_old[[a]]$Counts.QQnorm[ascl1.wins], na.rm = TRUE),
    gene == "NEUROD1" ~ sum(LL_rep_old[[a]]$Counts.QQnorm[neurod1.wins], na.rm = TRUE),
    gene == "POU2F3" ~ sum(LL_rep_old[[a]]$Counts.QQnorm[pou2f3.wins], na.rm = TRUE),
    gene == "ATOH1" ~ sum(LL_rep_old[[a]]$Counts.QQnorm[atoh1.wins], na.rm = TRUE), 
    .default = rep1
  ), 
  rep2 = case_when(
    gene == "ASCL1" ~ sum(LL_rep_old[[b]]$Counts.QQnorm[ascl1.wins], na.rm = TRUE),
    gene == "NEUROD1" ~ sum(LL_rep_old[[b]]$Counts.QQnorm[neurod1.wins], na.rm = TRUE),
    gene == "POU2F3" ~ sum(LL_rep_old[[b]]$Counts.QQnorm[pou2f3.wins], na.rm = TRUE),
    gene == "ATOH1" ~ sum(LL_rep_old[[b]]$Counts.QQnorm[atoh1.wins], na.rm = TRUE), 
    .default = rep2)
) -> data_reptf

data_reptf %>% 
  ggplot(aes(rep1,rep2)) +
  geom_point(size = .3, color = "gray") +
  geom_abline() +
  geom_density_2d() +
  stat_cor(size = base_size /.pt) +
  scale_x_continuous(transform = "log1p", breaks = c(0,3,10,30,100,300)) +
  scale_y_continuous(transform = "log1p", breaks = c(0,3,10,30,100,300)) +
  geom_point(data = data_reptf %>% 
               filter(grepl("POU2F3|ASCL1|NEUROD1|ATOH1", gene)), 
             color = "darkred", size = 1) +
  labs(x = "Repeat 1", y = "Repeat 2") +
  geom_label_repel(data = data_reptf %>% 
               filter(grepl("POU2F3|ASCL1|NEUROD1|ATOH1", gene)), 
               label.padding =  unit(1, "pt"),
             aes(label = gene), color = "darkred", size = base_size/.pt, 
             segment.size = .2) -> p
ggsaveG(p, paste0(figDirPaper, "JCI_resubmission/replicates_subtype_", outname), 50, 50)  
  


# RNA-ChIP p-value --------------------------------------------------------
r_matched <- unlist(cor.all["estimate",])
r_random  <- unlist(cor.all.rand["estimate",])
# Option A: Compare the two distributions (very common for manuscripts)
wilcox_result <- wilcox.test(r_matched, r_random, alternative = "greater")
final_p_value <- wilcox_result$p.value

# Option B: Empirical p-value (Fraction of random > matched mean)
actual_mean <- mean(r_matched, na.rm = TRUE)
empirical_p <- sum(r_random >= actual_mean, na.rm = T) / length(r_random)


# ATOH1 in healthy cohort -------------------------------------------------
atoh1.wins.new = c(557635:557637) # more or less the same as the windows in the old version

window2window_id %>% 
  filter(window_index %in% atoh1.wins.new) %>% 
  mutate(sig_name = "ATOH1") %>% 
  write.csv(paste0("~/Downloads/atoh1_wins.csv"), row.names = F)

# compile signature and evaluate on healthy reference 
system(paste("python3 ~/Documents/ModAnalysis/cli/compile_signatures.py", 
             "--winsig_input", "~/Downloads/atoh1_wins.csv",  
             "--winsig_output", "~/Downloads/atoh1_wins_filtered.csv",
             "--eval_healthy=True --eval_output=~/Downloads/atoh1_wins_healthy_eval.csv",
             "--ref_output", "~/Downloads/atoh1_wins_sig_compiled.csv"))


bind_rows(data_frame(sample_name_short = colnames(win_data_all), 
                     atoh1 = win_data_all[atoh1.wins,]/2.4) %>% #2.4 kb widtho of window
            inner_join(table1_info) %>%
            filter(grepl("^SCLC|NEC|Healthy", group)) %>% 
            select(atoh1, group), 
          data_frame(group = "Healthy", 
                     atoh1 = read.csv("~/Downloads/atoh1_wins_healthy_eval.csv") %>% 
                       filter(grepl("_H3K4me3", sample_name)) %>% 
                       pull(qbcw_norm_count))) -> data_atoh1

sample_counts <- data_atoh1 %>%
  group_by(group) %>%
  summarise(n = n()) %>%
  mutate(label = paste0(group, "\n(n = ", n, ")"))

# 2. Plot
ggplot(data_atoh1, aes(x = group, y = atoh1)) +
  geom_boxplot(outlier.shape = NA, lwd = base_line_size) +
  geom_jitter(data = data_atoh1 %>% 
                filter(group == "Healthy"), 
                       aes(color = group), 
                       width = 0.2, alpha = 0.5, size = .5) +
  geom_jitter(data = data_atoh1 %>% 
                filter(group != "Healthy"), 
              aes(color = group), 
              width = 0.2, size = 1) +
  # Use the custom labels for the x-axis
  scale_x_discrete(limits = sample_counts$group, labels = sample_counts$label) +
  # Professional styling
  scale_color_manual(values = c("Healthy" = "grey70", "SCLC" = "darkred", 
                                "NEC" = "darkred")) +
  labs(x = NULL, y = "normalized reads", title = "ATOH1 Signal Specificity") +
  theme(legend.position = "none") -> p
ggsaveG(p, paste0(figDirPaper, "JCI_resubmission/atoh1_specifict"), 50, 35)  
