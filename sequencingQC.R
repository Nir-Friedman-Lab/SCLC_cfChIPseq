# prepare QC table --------------------------------------------------------
sample_full = read.table("~/gavriel.fialkoff@mail.huji.ac.il - Google Drive/My Drive/CollatedFolders/Projects/NIH_SCLC/all_samples_March_2022.txt")$V1
sample_full = grep("SCLC", sample_full, invert = T, value = T)
sclc_combined = sub(".tagAlign.gz", "", list.files("~/Documents/SCLC_data/BED/H3K4me3/"))
sample_full = unique(c(sample_full, sclc_combined))
tibble(sample_name = sample_full, 
       sample_name_short = sub("_.*", "", sample_full)) %>% 
  filter(sample_name_short %in% a.samples) %>% 
  filter(!sample_name == "PS0011-1_NA_H3K4me3-39_26032021-13") ->
  data_long_short

sample.annotation %>% 
  filter(group %in% c("NEC", "NSCLC", "SCLC")) %>% 
  pull(Sample_id) -> samples_not_published
write.csv(data_long_short %>% 
            filter(sample_name_short %in% samples_not_published), "~/Documents/EGA_SCLC_JCI_data/manuscript_samples.csv", 
          row.names = F, quote = F)

qc.new = read.csv("~/Documents/SCLC_data/Output/H3K4me3/qc_all.csv", row.names = 1) 
qc.new %>% 
  rownames_to_column("sample_name") %>% 
  rename(total_unique = Total.uniq, 
         tss_pct_signal = X.Signal.at.TSS, 
         est_library_size = Total.uniq.est) %>% 
  filter(sample_name %in% data_long_short$sample_name) -> qc.new

library(RPostgres)
library(DBI)
dbcon <- dbConnect(RPostgres::Postgres(),
                   user = "gavriel.fialkoff@senseerahealth.com",
                   password = system(paste0("gcloud sql generate-login-token"), 
                                     intern = TRUE),
                   host = "34.165.66.84",
                   dbname = "senseera-prod",
                   sslmode='require')                     
dbGetQuery(dbcon, "SELECT sample_name, total_unique, tss_pct_signal, est_library_size
                    FROM sample_qc_view") %>% 
  mutate(sample_name_short = sub("_.*", "", sample_name)) -> qc.signal.yeild

bind_rows(
  qc.signal.yeild %>% 
    filter(sample_name %in% data_long_short$sample_name |
             sample_name == "SCLC0035-653_NA_H3K4me3-38_26032021-13") %>% 
    filter(!sample_name %in% qc.new$sample_name), 
  qc.new %>% 
    select(sample_name, total_unique, tss_pct_signal, est_library_size)
) %>% 
  mutate(sample_name = if_else(sample_name == "SCLC0035-653_NA_H3K4me3-38_26032021-13", 
                               "SCLC0037-653_NA_H3K4me3-38_26032021-13", 
                               sample_name), 
         sample_name_short = sub("_.*", "", sample_name)) -> 
  data_qc


# statistics --------------------------------------------------------------
data_qc %>% 
  pull(total_unique) %>% 
  mean() / 1e6

