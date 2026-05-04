flow_data_file = "~/gavriel.fialkoff@mail.huji.ac.il - Google Drive/Shared drives/Friedman Lab Shared Drive/BloodChIP/Analysis/Projects/NIH_SCLC/external_data/Flow Result Percentages for TMSCLC 20_18_9_7_3.xlsx"
readxl::excel_sheets(flow_data_file)


# I chose for every flow samples the corresonding patient with highest SCLC score
# metadata %>% 
# filter(grepl("SCLC0059|SCLC0061", Sample_id)) %>% 
#   select(Sample_id, SCLC.n) %>% 
#   mutate(ctdna = as.numeric(SCLC.n) * 100)
readxl::read_xlsx(flow_data_file, sheet = "Heatmap for all plates combined") %>% 
  dplyr::rename("gene" = "Description", 
                "SCLC0279-878" = "TMSCLC 18\r\n(% positive) - Jeff", 
                "SCLC0030-435" = "TMSCLC 3\r\n(% positive) - Brett", 
                "SCLC0059-681" = "TMSCLC 9\r\n(% positive) - Brett",
                "SCLC0061-686" = "TMSCLC 7\r\n(% positive) - Brett") %>% 
  select(c(gene, "SCLC0279-878", "SCLC0030-435", "SCLC0059-681", "SCLC0061-686")) ->
  flow_data


shared_genes = intersect(rownames(chip_data_all), flow_data$gene)
flow_data %>% 
  filter(gene %in% shared_genes) %>%
  column_to_rownames("gene") %>% 
  pheatmap::pheatmap()

log2(1+chip_data_all[shared_genes, colnames(flow_data)[2:ncol(flow_data)]]) %>% 
  pheatmap::pheatmap()
  
