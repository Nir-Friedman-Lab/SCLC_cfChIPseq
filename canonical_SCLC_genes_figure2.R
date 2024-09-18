# sclc vs. healthy canonical genes (genes shown in browser fig 2)
# canonical.genes = c("INSM1", "CHGA", "DLL3", "TAGLN3", "KIF5C", "CRMP1", "GAPDH")
canonical.genes = c("GAPDH", "DLL3", "INSM1", "CHGA", "CRMP1")
s = c(s.samples, h.samples)
data.frame(sample = s,
           group = sample.annotation[s,"group"], 
           t(chip_data_all[canonical.genes, s])) %>% 
  melt(id.vars = c("sample", "group") , variable.name = "gene") %>% 
  boxplotWOpoints(x = "gene", y = "value", fill = "group", 
                  ylab = "cfChIP reads") -> p
ggsave(paste0(figDirPaper, "figure2/boxplot_canonical_genes.pdf"), p, 
       height = 50, width = 50, units = "mm")
