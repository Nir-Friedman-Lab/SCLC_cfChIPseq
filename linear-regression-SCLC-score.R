# linear regression to estimate tumor load --------------------------------
diffGenesDir = paste0(baseDir, "Output/H3K4me3/DiffGenes/")
file.rename(paste0(diffGenesDir, list.files(diffGenesDir, "csv")),
            paste0(diffGenesDir, sub("_.*", ".csv", list.files(diffGenesDir, "csv"))))
diff.genes.files = list.files(diffGenesDir, ".csv")
diff.genes.files = diff.genes.files[sub(".csv", "", diff.genes.files) %in% s.samples]
diff.genes = data.frame(row.names = s.samples)
diff.genes$n.genes = sapply(s.samples, function(s) 
  sum(read.csv(paste0(baseDir, "Output/H3K4me3/DiffGenes/", s, 
                      ".csv"))$FoldChange..log2. > 0))
diff.genes.cutoff = round(quantile(diff.genes$n.genes, .90))
diff.genes.high.samples = rownames(diff.genes)[which(diff.genes$n.genes >= diff.genes.cutoff)]

sclc.arctype = rowMeans(chip_data_all[,diff.genes.high.samples])
reg.atlas = as.matrix(cbind(chip_data_all[,h.samples], sclc.arctype))
reg.mat = matrix(data = 0, nrow = length(a.samples), ncol = ncol(reg.atlas))
rownames(reg.mat) = a.samples; colnames(reg.mat) = colnames(reg.atlas)
for (s in a.samples) {
  print(grep(s, a.samples))
  res = round(nnls(A = reg.atlas[,grep(s, colnames(reg.atlas), invert = T)], b = chip_data_all[,s])$x, digits = 3)
  names(res) = colnames(reg.atlas)[grep(s, colnames(reg.atlas), invert = T)]
  reg.mat[s,names(res)] = res
}

data.frame(healthy = rowSums(reg.mat[,h.samples]), 
           SCLC = reg.mat[,"sclc.arctype"], 
           sample.annotation[rownames(reg.mat),]) %>% 
  mutate(SCLC.n = SCLC/(SCLC + healthy), 
         healthy.n = healthy/(SCLC + healthy)) -> estimated.tumor 
write.csv(estimated.tumor, paste0(baseDir, "tumor_estimation_diff_genes.csv"))

# currently, SCLC_score (which is called SCLC.n) appear in the metadata table
# metadata %>%
#   mutate(SCLC_score = estimated.tumor[Sample_id,"SCLC.n"]) -> metadata
# qplot(metadata$SCLC.n, metadata$SCLC_score)
# see comment at end of the SCLC-main script

# linear regression (tissue atlas) - CURRENTLY NOT USED ------------------------
reg.atlas = as.matrix(readRDS("~/gavriel.fialkoff@mail.huji.ac.il - Google Drive/My Drive/PHD_progress/data/gene.atals.clustered.rds"))

reg.mat = matrix(data = 0, nrow = length(a.samples), ncol = ncol(reg.atlas))
rownames(reg.mat) = a.samples; colnames(reg.mat) = colnames(reg.atlas)
for (s in a.samples) {
  print(grep(s, a.samples))
  res = round(nnls(A = reg.atlas, 
                   b = chip_data_all[Genes.notexcluded,s])$x, digits = 3)
  names(res) = colnames(reg.atlas)[grep(s, colnames(reg.atlas), invert = T)]
  reg.mat[s,names(res)] = res
}

high.tiss = names(which(colMaxs(reg.mat) > 0.05))
sort(reg.mat[s.samples,"erythroblast"])


reg.mat = data.frame(reg.mat)
reg.mat$group = sample.annotation[a.samples, "group"]
reg.mat$samp = rownames(reg.mat)
reg.mat.m = melt(reg.mat, id.vars = c("samp", "group"), variable.name = "tiss", 
                 value.name = "val")
d = reg.mat.m[reg.mat.m$tiss %in% high.tiss & 
                reg.mat.m$group %in% c("Healthy", "SCLC", "validation"),]
p = ggplot(d, aes(tiss, val, fill = group)) 
p = p + geom_boxplot()
p = p + facet_wrap(~tiss, scales = "free_x")
p = p + theme(axis.text.x = element_blank())
p = p + scale_fill_aaas()
ggsave(paste0(figDirPaper, "validation/tissue_composition.png"), p, width = 8, height = 6)
p = p + facet_wrap(~tiss, scales = "free")
ggsave(paste0(figDirPaper, "validation/tissue_composition_free.png"), p, width = 8, height = 6)
