library(edgeR)
library(Rsubread)
library(biomaRt)
library(edgeR)
library(tibble)
rna_data_old = readRDS(paste0(paperDir, "Figures/rna_data_matched_all.rds"))
# bam.file = "~/Documents/SCLC_data_from_NIH_box/RNA-seq/BAM/SB19_1645_RNA.bam"
# # gtf_file = paste0("~/gavriel.fialkoff@mail.huji.ac.il - Google Drive/",
# #                   "My Drive/cmv/2018-2019/ensembel_annotation.saf")
# # 
# gtf_file = paste0("~/Documents/hg19.refGene.gtf.gz")
# rna.counts = featureCounts(bam.file,annot.ext = gtf_file,
#                            isGTFAnnotationFile = T,
#                           #annot.inbuilt="hg19",
#                           useMetaFeatures=TRUE,
#                           GTF.attrType = "gene_name",
#                           # overlap between reads and features
#                           allowMultiOverlap=F, #check
#                           minOverlap=1,largestOverlap=T,readExtension5=0,
#                           # multi-mapping reads
#                           countMultiMappingReads=F,fraction=FALSE,
#                           # read filtering
#                           minMQS=10,
#                           ignoreDup=T,
#                           # strandness
#                           strandSpecific=0,
#                           # exon-exon junctions
#                           juncCounts=T,genome=NULL,
#                           # miscellaneous,
#                           isPairedEnd=T,
#                           nthreads=5)

# new_data = readRDS("~/Documents/SCLC_data_from_NIH_box/RNA-seq/SCLC_2024-04-04.rdata")
new_data = readRDS("~/Documents/SCLC_data_from_NIH_box/RNA-seq/SCLC_2024-07-30.rdata")
# new_data = readRDS("~/Documents/SCLC_data_from_NIH_box/RNA-seq/SCLC_2024-05-16.rdata")
exclude.genes = c() #c("IGLL5", "MIR3654") # genes that are very different between old and new analysis
rna_data_new = new_data$counts
notna = !is.na(rownames(rna_data_new)) & !rownames(rna_data_new) %in% exclude.genes
rna_data_new = rna_data_new[notna,]
new_data$annotation = new_data$annotation[notna,]
colnames(rna_data_new) = sub(".bam", "", colnames(rna_data_new))
colnames(rna_data_new) = sub(".Aligned.sortedByCoord.out", "", colnames(rna_data_new))
colnames(rna_data_new) = sub("_S24_part", "", colnames(rna_data_new))


metadata %>% 
  filter(`cfDNA Patient identical/duplicate` == "Identical") %>% 
  filter(Primary_tumor_sample %in% colnames(rna_data_new)) %>% 
  .$Primary_tumor_sample -> rna.chip.exist

setdiff(colnames(rna_data_new), rna.chip.exist) -> rna.chip.missing
# 2_WS_18_051587_S2 is a timepoint of SCLC0194-28
# 39_S18_11654_S37 (39_S18_11654) is a timepoint of SCLC0201-158
# 40_SF_20_216_B1_S39 is a timepoint of SCLC0275-602
# 48_NC_20_182_B1_S46 same date as 49_NC
# RS_19_5880 is a timepoint of SCLC0045-539_NA_H3K4me3-38_26032021-13 
# SB19_1389_RNA is a timepoint of SCLC0019-235_NA_H3K4me3-38_26032021-13
# ss_19_2398 is a timepoint of SCLC0037-478_NA_H3K4me3-39 
# SS19_1542_RNA is a timepoint of SCLC0030-435_NA_H3K4me3-39_26032021-13 

# missing RNA samples
metadata %>% 
  filter(!Primary_tumor_sample %in% (c("No tumor RNA", "No tumor RNAseq", NA))) %>% 
  filter(!Primary_tumor_sample %in% colnames(rna_data_new)) %>%
  filter(`cfDNA Patient identical/duplicate` == "Identical") %>% 
  .$Primary_tumor_sample %>% 
  unique() 

metadata %>% 
  filter(`cfDNA Patient identical/duplicate` == "Identical") %>% 
  select(Sample_id, Primary_tumor_sample) -> metadata.identical

# rna_data_new = rna_data_new[,rna.chip.exist]  

# sanity check
metadata.identical$Primary_tumor_sample[match(colnames(rna_data_new), 
                                   metadata.identical$Primary_tumor_sample)] == colnames(rna_data_new)
# colnames(rna_data_new) = metadata.identical$Sample_id[match(colnames(rna_data_new),
#                                                             metadata.identical$Primary_tumor_sample)]

# samples that were in initial rna_data but not in current
sub("_.*", "", colnames(rna_data_old))[!sub("_.*", "", colnames(rna_data_old)) %in% colnames(rna_data_new)]
# SCLC0046-572/SCLC0035-561 were the identical samples.  in the new validation 
# set there are better matching samples timewise 

rna.dge = DGEList(counts = rna_data_new, genes = rownames(rna_data_new))
# TMM-normalization
rna.dge = calcNormFactors(rna.dge, method = "TMM", )

rna.cpm = cpm(rna.dge)
rna.rpkm = rpkm(rna.dge, new_data$annotation$Length)
# rna.rpkm = log2(1+rna.rpkm)

hk.genes = CommonGenes[CommonGenes %in% rownames(rna_data_new)]
d = log10(1+rna.cpm[hk.genes,])
d = sweep(d, 1, rowMedians(d), "-")
data.frame(row.names = sub(".Aligned.sortedByCoord.out.bam", "", colnames(new_data$counts)), 
           depth = log10(colSums(new_data$counts))) -> col.a
pheatmap::pheatmap(d, show_rownames = F, fontsize = 7, annotation_col = col.a, 
                   filename = paste0(figDirPaper, "validation/rna_common_genes.png"), 
                   width = 12, height = 9)
# samples that seem very different then all other samples 
discard.samples = colnames(rna_data_new)[c(grep("24_AU|2387|CL0107|0422_T3R", colnames(rna_data_new)))]
d = d[,!colnames(d) %in% discard.samples]
pheatmap::pheatmap(d, show_rownames = F, fontsize = 5, annotation_col = col.a, 
                   filename = paste0(figDirPaper, "validation/rna_common_genes_filtered.png"), 
                   width = 12, height = 9)

rna_data_new = rna_data_new[,!colnames(rna_data_new) %in% discard.samples]
rna.dge = DGEList(counts = rna_data_new, genes = rownames(rna_data_new))
# TMM-normalization
rna.dge = calcNormFactors(rna.dge, method = "TMM", )

rna.cpm = cpm(rna.dge)
rna.rpkm = rpkm(rna.dge, new_data$annotation$Length)

saveRDS(rna.cpm, "~/gavriel.fialkoff@mail.huji.ac.il - Google Drive/Shared drives/Friedman Lab Shared Drive/BloodChIP/Analysis/Projects/NIH_SCLC/new_data/SCLC_cpm.rds")
saveRDS(rna.cpm, "~/Documents/SCLC_data_from_NIH_box/RNA-seq/SCLC_cpm.rds")
saveRDS(rna.rpkm, "~/Documents/SCLC_data_from_NIH_box/RNA-seq/SCLC_rpkm.rds")
saveRDS(rna.rpkm, "~/gavriel.fialkoff@mail.huji.ac.il - Google Drive/Shared drives/Friedman Lab Shared Drive/BloodChIP/Analysis/Projects/NIH_SCLC/new_data/SCLC_rpkm.rds")

new_data$stat %>% 
  column_to_rownames("Status") %>% 
  t() %>% 
  data.frame() %>% 
  select(Assigned) -> sample.cov
rownames(sample.cov) = sub(".bam", "", rownames(sample.cov)) 
sample.cov = sample.cov[exist.ind,, drop = F]
rownames(sample.cov) = metadata.identical$Sample_id[match(exist.ind, metadata.identical$Primary_tumor_sample)]  
  

matched.samples = intersect(colnames(rna_data), colnames(rna.rpkm))
matched.genes = intersect(rownames(rna_data),rownames(rna.rpkm))
rna_new_matched = as.matrix(rna.rpkm[matched.genes,matched.samples])
rna_old_matched = as.matrix(rna_data[matched.genes,matched.samples])


data.frame(gene = rep(rownames(rna_new_matched), ncol(rna_new_matched)),
           sample = rep(colnames(rna_new_matched), each = nrow(rna_new_matched)),
           new = as.numeric(rna_new_matched), 
           old = as.numeric(rna_old_matched)) %>% 
  left_join(new_data$annotation, join_by("gene" == "GeneID")) %>% 
  left_join(sample.cov %>% rownames_to_column("sample"), join_by(sample))-> df


DensityScatter(df$new,df$old, threshold = 1)
df %>% 
  sample_n(10000) %>% 
  ggplot(aes(new, old)) + 
  geom_point(size = .2) + 
  geom_abline() +
  labs(x = "new analysis", y = "old analysis") +
  geom_density_2d() +
  geom_point(data = df[df$gene == "NBPF14",], color = "red", size = 1) -> p
ggsave(paste0(figDirPaper, "old_vs_new_rna.pdf"), p, width = 60, height = 60, 
       units = "mm")

df %>% 
  ggplot(aes(new - old, Length)) + 
  geom_point(size = .2) -> p
ggsave(paste0(figDirPaper, "old_new_diff_vs_length.pdf"), p, width = 60, height = 60, 
       units = "mm")  

df %>% 
  ggplot(aes(new, old, color = Assigned)) + 
  geom_point(size = .2) -> p
ggsave(paste0(figDirPaper, "old_new_diff_vs_depth.pdf"), p, width = 60, height = 60, 
       units = "mm")  

x = c(subtypes, "ATOH1")
pheatmap::pheatmap(rna_new_matched[x,] - rna_old_matched[x,], 
                   annotation_col = sample.cov,
                   # filename = paste0(figDirPaper, "subytpe_difference_v1.png"),
                   fontsize = 7, width = 10, height = 8)

sort(table(rownames(which(rna_new_matched > 9 & rna_old_matched < 2, arr.ind = T))))
sort(table(rownames(which(rna_old_matched > 8 & rna_new_matched < 1, arr.ind = T))))
sort(rna_ol)

dim(rna_old_matched)
