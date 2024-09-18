# load libraries  ----------------------------------------
library(matrixStats)
library(reshape2)
library(ggplot2)
library(readxl)
library(ggrepel)
library(ggsci)
library(rtracklayer)
library(nnls)
library(plotly)
# source("~/gavriel.fialkoff@mail.huji.ac.il - Google Drive/My Drive/PHD_progress/deconvolution-functions.R")
source("generic-plots.R")
source("~/BloodChIP/git/Core/NMF-util.R")
# source("~/BloodChIP/git/Core/DensityScatter.R")
# source("~/BloodChIP/git/Apps/CompareSampleGroups.R")
library(ggthemes)
library(factoextra)
library(dendextend)
library(gridExtra)
library("stringr")
library(ggforce)
library(ggbreak)
library(ggpubr)
library(eulerr)
library(qusage)
library(cowplot)
library(grid)
library(gridExtra)
library(GGally)
library(zoo)
library(pheatmap)
library(ComplexHeatmap)
library(RColorBrewer)
library(pROC)
library(circlize)
library(ggrastr)
library(GSVA)
library(MASS)
library(conflicted)
library(reporter)
library(magrittr)
suppressPackageStartupMessages(library(ggstream))
suppressPackageStartupMessages(library(googlesheets4))
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
conflicts_prefer(base::as.matrix)
conflicts_prefer(base::setdiff)
conflicts_prefer(stats::dist)
conflicts_prefer(base::intersect)
conflicts_prefer(matrixStats::rowMedians)
conflicts_prefer(dplyr::count)
conflicts_prefer(dplyr::rename)
conflicts_prefer(stats::as.dist)
conflicts_prefer(matrixStats::rowMaxs)
# load setup files --------------------------------------------------------
LoadMeta = FALSE
LoadSig = TRUE
LoadQC = FALSE
LoadDescription = TRUE
LoadGenomicGC = FALSE
ReadOptions = FALSE

# use these to overwrite defaults
DefaultMod = "H3K4me3"
SourceDIR = "~/BloodChIP/cfChIP-Development/Pipeline/"
source(paste0(SourceDIR,"Core/cfChIP-Cmd-Common.R"))
source(paste0(SourceDIR,"Core/cfChIP-Functions.R"))
# set variables -----------------------------------------------------------
base_theme = theme_classic(base_size = base_size, base_family = "Helvetica", 
                           base_line_size = base_line_size)
theme_set(base_theme)

# variables
baseDir = "~/gavriel.fialkoff@mail.huji.ac.il - Google Drive/Shared drives/Friedman Lab Shared Drive/BloodChIP/Analysis/Projects/NIH_SCLC/new_data/"
# baseDirNew = "~/Documents/SCLC_data/"
paperDir = "~/gavriel.fialkoff@mail.huji.ac.il - Google Drive/My Drive/CollatedFolders/Projects/NIH_SCLC/cfChIP-paper/"
phdDir = "~/gavriel.fialkoff@mail.huji.ac.il - Google Drive/My Drive/PHD_progress/"
figDirPaper = paste0(paperDir, "Figures/")
tableDir = paste0(paperDir, "Tables/")
rawdataDir = paste0(baseDir, "Samples/H3K4me3/")
externalDataDir = "data/"
senseeraDir = "~/gfialk@gmail.com - Google Drive/Shared drives/Senseera-Nir-Friedman-Lab/"
senseeraDatadir = paste0(senseeraDir, "Data/Samples/H3K4me3/")

# sclc.histology = c("SCLC", "SCLC transformed from EGFR mut NSCLC",
#                    "EGFR-mutated transformed SCLC")
cutoff.yeild = 400000
cutoff.signal = 20
# cfm.color = colorRampPalette(c("blue", "black", "red"))(50)
cfm.color = colorRampPalette(c("#00A7FE", "black", "#E54C00"))(50)
cfm.color.rna = colorRampPalette(c("#42CE3F", "black", "#E20119"))(50)

cfm.breaks = seq(-3,3,length = 50)
sclc.breaks = c(0,.25,.5,.75,1); sclc.lab = c("Healthy\nlike","0.25", "0.5", 
                                              "0.75", "SCLC\nlike")
cfm.leg = expression(paste(Delta, " (mean log2 reads)"))
# colors from https://colorbrewer2.org/#type=diverging&scheme=BrBG&n=3
group.colors = list(Healthy = "#1b9e77", SCLC = "#d95f02", "NEC" = "#d95f02", 
                    NSCLC = "#7570b3", 
                    Covid = "#e7298a", roadmap = "#66a61e", AMI = "#e6ab02", 
                    Liver = "#a6761d", CRC = "#666666", low = "#fee8c8", high = "#990000", na.value = "#d9d9d9",
                    ASCL1 = "#e31a1c", NEUROD1 = "#ff7f00", "ASCL1-NEUROD1" = "#fb9a99", 
                    POU2F3 = "#1f78b4", YAP1 = "#33a02c", pou2f3.rna = "#1f78b4", 
                    ascl1.rna = "#e31a1c", neurod1.rna = "#ff7f00", yap1.rna = "#33a02c", 
                    resistant = "#e41a1c", sensitive = "#4daf4a", . = "#bdbdbd", 
                    Responder = "#984ea3", NonResponder = "#ff7f00", "NA" = "#d9d9d9", 
                    "SCLC\npre" = "#d95f02", "SCLC\npost" = "#d95f02", 
                    hm.low = "#00A7FE", hm.mid = "black", hm.high = "#E54C00")
group.colors.heatmap = list(NE = c(high = "#e31a1c", low = "#1f78b4"), ne.stat = c(NE = "#e31a1c", NonNE = "#1f78b4"), 
                            subtype = c(ASCL1 = "#e31a1c", NEUROD1 = "#ff7f00", "ASCL1-NEUROD1" = "#fb9a99", 
                                        POU2F3 = "#1f78b4", YAP1 = "#33a02c", "no-RNA"="#d9d9d9", "."="#d9d9d9"), 
                            ne.score = colorRampPalette(c("white", "#e31a1c"))(10), 
                            group = c(Healthy = "#1b9e77", SCLC = "#d95f02"), 
                            response = c(Responder = "#984ea3", NonResponder = "#ff7f00")) 
color.map = list(group = c(Healthy = "#1b9e77", SCLC = "#d95f02", NSCLC = "#7570b3", CRC = "#666666"))#, 
# roadmap = "#66a61e", AMI = "#e6ab02", Liver = "#a6761d"))

# for Heatmaps
gp = gpar(fontsize = base_size)
# SCLC score cutoffs
high.score.cutoff = .5
low.score.cutoff = .05

# SCLC cannonical genes
# sclc.genes = c("ASCL1", "NEUROD1", "POU2F3", "YAP1", "ATOH1", "BCL2", "NFIB", 
#                "RET", "SOX2", "DLL3", "FOXA2", "CHGA", "INSM1", "SYP")
sclc.subtypes = c("ASCL1", "NEUROD1", "YAP1", "POU2F3", "ATOH1")

# load data and metadata ---------------------------------------------------------------
gene.atlas = readRDS(paste0(phdDir, "data/gene.atlas.rds")) 
gene.atlas = gene.atlas$gene.atlas[Genes.notexcluded,]
win.atlas = readRDS("~/gavriel.fialkoff@mail.huji.ac.il - Google Drive/Shared drives/Friedman Lab Shared Drive/BloodChIP/Atlases/Win.Atlas.named.rdata")

# win_data_all = sapply(a.samples, function(s) readRDS(paste0(rawdataDir, s, ".rdata"))$Counts.QQnorm)
# rownames(chip_data_all) = rownames(gene.atlas)
# chip_data_new = sapply(new.samples.passQC, function(s) 
#   readRDS(paste0(1, s, ".rdata"))$GeneCounts.QQnorm)
# rownames(chip_data_new) = rownames(chip_data_all)
# win_data_new = sapply(new.samples.passQC, function(s) 
#   readRDS(paste0(senseeraDatadir, s, ".rdata"))$Counts.QQnorm)
# chip_data_all = cbind(chip_data_all, chip_data_new)
# win_data_all = cbind(win_data_all, win_data_new)
# saveRDS(win_data_all, paste0(baseDir, "win_data_all.rds"))
# saveRDS(chip_data_all, paste0(baseDir, "chip_data_all.rds"))

############## final version after sample correction ##########################
win_data_all = readRDS(paste0(baseDir, "win_data_all.rds"))
chip_data_all = readRDS(paste0(baseDir, "chip_data_all.rds"))
# SCLC0126-223 - gentic missmatch.  # change name of SCLC0035-653 to SCLC0037-653
chip_data_all = chip_data_all[,grep("SCLC0126-223", colnames(chip_data_all), invert = T)]
win_data_all = win_data_all[,grep("SCLC0126-223", colnames(win_data_all), invert = T)]
colnames(chip_data_all) = sub("35-653", "37-653", colnames(chip_data_all))
colnames(win_data_all) = sub("35-653", "37-653", colnames(win_data_all))
a.samples = colnames(chip_data_all)
chip.matched = readRDS(paste0(baseDir, "SCLC_ChIP_fixed.rds"))

#roadmap data 
log2(1+read.csv(paste0(externalDataDir, "57epigenomes.RPKM.gene-names.pc.gz"), 
                              sep = " ", header = 1, row.names = 1)) -> rna_roadmap
read.csv(paste0(externalDataDir, "EG.name.txt"), sep = "\t", 
                         header = 1, row.names = 1) -> roadmap.names 
colnames(rna_roadmap)[match(rownames(roadmap.names), colnames(rna_roadmap))] == rownames(roadmap.names)
colnames(rna_roadmap)[match(rownames(roadmap.names), colnames(rna_roadmap))] = roadmap.names$Universal_Human_Reference
rna_roadmap$E000 = NULL

# updated Aug 2024##
# samp.composition.clus = read.csv(paste0(paperDir, "Figures/sample_composition_clustered.csv"), row.names = 1)
# healthy.estimate = colSums(samp.composition.clus[rownames(samp.composition.clus) %in% c("megakaryocyte", "monocyte.macrophage", "Neutrophil", "eosinophil"),])
rna_data = readRDS("~/gavriel.fialkoff@mail.huji.ac.il - Google Drive/Shared drives/Friedman Lab Shared Drive/BloodChIP/Analysis/Projects/NIH_SCLC/new_data/SCLC_RNA_fixed.rds")
matched.genes = intersect(rownames(rna_data), rownames(chip_data_all))
url = "https://docs.google.com/spreadsheets/d/1hkxqRSBGzbQC-UhvHHFOB5qxd_SW3SjwMuMJj0BEwD0?usp=drive_fs"
sid = "1hkxqRSBGzbQC-UhvHHFOB5qxd_SW3SjwMuMJj0BEwD0"
gs4_deauth()
gs4_auth() 
googlesheets4::sheet_names(sid)
googlesheets4::read_sheet(sid, sheet = "Blood Samples") %>%
  as_tibble() -> metadata

googlesheets4::read_sheet(sid, sheet = "Timelines") %>%
  as_tibble() -> metadata.timeline

googlesheets4::read_sheet(sid, sheet = "Patients") %>%
  as_tibble() -> metadata.patients

googlesheets4::read_sheet(sid, sheet = "Biopsies") %>%
  as_tibble() -> metadata.biopsies

googlesheets4::read_sheet(sid, sheet = "RNA-ChIP pairing") %>%
  as_tibble() -> metadata.matching

metadata %>% 
  left_join(metadata.patients %>% select(SCLC.state, PatientID, Histology),
           by = "PatientID") %>%
  mutate(ctDNA = as.numeric(ctDNA), 
         CT_volume_ALL = as.numeric(CT_volume_ALL), 
         RECIST = as.numeric(RECIST), 
         cfDNA = as.numeric(cfDNA), 
         CTC = as.numeric(CTC))-> metadata


# qc = read.csv(paste0(baseDir, "Output/H3K4me3/May_2022_all_qc.csv"), row.names = 1)
# qc = read.csv("~/gavriel.fialkoff@mail.huji.ac.il - Google Drive/Shared drives/Friedman Lab Shared Drive/BloodChIP/Analysis/Projects/NIH_SCLC/new_data/Output/H3K4me3/qc_all.csv", 
#               row.names = 1)
#low.qc.samples = rownames(qc)[qc$TSS < cutoff.yeild | qc$X.Signal.at.TSS < cutoff.signal]
# low.qc.samples.short = sub("_.*", "", low.qc.samples)
# passed.qc.samples = a.samples[!a.samples %in% low.qc.samples.short]
# till here 

estimated.tumor = read.csv(paste0(baseDir, "tumor_estimation_diff_genes.csv"), row.names = 1)
estimated.tumor$group = factor(estimated.tumor$group)
healthy.ref = readRDS("/Users/gavrielfialkoff/Documents/cfChIP_Pipeline/SetupFiles/hg19/H3K4me3/HealthyMixModel.rds")
names(healthy.ref$HealthyGeneModel$mean.0) = healthy.ref$HealthyGeneModel$Gene
healthy.ref = healthy.ref$HealthyGeneModel$mean.0
healthy.ref = healthy.ref[Genes.notexcluded]

TSS.windows = readRDS("/Users/gavrielfialkoff/Documents/cfChIP_Pipeline/SetupFiles/hg19/H3K4me3/Windows.rds")

# define groups -----------------------------------------------------------
metadata %>% 
  filter(grepl("Pre", Timepoint), SCLC.state == "SCLC") %>% 
  pull(Sample_id) -> pre.samples

metadata %>% 
  filter(!grepl("Pre", Timepoint), SCLC.state == "SCLC") %>% 
  pull(Sample_id) -> post.samples

metadata %>% 
  filter(SCLC.state == "SCLC") %>% 
  pull(Sample_id) -> s.samples

metadata %>% 
  filter(SCLC.state == "NEC") %>% 
  pull(Sample_id) -> nec.samples

metadata %>% 
  filter(SCLC.state == "other", !grepl("^NSCLC", Histology)) %>% 
  pull(Sample_id) -> o.samples

a.samples[grep("PS", a.samples)] -> h.samples
a.samples[grep("CRC", a.samples)] -> c.samples
c(a.samples[grep("^LC", a.samples)],  metadata %>% 
    filter(Histology == "NSCLC") %>% 
    pull(Sample_id)) -> l.samples
sh.samples = c(s.samples, h.samples)

data.frame(row.names = colnames(chip_data_all), 
           sample = colnames(chip_data_all)) %>% 
  mutate(group = case_match(
    sample,
    c.samples ~ "CRC",
    h.samples ~ "Healthy",
    l.samples ~ "NSCLC",
    o.samples ~ "other", 
    nec.samples ~ "NEC",
    s.samples ~ "SCLC",
    .default = NA
  )) -> sample.annotation

healthy.mean = rowMeans(win_data_all[,h.samples])
diff.genes.common = read.table(paste0(figDirPaper, "common.diff.genes.txt"), header = F)$V1

# till here

# rna.samples = colnames(rna_data)
# # rna.samples.chip.passQC = rna.samples[rna.samples %in% sv.samples] # change May-2024
# rna.samples.chip.passQC = rna.samples[rna.samples %in% passed.qc.samples.short]
# # rna_samples not in sv.samples are low QC or NSCLC
# 
# metadata %>%
#   filter(Sample_id %in% rna.samples &
#            # Sample_id %in% sv.samples & # old version didn't have this filter
#            `Tumor_ctDNA_timing_match (determined based on treatment timing standpoint)` 
#          == "Yes") %>%
#   select(Sample_id) %>%
#   .$Sample_id -> rna.samples.matched.time
# 
# high.ctDNA.samples = s.samples[estimated.tumor[s.samples, "SCLC.n"] > high.score.cutoff]
# write.table(high.ctDNA.samples, paste0(baseDir, "Samples-High-SCLC.txt"), quote = F, col.names = F, row.names = F)
# high.ctDNA.samples.wRNA = high.ctDNA.samples[high.ctDNA.samples %in% rna.samples.chip.passQC]
# high.ctDNA.samples.wRNA.matched.time = high.ctDNA.samples[high.ctDNA.samples %in% rna.samples.matched.time]
# med.ctDNA.samples = sv.samples[estimated.tumor[sv.samples, "SCLC.n"] > med.score.cutoff]
# med.ctDNA.samples.wRNA = med.ctDNA.samples[med.ctDNA.samples %in% rna.samples]
# low.ctDNA.samples = sv.samples[estimated.tumor[sv.samples, "SCLC.n"] > low.score.cutoff]
# low.ctDNA.samples.wRNA = low.ctDNA.samples[low.ctDNA.samples %in% rna.samples]


# selected_groups = list("SCLC" = s.samples, "Healthy" = h.samples, "NSCLC" = l.samples, "CRC" = c.samples)
# write.csv(metadata[high.ctDNA.samples,c(6,7,8,11,12,53,19,20,36)], paste0(baseDir, "Samples-SCLC-attr.csv"))



# maybe move somewhere else
pou2f3.wins = c(1058780,1058781)
neurod1.wins = 243856
ascl1.wins = c(1127785, 1127786)
atoh1.wins = 436305



# "SCLC0176-120_NA_H3K4me3-39" has very high erythrocyte (which shares some genes with SCLC). this leads to overestimate of 
# its SCLC score and other discrepancies (e.g. between ASCL1 in RNA and ChIP).  same but less so for SCLC0267-407_NA_H3K4me3-38 and 
# SCLC0034-512_NA_H3K4me3-38

