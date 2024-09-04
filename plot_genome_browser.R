# prepare data ------------------------------------------------------------
trackDir = paste0(baseDir, "Tracks/H3K4me3/")
trackDirNew = paste0(baseDir, "new_data/Tracks/H3K4me3/")


# canonical genes ---------------------------------------------------------
# s.tracks = c("SCLC0146-1504_NA_H3K4me3-L", "SCLC0193-1753_NA_H3K4me3-L", 
#                   "SCLC0109-897_NA_H3K4me3-L", "SCLC0030-435_NA_H3K4me3-39_26032021-13", 
#                   "SCLC0174-261_NA_H3K4me3-38_26032021-13", "SCLC0176-120_NA_H3K4me3-39")
h.tracks = c("PS0001-2_NA_H3K4me3-39_26032021-13.bw", 
             "PS0002-1_NA_H3K4me3-12_05022021-6.bw",
             "PS0003-1_NA_H3K4me3-8_07012021-4.bw")
s.tracks = c("SCLC0146-1504_NA_H3K4me3-L_21062020-L.bw", 
             "SCLC0193-1753_NA_H3K4me3-L_21062020-L.bw", 
             "SCLC0109-897_NA_H3K4me3-L_21062020-L.bw", 
             "SCLC0030-435_NA_H3K4me3-39_26032021-13.bw", 
             "SCLC0174-261_NA_H3K4me3-38_26032021-13.bw", 
             "SCLC0176-120_NA_H3K4me3_NA.bw")
all.tracks = c(s.tracks, h.tracks)
all.path = paste0(trackDirNew, all.tracks)
s.color = rep(group.colors$SCLC, length(s.tracks))
h.color = rep(group.colors$Healthy, length(h.tracks))
color_code = c(s.color, h.color)
data.frame(sample_id = sub("_.*", "", all.tracks), color_code = color_code, 
           bw_path = all.path) %>%
  mutate(sample_id = sub("PS", "H", sample_id)) -> bw.data


gr = data.frame(name = c("GAPDH", "DLL3", "INSM1", "CHGA", "CRMP1"), 
                         chr = c("chr12", "chr19", "chr20", "chr14", "chr4"),
                         start = c(6635000,39985000, 20344000, 93375000, 5818000), 
                         end = c(6657000, 40000000,20356000, 93405000, 5900000))
gr$win = paste0(gr$chr, ":", gr$start, "-", gr$end)
gr$comma_win = paste0(gr$chr, ":", scales::comma(gr$start, scale = 1e-3),  "-", 
                      scales::comma(gr$end, scale = 1e-3), "(kb)")
p = plotBrowser(bw.data, gr)
ggsave(paste0(figDirPaper, "figure2/browser_canonical_genes.pdf"), p, 
       width = 120, height = 50, units = "mm")
ggsave(paste0(figDirPaper, "figure2/browser_canonical_genes.png"), p, width = 110, height = 125, units = "mm")
# ggplot custom made ------------------------------------------------------
bw.samples = c("SCLC0024-293_NA_H3K4me3-39", "SCLC0106-705_NA_H3K4me3-L", "SCLC0033-489_NA_H3K4me3-39_26032021-13", 
               "SCLC0030-435_NA_H3K4me3-39_26032021-13", "SCLC0152-94_NA_H3K4me3-38", "SCLC0147-1529_NA_H3K4me3-L_21062020-L", 
               "PS0001-2_NA_H3K4me3-39_26032021-13", "PS0002-1_NA_H3K4me3-12_05022021-6", 
               "PS0003-1_NA_H3K4me3-8_07012021-4")
bw.data.pou2f3 = data.frame(sample_id = gsub("_.*", "", bw.samples), 
                            bw_path = paste0(trackDir, bw.samples, ".bw"))
bw.data.pou2f3$color_code = group.colors$SCLC
bw.data.pou2f3$color_code[grep("PS", bw.data.pou2f3$sample_id)] = group.colors$Healthy
gr = data.frame(name = c("GAPDH", "POU2F3", "ATOH1"), 
                       chr = c("chr12", "chr11", "chr4"),
                       start = c(6635000, 120093000, 94740000), 
                       end = c(6657000, 120170000, 94760000))
gr$win = paste0(gr$chr, ":", gr$start, "-", gr$end)
gr$comma_win = paste0(gr$chr, ":", scales::comma(gr$start, scale = 1e-3),  "-", 
                      scales::comma(gr$end, scale = 1e-3), "(kb)")
p = plotBrowser(bw.data = bw.data.pou2f3, genomic_ranges = gr)
ggsave(paste0(figDirPaper, "figureS4/pou2f3_atoh1.pdf"), p, width = 90, height = 50, units = "mm")


# TF genes ----------------------------------------------------------------
gr = data.frame(name = c("GAPDH", "ASCL1", "NEUROD1", "YAP1"), 
                chr = c("chr12", "chr12", "chr2", "chr11"),
                start = c(6635000,103345000, 182525000, 101970000), 
                end = c(6657000,103365000, 182555000, 102070000))
gr$win = paste0(gr$chr, ":", gr$start, "-", gr$end)
gr$comma_win = paste0(gr$chr, ":", scales::comma(gr$start, scale = 1e-3),  "-", 
                      scales::comma(gr$end, scale = 1e-3), "(kb)")
p = plotBrowser(bw.data.pou2f3, gr)
ggsave(paste0(figDirPaper, "figureS4/TF_browser.pdf"), p, width = 150, height = 50, units = "mm")


# ATACseq differential peaks ----------------------------------------------
bw.data = read.csv(paste0(atacFigDir, "bw_data.csv"), comment.char = "$")
gr = read.csv(paste0(atacFigDir, "ATAC_diffpeaks.csv"), comment.char = "#")
gr$win = paste0(gr$chr, ":", gr$start, "-", gr$end)
gr$comma_win = paste0(gr$chr, ":", scales::comma(gr$start, scale = 1e-3),  "-", 
                      scales::comma(gr$end, scale = 1e-3), "(kb)")
p = plotBrowser(bw.data, gr)
p = p + theme(text = element_text(size = base_size/.pt))
ggsave(paste0(atacFigDir, "diff_peaks.pdf"), p, width = 170, height = 50, 
       units = "mm")

# old ---------------------------------------------------------------------




if (0) { 
  library(Homo.sapiens)
  library(TxDb.Hsapiens.UCSC.hg19.knownGene)
  library(BSgenome.Hsapiens.UCSC.hg19)
  library(ggbio)
  genes_all = AnnotationDbi::keys(Homo.sapiens, keytype='GENEID')
  generanges = AnnotationDbi::select(Homo.sapiens, 
                                     keys=genes_all, keytype='GENEID',
                                     columns=c('GENEID', 'SYMBOL', 'TXCHROM', 'TXSTART', 'TXEND'))
  # convert to GRanges object
  generanges = makeGRangesFromDataFrame(generanges, keep.extra.columns=T)
  # we'll use RSPO1 gene as an example
  generanges = generanges[!is.na(generanges$SYMBOL)]
  data(genesymbol, package = "biovizBase")
  p = autoplot(Homo.sapiens, which=generanges[generanges$SYMBOL=='RSPO1'], 
           gap.geom='chevron', col='blue', fill='blue') +  theme_bw()
  ll = sapply(1:nrow(gr), function(i) 
    autoplot(Homo.sapiens, which=GRanges(gr$win[i]), stat = "identity")) #, stat = "reduce")
  p1 = autoplot(Homo.sapiens, which=GRanges(gr$win[i]), stat = "identity")
  
  plot_grid(ll[[1]], ll[[2]], nrow = 1)
  
  x = TxDb.Hsapiens.UCSC.hg19.knownGene
  exons = exonsByOverlaps(x, GRanges(gr$win[1]))
  transcriptsByOverlaps(x, GRanges(gr$win[1]))  
  disjointExons(x)
  genes(Homo.sapiens)
  transcript(x)
  }





# ALPS package (very limited) ---------------------------------------------
library("ALPS")
gene_range = c(DLL3 = "chr19:39987586-40001109",
               INSM1 = "chr20:20343513-20356822")
p1 = plot_browser_tracks(data_table = data.browser, gene_range = unname(gene_range[1]), 
                        ref_gen = "hg19")
p2 = plot_browser_tracks(data_table = data.browser, gene_range = unname(gene_range[2]), 
                         ref_gen = "hg19")
plot_grid(p1, p2)
class(p1[[1]])

# signac package ----------------------------------------------------------
# based on https://satijalab.org/signac/articles/pbmc_vignette.html
library(Signac)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v75)
library(Seurat)
library(patchwork)
counts = Read10X_h5(filename = paste0(figDirPaper, "/signac_genome_browser/atac_v1_pbmc_10k_filtered_peak_bc_matrix.h5"))
metadata = read.csv(
  file = paste0(figDirPaper, "/signac_genome_browser/atac_v1_pbmc_10k_singlecell.csv"),
  header = TRUE,
  row.names = 1
)
chrom_assay = CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  genome = 'hg19',
  fragments = paste0(figDirPaper, 'signac_genome_browser/atac_v1_pbmc_10k_fragments.tsv.gz'),
  min.cells = 10,
  min.features = 200
)
pbmc = CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = metadata
)
class(pbmc)
# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v75)
# change to UCSC style since the data was mapped to hg19
seqlevelsStyle(annotations) <- 'UCSC'
x = CreateSeuratObject(counts = chip_data_all[,s.samples], assay = "chip")
granges(x)
Annotation(x) = annotations
bw_path = paste0(trackDir, all.tracks, ".bw")
data.browser = data.frame(sample_id = all.tracks, bw_path = bw_path)
region = TSS.windows[grep("DLL3", TSS.windows$name),]
region = region[5]
region = c("DLL3", "INSM1")

# add the gene information to the object
CoverageBrowser(x, region)
CoveragePlot(object = x, region = region, features = "gene names", 
             assay = "coverage", bigwig = bw_path, bigwig.type = "coverage", 
             annotation = T, extend.upstream = 1000, extend.downstream = 1000)
# next attempt karyoplotR -------------------------------------------------
# see https://bernatgel.github.io/karyoploter_tutorial//Examples/EncodeEpigenetics/EncodeEpigenetics.html
library(karyoploteR)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

pp <- getDefaultPlotParams(plot.type=1)
pp$leftmargin <- 0.15
pp$topmargin <- 15
pp$bottommargin <- 15
pp$ideogramheight <- .5
pp$data1inmargin <- 10
pp$data1outmargin <- 0
pp$fontSize = .7
pp$ymin = 0

region = toGRanges("chr19:39983586-40006109"); reg.name = "DLL3"
region = toGRanges("chr7:5563219-5573736"); reg.name = "ActB"

pdf(paste0(figDirPaper, "figure2/", reg.name, "_browser.pdf"))
kp <- plotKaryotype(zoom = region, cex=pp$fontSize, plot.params = pp)
genes.data = makeGenesDataFromTxDb(TxDb.Hsapiens.UCSC.hg19.knownGene,
                                    karyoplot=kp,
                                    plot.transcripts = TRUE, 
                                    plot.transcripts.structure = TRUE)
genes.data <- addGeneNames(genes.data)
genes.data <- mergeTranscripts(genes.data)
kpAddBaseNumbers(kp, tick.dist = 10000, minor.tick.dist = 10000,
                 add.units = TRUE, cex=.5, tick.len = 3)
# kpAddMainTitle(kp, region_name, cex=pp$fontSize)
kpPlotGenes(kp, data=genes.data, r0=0, r1=0.1, gene.name.cex = pp$fontSize, 
            transcript.name.position = "bottom", mark.height = .3)

# SCLC tracks
total.tracks <- length(s.tracks)+length(h.tracks)
out.at <- autotrack(1:length(s.tracks), total.tracks, margin = 0.3, r0=0.23)
kpAddLabels(kp, labels = "SCLC", r0 = out.at$r0, r1=out.at$r1, cex=pp$fontSize,
            srt=90, pos=1, label.margin = 0.14)
for(i in seq_len(length(s.tracks))) {
  bigwig.file <- data.browser$bw_path[i]
  at <- autotrack(i, length(s.tracks), r0=out.at$r0, r1=out.at$r1, margin = 0.1)
  kp <- kpPlotBigWig(kp, data=bigwig.file, ymin = pp$ymin, ymax="visible.region",
                     r0=at$r0, r1=at$r1, col = data.browser$color_code[i])
  # computed.ymax <- ceiling(kp$latest.plot$computed.values$ymax)
  # kpAxis(kp, ymin=0, ymax=computed.ymax, tick.pos = computed.ymax,
  #        r0=at$r0, r1=at$r1, cex=pp$fontSize)
  kpAxis(kp, ymin=1, r0=at$r0, r1=at$r1, cex=pp$fontSize)
  kpAddLabels(kp, labels = data.browser$sample_id[i], r0=at$r0, r1=at$r1, 
              cex=pp$fontSize, label.margin = 0)
}
# healthy tracks
out.at <- autotrack((length(s.tracks)+1):total.tracks, total.tracks, margin = 0.3, r0=0.23)
kpAddLabels(kp, labels = "Healthy", r0 = out.at$r0, r1=out.at$r1,
            cex=pp$fontSize, srt=90, pos=1, label.margin = 0.14)
for(i in seq_len(length(h.tracks))) {
  bigwig.file <- data.browser$bw_path[length(s.tracks) + i]
  at <- autotrack(i, length(h.tracks), r0=out.at$r0, r1=out.at$r1, margin = 0.1)
  kp <- kpPlotBigWig(kp, data=bigwig.file, ymin = pp$ymin, ymax = "visible.region",
                     r0=at$r0, r1=at$r1, col = data.browser$color_code[length(s.tracks) + i])
  # computed.ymax <- ceiling(kp$latest.plot$computed.values$ymax)
  # kpAxis(kp, ymin=0, ymax=computed.ymax, tick.pos = computed.ymax,
  #        r0=at$r0, r1=at$r1, cex=pp$fontSize)
  kpAddLabels(kp, labels = data.browser$sample_id[length(s.tracks) + i], r0=at$r0, r1=at$r1, 
              cex=pp$fontSize, label.margin = 0)
}
dev.off()



# gviz option -------------------------------------------------------------
library(Gviz)
from = 39983586
to = 40006109
chr = "chr19"
name = "DLL3"
gm = "hg19"
j = 9
dTrack2 = sapply(1:j, function(i) DataTrack(range = data.browser$bw_path[i], 
                                                             genome = gm, chromosome = chr, 
                                                             type = "h",
                                                             name = data.browser$sample_id[i]))
plotTracks(dTrack2, showTitle = T, from = from, to = to, windowSize = 100, 
           box.legend = F, col = data.browser$color_code[1:j], cex = 8,
           cex.legend = 8, cex.sampleNames = 8,fontsize.legend = 8,
           fill = data.browser$color_code[1:j], ylim = c(0,20))
refGenes = UcscTrack(genome = "hg19", chromosome = chr,
                      track = "xenoRefGene", from = from, to = to,
                      trackType = "GeneRegionTrack", 
                      rstarts = "exonStarts", rends = "exonEnds", 
                      gene = "name",  symbol = "name2", 
                      transcript = "name", strand = "strand",
                      fill = "#8282d2", stacking = "dense", 
                      name = "Other RefSeq")
dTrack2[[length(dTrack2) + 1]] = refGenes
axTrack = GenomeAxisTrack()
dTrack2[[length(dTrack2) + 1]] = axTrack
idxTrack = IdeogramTrack(genome="hg19", chromosome=chr)
dTrack2[[length(dTrack2) + 1]] = idxTrack
# pdf(paste0(figDirPaper, "temp/gviz.pdf"), width = 5, height = 5)
plotTracks(dTrack2, showTitle = T, from = from, to = to, windowSize = 100, 
           box.legend = F, col = data.browser$color_code[1:j], cex = 8,
           cex.legend = 8, cex.sampleNames = 8,fontsize.legend = 8,
           fill = data.browser$color_code[1:j], ylim = c(0,20))
# dev.off()

# trackplot option --------------------------------------------------------
# https://github.com/PoisonAlien/trackplot
source("https://github.com/PoisonAlien/trackplot/blob/master/R/trackplot.R?raw=true")
source("https://github.com/CRG-Barcelona/libbeato.git")
library("bwtool")
#Path to bigWig files
bigWigs = data.browser$bw_path

#Step-1. Extract the siganl for your loci of interst
loci = paste0(gene_ranges$chr[1], ":", gene_ranges$start[1], "-", gene_ranges$end[1])
track_data = track_extract(bigWigs = bigWigs, loci = loci)

#Step-1a (optional). Summarize trcks by condition
track_data = track_summarize(summary_list = track_data, condition = c("A", "B", "B", "C", "D"), stat = "mean")

#Step-2. 
#Basic Plot 
track_plot(summary_list = track_data)

#With gene models (by default autoamtically queries UCSC genome browser for hg19 transcripts)
track_plot(summary_list = track_data, draw_gene_track = TRUE, build = "hg38")

#With GTF file as source for gene models
track_plot(summary_list = track_data, draw_gene_track = TRUE, gene_model = "hg38_refseq.gtf.gz", isGTF = TRUE)

#Heighlight regions of interest

markregions = data.frame(
  chr = c("chr3", "chr3"),
  start = c(187743255, 187735888),
  end = c(187747473, 187736777),
  name = c("Promoter-1", "Promoter-2")
)

track_plot(
  summary_list = track_data,
  draw_gene_track = TRUE,
  show_ideogram = TRUE,
  build = "hg38",
  regions = markregions
)
