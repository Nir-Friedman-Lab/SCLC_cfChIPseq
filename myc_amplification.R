# amplification analysis
# prepare data ------------------------------------------------------------
amplifications = read.csv(paste0(baseDir, "cfChIP-paper/data/sclc_final_1_amps.csv"))
amplifications = amplifications[amplifications$group_name %in% s.samples,]
myc.wins = which(TSS.windows$name == "MYC;PVT1;MIR1204")
amplifications$start_win = as.numeric(str_split_fixed(amplifications$wins, ":", 2)[,1])
amplifications$end_win = as.numeric(str_split_fixed(amplifications$wins, ":", 2)[,2])
amplifications$start = as.numeric(str_split_fixed(amplifications$genomic_loc, "-", 2)[,1])
amplifications$end = as.numeric(str_split_fixed(amplifications$genomic_loc, "-", 2)[,2])
# remove extra 0 from chr name 
amplifications$chr[grep("[:digit:,0]", amplifications$chr)] = 
  sub("0", "", amplifications$chr[grep("[:digit:,0]", amplifications$chr)])
# inspect genes who are konw to have amplification in SCLC ----------------
# list based on https://www.nature.com/articles/s12276-019-0349-5/tables/1
amp.genes = c("MYC", "MYCL1", "MYCN", "CCNE1", "MET", "FGFR1", "IRS2", "NFIB", "SOX2", "SOX4", "INSM1")
sort(table(amplifications$gene))
g = "INSM1"; s = c(s.samples, h.samples)
df = data.frame(chip = t(chip_data_all[g,s]), col.annotation[s,]); names(df) = c("chip", "group")
high = quantile(df$chip,.9)
df$samp = substr(rownames(df),4,13)
p = ggplot(df, aes(group, chip, group = group, fill = group, label = samp )) + geom_boxplot() + geom_point(size = 2)
# p = p + scale_y_continuous(trans = "log10")
p = p + geom_text_repel(data = df[df$chip > high,], size = 2)
p
TSS.windows[grep(paste0(g,"$"), TSS.windows$name),]
amplifications[grep(g, amplifications$gene),c(1:5,7)]
amplifications[grep("509", amplifications$group_name),c(1:5,7)]
i = grep(g, amplifications$gene)
amplifications[i,1:7]
amp.genes %in% unique(TSS.windows$name[which(start(TSS.windows) > min(amplifications$start[i]) 
                  & end(TSS.windows) < max(amplifications$end[i])
                  & seqnames(TSS.windows) == unique(amplifications$chr[i]))])



#CCNE1, snap47? SSRP1 nicely seen in 37-478 which is low.qc.  in 37-677 there are high reads 


# heatmap of amplifications in samples ------------------------------------
amp.genes = c("MYC", "MYCL1", "MYCN", "FGFR1", "NFIB", "ZNF568", "IGFBP7")
data.amp = data.frame(matrix(data = 1, nrow = length(amp.genes), ncol = length(s.samples)))
rownames(data.amp) = amp.genes
colnames(data.amp) = s.samples
for (i in 1:nrow(amplifications)){
  print(i)
  S = amplifications$group_name[i]
  C = amplifications$cnv[i]
  G = unique(TSS.windows$name[which(start(TSS.windows) > min(amplifications$start[i]) 
                                    & end(TSS.windows) < max(amplifications$end[i])
                                    & seqnames(TSS.windows) == unique(amplifications$chr[i]))])
  G = unlist(strsplit(G, split = ";"))
  GG = amp.genes[amp.genes %in% G] 
  if (!isEmpty(GG)) {
    data.amp[GG,S] = C
    print(GG)
  }
}
data.amp = data.amp[,colSums(data.amp) > nrow(data.amp)]
ro = hclust(dist(data.amp))$order
co = hclust(dist(t(data.amp)))$order
data.amp = data.amp[ro,co]
data.amp[data.amp > 8] = 8
ca = metadata[colnames(data.amp),"Responder_NonResponder", drop = F]
pheatmap(data.amp, annotation_col = ca, width = 8, height = 5,
         filename = paste0(figDirPaper, "figure5/amplification_response.png"))
data.amp$gene = rownames(data.amp)
data.amp.m = reshape2::melt(data.amp, id.vars = "gene", variable.name = "samp")
data.amp.m$value[data.amp.m$value > 5] = 5
data.amp.m$samp = sub(".K4me3*", "", data.amp.m$samp)
data.amp.m$samp = sub("_.*", "", data.amp.m$samp)
data.amp.m$samp = sub("PL.", "SCLC", data.amp.m$samp)
data.amp.m$samp = gsub("\\.", "-", data.amp.m$samp)
data.amp.m$samp = sub("-milli.*", "", data.amp.m$samp)
p = heatmap(data.amp.m, "gene", "samp", "value", chigh = "red", clow = "black", 
        breaks = c(1,3,5), leg = c("X1", "X3", ">5"), lim = c(1,5), 
        xlab = "amplified region", leglab = "duplication") 
p = p + scale_y_discrete(position = "right")
ggsave(paste0(figDirPaper, "figure5/amplification_table.pdf"), p, width = 80, height = 90, units = "mm")
ggsave(paste0(figDirPaper, "figure5/amplification_table.png"), p, width = 80, height = 90, units = "mm")

s = c(s.samples, h.samples)
for (g in amp.genes) {
  print(g)
  data.gene = data.frame(chip = t(chip_data_all[g,s]), col.annotation[s,]); 
  names(data.gene) = c("chip", "group")
  data.gene$samp = rownames(data.gene)
  high = quantile(data.gene$chip,.98)
  data.gene$samp = sub(".K4me3*", "", data.gene$samp)
  data.gene$samp = sub("_.*", "", data.gene$samp)
  data.gene$samp = sub("PL.", "SCLC", data.gene$samp)
  data.gene$samp = gsub("\\.", "-", data.gene$samp)
  data.gene$samp = sub("-milli.*", "", data.gene$samp)
  p = boxplotWpoints(data.gene, "group", "chip", "group", "samp", title = g, ylab = "count/TSS")
  p = p + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  # p = ggplot(data.gene, aes(group, chip, group = group, fill = group, label = samp )) + geom_boxplot() + geom_point(size = 2)
  p = p + geom_text_repel(data = data.gene[data.gene$chip > high & rownames(data.gene) %in% colnames(data.amp),], size = base_size/.pt)
  ggsave(paste0(figDirPaper, "figure5/", g, "_boxplot.pdf"), p, width = 40, height = 50, units = "mm")
  p = qplot(t(rna_data[g,high.ctDNA.samples.wRNA]), t(log2(1+chip_data_all[g,high.ctDNA.samples.wRNA])), 
        xlab = "RNA", ylab = "ChIP", main = g)
  ggsave(paste0(figDirPaper, "figure5/", g, "chip_rna_scatter.png"), p, width = 40, height = 40, units = "mm")
}


# myc ---------------------------------------------------------------------
myc.samples = amplifications$group_name[amplifications$start_win < max(myc.wins) & amplifications$end_win > max(myc.wins)]
s = c(s.samples, h.samples)
data.myc = data.frame(samp = s, tumor = estimated.tumor[s, "SCLC.n"], 
                      chip = colSums(chip_data_all["MYC",s]), 
                      rna = rep(NA, length(s)), 
                      group = col.annotation[s,],
                      response = metadata[s, "Responder_NonResponder"],
                      timepoint = metadata[s, "Timepoint"], 
                      amp = rep(F, length(s)))
data.myc[myc.samples, "amp"] = T
data.myc[rna.samples.chip.passQC,"rna"] = t(rna_data["MYC",rna.samples.chip.passQC])
tail(data.myc)
# boxplot healhy SCLC - zoom in on high samples  -------------------------------------------------------
p = boxplotWpoints(data.myc, x = "group", y = "chip", fill = "group", size = 1, ylab = "MYC counts")
p = p + facet_zoom(ylim = c(80,250), zoom.size = 1)
# p = p + scale_y_break(c(250, 1000), scales = 1.5) #+ scale_y_break(c(2200, 4000), scales = 1.5)
ggsave(paste0(figDirPaper, "figure5/myc_boxplot.pdf"), p, width = 60, height = 60 , units = "mm")

d = data.myc[data.myc$timepoint == "Pre-treatment" & data.myc$response != "NA",]
p = boxplotWpoints(d, x = "response", y = "chip", fill = "response", ylab = "MYC counts")
p = p + scale_y_continuous(trans = "log2")
ggsave(paste0(figDirPaper, "figure5/myc_reponse.pdf"), width = 55, height = 55, units = "mm")
# myc vs sclc score -------------------------------------------------------
d = data.myc[data.myc$timepoint == "Pre-treatment" & !is.na(data.myc$response)
             & data.myc$response != "NA",]
p = scatter.plot(d, "tumor", "chip", color = "response", guieds = T, size = 1, cor = F, 
                 xlab = "SCLC score", ylab = "MYC counts", clab = "")
p = p + scale_y_continuous(trans = "log10")
p = p + scale_x_continuous(breaks = sclc.breaks, labels = sclc.lab)
p = p + theme(legend.position = c(.2,.95), legend.key.size = unit(2, 'mm'))
ggsave(paste0(figDirPaper, "figureS5/myc_vs_tumor.pdf"), p, width = 55, height = 55, units = "mm")

# rna vs. chip MYC --------------------------------------------------------
p = scatter.plot(data.myc, "rna", "chip", xlab = "RNA (RPKM)", ylab = "ChIP (reads/TSS)")
p = p + scale_y_continuous(trans = "log2")
ggsave(paste0(figDirPaper, "figure5/myc_rna_vs_chip.pdf"), p , width = 55, height = 55, units = "mm")
p = ggplot(data.myc, aes(rna, chip, label = substr(samp,4,13))) 
p = p + scale_y_continuous(trans = "log2")
p = p + geom_point(size = 1) + geom_text_repel(size = 1)
ggsave(paste0(figDirPaper, "figure5/myc_rna_vs_chip_named.pdf"), p , width = 55, height = 55, units = "mm")
estimated.tumor[s.samples[grep("33-489", s.samples)],]

# amplification summary ---------------------------------------------------
grep("MYC", amplifications$gene)
myc.samples


