## windows that differentiate between subtypes
ca = data.frame(row.names = metadata$Sample_id, metadata[,c("POU2F3", "ASCL1", "NEUROD1", "YAP1")], 
                sclc_score = estimated.tumor[metadata$Sample_id, "SCLC.n"])
ca[] = lapply(ca, as.numeric)

# pou2f3
data.pou2f3 = data.frame(rna = t(rna_data["POU2F3", rna.samples.chip.passQC]), 
                         chip = log2(1+colSums(win_data_all[pou2f3.wins, rna.samples.chip.passQC])), 
                         tumor = estimated.tumor[rna.samples.chip.passQC, "SCLC.n"], 
                         high = factor(1*(rna.samples.chip.passQC %in% high.ctDNA.samples.wRNA)))
reg = (lm(POU2F3 ~ chip*tumor, data = data.pou2f3))
summary(reg)
qplot(reg$fitted.values, data.pou2f3$POU2F3)
ggplot(data.pou2f3, aes(POU2F3, chip, size = tumor, color = high)) + geom_point()


##
hist(colSums(win_data_all[pou2f3.wins,]), breaks = "fd")
high.pou2f3 = names(which(colSums(win_data_all[pou2f3.wins,rna.samples.chip.passQC]) > 30))

# sanity check
high.pou2f3 %in% names(sort(rna_data["POU2F3", ], decreasing = T))[1:3]
low.pou2f3 = names(which(colSums(win_data_all[pou2f3.wins,high.ctDNA.samples]) < 5))
pou2f3.sig = which(rowMins(as.matrix(win_data_all[,high.pou2f3])) / rowMaxs(as.matrix(win_data_all[,low.pou2f3])) > 2 &
                     rowMins(as.matrix(win_data_all[,high.pou2f3])) > 10)
pou2f3.sig = setdiff(pou2f3.sig, pou2f3.wins)

data.pou2f3 = data.frame(gene = colSums(win_data_all[pou2f3.wins, s.samples]), 
                         sig = colSums(win_data_all[pou2f3.sig,s.samples]))
ggplot(data.pou2f3, aes(sig, gene)) + geom_point()
data.pou2f3[order(data.pou2f3$sig, decreasing = T)[3],]
d = log2(1+win_data_all[pou2f3.sig,s.samples])
co = hclust(dist(colSums(d)))$order
pheatmap(d[,co], show_colnames = F, show_rownames = T, cluster_cols = F, annotation_col = ca, treeheight_row = 0)#, 
         filename = paste0(figDirPaper, "figureS3/pou2f3.chip_high.png"), width = 8, height = 5, fontsize = 4)
d = sweep(d, 1, rowMeans(d), "-")
co = hclust(dist(colSums(d)))$order
pheatmap(d[,co], show_colnames = F, cluster_cols = F, annotation_col = ca, breaks = cfm.breaks, color = cfm.color)#,
         treeheight_row = 0, filename = paste0(figDirPaper, "figureS3/pou2f3.high.png"), width = 8, height = 5, fontsize = 4)

TSS.windows[pou2f3.sig[is.na(TSS.windows$name[pou2f3.sig])],]
pou2f3.genes = unique(TSS.windows$name[pou2f3.sig])
pou2f3.genes = pou2f3.genes[!is.na(pou2f3.genes)]
pou2f3.genes = pou2f3.genes[grep("^\\.", pou2f3.genes, invert = T)]
pou2f3.genes = c(pou2f3.genes, "GTF3C5","CELF4" )
pou2f3.genes = pou2f3.genes[pou2f3.genes %in% rownames(rna_data)]
d = rna_data[pou2f3.genes,]
pheatmap(d, show_colnames = F,annotation_col = ca)#, 
         width = 8, height = 5, filename = paste0(figDirPaper, "figureS3/pou2f3.high.rna.png"))
d = sweep(d, 1, rowMeans(d), "-")
pheatmap(d, show_colnames = F, annotation_col = ca, annotation_colors  = group.colors.heatmap, breaks = cfm.breaks, 
         color = cfm.color, filename = paste0(figDirPaper, "figureS3/pou2f3.high.rna.cfm.png"), 
         width = 12, height = 7)
d = rna_data[pou2f3.genes[rowMaxs(as.matrix(rna_data[pou2f3.genes,])) > 3],]; 
d = sweep(d, 1, rowMeans(d), "-")
pheatmap(d, show_colnames = F, annotation_col = ca, annotation_colors  = group.colors.heatmap, breaks = cfm.breaks, 
         color = cfm.color, filename = paste0(figDirPaper, "figureS3/pou2f3.high.rna.2.png"), 
         width = 12, height = 7)


# same for ascl1
# high.ascl1 = high.ctDNA.samples[which(chip_data_all["ASCL1",high.ctDNA.samples] > 90)]
# low.ascl1 = high.ctDNA.samples[which(chip_data_all["ASCL1",high.ctDNA.samples] < 0)]
high.ascl1 = s.samples[which(colSums(win_data_all[ascl1.promoter.ind,s.samples]) > 15)]
low.ascl1 = s.samples[which(colSums(win_data_all[ascl1.promoter.ind,s.samples]) < 1)]
high.ascl1 = colnames(rna_data)[which(rna_data["ASCL1",] > 4)]; high.ascl1 = high.ascl1[estimated.tumor[high.ascl1, "SCLC"] > .5]
low.ascl1 = colnames(rna_data)[which(rna_data["ASCL1",] < 2)]; low.ascl1 = low.ascl1[estimated.tumor[low.ascl1, "SCLC"] > .5]
ascl1.cor = which(rowMins(as.matrix(win_data_all[,high.ascl1])) / rowMedians(as.matrix(win_data_all[,low.ascl1])) > 1.5 &
                    rowMins(as.matrix(win_data_all[,high.ascl1])) > 10 & rowMeans(win_data_all[,h.samples]) < 3)
ascl1.cor = setdiff(ascl1.cor, ascl1.promoter.ind)
d = log2(1+win_data_all[ascl1.cor,c(s.samples)])
pheatmap(d[,high.ctDNA.samples], show_colnames = F, annotation_col = ca,
         filename = paste0(figDirPaper, "figureS3/ascl1.high.chip.raw.png"), width = 8, height = 5)
d = sweep(d, 1, rowMeans(d), "-")
pheatmap(d[,high.ctDNA.samples], show_colnames = F, annotation_col = ca,
         filename = paste0(figDirPaper, "figureS3/ascl1.high.chip.png"), width = 8, height = 5)
write.table(TSS.windows$name[ascl1.cor], paste0(figDirPaper, "ascl1.targets.csv"), row.names = F)
ascl1.genes = unique(TSS.windows$name[ascl1.cor])
ascl1.genes = ascl1.genes[!is.na(ascl1.genes)]
ascl1.genes = ascl1.genes[grep("^\\.", ascl1.genes, invert = T)]
ascl1.genes = ascl1.genes[ascl1.genes %in% rownames(rna_data)]
ascl1.genes = setdiff(ascl1.genes, "ASCL1")

d = rna_data[ascl1.genes,]
co = order(rna_data["ASCL1",])
pheatmap(d, show_colnames = F, annotation_col = ca, width = 10, height = 7, #breaks = cfm.breaks, color = cfm.color,
         filename = paste0(figDirPaper, "figureS3/ascl1.high.raw.png"))
# d = rna_data[ascl1.genes[rowMedians(as.matrix(rna_data[ascl1.genes,])) > 2],]
d = sweep(d, 1, rowMeans(d), "-")
pheatmap(d, show_colnames = F, annotation_col = ca, width = 10, height = 7, breaks = cfm.breaks, color = cfm.color,
         filename = paste0(figDirPaper, "figureS3/ascl1.high.rna.png"))

# same for neurod1
# high.neurod1 = colnames(rna_data)[which(rna_data["NEUROD1",] > 4)]; high.neurod1 = high.neurod1[estimated.tumor[high.neurod1,"SCLC"] > .5]
# low.neurod1 = colnames(rna_data)[which(rna_data["NEUROD1",] < 1)]; low.neurod1 = low.neurod1[estimated.tumor[low.neurod1,"SCLC"] > .5]
high.neurod1 = high.ctDNA.samples[which(win_data_all[neurod1.promoter.ind,high.ctDNA.samples] > 13)]
low.neurod1 = high.ctDNA.samples[which(win_data_all[neurod1.promoter.ind,high.ctDNA.samples] < 1)]
neurod1.cor = which(rowMins(as.matrix(win_data_all[,high.neurod1])) / rowQuantiles(as.matrix(win_data_all[,low.neurod1]), probs = .9) > 2 &
                      rowMins(as.matrix(win_data_all[,high.neurod1])) > 5 & rowMeans(win_data_all[,h.samples]) < 3 &
                      rowMeans(win_data_all[,h.samples]) < 3)
neurod1.cor = setdiff(neurod1.cor, neurod1.promoter.ind)
d = log2(1+win_data_all[neurod1.cor,c(s.samples)])
d = sweep(d, 1, rowMeans(d), "-")
pheatmap(d, show_colnames = F, show_rownames = F, annotation_col = ca,
         filename = paste0(figDirPaper, "figureS3/neurod1.high.png"), width = 8, height = 5)
write.table(TSS.windows$name[neurod1.cor], paste0(figDirPaper, "neurod1.targets.csv"), row.names = F)

neurod1.genes = unique(TSS.windows$name[neurod1.cor])
neurod1.genes = neurod1.genes[!is.na(neurod1.genes)]
neurod1.genes = neurod1.genes[grep("^\\.", neurod1.genes, invert = T)]
neurod1.genes = neurod1.genes[neurod1.genes %in% rownames(rna_data)]
neurod1.genes = setdiff(neurod1.genes, "NEUROD1")
d = rna_data[neurod1.genes,]
d = sweep(d, 1, rowMeans(d), "-")
pheatmap(d, show_colnames = F, annotation_col = ca, breaks = cfm.breaks, color = cfm.color,
         filename = paste0(figDirPaper, "figureS3/neurod1.high.rna.png"), width = 12, height = 9)

# combine ascl1 and neurod1
pheatmap(rna_data[c("NEUROD1", "ASCL1", "YAP1", "POU2F3"),])
high.ne = colnames(rna_data)[which(rna_data["NEUROD1",] > 4 | rna_data["ASCL1",] > 4)]
high.ne = high.ne[high.ne %in% high.ctDNA.samples]
low.ne = colnames(rna_data)[which(rna_data["NEUROD1",] < 3 & rna_data["ASCL1",] < 3)]
low.ne = low.ne[low.ne %in% high.ctDNA.samples]
ne.cor = which(rowMedians(as.matrix(win_data_all[,high.ne])) / rowMaxs(as.matrix(win_data_all[,low.ne])) > 2 &
                 rowMeans(as.matrix(win_data_all[,high.ne])) > 5 & rowMeans(win_data_all[,h.samples]) < 1)
d = log2(1+win_data_all[ne.cor,c(s.samples)])
d = sweep(d, 1, rowMeans(d), "-")
pheatmap(d, show_colnames = F, show_rownames = F, annotation_col = ca, color = cfm.color, breaks = cfm.breaks,
         filename = paste0(figDirPaper, "figureS3/ne.high.png"), width = 8, height = 5)
ne.genes = unique(TSS.windows$name[ne.cor])
ne.genes = ne.genes[grep("^\\.", ne.genes, invert = T)]
ne.genes = ne.genes[ne.genes %in% rownames(rna_data)]
ne.genes = ne.genes[!ne.genes %in% c("ASCL1", "NEUROD1")]
d = rna_data[ne.genes,]
pheatmap(d, show_colnames = F, annotation_col = ca, color = cfm.color, breaks = cfm.breaks,
         filename = paste0(figDirPaper, "figureS3/ne.high.png"), width = 8, height = 5)

# same for yap1
high.yap1 = metadata$Sample_id...1[which(as.numeric(metadata$YAP1) > 3)]; high.yap1 = high.yap1[estimated.tumor[high.yap1, "SCLC"] > .5]
high.yap1 = high.yap1[!is.na(high.yap1)]
low.yap1 = metadata$Sample_id...1[which(as.numeric(metadata$YAP1) < 1.5)]; low.yap1 = low.yap1[estimated.tumor[low.yap1, "SCLC"] > .5]
low.yap1 = low.yap1[!is.na(low.yap1)]
yap1.cor = which(rowMins(as.matrix(win_data_all[,high.yap1])) / rowMedians(as.matrix(win_data_all[,low.yap1])) > 2 &
                   rowMins(as.matrix(win_data_all[,high.yap1])) > 5)
d = log2(1+win_data_all[yap1.cor,c(s.samples)])
pheatmap(d, show_colnames = F, show_rownames = F, annotation_col = estimated.tumor)
#filename = paste0(figDirPaper, "figureS3/yap1.high.png"), width = 8, height = 5)
yap1.genes = unique(TSS.windows$name[yap1.cor])
yap1.genes = yap1.genes[!is.na(yap1.genes)]
yap1.genes = yap1.genes[grep("^\\.", yap1.genes, invert = T)]
yap1.genes = yap1.genes[yap1.genes %in% rownames(rna_data)]
yap1.genes = setdiff(yap1.genes, "YAP1")
d = rna_data[yap1.genes,]
d = sweep(d, 1, rowMeans(d), "-")
pheatmap(d, show_colnames = F, annotation_col = ca, breaks = cfm.breaks, color = cfm.color,
         filename = paste0(figDirPaper, "figureS3/yap1.high.rna.png"), width = 12, height = 9)

s = c(s.samples, h.samples)
s = rna.samples.matched.time[rna.samples.matched.time %in% s.samples]
data.cor.tf = data.frame(pou2f3.sig = colSums(win_data_all[pou2f3.sig,s]), ne.sig = colSums(win_data_all[ne.cor, s]),
                         ascl1.sig = colSums(win_data_all[ascl1.cor,s]), neurod1.sig = colSums(win_data_all[neurod1.cor,s]), 
                         yap1.sig = colSums(win_data_all[yap1.cor, s]), pou2f3.rna = as.numeric(metadata$POU2F3[match(s, metadata$Sample_id...1)]),
                         ascl1.rna = as.numeric(metadata$ASCL1[match(s, metadata$Sample_id...1)]),
                         neurod1.rna = as.numeric(metadata$NEUROD1[match(s, metadata$Sample_id...1)]),
                         yap1.rna = as.numeric(metadata$YAP1[match(s, metadata$Sample_id...1)]),
                         pou2f3.chip = colSums(win_data_all[pou2f3.genebody.ind, s]),
                         ascl1.chip = colSums(win_data_all[ascl1.promoter.ind, s]),
                         neurod1.chip = as.vector(t(win_data_all[neurod1.promoter.ind, s])),
                         yap1.chip = as.vector(t(chip_data_all["YAP1", s])))
data.cor.tf$all.rna = data.cor.tf$ascl1.rna + data.cor.tf$neurod1.rna + data.cor.tf$pou2f3.rna + data.cor.tf$yap1.rna
data.cor.tf$tumor = estimated.tumor$SCLC[match(rownames(data.cor.tf), rownames(estimated.tumor))]
data.cor.tf$subtype = metadata$NAPY_allcluster_expression[match(rownames(data.cor.tf), metadata$Sample_id...1)]
data.cor.tf[h.samples,"subtype"] = "healthy"

intersect(pou2f3.genes, lung_celltypes_genes$gene[lung_celltypes_genes$cell_type == "Ciliated"])
intersect(ascl1.genes, lung_celltypes_genes$gene[lung_celltypes_genes$cell_type == "Neuroendocrine"])

## TF signature vs. tumor load
data.cor.tf = data.cor.tf[data.cor.tf$tumor > 0,]
tf.sig = data.cor.tf$pou2f3.sig; tf.rna = data.cor.tf$pou2f3.rna; tf.name = "POU2F3"; 
tf.sig = data.cor.tf$ascl1.sig; tf.rna = data.cor.tf$ascl1.rna; tf.name = "ASCL1"; 
tf.sig = data.cor.tf$neurod1.sig; tf.rna = data.cor.tf$neurod1.rna; tf.name = "NEUROD1"; 
# tf.sig = data.cor.tf$neurod1.sig + data.cor.tf$ascl1.sig; tf.rna = data.cor.tf$neurod1.rna + data.cor.tf$ascl1.rna; tf.name = "NE"; 

p = ggplot(data.cor.tf, aes(tumor, tf.sig, color = tf.rna/all.rna)) + geom_point() 
p = p + scale_color_gradient(low = "#fdbb84", high = "#990000", na.value = "#d9d9d9")
p = p + labs(x = "tumor-load", y = "reads/signature", color = "relative RNA", title = tf.name) 
p = p + geom_vline(xintercept = .75, linetype = "dashed")
p = p + geom_smooth(method = "lm")
ggsave(paste0(figDirPaper, "figureS3/", tf.name, "_sig_vs_tumor.png"), width = 120, height = 80, units = "mm")

tf.fit = lm(tf.sig ~ tumor, data = data.cor.tf)$coefficients
p = ggplot(data.cor.tf, aes(tf.rna/all.rna, log2(tf.sig / (tf.fit[1] + tumor * tf.fit[2])), color = tumor)) 
p = p + geom_point(size = 1, shape = 16) 
p = p + scale_color_gradient(low = "#9ecae1", high = "#084594", na.value = "#d9d9d9")
p = p + labs(x = "relative RNA", y = "log2 (observed / expected)", alpha = "tumor load", title = tf.name)
p = p + geom_smooth(method = "lm") + stat_cor(size = base_size/.pt)
ggsave(paste0(figDirPaper, "figure4/", tf.name, "_sig_vs_exp_sig.pdf"), width = 60, height = 60, units = "mm")


# TSS TF vs. TF signature
p = ggplot(data.cor.tf, aes(ascl1.chip, ascl1.sig, size = tumor, color = ascl1.rna)) + geom_point() + scale_colour_viridis_c(option = "plasma")
p = p + labs(x = "ASCL1 (1+TSS reads)", y = "reads/signature", title = "ASCL1")
ggsave(paste0(figDirPaper, "figureS3/ascl1_vs_signature.png"), width = 8, height = 6)
p = ggplot(data.cor.tf, aes(pou2f3.chip, pou2f3.sig, size = tumor, color = pou2f3.rna)) + geom_point() + scale_colour_viridis_c(option = "plasma")
p = p + labs(x = "POU2F3 (1+TSS reads)", y = "reads/signature", title = "POU2F3")
ggsave(paste0(figDirPaper, "figureS3/pou2f3_vs_signature.png"), width = 8, height = 6)
p = ggplot(data.cor.tf, aes(neurod1.chip, neurod1.sig, size = tumor, color = neurod1.rna)) + geom_point() + scale_colour_viridis_c(option = "plasma")
p = p + labs(x = "NEUROD1 (1+TSS reads)", y = "reads/signature", title = "NEUROD1")
ggsave(paste0(figDirPaper, "figureS3/neurod1_vs_signature.png"), width = 8, height = 6)
p = ggplot(data.cor.tf, aes(yap1.chip, yap1.sig, size = tumor, color = yap1.rna)) + geom_point() + scale_colour_viridis_c(option = "plasma")
p = p + labs(x = "YAP1 (1+TSS reads)", y = "reads/signature", title = "YAP1")
ggsave(paste0(figDirPaper, "figureS3/yap1_vs_signature.png"), width = 8, height = 6)


