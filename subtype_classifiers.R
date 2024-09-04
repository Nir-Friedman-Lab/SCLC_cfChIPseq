# subtype annotations -----------------------------------------------------
liver.genes = read.table(paste0(figDirPaper, "liver_specific_genes.txt"))$V1
subtypes = c("ASCL1", "NEUROD1", "POU2F3", "YAP1")
annotation.subtype = metadata[,c("NE_SCORE", "NE_status", 
                                 "NAPY_allcluster_expression", "SCLC_score", 
                                 subtypes)]
names(annotation.subtype) = c("NE-score", "NE-status", "subtype", "SCLC_score", subtypes)
annotation.subtype[colnames(rna_data), "liver"] = colSums(rna_data[liver.genes,])
annotation.subtype[,subtypes] = lapply(annotation.subtype[,subtypes], as.numeric)
annotation.subtype = annotation.subtype[!rownames(annotation.subtype) %in% l.samples,] # remove NSCLC RNA samples 

annotation.color = group.colors.heatmap
annotation.color$subtype_cutoff = c(annotation.color$subtype, "MIXED" = "#d9d9d9")
annotation.color$max = annotation.color$subtype_cutoff
n = 8
annotation.color$ASCL1 = colorRamp2(seq(0,n,length = n), hcl.colors(n, "Reds", rev = T)) # brewer.pal(9, "Reds")
annotation.color$NEUROD1 = colorRamp2(seq(0,n,length = n), hcl.colors(n, "Oranges", rev = T)) # brewer.pal(9, "Oranges")
annotation.color$POU2F3 = colorRamp2(seq(0,n,length = n), hcl.colors(n, "Blues", rev = T)) #brewer.pal(9, "Blues")
annotation.color$YAP1 = colorRamp2(seq(0,n,length = n), hcl.colors(n, "Greens", rev = T)) #brewer.pal(9, "Greens")
annotation.color$SCLC_score = colorRamp2(seq(0,1,length = 10), hcl.colors(10, "Grays", rev = T))
annotation.color$treatment = c("pre-" = "#EE0000FF", "post-" = "#3B4992FF")
# label samples (based on RNA) --------------------------------------------

# new labeling based on RNA levels of TF above cutoff.  
sapply(subtypes, function(i) quantile(as.numeric(rna_data[i,]), seq(0,1,0.05)))
tf.cutoff = 3
rna.class = data.frame(1*(annotation.subtype[labeled.samples,subtypes] > tf.cutoff))
rna.class = rna.class[!is.na(rowSums(rna.class)),]
annotation.subtype$subtype_cutoff = NA
annotation.subtype[rownames(rna.class)[rna.class$ASCL1 & rowSums(rna.class[,1:4]) == 1], "subtype_cutoff"] = "ASCL1"
annotation.subtype[rownames(rna.class)[rna.class$NEUROD1 & rowSums(rna.class[,1:4]) == 1], "subtype_cutoff"] = "NEUROD1"
annotation.subtype[rownames(rna.class)[rna.class$ASCL1 & rna.class$NEUROD1], "subtype_cutoff"] = "ASCL1-NEUROD1"
annotation.subtype[rownames(rna.class)[rna.class$POU2F3 & rowSums(rna.class[,1:4]) == 1], "subtype_cutoff"] = "POU2F3"
annotation.subtype[rownames(rna.class)[rna.class$YAP1 & rowSums(rna.class[,1:4]) == 1], "subtype_cutoff"] = "YAP1"
annotation.subtype[rownames(rna.class)[rowSums(rna.class[,1:4]) == 0], "subtype_cutoff"] = NA
annotation.subtype[rownames(rna.class)[rowSums(rna.class[,1:4]) > 1 & !(rna.class$ASCL1 & rna.class$NEUROD1)], "subtype_cutoff"] = "MIXED"
annotation.subtype[grep("037", rownames(annotation.subtype)),"subtype_cutoff"] = "POU2F3" # it seems in all analysis that these
# samples are POU2F3.  also in RNA, they have high relative POU2F3. 
table(annotation.subtype$subtype_cutoff)

# labels based on maximum of TF
annotation.subtype$max = NA
annotation.subtype[labeled.samples, "max"] = apply(annotation.subtype[labeled.samples,subtypes], MARGIN = 1, which.max)
annotation.subtype[names(which(sapply(labeled.samples, function(s)
  max(annotation.subtype[s,subtypes]) < 2))), "max"] = NA
annotation.subtype$max = factor(annotation.subtype$max, levels = 1:4, labels = subtypes)
# 37 has a strong POU2F3 signal
annotation.subtype[grep("SCLC0037", rownames(annotation.subtype)), "max"] = "POU2F3"
sort(table(paste(annotation.subtype[,"subtype"], annotation.subtype[,"subtype_cutoff"], 
                 annotation.subtype[,"max"], sep = "_")))

annotation.par = c("subtype_cutoff", "subtype", "max", "SCLC_score", subtypes)
if (F) {
  pheatmap::pheatmap(t(rna.class), annotation_col = annotation.subtype[,annotation.par], 
                     # annotation_colors = annotation.color, 
                     fontsize = 6, 
                     treeheight_col = 0, main = "log2(expression) > 3", 
                     legend = F, fontsize_col = 5, cluster_rows = F,
                     # filename = paste0(figDirPaper, "subtypes/RNA_new_subtype.png"), 
                     width = 8, height = 5)
  
  pheatmap::pheatmap(t(annotation.subtype[labeled.samples,subtypes]), 
                     fontsize = 6, fontsize_col = 4,
                     annotation_col = annotation.subtype[,annotation.par],
                     # annotation_colors = annotation.color, 
                     # filename = paste0(figDirPaper, "subtypes/RNA.png"), 
                     width = 8, height = 5, treeheight_col = 3)
}

# there is a good agreement between original labeling and labeling by cutoff.  use cutoff to label samples
ascl1.samples = rownames(annotation.subtype)[which(annotation.subtype$subtype_cutoff == "ASCL1")]
neurod1.samples = rownames(annotation.subtype)[which(annotation.subtype$subtype_cutoff == "NEUROD1")]
ascl1.neurod1.samples = rownames(annotation.subtype)[which(annotation.subtype$subtype_cutoff == "ASCL1-NEUROD1")]
pou2f3.samples = rownames(annotation.subtype)[which(annotation.subtype$subtype_cutoff == "POU2F3")]
pou2f3.samples.low.qc = "SCLC0037-478_NA_H3K4me3-39" # low qc.  strong clear POU2F3 signal - (POU2F3 RNA is low)
yap1.samples = rownames(annotation.subtype)[which(annotation.subtype$subtype_cutoff == "YAP1")]

# cluster and learn classifier from RNA data ------------------------------
rna_roadmap[c("ASCL1", "NEUROD1", "POU2F3", "YAP1"), "Adult_Liver"]
hist(rowMeans(chip_data_all[,h.samples]), breaks = "fd", ylim = c(0,500), xlim = c(0,200))
low.healthy.genes = names(which(rowMeans(chip_data_all[,h.samples]) < 10))
low.healthy.genes = intersect(low.healthy.genes, rownames(rna_data))
length(low.healthy.genes)

rna.pca = prcomp(t(rna_data[low.healthy.genes, ]))
dim(rna.pca$x)
rna.pca = data.frame(PC1 = rna.pca$x[,1], PC2 = rna.pca$x[,2], 
                     subtype = annotation.subtype[colnames(rna_data), "subtype_cutoff"], 
                     liver = colSums(rna_data[liver.genes,]))
p = ggplot(rna.pca, aes(PC1, PC2, color = liver)) + geom_point()
ggsave(paste0(figDirPaper, "subtypes/RNA_pca_liver.png"), p, width = 6, height = 5)

# one of the POU2F3 does not have RNA. need to filter this one out before running the following lines 
ascl1.genes = names(which(rowMeans(rna_data[,ascl1.samples]) > 
                      rowMaxs(as.matrix(rna_data)[,c(neurod1.samples, pou2f3.samples, yap1.samples)]) & 
                      rowMeans(rna_data[,ascl1.samples]) > 3))
neurod1.genes = names(which(rowMeans(rna_data[,neurod1.samples]) > 
                            rowMaxs(as.matrix(rna_data)[,c(ascl1.samples, pou2f3.samples, yap1.samples)]) & 
                            rowMeans(rna_data[,neurod1.samples]) > 3))
pou2f3.genes = names(which(rowMeans(rna_data[,pou2f3.samples]) > 
                              rowMaxs(as.matrix(rna_data)[,c(ascl1.samples, neurod1.samples, yap1.samples)]) & 
                              rowMeans(rna_data[,pou2f3.samples]) > 3))
yap1.genes = names(which(rowMeans(rna_data[,yap1.samples]) > 
                             rowMaxs(as.matrix(rna_data)[,c(ascl1.samples, neurod1.samples, pou2f3.samples)]) & 
                             rowMeans(rna_data[,yap1.samples]) > 3))
class.genes = setdiff(c(ascl1.genes, neurod1.genes, pou2f3.genes, yap1.genes), liver.genes)
s = c(ascl1.samples, ascl1.neurod1.samples, neurod1.samples, pou2f3.samples, yap1.samples)
rna.pca = prcomp(t(rna_data[class.genes, s]))
rna.pca = data.frame(PC1 = rna.pca$x[,1], PC2 = rna.pca$x[,2], 
                     subtype = annotation.subtype[s, "subtype_cutoff"], 
                     liver = colSums(rna_data[liver.genes,s]))
p = ggplot(rna.pca, aes(PC1, PC2, color = liver)) + geom_point()
p = ggplot(rna.pca, aes(PC1, PC2, color = subtype)) + geom_point()
p = p + scale_color_manual(values = annotation.color$subtype_cutoff, limits = force)
ggsave(paste0(figDirPaper, "subtypes/RNA_pca_class_genes.png"), p, width = 6, height = 5)
d = rna_data[class.genes, ]
d = sweep(d, 1, rowMeans(d), "-")
d[d > 3] = 3
d[d < -3] = -3
pheatmap(d, show_rownames = F, show_colnames = T, annotation_col = annotation.subtype, 
         annotation_colors = annotation.color, fontsize_col = 4, width = 8, height = 6,
         filename = paste0(figDirPaper, "subtypes/rna_class_genes.png"))
# 37-677 indeed has high signal in ASCL1 and POU2F3 genes


# subtype signature windows ---------------------------------------------------------
table(annotation.subtype[high.ctDNA.labeled,"subtype_cutoff"])
table(annotation.subtype[med.ctDNA.labeled,"subtype_cutoff"])

if (F) {
  winMaxs = rowMaxs(win_data_all[,c(h.samples, sv.samples)])
  quantile(winMaxs)
  high.windows = which(winMaxs > 10 & healthy.mean < .3)
  
  glm.test = function(X,Class, Score) {
    df = data.frame(
      X = X,
      Class = Class * 1,
      Score = Score,
      C1 = Score * Class,
      C2 = Score * (1-Class)
    )
    # fit1 = lm(X ~ C1+Score, data = df)
    fit1 = lm(X ~ C1, data = df)
    fit2 = lm(X ~ Score, data = df)
    anov = anova(fit1)
    p.val1 = anov$`Pr(>F)`[1]
    p.val0 = anova(fit2)$`Pr(>F)`[1]
    rmse = sqrt(mean(fit1$residuals^2))
    c(coef(fit1)[2], rmse, coef(fit2)[2], p.val1, p.val0)
  }
  groups = list(ascl1 = list(pos = c(ascl1.samples, ascl1.neurod1.samples), 
                          neg = c(neurod1.samples, pou2f3.samples, yap1.samples)), 
                neurod1 = list(pos = c(neurod1.samples, ascl1.neurod1.samples), 
                            neg = c(ascl1.samples, pou2f3.samples, yap1.samples)), 
                pou2f3 = list(pos = c(pou2f3.samples, "SCLC0287-501_NA_H3K4me3-726_06072023-96"), 
                           neg = c(ascl1.samples, ascl1.neurod1.samples, neurod1.samples, 
                                   yap1.samples)), 
                yap1 = list(pos = yap1.samples, neg = c(ascl1.samples, 
                                                     ascl1.neurod1.samples, 
                                                     neurod1.samples, 
                                                     pou2f3.samples)))
  
  groups = list(ascl1 = list(pos = names(which(sample.label["ASCL1",])),  
                             neg = names(which(!sample.label["ASCL1",]))), 
                neurod1 = list(pos = names(which(sample.label["NEUROD1",])),  
                               neg = names(which(!sample.label["NEUROD1",]))), 
                pou2f3 = list(pos = names(which(sample.label["POU2F3",])),  
                              neg = names(which(!sample.label["POU2F3",]))), 
                yap1 = list(pos = names(which(sample.label["YAP1",])),  
                            neg = names(which(!sample.label["YAP1",]))))
  
  
  class.wins = c() # classifier windows
  for (subtype in names(groups)) {
    print(subtype)
    windows = high.windows
    pos = groups[[subtype]]$pos
    neg = groups[[subtype]]$neg
    sample.class = c(rep(T, length(pos)), rep(F, length(neg)))
    sample.score = estimated.tumor[c(pos, neg), "SCLC.n"]
    reg = sapply(windows, function(i) 
      glm.test(win_data_all[i, c(pos, neg)], sample.class, sample.score))
    reg = data.frame(t(reg))
    names(reg) = c("coef1", "rmse", "coef2", "p.val1", "p.val0")
    reg$q.val1 = p.adjust(reg$p.val1, method = "fdr")
    reg$q.val0 = p.adjust(reg$p.val0, method = "fdr")
    reg = data.frame(window = windows, reg, subtype = subtype)
    class.wins = rbind(class.wins, reg)
  }
  
  ggplot(class.wins[class.wins$q.val1 < q.val.cutoff,], aes(x = rmse, color = subtype)) + 
    stat_ecdf() + lims(x = c(0,10))
  rmse.cutoff = 3
  rmse.cutoff.ascl1 = 10
  
  ggplot(class.wins[class.wins$q.val1 < q.val.cutoff,], aes(x = coef1, color = subtype)) + 
    stat_ecdf() + xlim(c(-5,5))
  
  ggplot(class.wins[class.wins$q.val1 < q.val.cutoff,], aes(x = coef2, color = subtype)) + 
    stat_ecdf() + xlim(c(-5,5))
  high.coef = 2
  low.coef = 1
  low.coef.ascl1 = 2.5
  quantile(class.wins$coef2[class.wins$q.val1 < q.val.cutoff & class.wins$subtype == "ascl1"])
  quantile(class.wins$coef2[class.wins$q.val1 < q.val.cutoff & class.wins$subtype == "pou2f3"])
  q.val.cutoff = .01
  which(class.wins$subtype == "ascl1" & abs(class.wins$coef1) < 3 & class.wins$q.val1 < 1e-10)
  i = high.windows[164186]
  ggplot(class.wins[class.wins$subtype == "ascl1",], aes(coef1, -log10(q.val1))) + geom_point()
  ascl1.class.wins = class.wins$window[class.wins$subtype == "ascl1" & 
                                         class.wins$rmse < rmse.cutoff.ascl1 & 
                                         class.wins$q.val1 < q.val.cutoff & 
                                         # class.wins$q.val2 > .1 &
                                         abs(class.wins$coef2) < low.coef.ascl1 &
                                         class.wins$coef1 > high.coef]
  neurod1.class.wins = class.wins$window[class.wins$subtype == "neurod1" &
                                           class.wins$rmse < rmse.cutoff & 
                                           class.wins$q.val1 < q.val.cutoff & 
                                           # class.wins$q.val2 > .1 &
                                           # abs(class.wins$coef2) < low.coef &
                                           class.wins$coef1 > high.coef]
  pou2f3.class.wins = class.wins$window[class.wins$subtype == "pou2f3" & 
                                          class.wins$rmse < rmse.cutoff & 
                                          class.wins$q.val1 < q.val.cutoff & 
                                          # class.wins$q.val2 > .1 &
                                          abs(class.wins$coef2) < low.coef &
                                          class.wins$coef1 > high.coef]
  yap1.class.wins = class.wins$window[class.wins$coef1 == "yap1" & 
                                        class.wins$rmse < rmse.cutoff & 
                                        class.wins$q.val1 < q.val.cutoff & 
                                        # class.wins$q.val2 > .1 &
                                        abs(class.wins$coef2) < low.coef &
                                        class.wins$coef1 > high.coef]
  
  class.genes.chip = c(ascl1.class.wins, neurod1.class.wins, pou2f3.class.wins)
  data.sig.wins = data.frame(signature = c(rep("ASCL1", length(ascl1.class.wins)),
                                           rep("NEUROD1", length(neurod1.class.wins)), 
                                           rep("POU2F3", length(pou2f3.class.wins))),
                             window = class.genes.chip)
}


subtype = "pou2f3"; wins = pou2f3.class.wins
subtype = "neurod1"; wins = neurod1.class.wins
subtype = "ascl1"; wins = ascl1.class.wins
pos = groups[[subtype]]$pos; neg = groups[[subtype]]$neg
s = sv.samples
sample.class = rep("UK", length = length(sv.samples))
sample.class[s %in% pos] = "pos"
sample.class[s %in% neg] = "neg"
sample.sub = annotation.subtype[s, "subtype_cutoff"]
sample.score = estimated.tumor[s, "SCLC.n"]
X = colSums(win_data_all[wins, s])
# X = win_data_all[i, c(pos, neg)]
p = ggplot(data.frame(X = X, score = sample.score, sample = s, 
                      class = sample.class, subtype = sample.class), 
       aes(score, X)) + 
  geom_point(aes(color = subtype, text = sample)) + 
  scale_color_aaas() + 
  geom_smooth(method = "lm", se = F, aes(color = class))
  
plotly::ggplotly(p)

in.group = c(pou2f3.samples[pou2f3.samples %in% sv.samples])
in.group = c(in.group, "SCLC0287-515_NA_H3K4me3-726_06072023-96")
not.in.group = c(ascl1.samples, ascl1.neurod1.samples, yap1.samples)

in.group = c(ascl1.samples, ascl1.neurod1.samples)
in.group = in.group[in.group %in% sv.samples]
# in.group = c(in.group, "SCLC0287-515_NA_H3K4me3-726_06072023-96")
not.in.group = c(pou2f3.samples, neurod1.samples, yap1.samples)

in.group = c(neurod1.samples, ascl1.neurod1.samples)
in.group = in.group[in.group %in% sv.samples]
# in.group = c(in.group, "SCLC0287-515_NA_H3K4me3-726_06072023-96")
not.in.group = c(pou2f3.samples, ascl1.samples, yap1.samples)

i = pou2f3.class.wins[3]
X = win_data_all[i, c(in.group, not.in.group)]
Class = sample.class
Score = sample.score

round(quantile(win.q, seq(0,.3,.03)), digits = 3)
# pou2f3.class.wins = high.windows[which(win.q < 0.01)]
hist(win.p["rmse",], breaks = "fd", xlim = c(0,20))
hist(win.p["coef1",], breaks = "fd", xlim = c(-50,50), ylim = c(0,100))
hist(win.p["coef2",], breaks = "fd", xlim = c(-50,50), ylim = c(0,100))
neurod1.class.wins = high.windows[which(win.q < 0.01 & 
                                          win.p["coef1",] > 15 &
                                          # abs(win.p["coef2",]) < 5 &
                                          win.p["rmse",] < 10)]


quantile(win.p)


x = sapply(high.windows, function(i)
  t.test(log(1+win_data_all[i,in.group]), log(1+win_data_all[i,not.in.group]))$p.value)
q = p.adjust(x, method = "fdr")
quantile(q, na.rm = T)
glm
head(sort(q))

high.healthy = rowQuantiles(win_data_all[, h.samples], probs = .90)
hist(high.healthy, breaks = "fd", xlim = c(0,3), ylim = c(0,5e4))
low.healthy.ind = which(high.healthy < 0.5)
low.healthy.win.data = win_data_all[low.healthy.ind,sv.samples]
low.healthy.win.data = sweep(low.healthy.win.data, 2, 
                             rowMaxs(cbind(estimated.tumor[sv.samples, "SCLC.n"], rep(.1, length(sv.samples))))
                             , "/")
in.group = pou2f3.samples[pou2f3.samples %in% sv.samples]
not.in.group = c(ascl1.samples, ascl1.neurod1.samples, yap1.samples)
x = sapply(1:nrow(low.healthy.win.data), function(i)
  t.test(low.healthy.win.data[i,in.group], low.healthy.win.data[i,not.in.group])$p.value)
q = p.adjust(x, method = "fdr")
xx = low.healthy.ind[which(rowMeans(low.healthy.win.data[1:1e5, in.group]) > rowMeans(low.healthy.win.data[1:1e5, not.in.group]) &
                             q < .01)]
sort(colSums(win_data_all[xx, sv.samples]))
sort(colSums(win_data_all[xx, in.group]))
sort(colSums(win_data_all[xx, not.in.group]))
pheatmap::pheatmap(log2(1+win_data_all[xx, low.ctDNA.samples]), fontsize = 3)

if (F) {
  s = ascl1.samples
  in.group = s[s %in% med.ctDNA.samples.wRNA]; not.in.group = med.ctDNA.samples.wRNA[med.ctDNA.samples.wRNA
                                                                                     %in% c(neurod1.samples, pou2f3.samples, yap1.samples)]
  # in.group = s[s %in% high.ctDNA.labeled]; not.in.group = med.ctDNA.labeled[med.ctDNA.labeled
  #                                                                                    %in% c(neurod1.samples, pou2f3.samples, yap1.samples)]
  # ascl1.class.wins = which(rowQuantiles(win_data_all[,in.group], probs = .2) > 
  #                            rowQuantiles(win_data_all[,not.in.group], probs = .9) & 
  #                            rowQuantiles(win_data_all[,in.group], probs = .2) > 3 & healthy.mean < .3)
  ascl1.class.wins = which(rowQuantiles(win_data_all[,in.group], probs = 0) > 
                             rowQuantiles(win_data_all[,not.in.group], probs = 1) & 
                             rowQuantiles(win_data_all[,in.group], probs = .1) > 3)
  ascl1.q.val = p.adjust(sapply(ascl1.class.wins, function(w) 
    t.test(win_data_all[w, in.group], win_data_all[w, not.in.group])$p.value), 
    method = "fdr")
  s = yap1.samples
  in.group = s[s %in% med.ctDNA.labeled]; not.in.group = med.ctDNA.labeled[med.ctDNA.labeled 
                                                                                      %in% c(ascl1.samples, ascl1.neurod1.samples, neurod1.samples, pou2f3.samples)]
  yap1.class.wins = which(rowQuantiles(win_data_all[,in.group], probs = .2) > 
                            rowQuantiles(win_data_all[,not.in.group], probs = .9) & rowQuantiles(win_data_all[,in.group], probs = .2) > 3 & healthy.mean < .3)
  # no yap1 windows
  s = c(neurod1.samples)
  in.group = s[s %in% med.ctDNA.samples.wRNA]; not.in.group = med.ctDNA.samples.wRNA[med.ctDNA.samples.wRNA
                                                                                      %in% c(ascl1.samples, pou2f3.samples, yap1.samples)]
  # use only the highes samples of every patient
  metadata[in.group, "SCLC_score", drop = F]
  in.group = setdiff(in.group, "SCLC0126-45_NA_H3K4me3-39_26032021-13")
  # in.group = s[s %in% med.ctDNA.labeled]; not.in.group = high.ctDNA.labeled[high.ctDNA.labeled
  #                                                                                     %in% c(ascl1.samples, pou2f3.samples, yap1.samples)]
  neurod1.class.wins = which(rowQuantiles(win_data_all[,in.group], probs = .2) > rowQuantiles(win_data_all[,not.in.group], probs = .9) &
                               rowQuantiles(win_data_all[,in.group], probs = .2) > 3 & healthy.mean < .3)
  neurod1.q.val = p.adjust(sapply(neurod1.class.wins, function(w) 
    t.test(win_data_all[w, in.group], win_data_all[w, not.in.group])$p.value), 
    method = "fdr")
  s = c(ascl1.neurod1.samples)
  in.group = s[s %in% med.ctDNA.labeled]; not.in.group = med.ctDNA.labeled[med.ctDNA.labeled
                                                                            %in% c(pou2f3.samples, yap1.samples)]
  ascl1.neurod1.class.wins = which(rowQuantiles(win_data_all[,in.group], probs = .2) > rowQuantiles(win_data_all[,not.in.group], probs = .9) &
                               rowQuantiles(win_data_all[,in.group], probs = .2) > 3 & healthy.mean < .3)
  s = pou2f3.samples
  in.group = s[s %in% med.ctDNA.labeled]
  in.group = c(in.group, "SCLC0287-515_NA_H3K4me3-726_06072023-96")
  not.in.group = med.ctDNA.labeled[med.ctDNA.labeled %in% c(ascl1.samples, ascl1.neurod1.samples, yap1.samples)]
  # pou2f3.class.wins = which(rowMins(as.matrix(win_data_all[,in.group])) / rowQuantiles(win_data_all[,not.in.group], probs = .9) > 2 &
  #                             rowMins(as.matrix(win_data_all[,in.group])) > 10 & healthy.mean < .3)
  pou2f3.class.wins = which(rowMins(win_data_all[,in.group]) >
                              rowQuantiles(win_data_all[,not.in.group], probs = .9)  &
                              rowMins(win_data_all[,in.group]) > 3 & healthy.mean < .3)
  
  pou2f3.q.val = p.adjust(sapply(pou2f3.class.wins, function(w) 
    t.test(win_data_all[w, in.group], win_data_all[w, not.in.group])$p.value), 
    method = "fdr")
  # ascl1.class.wins = setdiff(ascl1.class.wins, ascl1.wins)
  # neurod1.class.wins = setdiff(neurod1.class.wins, pou2f3.class.wins)
  # pou2f3.class.wins = setdiff(pou2f3.class.wins, neurod1.class.wins)
  # ascl1.neurod1.class.wins = setdiff(ascl1.neurod1.class.wins, c(ascl1.class.wins, neurod1.class.wins))
  class.genes.chip = c(ascl1.class.wins, neurod1.class.wins, pou2f3.class.wins)
  data.sig.wins = data.frame(signature = c(rep("ASCL1", length(ascl1.class.wins)),
                                           # rep("ASCL1-NEUROD1", length(ascl1.neurod1.class.wins)),
                                           rep("NEUROD1", length(neurod1.class.wins)), 
                                           rep("POU2F3", length(pou2f3.class.wins))), 
                             # qval = c(ascl1.q.val, neurod1.q.val, pou2f3.q.val),
                             window = class.genes.chip)
  data.sig.wins$gene = TSS.windows[data.sig.wins$window,]$name
  tableS6 = data.frame(TSS.windows[data.sig.wins$window,], signature = data.sig.wins$signature)
  write.csv(tableS6, "~/BloodChIP/Analysis/Projects/NIH_SCLC/cfChIP-paper/Cancer Discover submission/TableS6.csv")
  write.csv(data.frame(tableS6, q.value = data.sig.wins$qval), 
            paste0(baseDir, "subtype_signature_windows.csv"), quote = F, row.names = F)
}

use_all = F # use all samples (demonstrate that choosing high score samples is 
            # important for finding diff windows)
if (use_all) {
  s = ascl1.samples
  in.group = s; not.in.group = setdiff(labeled.samples, s)
  ascl1.class.wins = which(rowQuantiles(win_data_all[,in.group], probs = .2) > 
                             rowQuantiles(win_data_all[,not.in.group], probs = .9) & 
                             rowQuantiles(win_data_all[,in.group], probs = .2) > 3 & healthy.mean < .3)
  s = yap1.samples
  in.group = s; not.in.group = setdiff(labeled.samples, s)
  yap1.class.wins = which(rowQuantiles(win_data_all[,in.group], probs = .2) > 
                            rowQuantiles(win_data_all[,not.in.group], probs = .9) & rowQuantiles(win_data_all[,in.group], probs = .2) > 3 & healthy.mean < .3)
  s = neurod1.samples
  in.group = s; not.in.group = setdiff(labeled.samples, s)
  neurod1.class.wins = which(rowQuantiles(win_data_all[,in.group], probs = .2) > rowQuantiles(win_data_all[,not.in.group], probs = .9) &
                               rowQuantiles(win_data_all[,in.group], probs = .2) > 3 & healthy.mean < .3)
  s = c(ascl1.neurod1.samples)
  in.group = s; not.in.group = setdiff(labeled.samples, s)
  ascl1.neurod1.class.wins = which(rowQuantiles(win_data_all[,in.group], probs = .2) > rowQuantiles(win_data_all[,not.in.group], probs = .9) &
                                     rowQuantiles(win_data_all[,in.group], probs = .2) > 3 & healthy.mean < .3)
  s = pou2f3.samples
  in.group = s; not.in.group = setdiff(labeled.samples, s)
  pou2f3.class.wins = which(rowMins(as.matrix(win_data_all[,in.group])) / rowQuantiles(win_data_all[,not.in.group], probs = .9) > 2 &
                              rowMins(as.matrix(win_data_all[,in.group])) > 10 & healthy.mean < .3)
  # ascl1.class.wins = setdiff(ascl1.class.wins, ascl1.wins)
  neurod1.class.wins = setdiff(neurod1.class.wins, pou2f3.class.wins)
  pou2f3.class.wins = setdiff(pou2f3.class.wins, neurod1.class.wins)
  # ascl1.neurod1.class.wins = setdiff(ascl1.neurod1.class.wins, c(ascl1.class.wins, neurod1.class.wins))
  class.genes.chip = c(ascl1.class.wins, neurod1.class.wins, pou2f3.class.wins)
  data.sig.wins = data.frame(signature = c(rep("ASCL1", length(ascl1.class.wins)),
                                           # rep("ASCL1-NEUROD1", length(ascl1.neurod1.class.wins)),
                                           rep("NEUROD1", length(neurod1.class.wins)), 
                                           rep("POU2F3", length(pou2f3.class.wins))), 
                             window = class.genes.chip)
  data.sig.wins$gene = TSS.windows[data.sig.wins$window,]$name
  write.csv(data.sig.wins, paste0(baseDir, "subtype_signature_windows.csv"), quote = F, row.names = F)
  tableS6 = data.frame(TSS.windows[data.sig.wins$window,], signature = data.sig.wins$signature)
  write.csv(tableS6, "~/BloodChIP/Analysis/Projects/NIH_SCLC/cfChIP-paper/Cancer Discover submission/TableS6.csv")
  write.csv(data.sig.wins, paste0(baseDir, "subtype_signature_windows.csv"), quote = F, row.names = F)
}
data.sig.wins = read.csv(paste0(baseDir, "subtype_signature_windows.csv"))
# Jan 2024
hits = findOverlaps(TSS.windows, GRanges(data.sig.wins))
windows = queryHits(hits)[subjectHits(hits)]
data.sig.wins$window = windows
# till here
# 
class.genes.chip = data.sig.wins$window
ascl1.class.wins = data.sig.wins$window[data.sig.wins$signature == "ASCL1"]
# ascl1.neurod1.class.wins = data.sig.wins$window[data.sig.wins$signature == "ASCL1-NEUROD1"]
neurod1.class.wins = data.sig.wins$window[data.sig.wins$signature == "NEUROD1"]
pou2f3.class.wins = data.sig.wins$window[data.sig.wins$signature == "POU2F3"]
# signature windowns in RNA samples ---------------------------------------
table(TSS.windows[data.sig.wins$window,]$type)
sig.win.genes = unique(data.sig.wins[, c("signature", "gene")])
for (i in 1:nrow(sig.win.genes)) {
  if (length(grep(";", sig.win.genes$gene[i])) > 0) {
    sig = sig.win.genes$signature[i]
    genes = unlist(strsplit(sig.win.genes$gene[i], ";"))
    sig.win.genes = sig.win.genes[-i,]
    for (j in 1:length(genes)) 
      sig.win.genes = rbind(sig.win.genes, data.frame(signature = sig, gene = genes[j]))
  }
}
sig.win.genes = sig.win.genes[sig.win.genes$gene %in% rownames(rna_data),]
sig.win.genes = sig.win.genes[order(sig.win.genes$signature),]
rownames(sig.win.genes) = sig.win.genes$gene
d = rna_data[sig.win.genes$gene,]
d = sweep(d, 1, rowMeans(d), "-")
d[d > 2] = 2
d[d < -2] = -2
hist(rowMeans(d))
ord = clusterSubgroups(d, annotation.subtype, "subtype_cutoff")
clr.an = group.colors.heatmap
clr.an$subtype_cutoff = clr.an$subtype
clr.an$subtype_cutoff = c(clr.an$subtype_cutoff, "MIXED" = "#d9d9d9", "NA" = "#d9d9d9")
pheatmap::pheatmap(d[,ord], show_rownames = F, show_colnames = F, 
                   cluster_rows = F, cluster_cols = F,
                   annotation_row = sig.win.genes[,"signature", drop = F],
                   annotation_col = annotation.subtype[,c("liver", "subtype_cutoff")], 
                   annotation_colors = clr.an, color = cfm.color,
                   fontsize_col = 4, width = 8, height = 6, 
                   filename = paste0(figDirPaper, "subtypes/signature_windows_RNA_samples.png"))

s = 
data.rna.sig = data.frame(row.names = colnames(rna_data), subtype = annotation.subtype[colnames(rna_data), "subtype_cutoff"],
                          ascl1 = colSums(rna_data[sig.win.genes$gene[sig.win.genes$signature == "ASCL1"],]), 
                          neurod1 = colSums(rna_data[sig.win.genes$gene[sig.win.genes$signature == "NEUROD1"],]),
                          pou2f3 = colSums(rna_data[sig.win.genes$gene[sig.win.genes$signature == "POU2F3"],]))
pheatmap::pheatmap(t(data.rna.sig[ord,c("ascl1", "neurod1", "pou2f3")]), 
                   cluster_rows = F, cluster_cols = F, fontsize_col = 4,
                   annotation_col = annotation.subtype[,c("liver", "subtype_cutoff")], 
                   annotation_colors = clr.an, 
                   width = 8, height = 6)
data.rna.sig.m = melt(data.rna.sig, id.vars = "subtype", variable.name = "signature", value.name = "score")
head(data.rna.sig.m)
p = ggplot(data.rna.sig.m, aes(subtype, score, group = subtype, fill = subtype)) 
p = p + geom_boxplot(alpha = .7, outlier.shape = NA)
p = p + geom_point(position = position_jitterdodge(), size = .5)
p = p + facet_wrap(~signature, scale = "free_y")
p = p + scale_color_manual(values = annotation.color$subtype_cutoff, limits = force)
p = p + scale_fill_manual(values = annotation.color$subtype_cutoff, limits = force)
ggsave(paste0(figDirPaper, "subtypes/signature_windows_RNA_sample_boxplots.png"), p, width = 10, height = 5)
# heatmap and PCA signature windows -----------------------------------------------
gp = gpar(fontsize = base_size)
# s = c(low.ctDNA.labeled, pou2f3.samples.low.qc); outname = "_low"
s = c(labeled.samples); outname = "_all"
s = s[annotation.subtype[s,"subtype_cutoff"] %in% c(subtypes, "ASCL1-NEUROD1")]
s = c(s, v.samples)
# s = labeled.samples
d = log2(1+win_data_all[class.genes.chip, s])
d = sweep(d, 1, rowMedians(d), "-")
quantile(d, seq(0,1,.1))
d[d > 2] = 2; 
scale_range = c(-2, 0, 2)
scale_label = c("x1/4", "median", "x4")
ord = rankSubgroups(annotation.subtype[s,], "subtype_cutoff", "SCLC_score", decreasing = T)
ord = clusterSubgroups(d, annotation.subtype[s,], group_by = "subtype_cutoff",
                       group_order = c("ASCL1", "ASCL1-NEUROD1", "NEUROD1", "POU2F3", "YAP1"))
# d = d[,ord]
gp = gpar(fontsize = base_size)
gap.r = cumsum(c(length(ascl1.class.wins), 
                 length(neurod1.class.wins), length(pou2f3.class.wins)))
annotation.r = data.sig.wins$signature
annotation.params  = list(title_position = "lefttop-rot", 
                          legend_height = unit(1, "cm"),
                          legend_width = unit(1, "pt"),
                          title_gp = gp, 
                          labels_gp = gp, 
                          legend_gp = gp)
heatmap_legend_param = c(annotation.params, 
                         list(title = "log2 (1 + reads)",
                              at = scale_range,
                              labels = scale_label))
tf.a = HeatmapAnnotation(df = annotation.subtype[colnames(d),subtypes], 
                          col = annotation.color,  
                          annotation_label = subtypes, 
                          show_legend = c(F,F,F,F), 
                          annotation_name_gp = gp, 
                          annotation_height = unit(c(1,1,1,1), "mm"), 
                          annotation_name_side = "left", 
                          simple_anno_size = unit(2, "mm"),
                          annotation_legend_param = annotation.params)
sta = annotation.subtype[colnames(d), "SCLC_score", drop = F]

subtype.a = HeatmapAnnotation(df = sta, 
                              annotation_label = c("SCLC score"),
                              show_legend = c(F),
                              simple_anno_size = unit(2, "mm"),
                              annotation_name_gp = gp,
                              col = annotation.color,
                              # annotation_height = unit(c(5), "mm"),
                              annotation_name_side = "left",
                              foo = anno_block(height = unit(15, "pt"), 
                                               labels = c("A", "A+N", "N", "P", "Y", "no RNA"), 
                                               gp = gpar(fill = c(group.colors.heatmap$subtype[1:5], "NA" = "grey"), lwd = 0), 
                                               labels_gp = gp))
row.lab = TSS.windows$name[class.genes.chip]
row.lab = row.lab[grep(";", row.lab, invert = T)]
row.ind = which(!is.na(row.lab) & row.lab != ".")
row.ind = sample(row.ind, 20)
gene.a = HeatmapAnnotation(foo = anno_mark(labels_gp = gp, 
                                           at = row.ind, 
                                           labels = row.lab[row.ind], 
                                           lines_gp = gpar(lwd = .2)), 
                           which = "row")
# pdf(file = paste0(figDirPaper, "subtypes/subtype_signature_heatmap_cfm", 
#                   outname, ".pdf"), width = 4.8, height = 2.7)
pdf(file = paste0(figDirPaper, "subtypes/subtype_signature_heatmap_cfm", 
                  outname, ".pdf"), width = 8, height = 2.7)
png(file = paste0(figDirPaper, "subtypes/subtype_signature_heatmap_cfm", 
                  outname, ".png"), width = 8, height = 2.7, units = "in", res = 800)
Heatmap(d, col = cfm.color, cluster_row_slices = T, 
        cluster_column_slices = T, show_column_names = F, 
        row_split = annotation.r, 
        row_title = "signature genomic regions",
        row_title_gp = gpar(fontsize = base_size), 
        column_split = annotation.subtype[colnames(d), "subtype_cutoff"], 
        column_gap = unit(1, "pt"), row_gap = unit(1, "pt"),
        column_title_gp = gpar(fontsize = base_size),
        column_title = " ",
        cluster_columns = F, cluster_rows = F, column_order = ord,
        top_annotation = subtype.a, use_raster = T, raster_quality = 10,
        # bottom_annotation = tf.a, 
        right_annotation = gene.a,
        heatmap_legend_param = heatmap_legend_param)
dev.off()

d = log2(1+win_data_all[class.genes.chip, s])
chip.pca = prcomp(t(d))
chip.pca = data.frame(PC1 = chip.pca$x[,1], PC2 = chip.pca$x[,2], PC3 = chip.pca$x[,3],
                      subtype = annotation.subtype[s, "subtype_cutoff"], 
                      SCLC_score = metadata[s, "SCLC_score"])
p = ggplot(chip.pca, aes(PC1, PC2, color = subtype)); outname = "pc1_pc2"; xlab = "PC1"; ylab = "PC2"
p = ggplot(chip.pca, aes(PC1, PC3, color = subtype)); outname = "pc1_pc3"; xlab = "PC1"; ylab = "PC3"
p = ggplot(chip.pca, aes(PC2, PC3, color = subtype)); outname = "pc2_pc3"; xlab = "PC2"; ylab = "PC3"
p = p + geom_point(size = 1, aes(alpha = SCLC_score)) 
p = p + geom_point(data = chip.pca[!is.na(chip.pca$subtype),], size = 1, shape = 1, stroke = .3)
p = p + scale_color_manual(values = annotation.color$subtype_cutoff, limits = force)
p = p + guides(color = "none") + labs(x = xlab, y = ylab, alpha = "SCLC score") #+ guides(alpha = "none")
p = p + theme(legend.key.size = unit(2, "mm"), aspect.ratio = 1, legend.position = c(.8,0.8))
ggsave(paste0(figDirPaper, "subtypes/signature", outname, ".pdf"), p, width = 55, height = 55, units = "mm")

# signature vs. SCLC score ----------------------------------------
sWh = c(sv.samples, h.samples)
data.subtypes.sig = data.frame(row.names = sWh, sample = sWh,
                               ASCL1 = colSums(win_data_all[ascl1.class.wins,sWh]), 
                               # YAP1 = colSums(win_data_all[yap1.class.wins,sWh]),
                               NEUROD1 = colSums(win_data_all[neurod1.class.wins,sWh]), 
                               POU2F3 = colSums(win_data_all[pou2f3.class.wins,sWh]),
                               SCLC_score = estimated.tumor[sWh, "SCLC.n"],
                               ASCL1.score = colSums(win_data_all[ascl1.class.wins, sWh]) / estimated.tumor[sWh, "SCLC.n"], 
                               # YAP1.score = colSums(win_data_all[yap1.class.wins, sWh]) / estimated.tumor[sWh, "SCLC.n"], 
                               NEUROD1.score = colSums(win_data_all[neurod1.class.wins, sWh]) / estimated.tumor[sWh, "SCLC.n"], 
                               POU2F3.score = colSums(win_data_all[pou2f3.class.wins, sWh]) / estimated.tumor[sWh, "SCLC.n"],
                               subtype = annotation.subtype[sWh, "subtype_cutoff"])

data.m = melt(data.subtypes.sig, id.vars = c("sample", "subtype", "SCLC_score"), 
              measure.vars = c("ASCL1", "NEUROD1", "POU2F3"), variable.name = "signature", 
              value.name = "reads")
data.m$subtype = as.factor(data.m$subtype)
# known.subytpe = c(rownames(annotation.subtype)[!is.na(annotation.subtype$subtype_cutoff)], h.samples)
# data.m = data.m[data.m$sample %in% known.subytpe,]
p = ggplot(data.m, aes(SCLC_score, reads, color = subtype))
p = p + geom_point(size = .5, aes(text = sample))
p = p + scale_color_manual(values = annotation.color$subtype_cutoff, limits = force)
p = p + facet_wrap(~signature, scales = "free") + labs(y = "reads/signature")
p = p + theme(strip.text = element_text(face = "bold"), 
              strip.background = element_blank(), 
              legend.key.size = unit(2, "mm"))
p = p + guides(color = "none")
p = p + geom_vline(xintercept = low.score.cutoff, linetype = "dashed", size = base_line_size)
ggsave(paste0(figDirPaper, "subtypes/signature_vs_score.pdf"), p, width = 110, height = 50, units = "mm")
p = plotly::ggplotly(p)
# plotly::or(p = p, file = paste0(figDirPaper, "subtype/signature_vs_score.html"))
p = p + geom_smooth(data = data.m[which(data.m$subtype %in% subtypes),], method = "rlm", se = F, linewidth = .3)
ggsave(paste0(figDirPaper, "subtypes/signature_vs_score_regline.pdf"), p, width = 110, height = 50, units = "mm")
ggsave(paste0(figDirPaper, "subtypes/signature_vs_score_regline.png"), p, width = 110, height = 50, units = "mm")


# signature score (heatmap, boxplot and scatter comparing scores) -----------------------------------------------
# ord = rankSubgroups(data.subtypes.sig[s,], "subtype", "SCLC_score", decreasing = T)
d = log2(1+t(data.subtypes.sig[s, c("ASCL1", "NEUROD1", "POU2F3")])); 
d = log2(1+t(data.subtypes.sig[low.ctDNA.samples, c("ASCL1.score", "NEUROD1.score", "POU2F3.score")])); 
outname = "signature_sum"; 
min.c = quantile(d, 0.05)
max.c = quantile(d, 0.95)
col = colorRamp2(seq(min.c, max.c,length = 10), hcl.colors(10, "Reds", rev = T))
title = "log2(reads)"
scale_label = c("x1/8", "median", "x8") # change
cfm = T
if (cfm) {
  d = sweep(d, 1, rowMedians(d, na.rm = T), "-")
  quantile(d,na.rm = T)
  d[d > 3] = 3
  d[d < -3] = -3
  scale_range = c(-3, 0, 3);
  outname = "signature_sum_cfm"; 
  title = cfm.leg
  col = cfm.color
}
ord = clusterSubgroups(d, annotation.subtype[colnames(d),], group_by = "subtype_cutoff",
                       group_order = c("ASCL1", "ASCL1-NEUROD1", "NEUROD1", "POU2F3", "YAP1"))
colnames(d) = sub("_.*", "", colnames(d))
hlp = heatmap_legend_param
hlp$title = ""
hlp$at = scale_range
hlp$labels = scale_label
pdf(paste0(figDirPaper, "subtypes/", outname, ".pdf"), height = 2, width = 6)
png(paste0(figDirPaper, "subtypes/", outname, ".png"), width = 120, height = 40, 
    units = "mm", res = 500)
Heatmap(d, show_column_names = F, 
        column_title = " ",
        column_gap = unit(1, "pt"), 
        # column_order = ord, 
        column_names_gp = gpar(fontsize = 3),
        cluster_columns = T,  
        cluster_rows = F,
        cluster_column_slices = F, 
        row_title = "signature", 
        row_title_side = "left", 
        row_names_gp = gp,
        row_title_gp = gp, 
        row_names_side = "left",
        heatmap_legend_param = hlp,
        top_annotation = subtype.a, 
        col = col, 
        column_split = annotation.subtype[colnames(d), "subtype_cutoff"], 
        heatmap_height = unit(5, "cm"))
dev.off()

# data.subtypes.sig$total_sig = apply(data.subtypes.sig[, c("ASCL1", "NEUROD1", "POU2F3")], 1, "sum")
data.subtypes.sig$subtype[data.subtypes.sig$subtype %in% c("MIXED", "LOW-TF")] = NA
data.m = melt(data.subtypes.sig, id.vars = c("sample", "subtype", "SCLC_score"), 
              measure.vars = c("ASCL1.score", "NEUROD1.score", "POU2F3.score"), variable.name = "signature", 
              value.name = "score")
data.m$subtype = as.factor(data.m$subtype)
head(data.m)
ascl1.cutoff = 1050
neurod1.cutoff = 1200
pou2f3.cutoff = 1200
data.m$
p = ggplot(data.m[!is.na(data.m$subtype) & data.m$sample %in% low.ctDNA.labeled,], 
           aes(subtype, score, fill = subtype)); outname = "signature_boxplot_norm"; ylab = "signature score"
p = ggplot(data.m[data.m$sample %in% low.ctDNA.samples,], 
           aes(subtype, score, fill = subtype, group = subtype, alpha = SCLC_score, text = sample)); outname = "signature_boxplot_val"; ylab = "signature score"
# p = ggplot(data.m, aes(subtype, score/total_sig, fill = subtype, label = substr(sample, 5,11))); outname = "signature_boxplot_norm2"; ylab = "norm. reads/signature"
# p = ggplot(data.m, aes(subtype, score, fill = subtype, label = substr(sample, 5,11))); outname = "signature_boxplot"
# p = boxplotWpoints(data.m, "signature", "score", "subtype", plot_stat = F, jitter_w = .1, dodge_w = .85)
p = p + geom_point(position = position_jitterdodge(jitter.width = .5, dodge.width = .95), size = .5, aes(color = subtype))
p = p  + geom_boxplot(size = base_line_size, outlier.shape = NA, coef = 0, alpha = .7, position = position_dodge(width = .95))
p = p + scale_color_manual(values = annotation.color$subtype_cutoff, limits = force)
p = p + scale_fill_manual(values = annotation.color$subtype_cutoff, limits = force)
p = p + theme(strip.text = element_text(face = "bold"), 
              strip.background = element_blank(), 
              legend.key.size = unit(3, "mm"))
p = p + guides(fill = "none", color = "none") + labs(y = ylab) 
p = p + facet_wrap(~signature, drop = T, scales = "free") 
p = p + scale_x_discrete(breaks = c("ASCL1", "ASCL1-NEUROD1", "NEUROD1", "POU2F3", "YAP1"), 
                           labels = c("A", "A+N", "N", "P", "Y"))
# p = p + geom_vline(xintercept = which(levels(data.m$subtype) == "YAP1") + .5, 
                   # size = 0.3, linetype = "dashed")
plotly::ggplotly(p)
p = p + geom_hline(data = data.m[data.m$signature == "ASCL1.score",], aes(yintercept = ascl1.cutoff), 
                   linetype = "dashed", linewidth = .3)
p = p + geom_hline(data = data.m[data.m$signature == "NEUROD1.score",], aes(yintercept = neurod1.cutoff), 
                   linetype = "dashed", linewidth = .3)
p = p + geom_hline(data = data.m[data.m$signature == "POU2F3.score",], aes(yintercept = pou2f3.cutoff), 
                   linetype = "dashed", linewidth = .3)
ggsave(paste0(figDirPaper, "subtypes/", outname, ".pdf"), p, height = 55, width = 110, units = "mm")
ggsave(paste0(figDirPaper, "subtypes/", outname, ".png"), p, height = 55, width = 110, units = "mm")


a = "ASCL1"; b = "NEUROD1"; x.cutoff = ascl1.cutoff; y.cutoff = neurod1.cutoff; vline = T
a = "ASCL1"; b = "POU2F3"; x.cutoff = ascl1.cutoff; y.cutoff = pou2f3.cutoff; vline = F
a = "NEUROD1"; b = "POU2F3"; x.cutoff = neurod1.cutoff; y.cutoff = pou2f3.cutoff; vline =F

outname = paste0(a, "_", b, "_score.pdf")
p = ggplot(data.subtypes.sig[data.subtypes.sig$sample %in% low.ctDNA.labeled & 
                               !is.na(data.subtypes.sig$subtype),], 
           aes_string(paste0(a, ".score"), paste0(b, ".score"), color = "subtype")) + geom_point(size = .5)
p = p + scale_color_manual(values = annotation.color$subtype, limits = force)
p = p + labs(x = paste0(a, " score"), y = paste0(b, " score"))
p = p + theme(aspect.ratio = 1) + guides(color = "none")
p = p + geom_segment(aes(x.cutoff, -1, xend = x.cutoff, yend = y.cutoff), linetype = "dashed",
                     color = "black", linewidth = .3)
if (vline)
  p = p + geom_vline(xintercept = x.cutoff, linetype = "dashed",size = .3)
p = p + geom_hline(yintercept = y.cutoff, linetype = "dashed",linewidth = .3)
ggsave(paste0(figDirPaper, "subtypes/", outname), p, width = 45, height = 45, units = "mm")

# ROC signature score -----------------------------------------------------
data.subtype.roc = data.subtypes.sig[low.ctDNA.labeled, ]
data.subtype.roc = data.subtype.roc[!is.na(data.subtype.roc$subtype),]
roc.ascl1 = roc(data.subtype.roc$subtype %in% c("ASCL1", "ASCL1-NEUROD1"), 
                data.subtype.roc$ASCL1.score)
roc.neurod1 = roc(data.subtype.roc$subtype %in% c("NEUROD1", "ASCL1-NEUROD1"), 
                  data.subtype.roc$NEUROD1.score)
roc.pou2f3 = roc(data.subtype.roc$subtype == "POU2F3", 
                 data.subtype.roc$POU2F3.score)
p = ggroc(list(ASCL1 = roc.ascl1, NEUROD1 = roc.neurod1, POU2F3 = roc.pou2f3), size = .5)
# p = p + ggtitle(paste0('ROC Curve ', '(AUC = ', auc, ')'))
p = p + geom_abline(intercept = 1, slope = 1, size = base_line_size, linetype = "dashed")
p = p + scale_color_manual(values = annotation.color$subtype_cutoff, limits = force, 
                           labels = c(paste0("ASCL1 (auc = ", round(roc.ascl1$auc, 3), ")"),
                                      paste0("NEUROD1 (auc = ", round(roc.neurod1$auc, 3), ")"),
                                      paste0("POU2F3 (auc = ", round(roc.pou2f3$auc, 3), ")")))
p = p + labs(color = "") + theme(legend.position = c(.75,.25))
ggsave(paste0(figDirPaper, "subtypes/ROC_subtypes.pdf"), p, width = 55, height = 55, units = "mm")

# correlation of signature score and RNA TF -------------------------------
pou2f3.genes.chip = intersect(rownames(rna_data), unique(TSS.windows$name[pou2f3.class.wins]))
ascl1.genes.chip = intersect(rownames(rna_data), unique(TSS.windows$name[ascl1.class.wins]))
neurod1.genes.chip = intersect(rownames(rna_data), unique(TSS.windows$name[neurod1.class.wins]))
p = qplot(t(rna_data["ASCL1", ]), colSums(rna_data[ascl1.genes.chip, ]), xlab = "lo2(RPKM) - RNA", 
          ylab = "reads/siganture - ChIP", main = "ASCL1") + stat_cor()
ggsave(paste0(figDirPaper, "subtypes/ASCL1_chip_sig_vs_RNA.png"), p, width = 4, height = 4)
p = qplot(t(rna_data["NEUROD1", ]), colSums(rna_data[neurod1.genes.chip, ]), xlab = "lo2(RPKM) - RNA", 
          ylab = "reads/siganture - ChIP", main = "NEUROD1") + stat_cor()
ggsave(paste0(figDirPaper, "subtypes/NEUROD1_chip_sig_vs_RNA.png"), p, width = 4, height = 4)
p = qplot(t(rna_data["POU2F3", high.ctDNA.samples.wRNA]), colSums(rna_data[pou2f3.genes.chip, high.ctDNA.samples.wRNA]), xlab = "lo2(RPKM) - RNA", 
          ylab = "reads/siganture - ChIP", main = "POU2F3") + stat_cor()
ggsave(paste0(figDirPaper, "subtypes/POU2F3_chip_sig_vs_RNA.png"), p, width = 4, height = 4)
# classify samples -------------------------------------
# preparation: discard samples with very low ctDNA
annotation.subtype$subtype_chip = "YAP1?"
no.ctDNA.samples = sv.samples[estimated.tumor[s.samples, "SCLC.n"] < low.score.cutoff]
ctDNA.samples = sv.samples[estimated.tumor[s.samples, "SCLC.n"] >= low.score.cutoff]
annotation.subtype[no.ctDNA.samples, "subtype_chip"] = "Unknown"
# s = c(ctDNA.samples, pou2f3.samples.low.qc)
s = ctDNA.samples

# step1: evaluate pou2f3 score
s[data.subtypes.sig[s, "POU2F3.score"] > pou2f3.cutoff]
annotation.subtype[s[which(data.subtypes.sig[s, "POU2F3.score"] > pou2f3.cutoff)], "subtype_chip"] = "POU2F3"

# step2: evaluate ascl1 and neurod1 score
ascl1.high.score = s[data.subtypes.sig[s, "ASCL1.score"] > ascl1.cutoff]
neurod1.high.score = s[data.subtypes.sig[s, "NEUROD1.score"] > neurod1.cutoff]
ascl1.neurod1.high.score = intersect(ascl1.high.score, neurod1.high.score)
annotation.subtype[ascl1.high.score, "subtype_chip"] = "ASCL1" 
annotation.subtype[neurod1.high.score, "subtype_chip"] = "NEUROD1" 
annotation.subtype[ascl1.neurod1.high.score, "subtype_chip"] = "ASCL1-NEUROD1" 

table(annotation.subtype$subtype_chip) / sum(table(annotation.subtype[ctDNA.samples, "subtype_chip"]))
table(annotation.subtype$subtype_cutoff) / sum(table(annotation.subtype$subtype_cutoff))
write.csv(annotation.subtype, paste0(figDirPaper, "subtypes/annotations_subtype.csv"))

if (F) {  # old version
  # step1: identify POU2f3 samples
  obs = colSums(win_data_all[pou2f3.class.wins,s])
  pred.pos = predict(lm.pou2f3, metadata[s,"SCLC_score", drop = F])
  pred.neg = predict(lm.Npou2f3, metadata[s,"SCLC_score", drop = F])
  qplot(abs(pred.pos - obs), abs(pred.neg - obs))
  pou2f3.pos = names(which(abs(pred.neg - obs) / abs(pred.pos - obs) > 1))
  annotation.subtype[pou2f3.pos, "subtype_chip"] = "POU2F3"
  s = setdiff(s, pou2f3.pos)
  
  # step2: identify NEUROD1 samples
  obs = colSums(win_data_all[neurod1.class.wins,s])
  pred.pos = predict(lm.neurod1, metadata[s,"SCLC_score", drop = F])
  pred.neg = predict(lm.Nneurod1 , metadata[s,"SCLC_score", drop = F])
  qplot(abs(pred.pos - obs), abs(pred.neg - obs)) + geom_abline()
  neurod1.pos = names(which(abs(pred.neg - obs) / abs(pred.pos - obs) > 1))
  annotation.subtype[neurod1.pos, "subtype_chip"] = "NEUROD1"
  s = setdiff(s, neurod1.pos)
  
  # step3: identify ASCL1 samples
  obs = colSums(win_data_all[ascl1.class.wins,s])
  pred.pos = predict(lm.ascl1, metadata[s,"SCLC_score", drop = F])
  pred.neg = predict(lm.Nascl1 , metadata[s,"SCLC_score", drop = F])
  qplot(abs(pred.pos - obs), abs(pred.neg - obs)) + geom_abline()
  ascl1.pos = names(which(abs(pred.neg - obs) / abs(pred.pos - obs) > 1 & abs(pred.neg - obs) > 50))
  annotation.subtype[ascl1.pos, "subtype_chip"] = "ASCL1"
  
  table(paste(annotation.subtype$subtype_cutoff, annotation.subtype$subtype_chip, sep = "-"))
}
# compare RNA vs. ChIP classification -------------------------------------
x = rownames(annotation.subtype)[annotation.subtype$subtype_cutoff %in% c(subtypes, "ASCL1-NEUROD1")]
x = x[annotation.subtype[x, "subtype_chip"] != "Unknown"]
# x = x[grep("Unknown", annotation.subtype[x, "subtype_chip"], invert = T)]
table(paste(annotation.subtype[x, "subtype_cutoff"], annotation.subtype[x, "subtype_chip"], sep = "-------"))
sum(x %in% med.ctDNA.labeled)
sum(!x %in% med.ctDNA.labeled)

comp.subtype = annotation.subtype
comp.subtype$sample = rownames(comp.subtype)

st.lev = c("ASCL1", "ASCL1-NEUROD1", "NEUROD1", "POU2F3", "YAP1", "Unknown_(YAP1?)", "Unknown")
st.lab = c("1", "2", "3", "4", "5", "6", "7")
comp.subtype$subtype_chip = factor(comp.subtype$subtype_chip, levels = st.lev, labels = st.lab)
comp.subtype$subtype_cutoff = factor(comp.subtype$subtype_cutoff, levels = st.lev, labels = st.lab)
d = apply(comp.subtype[,c("subtype_chip", "subtype_cutoff")], 2, as.numeric)
d[is.na(d)] = 0
ro = hclust(dist(d))$order
ro = rownames(annotation.subtype)[ro]
comp.subtype.m = melt(comp.subtype[ro,],  id.vars = "sample",
                            measure.vars = c("subtype_cutoff", "subtype_chip"),
                            variable.name = "assay", value.name = "subtype")
tail(comp.subtype.m)
comp.subtype.m$sample = factor(comp.subtype.m$sample, levels = ro)
s = s.samples
s = unique(comp.subtype.m$sample[!is.na(comp.subtype.m$subtype)]); outname = "subtype_RNA_vs_ChIP"
s = s[s %in% rna.samples]; outname = "subtype_RNA_vs_ChIP_wRNA"
s = s[s %in% med.ctDNA.samples.wRNA]; outname = "subtype_RNA_vs_ChIP_med_wRNA" 
p = ggplot(comp.subtype.m[comp.subtype.m$sample %in% s,], aes(assay, sample, fill = subtype, color = "black")) 
p = p + geom_tile() + scale_fill_manual(values = c("1" = unname(annotation.color$subtype_cutoff["ASCL1"]), 
                                                     "2" = unname(annotation.color$subtype_cutoff["ASCL1-NEUROD1"]), 
                                                     "3" = unname(annotation.color$subtype_cutoff["NEUROD1"]),
                                                     "4" = unname(annotation.color$subtype_cutoff["POU2F3"]), 
                                                     "5" = unname(annotation.color$subtype_cutoff["YAP1"]),
                                                     "6" = unname(annotation.color$subtype_cutoff["YAP1"]), 
                                                     "7" = "grey",
                                                     "0" = "grey"), labels = st.lev)
p = p + guides(color = "none", fill = "none")
p = p + scale_x_discrete(breaks=c("subtype_cutoff","subtype_chip"), labels=c("RNA", "ChIP"))
p = p + scale_y_discrete(breaks=s, labels=sub("_.*", "", s))
p = p + labs(x = "Assay", y = "Sample") 
p = p + theme(axis.ticks = element_blank(), line = element_blank(), 
              axis.text.y = element_text(hjust = 0, size = 4), legend.key.size = unit(2, "mm"), 
              axis.text.x = element_text(angle = 90))
ggsave(paste0(figDirPaper, "subtypes/", outname, ".pdf"), p, width = 35, height = 100, units = "mm")
# subtype during treatment (subtype switch) ------------------------------------------------
s = s.samples
mult.timepoint = metadata %>% count(Patient_id, Primary_tumor_sample) %>% 
  count(Patient_id) %>% filter(n > 1) %>% select(Patient_id)
data.sub.time = cbind(sample = rownames(annotation.subtype), 
                      annotation.subtype, metadata[rownames(annotation.subtype), 
                      c("Timepoint", "Patient_id", "Primary_tumor_sample")])

data.sub.time$Timepoint = factor(data.sub.time$Timepoint, levels = c("Pre-treatment", "Post-treatment", 
                                                                     "Disease-progression", "Ongoing-response", 
                                                                     "Best-response"))
data.sub.time.m = melt(data.sub.time, id.vars = c("Patient_id", "Primary_tumor_sample", 
                                                  "Timepoint", "sample"), 
                       measure.vars = c("subtype_cutoff", "subtype_chip"),
                       variable.name = "assay", value.name = "subtype")
data.sub.time.m$subtype = as.factor(data.sub.time.m$subtype)
data = data.sub.time.m[data.sub.time.m$assay == "subtype_cutoff" &
                         data.sub.time.m$subtype %in% c(subtypes, "ASCL1-NEUROD1") &
                         data.sub.time.m$Patient_id %in% mult.timepoint$Patient_id,]
data = data %>% complete(Patient_id, Timepoint)
p = ggplot(data, aes(Timepoint, Patient_id, fill = subtype)) 
p = p + geom_tile(color = "black") 
p = p + scale_fill_manual(values = group.colors, na.value = "#f0f0f0")
p = p + labs(x = "", y = "Patient", title = "tumor RNA-seq")
p = p + theme(axis.line = element_blank(), #axis.text.y = element_blank(), 
              axis.ticks = element_blank(), legend.key.size = unit(3,"mm"), 
              axis.text.x = element_text(angle = 90))
p
ggsave(paste0(figDirPaper, "subtypes/subtype_during_treatment.pdf"), p, 
       width = 60, height = 55, units = "mm")

data = data.sub.time.m[data.sub.time.m$assay == "subtype_chip" &
                         data.sub.time.m$subtype %in% c(subtypes, "ASCL1-NEUROD1"),]
# EMT score ---------------------------------------------------------------
# Byers (cancer cell 2021) https://www.sciencedirect.com/science/article/pii/S1535610820306620?via%3Dihub
# defines a Inflammatory subtype (instead of YAP1), which has a high EMT score
emt.genes = read.csv(paste0(baseDir, "data_from_NIH/EMT-signature-Byers_Clin_Cancer_Res-2013.csv"), header = F)$V1
emt.genes.rna = emt.genes[emt.genes %in% rownames(rna_data)]
d = rna_data[emt.genes.rna,]
d = sweep(d, 1, rowMedians(as.matrix(d)), "-")
d[d > 3] = 3
d[d < -3] = -3
pheatmap(d, show_rownames = T, show_colnames = F, fontsize_row = 4,
         annotation_col = annotation.subtype[,"subtype", drop = F])

emt.score = data.frame(subtype = annotation.subtype[rownames(rna.subtypes),"subtype"], 
                       emt.score = colSums(rna_data[emt.genes.rna,rownames(rna.subtypes)]))
emt.score$emt.score =  emt.score$emt.score / mean(emt.score$emt.score)
p = ggplot(emt.score, aes(subtype, emt.score, group = subtype)) + geom_boxplot(outlier.shape = NA) + 
  geom_point(position = position_jitter(width = .2)) + labs(y = "EMT score", x = "")
ggsave(paste0(figDirPaper, "subtypes/emt_score_rna.png"), p, width = 6, height = 5)
# this does not seem like figure 2 in the Bayer paper

g = "SYP"#  "CHGA""REST" 
data.gene = data.frame(subtype = annotation.subtype[rownames(rna.subtypes),"subtype"], 
                       t(rna_data[g,rownames(rna.subtypes)]))
names(data.gene) = c("subtype", "counts")
p = ggplot(data.gene, aes(subtype, counts, group = subtype, label = sub("_.*", "", rownames(data.gene)))) 
p = p + geom_boxplot(outlier.shape = NA) +  geom_point(position = position_jitter(width = .2)) 
p = p + labs(y = "log2(RPKM)", x = "", title = g) #+ geom_text_repel(size = 2)
ggsave(paste0(figDirPaper, "subtypes/", g, "_subtypes.png"), p, width = 6, height = 5)


d = sweep(emt.score, 2, colMeans(emt.score), "-")
d[d > 3] = 3
d[d < -3] = -3
pheatmap(d, fontsize_row = 4)

head(emt.score)
emt.genes = emt.genes[emt.genes %in% rownames(chip_data_all)]
hist(rowMeans(chip_data_all[emt.genes, h.samples]), breaks = "fd")
emt.genes = names(which(rowMeans(chip_data_all[emt.genes, h.samples]) < 10))
d = log2(1+chip_data_all[emt.genes,high.ctDNA.samples])
d = sweep(d, 1, rowMeans(d), "-")
d[d > 3] = 3
d[d < -3] = -3
pheatmap(d, show_rownames = T, show_colnames = F, fontsize_row = 4,
         annotation_col = annotation.subtype[,c("subtype", "SCLC_score")])

inflam.dir = paste0(baseDir, "cfChIP-paper/data/enrichR_gene_sets/")
for (sig in list.files(inflam.dir, "*.gmt")) { 
  print(sig)
  inflamation.sig = read.gmt(paste0(inflam.dir, sig))[[1]]
  inflamation.sig = inflamation.sig[inflamation.sig %in% matched.genes]
  hist(rowMeans(chip_data_all[inflamation.sig, h.samples]), breaks = "fd")
  inflamation.sig = names(which(rowMeans(chip_data_all[inflamation.sig, h.samples]) < 10))
  inf.score = data.frame(subtype = annotation.subtype[med.ctDNA.samples,"subtype"], 
                         inf.score = colSums(chip_data_all[inflamation.sig,med.ctDNA.samples]))
  inf.score$inf.score =  inf.score$inf.score
  p = ggplot(inf.score, aes(subtype, inf.score, group = subtype, color = subtype)) + geom_boxplot() + 
    geom_point(position = position_jitter(width = .2)) + labs(y = "reads/inflmation signature")
  p = p + scale_color_manual(values = annotation.color$subtype_cutoff, limits = force)
  ggsave(paste0(inflam.dir, "figures/", gsub("%", "-", sig), ".png"), p, width = 8, height = 6)
}

byers.genes = data.frame(read_xlsx("~/BloodChIP/Analysis/Projects/NIH_SCLC/data_from_NIH/1-s2.0-S1535610820306620-mmc2.xlsx", 
                        skip = 1, col_names = T))
byers.genes = byers.genes$...1
byers.genes = byers.genes[byers.genes %in% rownames(rna_data)]
d = rna_data[byers.genes,rownames(rna.subtypes)]
d = sweep(d, 1, rowMeans(d), "-")
d[d > 3] = 3
d[d < -3] = -3
pheatmap(d, show_rownames = F, show_colnames = T, annotation_col = annotation.subtype[,"subtype", drop = F])
pheatmap(d, show_rownames = F, show_colnames = F, annotation_col = data.frame(rna.subtypes), 
         filename = paste0(figDirPaper, "subtypes/bayer_signature_rna.png"), width = 8, height = 6)




# PCA variable genes ------------------------------------------------------
X = healthy.ref$Gene.var < .3
sum(X)
s = med.ctDNA.labeled
XX = rowSds(chip_data_all[X,s]) > 4
sum(XX)
d = win_data_all[XX,s]
chip.pca = prcomp(t(d))
chip.pca = data.frame(PC1 = chip.pca$x[,1], PC2 = chip.pca$x[,2], PC3 = chip.pca$x[,3],
                      subtype = annotation.subtype[s, "subtype_cutoff"], 
                      SCLC_score = metadata[s, "SCLC_score"])
p = ggplot(chip.pca, aes(PC1, PC2, color = subtype)); outname = "pc1_pc2"; xlab = "PC1"; ylab = "PC2"
p = p + geom_point(size = 1) 
p = p + scale_color_manual(values = annotation.color$subtype_cutoff, limits = force)
p
# compare gene groups -----------------------------------------------------
diffgeneDir = paste0(figDirPaper, "subtypes/groups/")
diffgene.a = read.csv(paste0(diffgeneDir, "ascl1_siggenes.csv"))
diffgene.a = diffgene.a$X[diffgene.a$Q.Value...log10. > 2 & diffgene.a$FoldChange..log2. < 0]
diffgene.y = read.csv(paste0(diffgeneDir, "yap1_siggenes.csv"))
diffgene.y = diffgene.y$X[diffgene.y$Q.Value...log10. > 2 & diffgene.y$FoldChange..log2. < 0]
diffgene.p = read.csv(paste0(diffgeneDir, "pou2f3_siggenes.csv"))
diffgene.p = diffgene.p$X[diffgene.p$Q.Value...log10. > 2 & diffgene.p$FoldChange..log2. < 0]
diffgene.all = c(diffgene.a, diffgene.n, diffgene.p, diffgene.y)

diffwin.n = read.csv(paste0(diffgeneDir, "neurod1_sigwins_significant_wins.csv"))
diffwin.n = diffwin.n$X.1[diffwin.n$qvalue > 2 & diffwin.n$Significant.up]
diffwin.y = read.csv(paste0(diffgeneDir, "yap1_sigwins_significant_wins.csv"))
diffwin.y = diffwin.y$X.1[diffwin.y$qvalue > 2 & diffwin.y$Significant.up]
diffwin.p = read.csv(paste0(diffgeneDir, "pou2f3_sigwins_significant_wins.csv"))
diffwin.p = diffwin.p$X.1[diffwin.p$qvalue > 2 & diffwin.p$Significant.up]
diffwin.a = read.csv(paste0(diffgeneDir, "ascl1_sigwins_significant_wins.csv"))
diffwin.a = diffwin.a$X.1[diffwin.a$qvalue > 2 & diffwin.a$Significant.up]
diffwin.a = setdiff(diffwin.a, c(diffwin.n, diffwin.y, diffwin.p))
diffwin.n = setdiff(diffwin.n, c(diffwin.a, diffwin.y, diffwin.p))
diffwin.y = setdiff(diffwin.y, c(diffwin.n, diffwin.a, diffwin.p))
diffwin.p = setdiff(diffwin.p, c(diffwin.n, diffwin.y, diffwin.a))
diffwin.all = c(diffwin.a, diffwin.n, diffwin.p, diffwin.y)
# win.annotation = data.frame(row.names = diffwin.all, 
#                             subtype = c(rep("ASCL1", length(diffwin.a)), 
#                             rep("NEUROD1", length(diffwin.n)), 
#                             rep("POU2F3", length(diffwin.p)),
#                             rep("YAP1", length(diffwin.y))))
d = log2(1+win_data_all[diffwin.all, s])
# d = log2(1+chip_data_all[diffgene.y, med.ctDNA.samples])
d = sweep(d, 1, rowMeans(d), "-")
d[d > 3] = 3
d[d < -3] = -3
col.a = annotation.subtype[, c("subtype_cutoff", "SCLC_score")]
col.a$subtype_cutoff[col.a$subtype_cutoff %in% c("MIXED", "LOW-TF")] = NA
names(col.a) = c("subtype (RNA)", "SCLC_score")
annotation.color[["subtype (RNA)"]] = annotation.color$subtype[!names(annotation.color$subtype) %in% c("no-RNA", ".")]
# rownames(d) = diffwin.all
gaps.r = cumsum(c(length(diffwin.a), length(diffwin.n), 
                   length(diffwin.p)))
pheatmap::pheatmap(d[,ord], show_rownames = F, cluster_rows = F, cluster_cols = F,
         fontsize_col = 4, show_colnames = F, width = 10, height = 8,
         annotation_col = col.a, #annotation_row = win.annotation, 
         annotation_colors = annotation.color, drop_levels = T, 
         breaks = cfm.breaks, color = cfm.color, gaps_row = gaps.r)




