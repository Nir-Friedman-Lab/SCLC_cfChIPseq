# define training and test groups - healthy vs. SCLC -----------------------------------------
high.cutoff = .7
high.score = quantile(estimated.tumor$SCLC.n, high.cutoff)
high.sclc = metadata$Sample_id[which(metadata$SCLC_score > high.score)]
low.sclc = metadata$Sample_id[which(metadata$SCLC_score <= high.score)]

s.train = c(sample(high.sclc, size = ceiling(length(high.sclc) * .8), replace = F), 
                 sample(low.sclc, size = ceiling(length(low.sclc) * .8), replace = F))
h.train = sample(h.samples, size = ceiling(length(h.samples) * .8), replace = F)
s.train = setdiff(s.train, l.samples)
all.train = c(s.train, h.train)

s.test = setdiff(s.samples, s.train)
h.test = setdiff(h.samples, h.train)
all.test = c(s.test, h.test)
classifier.samples = c(all.train, all.test)

# healthy vs.  SCLC calssifier --------------------------------------------
s.high.wins = rowQuantiles(win_data_all[,s.train], probs = .3)
hist(s.high.wins, breaks = "fd", xlim = c(0,20), ylim = c(0,500))
h.low.wins = rowQuantiles(win_data_all[,h.train], probs = .9)
sclc.wins = which(s.high.wins >  h.low.wins & s.high.wins > 3)
qplot(norm.f, metadata[s.samples, "SCLC_score"])
data.classifier = data.frame(row.names = classifier.samples, group = c(rep("SCLC.train", length(s.train)), rep("H.train", length(h.train)), 
                                       rep("SCLC.test", length(s.test)), rep("H.test", length(h.test))), 
                             counts = colSums(win_data_all[sclc.wins,c(classifier.samples)]))
data.classifier$group = factor(data.classifier$group, levels = c("H.train", "SCLC.train", "H.test", "SCLC.test"))
head(data.classifier$group)
p = boxplotWpoints(data.classifier, "group", "counts", "group", ylab = "counts/signature", plot_stat = F, size = 2) + scale_color_aaas() + scale_fill_aaas()
ggsave(paste0(figDirPaper, "healthy_vs_SCLC.png"), p, height = 5, width = 5)
# p = ggplot(data.classifier, aes(group, counts, color = group)) + geom_jitter() + scale_color_aaas()


# ROC curve ---------------------------------------------------------------
data.train = data.classifier[all.train,]
data.test = data.classifier[all.test,]
data.train$group = droplevels(data.train$group)
data.test$group = droplevels(data.test$group)
roc.train = roc(data.train, "group", "counts")
roc.test = roc(data.test, "group", "counts")
p = ggroc(list(train = roc.train, test = roc.test), size = base_line_size)
# p = p + ggtitle(paste0('ROC Curve ', '(AUC = ', auc, ')'))
p = p + geom_abline(intercept = 1, slope = 1, size = base_line_size, linetype = "dashed")
p = p + scale_color_aaas(labels = c(paste0("train (auc = ", round(roc.train$auc, 3), ")"),
                                    paste0("test (auc = ", round(roc.test$auc, 3), ")")))
p = p + labs(color = "") + theme(legend.position = c(.75,.25))
ggsave(paste0(figDirPaper, "ROC_SCLC_vs_Healthy.png"), p, width = 4, height = 4)


# define training and test groups - NE vs. nonNE -----------------------------------------
# trying to enlarge the group by using calssification based on ChIP. 

# attempt to use ChIP for definding YAP1 - failed 
hist(yap1, breaks = "fd")
qplot(yap1, metadata[s.samples, "SCLC_score"])
head(sort(chip_data_all["YAP1",h.samples]), 10)
not.yap = s.samples[chip_data_all["YAP1",s.samples] == 0]
metadata[not.yap, c("NE_SCORE", "NE_status", "NAPY_allcluster_expression", "ASCL1", "NEUROD1", "POU2F3", "YAP1", "SCLC_score")]
# 033-489 seems wierd.  high yap1 in RNA but not 0 in ChIP normalized gene counts. can't trust ChIP to include or exclude yap1

pheatmap(data.subtypes[metadata$Sample_id[which(metadata$SCLC_score > .5)],], annotation_row = subtype.max)

qplot(colSums(win_data_all[pou2f3.wins,s.samples]), metadata[s.samples,"SCLC_score"])
# there are 4 samples that seem to be POU2F3.  too few to use for training a classifier. 
ne.group = c("ASCL1", "NEUROD1", "ASCL1-NEUROD1")
is.ne = metadata$NAPY_allcluster_expression %in% ne.group
nonne.group = c("POU2F3", "YAP1")
is.nonne = metadata$NAPY_allcluster_expression %in% nonne.group

high.cutoff = .7
high.score = quantile(estimated.tumor$SCLC.n, high.cutoff)
# high.ne = metadata$Sample_id[which(metadata$SCLC_score > high.score & metadata$NE_status == "NE")]
# low.ne = metadata$Sample_id[which(metadata$SCLC_score <= high.score & metadata$NE_status == "NE")]
# high.nonne = metadata$Sample_id[which(metadata$SCLC_score > high.score & metadata$NE_status == "NonNE")]
# low.nonne = metadata$Sample_id[which(metadata$SCLC_score <= high.score & metadata$NE_status == "NonNE")]
high.ne = metadata$Sample_id[which(metadata$SCLC_score > high.score & is.ne)]
low.ne = metadata$Sample_id[which(metadata$SCLC_score <= high.score & is.ne)]
high.nonne = metadata$Sample_id[which(metadata$SCLC_score > high.score & is.nonne)]
low.nonne = metadata$Sample_id[which(metadata$SCLC_score <= high.score & is.nonne)]

ne.train = c(sample(high.ne, size = ceiling(length(high.ne) * .8), replace = F), 
            sample(low.ne, size = ceiling(length(low.ne) * .8), replace = F))
nonne.train = c(sample(high.nonne, size = ceiling(length(high.nonne) * .8), replace = F), 
                sample(low.nonne, size = ceiling(length(low.nonne) * .8), replace = F))
subtype.train = c(ne.train, nonne.train)

# ne.test = setdiff(metadata$Sample_id[metadata$NE_status == "NE"], ne.train)
# nonne.test = setdiff(metadata$Sample_id[metadata$NE_status == "NonNE"], nonne.train)
ne.test = setdiff(metadata$Sample_id[metadata$NAPY_allcluster_expression %in% c("ASCL1", "NEUROD1", "ASCL1-NEUROD1")], ne.train)

# exclude sample with very high erythroblast that causes artifacts
ne.train = setdiff(ne.train, "SCLC0176-120_NA_H3K4me3-39")
ne.test = setdiff(ne.test, "SCLC0176-120_NA_H3K4me3-39")

nonne.test = setdiff(metadata$Sample_id[metadata$NAPY_allcluster_expression %in% c("POU2F3", "YAP1")], nonne.train)
subtype.test = c(ne.test, nonne.test)
subtype.samples = c(subtype.train, subtype.test)



# NE vs. nonNE calssifier --------------------------------------------
# ne.high.wins = rowQuantiles(win_data_all[,ne.train], probs = .6)
# nonne.low.wins = rowQuantiles(win_data_all[,nonne.train], probs = .8)
# nonne.high.wins = rowQuantiles(win_data_all[,nonne.train], probs = .6)
# ne.low.wins = rowQuantiles(win_data_all[,ne.train], probs = .8)
# ne.wins = which(ne.high.wins / nonne.low.wins > 3 & ne.high.wins > 3)
# nonne.wins = which(nonne.high.wins / ne.low.wins > 3 & nonne.high.wins > 3)

# use only high samples to train classifier 
ne.med.wins = rowMedians(win_data_all[,ne.train[ne.train %in% high.ne]])
nonne.med.wins = rowMedians(win_data_all[,nonne.train[nonne.train %in% high.nonne]])
healthy.med.wins = rowMedians(win_data_all[,h.samples])

ne.wins = which(ne.med.wins / nonne.med.wins > 2 & ne.med.wins > 4 & healthy.med.wins < .1)
nonne.wins = which(nonne.med.wins / ne.med.wins > 2 & nonne.med.wins > 4 & healthy.med.wins < .1)

data.classifier = data.frame(samp = subtype.samples, group = c(rep("NE.train", length(ne.train)), rep("nonNE.train", length(nonne.train)), 
                                                                       rep("NE.test", length(ne.test)), rep("nonNE.test", length(nonne.test))), 
                             ne = colSums(win_data_all[ne.wins, subtype.samples]), 
                               nonne = colSums(win_data_all[nonne.wins, subtype.samples]))
data.classifier$ratio = data.classifier$ne / data.classifier$nonne
data.classifier$group = factor(data.classifier$group, levels = c("NE.train", "nonNE.train", "NE.test", "nonNE.test"))
head(data.classifier)
p = boxplotWpoints(data.classifier, "group", "ratio", "group", ylab = "counts/signature", plot_stat = F) + scale_color_aaas() + scale_fill_aaas()
p = p + geom_text_repel(data = data.classifier, aes(group, ratio, label = substr(samp, 5,12)), size = 2)
ggsave(paste0(figDirPaper, "NE_vs_nonNE.png"), p, height = 3, width = 3)

p = ggplot(data.classifier, aes(ne, nonne, color = group, label = substr(samp, 5,12))) + geom_point() 
p = p + geom_text(data = data.classifier[data.classifier$group == "NE.test",], size = 2)
p
# ROC curve NE vs. nonNE---------------------------------------------------------------
data.train = data.classifier[subtype.train,]
data.test = data.classifier[subtype.test,]
data.train$group = droplevels(data.train$group)
data.test$group = droplevels(data.test$group)
roc.train = roc(data.train, "group", "ratio")
roc.test = roc(data.test, "group", "ratio")
p = ggroc(list(train = roc.train, test = roc.test), size = base_line_size)
# p = p + ggtitle(paste0('ROC Curve ', '(AUC = ', auc, ')'))
p = p + geom_abline(intercept = 1, slope = 1, size = base_line_size, linetype = "dashed")
p = p + scale_color_aaas(labels = c(paste0("train (auc = ", round(roc.train$auc, 3), ")"),
                                    paste0("test (auc = ", round(roc.test$auc, 3), ")")))
p = p + labs(color = "") + theme(legend.position = c(.75,.25))
ggsave(paste0(figDirPaper, "ROC_NE_vs_nonNE.png"), p, width = 3, height = 3)
