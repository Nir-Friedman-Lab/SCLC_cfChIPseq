# heatmap -----------------------------------------------------------------
# diff.genes are defined in linear-regression-SCLC-score.R
s = c(s.samples, nec.samples, h.samples, l.samples, c.samples)
log2(1+chip_data_all[diff.genes.common, s]) %>% 
  as.data.frame() %>%
  sweep(1, rowMeans(.), "-") %>% 
  mutate(across(everything(), ~ ifelse(. > 3, 3, .))) %>% 
  mutate(across(everything(), ~ ifelse(. < -3, -3, .))) ->
  data.signature

test.groups = c("SCLC", "NEC", "Healthy", "NSCLC", "CRC")

colnames(data.signature)[rankSubgroups(estimated.tumor[colnames(data.signature),],
                        group_by = "group", rank_by = "SCLC.n",
                        decreasing = T, group_order = test.groups)] ->
  col.ord

hc = PClusterMatrixRows(as.matrix(data.signature))
hc.c = hclust(dist(t(data.signature)))
hc.r = hclust(dist(data.signature))

data.signature = data.signature[hc$labels[hc$order], col.ord]
heatmap.groups = sample.annotation[col.ord,"group"]
gap_cols = cumsum(table(factor(heatmap.groups, levels = unique(heatmap.groups))))
breaks = colnames(data.signature)[ceiling(rollmean(c(0,gap_cols),2))]
lab = c(paste0("SCLC\n(n=", length(s.samples), ")"), 
        paste0("NEC\n(n=", length(nec.samples), ")"), 
        paste0("Healthy\n(n=", length(h.samples), ")"), 
        paste0("NSCLC\n(n=", length(l.samples), ")"), 
        paste0("CRC\n(n=", length(c.samples), ")"))



# data.signature %>% 
#   rownames_to_column("gene") %>% 
#   melt(id.vars = "gene", variable.name = "samp") %>%
#   heatmap(x = "samp", y = "gene", fill = "value", 
#         chigh = "#E54C00", clow = "#00A7FE", cmed = "black", lim = c(-3,3), 
#         breaks = c(-3,0,3), leglab = "counts", leg = c("X1/8", "mean", "X8")) +
#   geom_vline(xintercept = gap_cols+.5, color = "white") + 
#   scale_x_discrete(breaks = breaks, labels = lab, position = "top") +
#   theme(axis.text.y = element_blank(), axis.ticks = element_blank(),
#         axis.text.x = element_text(face = "bold"), 
#         legend.key.height = unit(1, 'mm'), legend.key.width = unit(2,"mm"), 
#         legend.title = element_text(size=5), legend.text = element_text(size=3)) +
#   labs(y = paste(length(diff.genes.common), "genes")) -> p
# p = rasterise(p, dpi = 500)
# ggsave(paste0(figDirPaper, "figure2/signature_heatmap.pdf"), p, width = 120, 
#        height = 65, units = "mm")
# saveRDS(data.signature, paste0(baseDir, "signature_matrix.rds"))


if (T) {
  d = data.signature[,col.ord]
  scale_range = c(-3,0,3)
  scale_label = c("x1/8", "med.", "x8")
  annotation.params  = list(title_position = "lefttop-rot", 
                            legend_height = unit(1, "cm"),
                            legend_width = unit(1, "pt"),
                            title_gp = gp, 
                            labels_gp = gp, 
                            legend_gp = gp)
  heatmap_legend_param = c(annotation.params, 
                           list(title = "",
                                at = scale_range,
                                labels = scale_label))
  rbind(data.frame(ind = 1:nrow(d), 
                   gene = rownames(d)) %>% 
          filter(gene %in% c("DLL3", "INSM1", "CHGA", "CRMP1")), 
        data.frame(ind = 1:nrow(d), 
                   gene = rownames(d)) %>% 
          slice_sample(n = 15)) -> data.gene 
  HeatmapAnnotation(foo = anno_mark(labels_gp = gp, 
                                    at = data.gene$ind, 
                                    labels = data.gene$gene,
                                    lines_gp = gpar(lwd = .2), 
                                    link_width = unit(2, "mm")), which = "row", 
                    annotation_name_gp = gp, 
                    gp = gp, 
                    show_legend = T) -> gene.a
  
  subtype.a = HeatmapAnnotation(df = sample.annotation[colnames(d), "group", drop = F],
                                annotation_label = "",
                                show_legend = F,
                                simple_anno_size = unit(2, "mm"),
                                annotation_name_gp = gp,
                                foo = anno_block(height = unit(15, "pt"),
                                                 # labels = c("SCLC", "NEC", 
                                                 #            "Healthy", "NSCLC", "CRC"),
                                                 labels = lab,
                                                 gp = gpar(lwd = 0, fill = 0),
                                                 # gp = gpar(fill = c(group.colors$$subtype[1:5], "NA" = "grey"), lwd = 0),
                                                 labels_gp = gp))
  pdf(file = paste0(figDirPaper, "figure2/signature_heatmap_V1.pdf"), width = 5, 
      height = 7)
  png(file = paste0(figDirPaper, "figure2/signature_heatmap_V1.png"))
  Heatmap(d, 
          heatmap_width = unit(115, "mm"),
          heatmap_height = unit(60, "mm"),
          col = cfm.color, 
          right_annotation = gene.a,
          cluster_column_slices = F,
          show_column_names = F,
          show_row_names = F,
          show_row_dend = F,
          row_title = paste(nrow(d), "genes"),
          row_title_gp = gpar(fontsize = base_size),
          column_split = factor(sample.annotation[colnames(d), "group"],
                                levels = c("SCLC", "NEC", "Healthy", "NSCLC", "CRC")),
          column_gap = unit(1, "pt"),
          column_title_gp = gpar(fontsize = base_size),
          column_title = " ",
          cluster_columns = F,
          top_annotation = subtype.a,
          use_raster = T, raster_quality = 10,
          heatmap_legend_param = heatmap_legend_param)
  dev.off()
}
# write.table(data.signature[hc.r$labels[hc.r$order], hc.c$labels[hc.c$order]],file = paste0("~/treeview_files/SCLC/heatmap_sclc_clustered.cdt"), 
#             sep = '\t', row.names = T, col.names = TRUE, na="", quote=FALSE)
# r2cdt does not write gene names for some reason. I ran the raw code https://rdrr.io/bioc/ctc/src/R/r2cdt.R
# and before saving the table, added the following line: data$NAME = rownames(data)
# r2cdt(hc = hc.c, hr = hc.r, data = data.signature, file = paste0("~/treeview_files/SCLC/heatmap_sclc_clustered.cdt"))
# r2atr(hc = hc.c, distance = "euclidean", file = paste0("~/treeview_files/SCLC/heatmap_sclc_clustered.atr"))
# r2gtr(hr = hc.r, distance = "maximum", file = paste0("~/treeview_files/SCLC/heatmap_sclc_clustered.gtr"))

# save gene signature (for enrichment tests etc.)
# boxplot of reads/SCLC-signature -----------------------------------------
comparisons = list(c("SCLC", "NEC"), c("SCLC", "Healthy"), 
                   c("SCLC", "NSCLC"), c("SCLC", "CRC"))
data.frame(group = factor(sample.annotation[s,"group"], 
                                       levels = test.groups), 
         counts = colSums(chip_data_all[diff.genes.common, s]), 
         tumor = estimated.tumor[s, "SCLC.n"])  -> data.sig.sum
data.sig.sum %>% 
  boxplotWpoints(x = "group", y = "counts", fill = "group", 
                 ylab = "reads/signature", comparisons = comparisons) -> p 
ggsave(paste0(figDirPaper, "figure2/reads_per_gene_signature_boxplot.pdf"), p, 
       width = 50, height = 55, units = "mm")

# ROC - SCLC score - healthy vs. SCLC (not used) -------------------------------
data.roc = data.sig.sum[c(s.samples, h.samples),]
data.roc$group = droplevels(data.roc$group) 
roc.pre = roc(data.roc[c(pre.samples, h.samples), "group"], 
              data.roc[c(pre.samples, h.samples),"counts"])
roc.post = roc(data.roc[c(post.samples, h.samples), "group"], 
               data.roc[c(post.samples, h.samples),"counts"])
roc.both = roc(data.roc$group, data.roc$counts)
ggroc(list(pretreatment = roc.pre, posttreatment = roc.post, all = roc.both), 
      size = base_line_size) + 
  geom_abline(intercept = 1, slope = 1, size = base_line_size, 
              linetype = "dashed") + 
  scale_color_aaas(labels = c(paste0("pre (auc = ", round(roc.pre$auc, 3), ")"),
                              paste0("post (auc = ", round(roc.post$auc, 3), ")"),
                              paste0("all (auc = ", round(roc.both$auc, 3), ")"))) + 
  labs(color = "") + theme(legend.position = c(.75,.25)) -> p
ggsave(paste0(figDirPaper, "figureS2/ROC_sig_genes.pdf"), p, width = 55, height = 55, units = "mm")
