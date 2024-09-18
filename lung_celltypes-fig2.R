# lung cell markers in samples (from single-cell lung atlas) --------------
# prepare data ------------------------------------------------------------
lung.ind = colnames(gene.atlas)[grep("lung", colnames(gene.atlas))]
ReCompute = F
if (ReCompute) { 
  sc.lung.file = "~/gavriel.fialkoff@mail.huji.ac.il - Google Drive/Shared drives/Friedman Lab Shared Drive/BloodChIP/Analysis/Projects/NIH_SCLC/sc_Lung_atlas_data/41586_2020_2922_MOESM6_ESM.xlsx"
  sc.clusters = excel_sheets(sc.lung.file)
  sc.clusters = sc.clusters[grep("SS2", sc.clusters)]
  # lung.cell.types = read.csv(paste0(baseDir, "sc_Lung_atlas_data/lung_cell_types.csv"))
  
  data.lung = data.frame()
  genes.per.clus = data.frame()
  cluster.genes = list()
  min.genes = 10
  for (i in 1:length(sc.clusters)) {
    clus = sc.clusters[i]
    print(i)
    clus.name = names(read_xlsx(sc.lung.file, sheet = clus, range = "A1:A1"))
    read_xlsx(sc.lung.file, sheet = clus, skip = 1) %>%
      filter(pct_out_cluster < .1 & avg_logFC > 2 & p_val_adj < .1 &
               Gene %in% rownames(chip_data_all) & 
               healthy.ref[Gene] < 3) %>%
      select(Gene) %>% .$Gene -> sc.lung.markers
    
    cluster.genes[[clus.name]] = sc.lung.markers
    if (length(sc.lung.markers) < min.genes) { next}
    genes.per.clus = rbind(genes.per.clus, 
                           data.frame(cluster = clus.name, 
                                      n.genes = length(sc.lung.markers)))
    cat(clus.name, " - ", length(sc.lung.markers), "\n")
    df = data.frame(chip_data_all[sc.lung.markers,],
                    cluster = rep(clus.name, length(sc.lung.markers)), check.names = F)
    # roadmap.lung = data.frame(LUNG = gene.atlas[sc.lung.markers,"LNG"])
    roadmap.lung = data.frame(gene.atlas[sc.lung.markers,lung.ind])
    # roadmap.lung = data.frame(lung = mean(gene.atlas[sc.lung.markers,lung.ind]))
    data.lung = rbind(data.lung, cbind(df,roadmap.lung))
  }
  
  data.cluster.genes = melt(cluster.genes)
  names(data.cluster.genes) = c("gene", "cell_type")
  write.csv(data.cluster.genes, paste0(figDirPaper, 'lung_cell_type_markers.csv', 
                                       row.names = F))
  write.csv(genes.per.clus, paste0(figDirPaper, '/lung_cell_type_n_markers.csv', 
                                   row.names = F))

  data.lung %>% 
    group_by(cluster) %>% 
    summarise(across(everything(), sum)) %>% 
    pivot_longer(cols = c(-cluster), names_to = "sample")%>%
    pivot_wider(names_from = c(cluster)) %>% 
    mutate(group = sample.annotation[sample,"group"]) ->
    data.lung_celltypes
    
  write.csv(data.lung_celltypes, 
            paste0(figDirPaper, "reads_per_lung_cell_type.csv"))
  }

data.lung_celltypes = read.csv(paste0(figDirPaper, "reads_per_lung_cell_type.csv"), 
                               row.names = 1, check.names = F)

# heatmap of lung cell-types over all samples  ----------------------------
data.lung_celltypes %>% 
  filter(sample %in% sh.samples) %>% 
  column_to_rownames("sample") %>% 
  select(-c(group)) -> d
d.log = t(log2(1+d))
data.frame(row.names = rownames(estimated.tumor), 
           group = estimated.tumor$group, 
           SCLC.score = estimated.tumor$SCLC.n) -> ca
col.ord = clusterSubgroups(d.log, ca, "group", group_order = c("SCLC", "Healthy"))
# col.ord = rankSubgroups(estimated.tumor[colnames(d.log),], "group", "SCLC.n", 
#                         decreasing = T, group_order = c("SCLC", "Healthy"))
pheatmap::pheatmap(d.log[,col.ord], 
                   show_colnames = F, annotation_col = ca, 
                   cluster_cols = F,
                   color = cfm.color,
                   fontsize = 6, treeheight_row = 0, height = 4,
                   annotation_colors = group.colors.heatmap, width = 3,
                   # filename = paste0(figDirPaper, "figureS2/lung_cell_type_heatmap.pdf"),
                   annotation_legend = F,) 
dev.off()


# boxplot lung cell type distribution  ------------------------------------
read.csv(paste0(figDirPaper, "lung_cell_type_markers.csv")) -> 
  lung_celltypes_genes

lung_celltypes_genes %>% 
  count(cell_type) -> lung_celltypes_genes_n
rownames(lung_celltypes_genes_n) = lung_celltypes_genes_n$cell_type

data.lung_celltypes %>%
  filter(sample %in% sh.samples) %>%
  melt(id.vars = c("sample", "group"), 
       variable.name = "cell_type", value.name = "val") %>%
  mutate(val.per.gene = val/lung_celltypes_genes_n[cell_type, "n"]) %>%
  mutate(cell_type = reorder(cell_type, val.per.gene, mean)) -> data.lc.m

boxplotWOpoints(data.lc.m, "cell_type", "val.per.gene", "group", plot_stat = F) +
  # ylim(c(NA,6000)) +
  # scale_y_break(c(1000,3000)) +
  coord_flip() + 
  labs(x = "lung cell-type", y = "reads/gene") -> p.box

lung_celltypes_genes_n %>%
  filter(cell_type %in% levels(data.lc.m$cell_type)) %>%
  mutate(cell_type = factor(cell_type, levels(data.lc.m$cell_type))) %>%
  ggplot(aes(cell_type, n)) + 
  geom_bar(stat="identity") +
  labs(y = "genes/marker", x = "") +
  # scale_y_continuous(trans = "reverse", position = "left") +
  # scale_x_discrete(position = "top") + 
  coord_flip() +
  theme(axis.text.y = element_blank()) -> p.bar
p = plot_grid(p.box, p.bar, align = "h", nrow = 1, rel_widths = c(7/10, 3/10))
ggsave(paste0(figDirPaper, "figureS2/lung_cell_types_sum.pdf"), p, width = 90, 
       height = 85, units = "mm")


# zoom in on some of the cell types -------------------------------------
selected.lung.tiss = c("Alveolar Epithelial Type 1", "B", "Lymphatic", 
                       "Ciliated", "Neuroendocrine")
data.lc.m %>%
  filter(cell_type %in% selected.lung.tiss) %>%
  mutate(cell_type = sub("B", "B cells", cell_type)) %>%
  mutate(cell_type = sub("Alveolar Epithelial Type 1", 
                         "Alveolar Epithelial\nType 1", cell_type)) %>%
  mutate(cell_type = reorder(cell_type, val.per.gene, mean)) %>%
  boxplotWOpoints("cell_type", "val.per.gene", "group") +
  coord_flip() +
  labs(x = "lung cell-type", y = "reads/gene") -> p
ggsave(paste0(figDirPaper, "figure2/lung_cell_types_sum_selected.pdf"), p, 
       width = 50, height = 55, units = "mm")

# signal in lung cell types relative to roadmap lung ("absolute") ---------
# num.ind = !names(data.lung_celltypes) %in% c("sample", "group")
# lung.signal = colMeans(data.lung_celltypes[lung.ind, num.ind])
data.lung_celltypes %>% 
  filter(sample %in% lung.ind) %>% 
  select(-c("sample", "group")) %>% 
  colMeans() -> lung.signal

data.lung_celltypes %>%
  column_to_rownames("sample") %>% 
  select(!group) %>%
  sweep(2,lung.signal, "/") %>% 
  mutate(sample = rownames(.)) %>%
  left_join(estimated.tumor %>% select(sample, SCLC.n, group), by = "sample") %>% 
  melt(id.vars = c("sample", "group", "SCLC.n"), variable.name = "cell_type", 
       value.name = "val") %>% 
  filter(sample %in% c(s.samples), SCLC.n > high.score.cutoff) %>%
  mutate(val = log2(val), cell_type = reorder(cell_type, val, median)) -> data.lct

data.lct %>%
  boxplotWOpoints("cell_type", "val", "group", plot_stat = F) + 
  geom_hline(yintercept = 0, linetype="dashed", linewidth = base_line_size) +
  coord_flip() + 
  labs(x = "lung cell-type", y = "log2(sample/Roadmap lung)") -> p
 
ggsave(paste0(figDirPaper, "figure2/cell_types_relative_sum_high.pdf"), 
       p, width = 50, height = 55, units = "mm") 

data.lct %>%
  filter(cell_type %in% selected.lung.tiss) %>%
  boxplotWOpoints("cell_type", "val", "group", plot_stat = F) + 
  geom_hline(yintercept = 0, linetype="dashed", linewidth = base_line_size) +
  coord_flip() + 
  labs(x = "lung cell-type", y = "log2(sample/Roadmap lung)") -> p
ggsave(paste0(figDirPaper, "figure2/selected_cell_types_relative_sum_high.pdf"), 
       p, width = 50, height = 55, units = "mm") 

# roadmap.lung.dist = as.numeric(data.lung_celltypes[rownames(data.lung_celltypes) == "LUNG",num.ind])
# roadmap.lung.dist = as.numeric(colMedians(as.matrix(data.lung_celltypes[rownames(data.lung_celltypes) %in% lung.ind, num.ind])))
# data.lung_celltypes.relative = data.lung_celltypes
# data.lung_celltypes.relative[,num.ind] = sweep(data.lung_celltypes.relative[,num.ind], 2, roadmap.lung.dist, "/")
# data.lcr = melt(data.lung_celltypes.relative, id.vars = c("samp", "group"), variable.name = "cell_type", value.name = "val")
# d = data.lcr[data.lcr$group %in% c("SCLC"),]; outname = "lung_cell_types_relative_sum"
# d = data.lcr[data.lcr$samp %in% high.ctDNA.samples,]; outname = "lung_cell_types_relative_sum_high"
# d = data.lcr[data.lcr$group %in% c("SCLC") & data.lcr$cell_type %in% selected.lung.tiss,]; outname = "lung_selected_cell_types_relative_sum"
# d = data.lcr[data.lcr$samp %in% high.ctDNA.samples & data.lcr$cell_type %in% selected.lung.tiss,]; outname = "lung_selected_cell_types_relative_sum_high"
# d$cell_type = with(d, reorder(cell_type, val, median))
# d$val = log2(d$val)
# p = boxplotWOpoints(d, "cell_type", "val", "group", plot_stat = F)
# p = p + geom_hline(yintercept = 0,  linetype="dashed", size = base_line_size) 
# # p = p + geom_bar(stat = "summary", fun = "mean") + geom_hline(yintercept = 0,  linetype="dashed", size = 1); outname = paste0(outname, "_bar")
# p = p + coord_flip() + labs(x = "lung cell-type", y = "log2(sample/Roadmap lung)") 
# ggsave(paste0(figDirPaper, "figure2/", outname, ".pdf"), p, width = 50, 
#        height = 55, units = "mm")


# cell types vs. NE score (not used)-------------------------------------------------
names(metadata)
ne.params = c("ASCL1", "NEUROD1", "YAP1", "POU2F3", "NE_SCORE")
data.ct.ne = data.frame(metadata[s.samples, ne.params], 
                        data.lung_celltypes[s.samples,selected.lung.tiss])
data.ct.ne$ASCL1 = as.numeric(data.ct.ne$ASCL1)
data.ct.ne$NEUROD1 = as.numeric(data.ct.ne$NEUROD1)
data.ct.ne$YAP1 = as.numeric(data.ct.ne$YAP1)
data.ct.ne$POU2F3 = as.numeric(data.ct.ne$POU2F3)
data.ct.ne$NE_SCORE = as.numeric(data.ct.ne$NE_SCORE)
ggpairs(data.ct.ne)

data.ct.ne$all = data.ct.ne$ASCL1 + data.ct.ne$NEUROD1 + data.ct.ne$YAP1 + data.ct.ne$POU2F3
data.ct.ne$samp = rownames(data.ct.ne)

ggplot(data.ct.ne[rna.samples.chip.passQC,], aes(log(Neuroendocrine), log(Ciliated), alpha = POU2F3/all)) + geom_point()
ggplot(data.ct.ne[rna.samples.chip.passQC,], aes(log(Neuroendocrine), log(Ciliated), alpha = NE_SCORE)) + geom_point()
p = ggplot(data.ct.ne[high.ctDNA.samples.wRNA.matched.time,], aes(Neuroendocrine/Ciliated, NE_SCORE, label = substr(samp, 5,12))) + 
  geom_point() + stat_cor(size = base_size/.pt) + geom_text_repel(size = 2) 
ggsave(paste0(figDirPaper, "figure4/ne_vs_ciliated.pdf"), p, width = 3, height = 3)
###sweep(data.ct.ne[,c("ASCL1", "NEUROD1", "YAP1", "POU2F3")], 1, data.ct.ne$all, "/")

# ciliated/neuroendocrine ratio -------------------------------------------
# it seems that high pou2f3 samples have a high ciliated/neuroendocrine ratio

s = rna.samples.matched.time[rna.samples.matched.time %in% s.samples]
data.frame(s = s, 
           estimated.tumor[s,], 
           c = data.lung_celltypes[s, "Ciliated"],
           n = data.lung_celltypes[s, "Neuroendocrine"], 
           t(rna_data["POU2F3",s]), 
           t(rna_data["ASCL1",s]),
           t(rna_data["NEUROD1",s]), 
           t(rna_data["YAP1",s])) %>%
  mutate(sclc = estimated.tumor[s,"SCLC.n"], 
         all.rna = ASCL1 + NEUROD1 + YAP1 + POU2F3) %>%
  filter(estimated.tumor[s, "SCLC.n",] > low.score.cutoff, 
         s != "SCLC0176-120") %>% 
  # 176-120 has high erythrocyte which causes to overestimate of the sclc score
  ggplot(aes(POU2F3/all.rna, log2(c/n), color = sclc, label = s)) +
  geom_point(shape = 16, size = 1, alpha = .8) +
  labs(x = "relative POU2F3 - RNA", y = "log2 (ciliated/neuroendocrine) - ChIP", 
       color = "SCLC score") +
  geom_vline(xintercept = .2, linetype="dashed", color = "red", linewidth = .5) +
  # geom_text_repel(size = 2) +
  theme(legend.position = c(.7,.3), legend.key.size = unit(2, 'mm')) -> p
ggsave(paste0(figDirPaper, "figureS4/ciliated_ne_ratio.pdf"), p, 
       width = 60, height = 60, units = "mm")
