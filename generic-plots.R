library(dplyr)
library(cowplot)
library(enrichR,quietly = T)
read.table(paste0(externalDataDir, "ncbiRefSeqSelect.txt"), 
           header = T, sep = "\t", check.names = F) %>% 
  # remove genes in non classical chromosomes
  filter(!grepl("_", chrom)) -> genes

base_size = 6
base_line_size = .1

catn = function(...) { cat(...,"\n") }

# generic scatter plot 
scatter.plot = function(df, x, y, xlab = "", ylab = "", title = "", flab = "", clab = "", shape = 16, 
                        size = 1, color = NULL, alpha = NULL, guieds = F, cor = T,
                        trans.x = NULL, trans.y = NULL, xlim = NULL, ylim = NULL, text.x = NULL) {
  pval = format(max(2e-15, cor.test(df[[x]], df[[y]])$p.value), digits = 2)
  r = round(cor.test(df[[x]], df[[y]])$estimate, digits = 2)
  maxY = max(df[[y]], na.rm = T)
  minX = min(df[[x]], na.rm = T)
  p = ggplot(df, aes(.data[[x]],.data[[y]])) 
  if (!is.null(color)) {
    p = ggplot(df, aes(.data[[x]],.data[[y]], color = .data[[color]])) 
  }
  if (!is.null(alpha)) {
    p = ggplot(df, aes(.data[[x]],.data[[y]], color = .data[[color]]), alpha = .data[[alpha]]) 
  }
  p = p + geom_point(shape = shape, size = size)
  # https://stackoverflow.com/questions/65076492/ggplot-size-of-annotate-vs-size-of-element-text
  if (cor) {
    p = p + annotate(geom = "text", x=minX, y=1.2*maxY, label = paste0("R = ", r, ", p < ", pval),
                     vjust="inward",hjust="inward", color="black", fontface = 'italic', size = base_size/.pt)
    # p = p + stat_cor(label.y = 1, r.digits = 2, p.digits = 4)
  }
  
  if (!is.null(trans.x)) {
    p = p + scale_x_continuous(trans = trans.x)
  }
  if (!is.null(trans.y)) {
    p = p + scale_y_continuous(trans = trans.y)
  }
  if (!is.null(color)) {
    p = p + scale_color_manual(values = group.colors, limits = force)
  }
  if (!guieds) {
    p = p + guides(color = "none")
  }
  p = p + coord_fixed(ratio = 1, xlim = xlim, ylim = ylim)
  p = p + theme(aspect.ratio = 1, axis.text.x = text.x, axis.text.y = element_text(hjust = .5))
  p = p + labs(x = xlab, y = ylab, title = title, color = clab, fill = flab)
}

# boxplot with points
boxplotWpoints = function(df, x, y, fill, label = NULL, xlab = "", ylab = "", title = "", size = .5,
                          alpha = .5, shape = 16, guides = F, stat_test = "t.test", plot_stat = T, 
                          comparisons = NULL, jitter_w = .4, dodge_w = .75, dodge) {
  p = ggplot(df, aes(.data[[x]],.data[[y]], fill = .data[[fill]]))
  if (!is.null(label))
    p = ggplot(df, aes(.data[[x]],.data[[y]], fill = .data[[fill]], label = .data[[label]]))
  # if (jitter) 
  #   p = p + geom_jitter(width = .2, height = 0, alpha = alpha, size = size, shape = shape)
  # else 
  p = p + geom_point(aes(color = .data[[fill]]), alpha = alpha, size = size, shape = shape, 
                       position = position_jitterdodge(jitter.width = jitter_w, dodge.width = dodge_w))
  p = p + geom_boxplot(size = base_line_size, outlier.shape = NA, coef = 0, alpha = .7, position = position_dodge(width = dodge_w)) 
  p = p + labs(x = xlab, y = ylab, title = title) 
  if (!guides) {
    p = p + guides(fill = "none", color = "none")
  }
  if (all(df[[fill]] %in% names(group.colors))) {
    p = p + scale_fill_manual(values = group.colors, limits = force)
    p = p + scale_color_manual(values = group.colors, limits = force)
  }
  if (plot_stat)
    p = p + stat_compare_means(label = "p.signif", method = stat_test, 
                               bracket.size = base_line_size,tip.length = 0.01,
                               step.increase = 0.07,
                               size = base_size/.pt, comparisons = comparisons)
  p = p + theme(axis.text.y = element_text(hjust = .5))
  p
}
#boxplot without points
boxplotWOpoints = function(df, x, y, fill, xlab = "", ylab = "", title = "", guides = F, 
                           stat_test = "t.test", plot_stat = T, comparisons = NULL) {
  p = ggplot(df, aes(.data[[x]],.data[[y]], fill = .data[[fill]]))
  p = p + geom_boxplot(size = base_line_size, outlier.size = base_line_size, outlier.shape = 16, outlier.alpha = .7)
  p = p + labs(x = xlab, y = ylab, title = title)
  if (!guides) {
    p = p + guides(fill = "none")
  }
  if (all(df[[fill]] %in% names(group.colors))) {
    p = p + scale_fill_manual(values = group.colors, limits = force)
  }
  if (plot_stat)
    p = p + stat_compare_means(label = "p.signif", method = stat_test, 
                               bracket.size = base_line_size,tip.length = 0.01,
                               size = base_size/.pt, comparisons = comparisons)
  p = p + theme(axis.text.y = element_text(hjust = .5))
  p
}

#heatmap
heatmap = function(df, x, y, fill, chigh = "red", clow = "blue", cmed = "black", lim = c(-3,3), 
                   breaks = c(-3,0,3), xlab = "", ylab = "", title = "", divergent = T,
                   leglab ="", leg = c("X1/4", "Median", "X4")) {
  df[,x] = factor(df[,x], levels = unique(df[,x]))
  df[,y] = factor(df[,y], levels = unique(df[,y]))
  p = ggplot(df, aes(.data[[x]],.data[[y]], fill = .data[[fill]]))
  p = p  + geom_raster() 
  if (divergent)
    p = p + scale_fill_gradient2(limits = lim, low = clow, high = chigh, mid = cmed, 
                               breaks = breaks, labels = leg)
  else
    p = p + scale_fill_gradient(limits = lim, low = clow, high = chigh, breaks = breaks, labels = leg)
  p = p + theme(axis.line = element_blank())
  p = p + labs(x = xlab, y = ylab, fill = leglab, title = title)
  p
}

# input: bigwig track.  output: data.frame of binned signal 
bin_signal = function(bw.path, genomic_range, chr, position, gene, bin.size = 50) {
  bw = import(bw.path, selection = genomic_range, as="RleList")
  data.bw = as.data.frame(bw[[chr]][position])
  data.bw$pos = position
  n.bins = ceiling(nrow(data.bw)/bin.size) 
  cov.b = aggregate(data.bw$value, list(cut(data.bw$pos, n.bins, include.lowest = T)), FUN = mean)
  pos.b = aggregate(data.bw$pos, by = list(cut(data.bw$pos, n.bins, include.lowest = T)), 
                    FUN = function(x) median(x))
  df = data.frame(pos = pos.b, cov = cov.b)
  df$gene = gene
  df
}

# input: genomic ranges of exons
plotTranscript = function(gene_name ,coord_start, coord_end, gene_win, 
                          color = "black", strand_color = "grey") {
  df = data.frame()
  for (i in 1:length(gene_name)) { 
    g = gene_name[i]
    gw = gene_win[i]
    xmin = coord_start[i]
    xmax = coord_end[i]
    gene_ind = grep(paste0("^", g, "$"), genes$name2)
    if (length(gene_ind) > 1) {
      catn(g, "appears multipul times in annotations. skipping")
      next
    }
    exon_starts = as.numeric(unlist(strsplit(genes[gene_ind, "exonStarts"], ",")))
    exon_ends = as.numeric(unlist(strsplit(genes[gene_ind, "exonEnds"], ",")))
    strand = genes[gene_ind, "strand"]
    if (strand == "-") {
      gene_s = max(exon_ends)
      gene_e = min(exon_starts)
    } else if (strand == "+") {
      gene_s = min(exon_starts)
      gene_e = max(exon_ends)
    }
      
    chr = genes[gene_ind, "chrom"]
    df = rbind(df, data.frame(s = exon_starts, e = exon_ends, chr = chr, gene = g, 
                              xmin = xmin, xmax = xmax, gene_s = gene_s, gene_w = gw, 
                              gene_e = gene_e, starnd = strand))
  }
  df$gene = factor(df$gene, levels = gene_name)
  p = ggplot(df, aes(x = s, y = 1, xend = e, yend = 1, label = gene_w)) 
  p = p + geom_segment(aes(x = gene_s, xend = gene_e), size = 3, color = strand_color, 
                       arrow = arrow(angle = 20, type = "closed", length = unit(2, "mm")),
                       lineend = "butt", linejoin = "mitre")
  p = p + geom_segment(aes(x = s, xend = e), size = 2.5, color = color)
  # p = p + geom_rect(aes(xmin = s, xmax = e), fill = "black")
  
  # add empty white dots to set x scale 
  p = p + geom_point(aes(x = xmin), color = "transparent")
  p = p + geom_point(aes(x = xmax), color = "transparent")
  p = p + geom_line(linewidth = base_line_size, color = color) 
  p = p + labs(x = "", y = "") 
  p = p + scale_x_continuous(breaks = df$xmax - (df$xmax - df$xmin)/2, 
                             labels = df$gene_w)
  # p + geom_text(data = df[as.numeric(rownames(unique(df[,c("gene", "xmin", "xmax", "gene_w")]))),], 
  #                   mapping = aes(x = s, y = .9), size = 5/.pt, hjust = .5, inherit.aes = T) 
  p = p + facet_wrap(~ gene, scale = "free_x", nrow = 1)
  p = p + theme(axis.line = element_blank(), axis.ticks = element_blank(), 
                axis.text.y = element_blank(), strip.background = element_blank(), 
                panel.spacing.x = unit(4, "mm"))
  p = p + lims(y = c(0.8,1.1)) 
  p 
  
}

# input: genomic ranges - data frame with fields: name, chr, start, end, win
#        bw.data - data frame with fields: sample_id, color_code, bw_path
# https://rockefelleruniversity.github.io/Bioconductor_Introduction/exercises/answers/GS_answers.html
plotBrowser = function(bw.data, genomic_ranges, scale_tracks = F, scale_by = 1) {
  # 'autoscale-like'.  set y limit to be the maximum of all samples. scale_by - index of gene to scale y axis by.
  maxY = sapply(1:nrow(bw.data), function(i) 
    max(max(import(bw.data$bw_path[i], selection = GRanges(genomic_ranges$win[scale_by]), as="RleList"))))
  pl = list()
  p_gene = plotTranscript(genomic_ranges$name, genomic_ranges$start, genomic_ranges$end, 
                          genomic_ranges$win, color = "black", strand_color = "gray") 
  p_gene = p_gene + theme(plot.margin = unit(c(u=0,d=-8,r=0,l=0), "mm"))
  
  pl[[1]] = p_gene
  for (i in 1:nrow(bw.data)){
    bw.path = bw.data$bw_path[i]
    cl = bw.data$color_code[i]
    samp = bw.data$sample_id[i]
    cat("creating plot of ", samp, "\n")
    df = sapply(1:nrow(genomic_ranges), function(i) bin_signal(bw.path, GRanges(genomic_ranges$win[i]), genomic_ranges$chr[i], 
                                                   as.numeric(seq(genomic_ranges$start[i], genomic_ranges$end[i])), 
                                                   genomic_ranges$name[i]))
    data.cov = data.frame(pos = unlist(sapply(1:ncol(df), function(i) rbind(df[2,i]))),
                          cov = unlist(sapply(1:ncol(df), function(i) rbind(df[4,i]))),
                          gene = unlist(sapply(1:ncol(df), function(i) rbind(df[5,i]))))
    data.cov$gene = factor(data.cov$gene,levels = unique(data.cov$gene))
    p = ggplot(data.cov, aes(pos, cov)) 
    p = p + labs(x = "", y = samp) + guides(fill = "none")
    p = p + theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), 
                  strip.text.x = element_blank(), axis.ticks.y = element_blank(), 
                  axis.text.y = element_blank(), axis.line.y = element_blank(), 
                  axis.title.y = element_text(angle = 0, vjust = .5, size = 5), 
                  strip.background = element_blank(), strip.placement = "outside", 
                  plot.margin = unit(c(u=-2,d=-2,r=0,l=0), "mm"), 
                  panel.spacing.x = unit(4, "mm"))
    p = p + facet_wrap(~gene, scale = "free_x", strip.position = "bottom", nrow = 1)
    # if (i == 1) # empty plot for the purpose of padding top of plot  
    #   pl[[1]] = p + geom_area(fill = "white", color = "white") + ylab("") + 
    #               theme(axis.line.x = element_line(color = "white"))
    p = p + geom_area(fill = cl) # stat_smooth(geom = 'area')
    if (scale_tracks) 
      p = p + ylim(NA, maxY)
    if (i == nrow(bw.data)) 
      p = p + theme(plot.margin = unit(c(u=-2,d=0,r=0,l=0), "mm"))
    pl[[i+1]] = p
  }
  p = plot_grid(plotlist = pl, ncol = 1, align = "v", axis = "rl", rel_heights = c(4, rep(1, nrow(bw.data))))
  p 
}

# data - data frame with fields: group_by, rank_by
# group_by - names of field to group by 
# rank_by - name of field to rank by
rankSubgroups = function(data, group_by, rank_by, 
                         decreasing = F, group_order = NULL, group_by_sec = NULL) {
  ord = order(data[[rank_by]], decreasing = decreasing)
  data = data[ord, ]
  if (is.null(group_order))
    group_order = levels(factor(data[[group_by]], exclude = NULL))
  rank = unlist(sapply(group_order, function(g) which(mapply(`%in%`, data[[group_by]],g))), 
                use.names = F)
  if (!is.null(group_by_sec)) {
    group_order_sec = levels(factor(data[[group_by_sec]]))
    group.pairs = expand.grid(group_order_sec,group_order)
    group.pairs$Var1 = factor(group.pairs$Var1, exclude = F)
    group.pairs$Var2 = factor(group.pairs$Var2, exclude = F)
    rank = unlist(sapply(1:nrow(group.pairs), function(g) { 
        # comparison that works also for NA (based on https://stackoverflow.com/questions/37610056/how-to-treat-nas-like-values-when-comparing-elementwise-in-r)
        which(mapply(`%in%`, data[[group_by_sec]], group.pairs[g,1]) & 
                mapply(`%in%`, data[[group_by]], group.pairs[g,2]))
      }), use.names = F)
    }
  # na.ind = which(is.na(data[[group_by]]))
  # return(ord[c(rank, na.ind)])
  ord[rank]
}

# data - data frame of data to display as heatmap
# groups - data frame with field 'group'. rownames should be same as colnames of data. 
# group_by - group in groups to group by
# group_order - optional order of group levels 
clusterSubgroups = function(data, groups, group_by, group_order = NULL) {
  groups = groups[colnames(data),, drop = F]
  if (any(is.null(group_order))) # patch
    group_order = levels(factor(groups[[group_by]]))
  ord = c()
  for (g in group_order) {
    group.samp = rownames(groups)[which(groups[[group_by]] == g)]
    group.ind = match(group.samp, colnames(data))
    if (length(group.samp) > 1){
      group.ord = hclust(dist(t(data[,group.samp])))$order
      ord = c(ord, group.ind[group.ord])
    } else {
      ord = c(ord, group.ind)
    }
  }
  ind.na = which(is.na(groups[[group_by]]))
  if (length(ind.na) > 0){
    if (length(ind.na) == 1) {
      ord = c(ord, ind.na)
    }
    else {
      ord = c(ord, ind.na[hclust(dist(t(data[,ind.na])))$order])
    }
  }
  return(ord)
}



plotTranscript_progress = function(gene_name ,genomic_windows, coord_start, coord_end, color = "black", 
                          strand_color = "grey", chromosomes) {
  df = data.frame()
  for (i in 1:length(chromosomes)) { 
    print(i)
    # gene_ind = grep(paste0("^", g, "$"), genes$name2)
    # g = gene_name[i]
    xmin = coord_start[i]
    xmax = coord_end[i]
    gene_ind = which(genes$chrom == chromosomes[i] & genes$txStart > xmin & genes$txEnd < xmax)
    gene_names = genes$name2[gene_ind]
    print(gene_names)
    # exon_starts = as.numeric(unlist(strsplit(genes[gene_ind, "exonStarts"], ",")))
    # exon_ends = as.numeric(unlist(strsplit(genes[gene_ind, "exonEnds"], ",")))
    exon_starts = strsplit(genes[gene_ind, "exonStarts"], ",")
    exon_ends = strsplit(genes[gene_ind, "exonEnds"], ",")
    names(exon_starts) = gene_names
    names(exon_ends) = gene_names
    strands = genes[gene_ind, "strand"]
    gene_s = list()
    gene_e = list()
    for (s in 1:length(strands)) {
      gene.name = gene_names[s]
      strand = strands[s]
      if (strand == "-") {
        gene_s[[gene.name]] = max(as.numeric(exon_starts[[gene.name]]))
        gene_e[[gene.name]] = min(as.numeric(exon_ends[[gene.name]]))
      } else if (strand == "+") {
        gene_s[[gene.name]] = max(as.numeric(exon_ends[[gene.name]]))
        gene_e[[gene.name]] = min(as.numeric(exon_starts[[gene.name]]))
      } 
    }
    # chr = genes[gene_ind, "chrom"]
    
    df = rbind(df, data.frame(s = as.numeric(unlist(exon_starts)), 
                              e = as.numeric(unlist(exon_ends)), 
                              gene.name = gene.names[i]))
  }
  
  p = ggplot(df, aes(x = s, y = 1, xend = e, yend = 1)) 
  p = p + geom_segment(linewidth = 2.5, color = color)
  
  p = p + geom_segment(aes(x = gene_s, xend = gene_e), size = 3, color = strand_color, 
                       arrow = arrow(angle = 20, type = "closed", length = unit(2, "mm")),
                       lineend = "butt", linejoin = "mitre")
  
  # xmin = xmin, xmax = xmax, gene_s = gene_s, 
  # gene_e = gene_e, starnd = strand
  
  
  
  df$gene = factor(df$gene, levels = gene_name)
  
  # loop over genes 
  p = p + geom_segment(aes(x = s, xend = e), size = 2.5, color = color)
  # p = p + geom_rect(aes(xmin = s, xmax = e), fill = "black")
  
  # add empty white dots to set x scale 
  p = p + geom_point(aes(x = xmin), color = "transparent")
  p = p + geom_point(aes(x = xmax), color = "transparent")
  p = p + geom_line(size = base_line_size, color = color) 
  p = p + labs(x = "", y = "") + ylim(c(0.5,1.3))
  p = p + facet_wrap(~ gene, scale = "free_x", nrow = 1)
  p = p + theme(axis.line = element_blank(), axis.ticks = element_blank(), 
                axis.text = element_blank(), strip.background = element_blank(), 
                panel.spacing.x = unit(4, "mm"))
  p
  
}



# x is data.frame / matrix with rownames = genes, colnames = samples
# you don't have to cluster rows / columns (see args)
# keep.hclust returns the hclust object for row and column
hclust2treeviewML <- function (x, file = "cluster.cdt", method = "euclidean", link = "complete", 
                               keep.hclust = FALSE, cluster_row = T,
                               cluster_col =T) 
{
  if(is.matrix(x)){x <- as.data.frame(x)}
  pkg <- loadedNamespaces()   
  if (!"ctc" %in% pkg){library(ctc)} # calls the package if it's missing
  if ("amap" %in% .packages(all = TRUE)) { # beter hclust
    hr <- hcluster(x, method = method, link = link)
    if (!cluster_row){ hr$order = 1:length(hr$labels)}
    hc <- hcluster(t(x), method = method, link = link)
    if (!cluster_col){ hc$order = 1:length(hc$labels)}
  }
  else { # if dependency is missing use base hclust
    hr <- hclust(dist(x, method = method), method = link)
    if (!cluster_row){ hr$order = 1:length(hr$labels)}
    hc <- hclust(dist(t(x), method = method), method = link)
    if (!cluster_col){ hc$order = 1:length(hc$labels)}
  }
  basefile = strsplit(file, "\\.")[[1]]
  if (length(basefile > 1)) 
    basefile = paste(basefile[-length(basefile)], collapse = ".")
  # this is where the magic happens
  r2atr(hc, file = paste(basefile, ".atr", sep = ""))
  r2gtr(hr, file = paste(basefile, ".gtr", sep = ""))
  r2cdt(hr, hc, cbind(rownames(x),x),
        file = paste(basefile, ".cdt", sep = ""),labels = T)
  if (keep.hclust) 
    return(list(hr, hc))
  else return(1)
}


enricher <- function(genes,
                     output.path="./",
                     outputfilename="enrichr.output.csv",
                     threshold=0.05){
  setEnrichrSite("Enrichr")
  dbs.we.want <-c("Reactome_2022" ,"BioPlanet_2019","Cancer_Cell_Line_Encyclopedia","Human_Gene_Atlas",
                  "ARCHS4_Tissues", "ARCHS4_Cell-lines","PanglaoDB_Augmented_2021", "HuBMAP_ASCTplusB_augmented_2022")
  en <- do.call(rbind.data.frame,enrichr(genes,dbs.we.want))
  print(dim(en))
  print(threshold)
  en <- en[en$Adjusted.P.value < threshold, ]
  print(dim(en))
  cat("found", nrow(en), "significant terms \n")
  cat("writing results to: ", paste0(output.path, outputfilename), "\n")
  write.csv(en, paste0(output.path,outputfilename))
}
