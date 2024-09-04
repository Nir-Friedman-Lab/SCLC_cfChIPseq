# Kaplan Meier curve
library(survival)
library(survminer)
# library(ggfortify)

metadata %>% 
  filter(Sample_id %in% pre.samples) %>% 
  group_by(PatientID) %>% 
  slice_min(order_by = sampling_date) %>% 
  ungroup() %>%
  left_join(metadata.patients %>% select(PatientID, `Mortality Status`, 
                                         `Date of Death`), 
            by = "PatientID") %>% 
  rename(survival = sampling_death_duration, 
         mortal = `Mortality Status`)  %>% 
    filter(!is.na(mortal)) -> data.surv 

data.surv %>% 
  # filter(!is.na(mortal)) %>% 
  pull(SCLC.n) %>% 
  median() -> surv.cutoff

data.surv %>% 
  # mutate(mortal = if_else(!is.na(`Date of Death`), 1, mortal)) %>% # TODO patch till this is fixed in metadata
  mutate(above_median = SCLC.n > surv.cutoff, 
         survival = as.integer(survival)) -> data.surv
  
data.surv %>% 
  count(above_median)
  # filter(`cfDNA Patient identical/duplicate` == "Identical" & # pre-treatment samples
  #          Timepoint == "Pre-treatment" & Sample_id %in% s.samples) %>% 
  # select(c(Sample_id, SCLC_score, sampling_death_duration, 
  #          `Mortality Status`)) %>%
  # mutate(score = SCLC_score, survival = sampling_death_duration, 
  #        mortal = `Mortality Status`) -> data.surv 
# cutoff = median(data.surv$SCLC.n)

# data.surv = metadata[,c("SCLC_score", "sampling_death_duration", "Mortality Status")]
# names(data.surv) = c("score", "survival", "mortal")
# data.surv$survival = as.numeric(data.surv$survival)
# data.surv = data.surv[!is.na(data.surv$score) & !(is.na(data.surv$survival)),]
# data.surv = data.surv[pre.samples,]
# cutoff = median(metadata[rownames(data.surv), "SCLC_score"], na.rm = T)
# data.surv$above_median = data.surv$score > cutoff
y = Surv(time = data.surv$survival
         , event = data.surv$mortal
         )
# y = Surv(time = data.surv$survival)
kmfit = survfit(y ~ data.surv$above_median, data = data.surv)
pval = round(survdiff(y ~ data.surv$above_median)$pvalue, digits = 4)
km = summary(kmfit) # , scale = 30 #survival in month 
km.data = data.frame(time = km$time, survival = km$surv, 
                     group = km$strata, n.risk = km$strata)
tit = paste0("p < ", pval)
levels(km.data$group) = c("< median", "> median")
# p = ggplot(km.data, aes(time, 100*survival))
# p = p + geom_step(aes(color = group)) 
# p = p + geom_point(aes(color = group), size = .2, show.legend = F)
# p = p + scale_color_aaas()
# p = p + labs(x = "survival from diagnosis (months)", y = "Percent survival", 
#              color = "SCLC-score", title = tit)
# p = p + theme(legend.position = c(.7,.7), legend.key.size = unit(2, 'mm'),  
#           plot.title = element_text(hjust = 0.5))
# p = p + geom_text(data = data.frame(x = 15, y = 1, text = paste("p < ", pval)), 
#                   aes(x,y,label = text), size = base_size/.pt, fontface = "italic")
# ggsave(paste0(figDirPaper, "figure1/km_curve.pdf"), p, width = 55, height = 55, units = "mm")
# from https://rpkgs.datanovia.com/survminer/

pdf(paste0(figDirPaper, "figure1/km_curve_v1.pdf"), width = 2.1, height = 2.1, 
    pointsize = base_size)
p = ggsurvplot(fit = kmfit, data = data.surv, pval = T, risk.table = TRUE, 
               ggtheme = base_theme, tables.theme = base_theme, 
               legend.labs = c("< median", "> median"), size = .5,
               xlab = "survival from diagnosis (months)", censor.shape = NA,
               palette = c(group.colors$SCLC, "#fdae6b"), 
               pval.coord = c(20, .9),
               # surv.median.line = "hv",     
               risk.table.y.text = FALSE, tables.height = .3,
               fontsize = base_size/.pt, pval.size = base_size/.pt)
p
print(p, newpage = F)
dev.off()
