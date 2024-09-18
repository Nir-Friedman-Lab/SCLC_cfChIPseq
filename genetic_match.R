read.csv("~/Downloads/rLODs.csv") %>% 
  filter(Count > 50) -> rlod.table

# first in row is genetically identical to seconed in row
false.negative = c("SCLC0315", "SCLC0021",
                   "SCLC0205", "SCLC0061", 
                   "SCLC0312", "SCLC0026", 
                   "SCLC0201", "SCLC0040", 
                   "SCLC0327", "SCLC0010", 
                   "SCLC0351", "SCLC0062", 
                   "SCLC0351", "SCLC0058", 
                   "SCLC0353", "SCLC0256", 
                   "SCLC0319-173", "SCLC0014", 
                   "SCLC0035-653", "SCLC0037", 
                   "SCLC00311", "SCLC0027", 
                   "SCLC0344", "SCLC0251", 
                   "SCLC0194", "SCLC0024", 
                   "SCLC0346", "SCLC0342", 
                   "SCLC0289", "SCLC0066")

# false negative 
rlod.table %>% 
  # filter(rLOD > .75, SameID == F, Type == "mixed") %>%
  filter(rLOD > .75, SameID == F, Type == "ChIP") %>%
  filter(Count > 200) %>%
  filter(!grepl(paste(false.negative, collapse="|"), row),
         !grepl(paste(false.negative, collapse="|"), col)) %>%
  arrange(-rLOD) %>% 
  head(50)

s = "SCLC0261"
s = "SCLC0005-24"
rlod.table %>% 
  filter(grepl(s, row) | grepl(s, col)) %>% 
  filter(Count > 200) %>% 
  arrange(-rLOD) %>% 
  head(20)

false.positive = c("SCLC0126-223", "SCLC0035-653", "sb_19_6623_R1", "SCLC0114-351", 
                   "SCLC0150-669", "19_2158", "19_5931")
ignore = c(false.positive, "2698_S13", "CL0196", "CL0124", "CL0186")
# SCLC0126-223 - Genetic dissimilar to other samples from the same patient (and any other sample)
# SCLC0035-653 - Genetic match to samples from NCI037, mismatch with other samples from NCI035
# sb_19_6623_R1 - Genetic match to samples from NCI037, mismatch with other samples from NCI035 ?
# SCLC0114-351_NA_H3K4me3_NA - Genetic dissimilar to CL0114_T1R_T_R1 (and any other sample)
# SCLC0150-669_NA_H3K4me3-39_26032021-13 - Weak genetic similarity to other samples from the same patient
# SS_19_2158_R1 - Genetic dissimilar to other samples from the same patient (and any other sample)
# CL0114_T1R_T_R - same
# sb_19_5931_R1 - Genetic match to samples (RNA and ChIP) of NCI035
# CL0124_T1R_T_R1 - week genetic match to matching ChIP but seems real
# CL0196 - week geneitc match to some of the corresponding ChIP but seems real
# CL0124/CL0186/NCI0422 - same
# 
 
# false positive
rlod.table %>% 
  filter(rLOD < .75, SameID == T) %>% 
  filter(!grepl(paste(ignore, collapse="|"), row),
         !grepl(paste(ignore, collapse="|"), col)) %>%
  arrange(rLOD) %>% 
  head(20)






rlod.table %>% 
  ggplot(aes(Count, rLOD, color = SameID)) +
  geom_point(size = .3) +
  geom_point(data = rlod.table %>% filter((SameID == F & rLOD > .75) |
                                            (SameID == T & rLOD < .75)), size = 1) +
  scale_color_aaas() +
  facet_wrap(~Type) +
  scale_x_log10()


discarde.samples = c("SCLC0035-653_NA_H3K4me3-38_26032021-13", 
                     "SCLC0126-223_NAH3K4me3-38_26032021-13", 
                     "SCLC0150-669_NA_H3K4me3-39_26032021-13",
                     "SS_19_2158_R1")

# SS_19_2158_R1 - Genetic dissimilar to other samples from the same patient (and any other sample)