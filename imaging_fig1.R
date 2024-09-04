# prepare data ------------------------------------------------------------
patient.id = "CL0191"
patient.ind = which(metadata$Patient_id == patient.id)
plasma.samples = metadata$Sample_id[patient.ind]
plasme.samp.dates = as.Date(metadata$sampling_date[patient.ind])
diagnosis.date = as.Date("2017-10-12")
time.from.diag = round(as.numeric(difftime(plasme.samp.dates, diagnosis.date, units = "days") / 30))
death.date = as.Date("2019-10-18")


# treatment and response --------------------------------------------------
data.image = data.frame(time = time.from.diag , score = metadata$SCLC_score[patient.ind])
data.treatment = data.frame(time = c(7,10,11,18), 
                            treatment = c("Tr. 1", "Tr. 2", "RT", "Tr. 3"), 
                            response = c("PD", "SD -> PD", NA,
                                         "CR"))

p = ggplot(data.image, aes(time, score)) + geom_line() + geom_point()
p = p + labs(x = "time from diagnosis (months)", y = "SCLC score")
p = p + scale_x_continuous(breaks = c(time.from.diag), 
                           labels = c(as.character(time.from.diag)), 
                           limits = c(6,18), position = "top")

p = p + geom_text(data = data.treatment, aes(x = rollmean(c(6,time),2), y = .9, label = treatment), 
                  size = base_size/.pt)
p = p + geom_text(data = data.treatment, aes(x = rollmean(c(6,time),2), y = 0.7, label = response), 
                  size = base_size/.pt)
p = p + geom_vline(data = data.treatment, aes(xintercept = time), size = base_line_size, linetype = "dashed") 
p = p + geom_hline(aes(yintercept = 0.6), size = base_line_size) 
p = p + geom_hline(aes(yintercept = 0.8), size = base_line_size) 
p = p + scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.7, .9), 
                                         labels = c("0", "0.25", "0.5", "response", "treatment"), 
                                         limits = c(0, 1))
p = p + theme(axis.title.y = element_text(margin = margin(r = -5)))
p
ggsave(paste0(figDirPaper, "figure1/imaging/imaging.pdf"), p, width = 60, height = 35, units = "mm")
  
