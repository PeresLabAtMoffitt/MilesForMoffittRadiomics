################################################################################################# I ### Survivals----
radiomics_surv <- radiomics %>% distinct(patient_id, .keep_all = TRUE)
# From Dx
# overall survival by 4 trt groups
ggsurvplot(survfit(Surv(months_at_dx_followup, os_event) ~ treatment_type, data=radiomics_surv),
           title = "overall survival from date of diagnosis", 
           font.main = c(24, "bold", "black"), 
           font.x = c(20, "bold", "black"), font.y = c(20, "bold", "black"), 
           font.legend = c(20, "bold", "black"), font.tickslab = c(18, "bold", "black"),
           size = 1.5,
           
           pval=TRUE, palette = c("black", "green","blue", "red"),
           # Add risk table
           risk.table = "abs_pct", risk.table.title = "Risk table (number(%))",
           risk.table.y.text = FALSE, risk.table.fontsize = 6,tables.height = 0.3,
           tables.theme = theme_survminer(base_size = 5,
                                          font.main = c(16, "bold", "black"),
                                          font.x = c(16, "bold", "black"),
                                          font.y = c(16, "bold", "transparent"),
                                          font.tickslab = c(19, "bold", "black")
           ),
           
           legend.title = "", xlab = "Time (in months)", ylab="Survival Probability"  
           # # legend.labs=c("chemo only", "surgery only", "upfront neoadjuvant", "upfront surgery")
           )


# color = "HCT_at_all_time",
# linetype = "CH_status",
# pval = TRUE,
# conf.int = FALSE,


# recurrence free survival by 4 trt groups
ggsurvplot(survfit(Surv(months_of_dx_rec_free, rec_event) ~ treatment_type, data=radiomics_surv),
           title = "recurrence free survival from date of diagnosis", 
           font.main = c(24, "bold", "black"), 
           font.x = c(20, "bold", "black"), font.y = c(20, "bold", "black"), 
           font.legend = c(20, "bold", "black"), font.tickslab = c(18, "bold", "black"), 
           size = 1.5,
           
           pval=TRUE, palette = c("black", "green","blue", "red"), 
           # Add risk table
           risk.table = "abs_pct", risk.table.title = "Risk table (number(%))",
           risk.table.y.text = FALSE, risk.table.fontsize = 6,tables.height = 0.3,
           tables.theme = theme_survminer(base_size = 5,
                                          font.main = c(16, "bold", "black"),
                                          font.x = c(16, "bold", "black"),
                                          font.y = c(16, "bold", "transparent"),
                                          font.tickslab = c(19, "bold", "black")
           ),
           
           legend.title = "", xlab = "Time (in months)", ylab="Recurrence Probability"  
           # legend.labs=c("chemo only", "surgery only", "upfront neoadjuvant", "upfront surgery")
           )

# From surgery
# overall survival calculated from date of upfront surgery 
ggsurvplot(survfit(Surv(months_at_surg_followup, os_event) ~ treatment_type, data=radiomics_surv), 
           title = "overall survival from date of upfront surgery", 
           font.main = c(24, "bold", "black"), 
           font.x = c(20, "bold", "black"), font.y = c(20, "bold", "black"), 
           font.legend = c(20, "bold", "black"), font.tickslab = c(18, "bold", "black"), 
           size = 1.5,
           
           pval=TRUE, palette = c("black", "green","blue", "red"), 
           # Add risk table
           risk.table = "abs_pct", risk.table.title = "Risk table (number(%))",
           risk.table.y.text = FALSE, risk.table.fontsize = 6,tables.height = 0.3,
           tables.theme = theme_survminer(base_size = 5,
                                          font.main = c(16, "bold", "black"),
                                          font.x = c(16, "bold", "black"),
                                          font.y = c(16, "bold", "transparent"),
                                          font.tickslab = c(19, "bold", "black")
           ),
           
           legend.title = "", xlab = "Time (in months)", ylab="Survival Probability"#,  
           # # legend.labs=c("chemo only", "surgery only", "upfront neoadjuvant", "upfront surgery")
           )

# recurrence free survival calculated from date of upfront surgery 
ggsurvplot(survfit(Surv(months_of_surg_rec_free, rec_event) ~ treatment_type, data=radiomics_surv),
           title = "recurrence free survival from date of upfront surgery ", 
           font.main = c(24, "bold", "black"), 
           font.x = c(20, "bold", "black"), font.y = c(20, "bold", "black"), 
           font.legend = c(20, "bold", "black"), font.tickslab = c(18, "bold", "black"), 
           size = 1.5,
           
           pval=TRUE, palette = c("black", "green","blue", "red"), 
           # Add risk table
           risk.table = "abs_pct", risk.table.title = "Risk table (number(%))",
           risk.table.y.text = FALSE, risk.table.fontsize = 6,tables.height = 0.3,
           tables.theme = theme_survminer(base_size = 5,
                                          font.main = c(16, "bold", "black"),
                                          font.x = c(16, "bold", "black"),
                                          font.y = c(16, "bold", "transparent"),
                                          font.tickslab = c(19, "bold", "black")
           ),
           
           legend.title = "", xlab = "Time (in months)", ylab="Recurrence Probability"#,  
           # # legend.labs=c("chemo only", "surgery only", "upfront neoadjuvant", "upfront surgery")
)

# From neoadjuvant
# overall survival calculated from date of neoadjuvant 
ggsurvplot(survfit(Surv(months_at_neo_followup, os_event) ~ treatment_type, data=radiomics_surv),
           title = "overall survival from date of neoadjuvant", 
           font.main = c(24, "bold", "black"), 
           font.x = c(20, "bold", "black"), font.y = c(20, "bold", "black"), 
           font.legend = c(20, "bold", "black"), font.tickslab = c(18, "bold", "black"), 
           size = 1.5,
           
           pval=TRUE, palette = c("black", "green","blue", "red"), 
           # Add risk table
           risk.table = "abs_pct", risk.table.title = "Risk table (number(%))",
           risk.table.y.text = FALSE, risk.table.fontsize = 6,tables.height = 0.3,
           tables.theme = theme_survminer(base_size = 5,
                                          font.main = c(16, "bold", "black"),
                                          font.x = c(16, "bold", "black"),
                                          font.y = c(16, "bold", "transparent"),
                                          font.tickslab = c(19, "bold", "black")
           ),
           
           legend.title = "", xlab = "Time (in months)", ylab="Survival Probability"#,  
           # # legend.labs=c("chemo only", "surgery only", "upfront neoadjuvant", "upfront surgery")
)

# recurrence free survival calculated from date of neoadjuvant 
ggsurvplot(survfit(Surv(months_of_neo_rec_free, rec_event) ~ treatment_type, data=radiomics_surv),
           title = "recurrence free survival from date of neoadjuvant", 
           font.main = c(24, "bold", "black"), 
           font.x = c(20, "bold", "black"), font.y = c(20, "bold", "black"), 
           font.legend = c(20, "bold", "black"), font.tickslab = c(18, "bold", "black"),
           size = 1.5,
           
           pval=TRUE, palette = c("black", "green","blue", "red"), 
           # Add risk table
           risk.table = "abs_pct", risk.table.title = "Risk table (number(%))",
           risk.table.y.text = FALSE, risk.table.fontsize = 6,tables.height = 0.3,
           tables.theme = theme_survminer(base_size = 5,
                                          font.main = c(16, "bold", "black"),
                                          font.x = c(16, "bold", "black"),
                                          font.y = c(16, "bold", "transparent"),
                                          font.tickslab = c(19, "bold", "black")
           ),
           
           legend.title = "", xlab = "Time (in months)", ylab="Recurrence Probability"#,  
           # # legend.labs=c("chemo only", "surgery only", "upfront neoadjuvant", "upfront surgery")
)

# From treatment
# overall survival calculated from date of treatment 
ggsurvplot(survfit(Surv(months_at_treat_followup, os_event) ~ treatment_type, data=radiomics_surv),
           title = "overall survival from date of treatment", 
           font.main = c(24, "bold", "black"), 
           font.x = c(20, "bold", "black"), font.y = c(20, "bold", "black"), 
           font.legend = c(20, "bold", "black"), font.tickslab = c(18, "bold", "black"), 
           size = 1.5,
           
           pval=TRUE, palette = c("black", "green","blue", "red"), 
           # Add risk table
           risk.table = "abs_pct", risk.table.title = "Risk table (number(%))",
           risk.table.y.text = FALSE, risk.table.fontsize = 6,tables.height = 0.3,
           tables.theme = theme_survminer(base_size = 5,
                                          font.main = c(16, "bold", "black"),
                                          font.x = c(16, "bold", "black"),
                                          font.y = c(16, "bold", "transparent"),
                                          font.tickslab = c(19, "bold", "black")
           ),
           
           legend.title = "", xlab = "Time (in months)", ylab="Survival Probability"#,  
           # # legend.labs=c("chemo only", "surgery only", "upfront neoadjuvant", "upfront surgery")
)

# recurrence free survival calculated from date of treatment 
ggsurvplot(survfit(Surv(months_of_treat_rec_free, rec_event) ~ treatment_type, data=radiomics_surv),
           title = "recurrence free survival from date of treatment", 
           font.main = c(24, "bold", "black"), 
           font.x = c(20, "bold", "black"), font.y = c(20, "bold", "black"), 
           font.legend = c(20, "bold", "black"), font.tickslab = c(18, "bold", "black"), 
           size = 1.5,
           
           pval=TRUE, palette = c("black", "green","blue", "red"), 
           # Add risk table
           risk.table = "abs_pct", risk.table.title = "Risk table (number(%))",
           risk.table.y.text = FALSE, risk.table.fontsize = 6,tables.height = 0.3,
           tables.theme = theme_survminer(base_size = 5,
                                          font.main = c(16, "bold", "black"),
                                          font.x = c(16, "bold", "black"),
                                          font.y = c(16, "bold", "transparent"),
                                          font.tickslab = c(19, "bold", "black")
           ),
           
           legend.title = "", xlab = "Time (in months)", ylab="Recurrence Probability"#,  
           # # legend.labs=c("chemo only", "surgery only", "upfront neoadjuvant", "upfront surgery")
)



# In stage3/4 group
restricted_stage_trt_type_ <- radiomics_surv %>% filter(str_detect(tnm_stage, "3|4"))
# overall survival in stage3/4 group of patients by upfront noeadj vs. upfront surgery groups
ggsurvplot(survfit(Surv(months_at_dx_followup, os_event) ~ treatment_type, data=restricted_stage_trt_type_),
           title = "overall survival from date of diagnosis",
           font.main = c(24, "bold", "black"),
           font.x = c(20, "bold", "black"), font.y = c(20, "bold", "black"),
           font.legend = c(20, "bold", "black"), font.tickslab = c(18, "bold", "black"),
           size = 1.5,
           
           pval=TRUE, palette = c("black", "green","blue", "red"), 
           # Add risk table
           risk.table = "abs_pct", risk.table.title = "Risk table (number(%))",
           risk.table.y.text = FALSE, risk.table.fontsize = 6,tables.height = 0.3,
           tables.theme = theme_survminer(base_size = 5,
                                          font.main = c(16, "bold", "black"),
                                          font.x = c(16, "bold", "black"),
                                          font.y = c(16, "bold", "transparent"),
                                          font.tickslab = c(19, "bold", "black")
           ),
           
           legend.title = "", xlab = "Time (in months)", ylab="Survival Probability"#,  
           # # legend.labs=c("upfront neoadjuvant", "upfront surgery")
           )

# recurrence free survival in stage3/4 group of patients by upfront noeadj vs. upfront surgery groups
ggsurvplot(survfit(Surv(months_of_dx_rec_free, rec_event) ~ treatment_type, data=restricted_stage_trt_type_),
           title = "recurrence free survival from date of diagnosis",
           font.main = c(24, "bold", "black"),
           font.x = c(20, "bold", "black"), font.y = c(20, "bold", "black"),
           font.legend = c(20, "bold", "black"), font.tickslab = c(18, "bold", "black"),
           size = 1.5,
           
           pval=TRUE, palette = c("black", "green","blue", "red"), 
           # Add risk table
           risk.table = "abs_pct", risk.table.title = "Risk table (number(%))",
           risk.table.y.text = FALSE, risk.table.fontsize = 6,tables.height = 0.3,
           tables.theme = theme_survminer(base_size = 5,
                                          font.main = c(16, "bold", "black"),
                                          font.x = c(16, "bold", "black"),
                                          font.y = c(16, "bold", "transparent"),
                                          font.tickslab = c(19, "bold", "black")
           ),
           
           legend.title = "", xlab = "Time (in months)", ylab="Recurrence Probability"#,  
           # # legend.labs=c("upfront neoadjuvant", "upfront surgery")
           )

