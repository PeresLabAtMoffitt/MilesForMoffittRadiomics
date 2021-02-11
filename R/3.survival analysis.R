################################################################################################# I ### Survivals----
radiomics_surv <- radiomics %>% distinct(mrn, .keep_all = TRUE)
# From Dx
# overall survival by 4 trt groups
ggsurvplot(survfit(Surv(months_at_dx_followup, os_event) ~ treatment_type, data=radiomics_surv),
           palette = c("black", "green","blue", "red"),
           pval=TRUE, risk.table=TRUE, font.x=c("bold"),font.y=c("bold"),
           legend.title = "", xlab = "Time (in months)", ylab="Survival Probability",  
           legend.labs=c("chemo only", "surgery only", "upfront neoadjuvant", "upfront surgery"))

# recurrence free survival by 4 trt groups
ggsurvplot(survfit(Surv(months_of_dx_rec_free, rec_event) ~ treatment_type, data=radiomics_surv),
           palette = c("black", "green","blue", "red"),
           pval=TRUE, risk.table=TRUE, font.x=c("bold"),font.y=c("bold"),
           legend.title = "", xlab = "Time (in months)", ylab="Recurrence Probability",  
           legend.labs=c("chemo only", "surgery only", "upfront neoadjuvant", "upfront surgery"))

# From surgery
# overall survival calculated from date of upfront surgery 
ggsurvplot(survfit(Surv(months_at_surg_followup, os_event) ~ treatment_type, data=radiomics_surv), 
           palette = c("black", "green","blue", "red"),
           pval=TRUE, risk.table=TRUE, font.x=c("bold"),font.y=c("bold"),
           legend.title = "", xlab = "Time (in months)", ylab="Survival Probability"#,  
           # legend.labs=c("chemo only", "surgery only", "upfront neoadjuvant", "upfront surgery")
           )

# recurrence free survival calculated from date of upfront surgery 
ggsurvplot(survfit(Surv(months_of_surg_rec_free, rec_event) ~ treatment_type, data=radiomics_surv),
           palette = c("black", "green","blue", "red"),
           pval=TRUE, risk.table=TRUE, font.x=c("bold"),font.y=c("bold"),
           legend.title = "", xlab = "Time (in months)", ylab="Recurrence Probability"#,  
           # legend.labs=c("chemo only", "surgery only", "upfront neoadjuvant", "upfront surgery")
)

# From neoadjuvant
# overall survival calculated from date of neoadjuvant 
ggsurvplot(survfit(Surv(months_at_neo_followup, os_event) ~ treatment_type, data=radiomics_surv),
           palette = c("black", "green","blue", "red"),
           pval=TRUE, risk.table=TRUE, font.x=c("bold"),font.y=c("bold"),
           legend.title = "", xlab = "Time (in months)", ylab="Survival Probability"#,  
           # legend.labs=c("chemo only", "surgery only", "upfront neoadjuvant", "upfront surgery")
)

# recurrence free survival calculated from date of neoadjuvant 
ggsurvplot(survfit(Surv(months_of_neo_rec_free, rec_event) ~ treatment_type, data=radiomics_surv),
           palette = c("black", "green","blue", "red"),
           pval=TRUE, risk.table=TRUE, font.x=c("bold"),font.y=c("bold"),
           legend.title = "", xlab = "Time (in months)", ylab="Recurrence Probability"#,  
           # legend.labs=c("chemo only", "surgery only", "upfront neoadjuvant", "upfront surgery")
)

# From treatment
# overall survival calculated from date of treatment 
ggsurvplot(survfit(Surv(months_at_treat_followup, os_event) ~ treatment_type, data=radiomics_surv),
           palette = c("black", "green","blue", "red"),
           pval=TRUE, risk.table=TRUE, font.x=c("bold"),font.y=c("bold"),
           legend.title = "", xlab = "Time (in months)", ylab="Survival Probability"#,  
           # legend.labs=c("chemo only", "surgery only", "upfront neoadjuvant", "upfront surgery")
)

# recurrence free survival calculated from date of treatment 
ggsurvplot(survfit(Surv(months_of_treat_rec_free, rec_event) ~ treatment_type, data=radiomics_surv),
           palette = c("black", "green","blue", "red"),
           pval=TRUE, risk.table=TRUE, font.x=c("bold"),font.y=c("bold"),
           legend.title = "", xlab = "Time (in months)", ylab="Recurrence Probability"#,  
           # legend.labs=c("chemo only", "surgery only", "upfront neoadjuvant", "upfront surgery")
)



# In stage3/4 group
restricted_stage_trt_type_ <- radiomics_surv %>% filter(str_detect(tnm_stage, "3|4"))
# overall survival in stage3/4 group of patients by upfront noeadj vs. upfront surgery groups
ggsurvplot(survfit(Surv(months_at_dx_followup, os_event) ~ treatment_type, data=restricted_stage_trt_type_),
           # palette = c("blue", "red"),
           pval=TRUE, risk.table=TRUE, font.x=c("bold"),font.y=c("bold"),
           legend.title = "", xlab = "Time (in months)", ylab="Survival Probability"#,  
           # legend.labs=c("upfront neoadjuvant", "upfront surgery")
           )

# recurrence free survival in stage3/4 group of patients by upfront noeadj vs. upfront surgery groups
ggsurvplot(survfit(Surv(months_of_dx_rec_free, rec_event) ~ treatment_type, data=restricted_stage_trt_type_),
           # palette = c("blue", "red"),
           pval=TRUE, risk.table=TRUE, font.x=c("bold"),font.y=c("bold"),
           legend.title = "", xlab = "Time (in months)", ylab="Recurrence Probability"#,  
           # legend.labs=c("upfront neoadjuvant", "upfront surgery")
           )

