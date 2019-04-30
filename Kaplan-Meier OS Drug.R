####Kaplan-Meier for OS By Drug
os_pemetrexed1 <- merge(pemetrexed, os, by = "sample")
os_pemetrexed <- na.omit(os_pemetrexed1)

os_pemetrexed_survival <- survfit(Surv(os_pemetrexed$X_TIME_TO_EVENT, os_pemetrexed$X_EVENT) ~ 1)
summary(os_pemetrexed_survival)                  
plot(os_pemetrexed_survival, main="OS on Pemetrexed", xlab="Months", ylab="Overall Survival")


os_paclitaxel1 <- merge(paclitaxel, os, by = "sample")
os_paclitaxel <- na.omit(os_paclitaxel1)

os_paclitaxel_survival <- survfit(Surv(os_paclitaxel$X_TIME_TO_EVENT, os_paclitaxel$X_EVENT) ~ 1)
summary(os_paclitaxel_survival)
plot(os_paclitaxel_survival, main="OS on Paclitaxel", xlab="Months", ylab="Overall Survival")


os_cisplatin1 <- merge(cisplatin, os, by = "sample")
os_cisplatin <- na.omit(os_cisplatin1)

os_cisplatin_survival <- survfit(Surv(os_cisplatin$X_TIME_TO_EVENT, os_cisplatin$X_EVENT) ~ 1)
summary(os_cisplatin_survival)
plot(os_cisplatin_survival, main="OS on Cisplatin", xlab="Months", ylab="Overall Survival")


os_carboplatin1 <- merge(carboplatin, os, by = "sample")
os_carboplatin <- na.omit(os_carboplatin1)

os_carboplatin_survival <- survfit(Surv(os_carboplatin$X_TIME_TO_EVENT, os_carboplatin$X_EVENT) ~ 1)
summary(os_carboplatin_survival)
plot(os_carboplatin_survival, main="OS on Carboplatin", xlab="Months", ylab="Overall Survival")


os_erlotinib1 <- merge(erlotinib, os, by = "sample")
os_erlotinib <- na.omit(os_erlotinib1)

os_erlotinib_survival <- survfit(Surv(os_erlotinib$X_TIME_TO_EVENT, os_erlotinib$X_EVENT) ~ 1)
summary(os_erlotinib_survival)
plot(os_erlotinib_survival, main="OS on Erlotinib", xlab="Months", ylab="Overall Survival")

####Kaplan-Meiers
stagei <- os_stage[c(os_stage$tumor_stage.diagnoses == "stage i"), ]
stageii <- os_stage[c(os_stage$tumor_stage.diagnoses == "stage ii"), ]
stageiii <- os_stage[c(os_stage$tumor_stage.diagnoses == "stage iii"), ]
stageiv <- os_stage[c(os_stage$tumor_stage.diagnoses == "stage iv"), ]

stagei_survival <- survfit(Surv(stagei$X_TIME_TO_EVENT, stagei$X_EVENT) ~ 1)
summary(stagei_survival)                    
plot(stagei_survival, main="OS for Stage i ", xlab="Months", ylab="Overall Survival")

stageii_survival <- survfit(Surv(stageii$X_TIME_TO_EVENT, stageii$X_EVENT) ~ 1)
summary(stageii_survival)                        
plot(stageii_survival, main="OS for Stage ii ", xlab="Months", ylab="Overall Survival")

stageiii_survival <- survfit(Surv(stageiii$X_TIME_TO_EVENT, stageiii$X_EVENT) ~ 1)
summary(stageiii_survival)                        
plot(stageiii_survival, main="OS for Stage iii ", xlab="Months", ylab="Overall Survival")

stageiv_survival <- survfit(Surv(stageiv$X_TIME_TO_EVENT, stageiv$X_EVENT) ~ 1)
summary(stageiv_survival)                        
plot(stageiv_survival, main="OS for Stage iv ", xlab="Months", ylab="Overall Survival")