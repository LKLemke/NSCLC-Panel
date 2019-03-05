## Uploading Data
# Clinical Data
setwd("~/Desktop/nsclc/")
read.csv(file = "xena_gdctcga_adenocarcinoma_phenotype.csv", header = TRUE)
phen.dat <- read.csv(file = "xena_gdctcga_adenocarcinoma_phenotype.csv", header = TRUE)

phen.dat[phen.dat==""] <- NA

# Survival Data
setwd("~/Desktop/nsclc/")
read.csv(file = "xena_gdctcga_adenocarcinoma_survival.csv", header = TRUE)
os1 <- read.csv(file = "xena_gdctcga_adenocarcinoma_survival.csv", header = TRUE)
os <- os1[ , c(1, 2, 4)]

# Merge Phenotype and Survivial Data
names(phen.dat)[1]<-"sample"
all.dat <- merge(phen.dat, os,
                 by = "sample")

# HTSeq-FPKM Data -------------------------to do
read.csv(file = "xena_gdctcga_adenocarcinoma_htseqfpkm.csv", header = TRUE)
exp.dat1 <- read.csv(file = "xena_gdctcga_adenocarcinoma_htseqfpkm.csv", header = TRUE)
exp.dat <- t(exp.dat1)
maybe <- as.matrix(exp.dat)
maybe5 <- maybe[2:nrow(maybe), 2:ncol(maybe)]
rownames(maybe5) <- maybe5[2:nrow(maybe5), 1]

maybe1 <- maybe5[1,1]
View(maybe1)
#--------------------------------------------

# Stages
stage1 <- all.dat[ , c(1,123)]
stage1$tumor_stage.diagnoses <- gsub("not reported", "", stage1$tumor_stage.diagnoses)
stage1[stage1==""] <- NA
stage <- na.omit(stage1)

#Compiled Stages (stage i = stage i, stage ia, stage ib, etc)

stage$tumor_stage.diagnoses <- gsub("stage ia", "stage i", stage$tumor_stage.diagnoses)
stage$tumor_stage.diagnoses <- gsub("stage ib", "stage i", stage$tumor_stage.diagnoses)

stage$tumor_stage.diagnoses <- gsub("stage iia", "stage ii", stage$tumor_stage.diagnoses)
stage$tumor_stage.diagnoses <- gsub("stage iib", "stage ii", stage$tumor_stage.diagnoses)

stage$tumor_stage.diagnoses <- gsub("stage iiia", "stage iii", stage$tumor_stage.diagnoses)
stage$tumor_stage.diagnoses <- gsub("stage iiib", "stage iii", stage$tumor_stage.diagnoses)

counts <- table(stage)

barplot(counts, 
        main="Staging Data",
        xlab="Tumor Stage", 
        ylab="Frequency")

# Drug Therapy Frequency
drugs <- na.omit(all.dat$drug_name)
drugs1 <-table(drugs)

# Drug Therapy Data
drug2 <- all.dat[ , c(1, 33)]
drug1 <- na.omit(drug2)

pemetrexed <- drug1[c(drug1$drug_name == "Pemetrexed"), ]
paclitaxel <- drug1[c(drug1$drug_name == "Paclitaxel"), ]
cisplatin <- drug1[c(drug1$drug_name == "Cisplatin"), ]
carboplatin <- drug1[c(drug1$drug_name == "Carboplatin"), ]
erlotinib  <- drug1[c(drug1$drug_name == "Erlotinib"), ]

drug <- rbind(pemetrexed, paclitaxel, cisplatin, carboplatin, erlotinib)

# Race
race1 <- all.dat[, c(1,107)]
race1$race.demographic <- gsub("not reported", "", race1$race.demographic)
race <- na.omit(race1)

# Gender
gender1 <- all.dat[ , c(1,106)]
gender1$gender.demographic <- gsub("not reported", "", gender1$gender.demographic)
gender <- na.omit(gender1)

# Days to Start Drug Therapy
ddt1 <- all.dat[ , c(1,25)]
ddt <- na.omit(ddt1)

range(ddt$days_to_drug_therapy_start)

# Cigarettes/Day
cpd1 <- all.dat[ , c(1,126)]
cpd <- na.omit(cpd1)

range(cpd$cigarettes_per_day.exposures)

# Years Smoked
ys1 <- all.dat[ , c(1,127)]
ys <- na.omit(ys1)

range(ys$years_smoked.exposures)

#SURVIVAL ANALYSIS
library(survival)

##OVERALL SURVIVAL

###Drug/OS
os_drug1 <- merge(drug, os,
                  by = "sample")
os_drug <- na.omit(os_drug1)
dim(os_drug)

os_drug_survival <- survdiff(Surv(os_drug$X_TIME_TO_EVENT, os_drug$X_EVENT) ~ os_drug$drug_name)
os_drug_survival

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

###Stage/OS
os_stage1 <- merge(stage, os,
                   by = "sample")
os_stage <- na.omit(os_stage1)

os_stage_survival <- survdiff(Surv(os_stage$X_TIME_TO_EVENT, os_stage$X_EVENT)
                              ~ os_stage$tumor_stage.diagnoses)
os_stage_survival

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

###Stage+Drug/OS
os_stage_drug1 <- merge(os_stage, drug,
                        by = "sample")
os_stage_drug <- na.omit(os_stage_drug1)


os_stage_drug_survival <- survdiff(Surv(os_stage_drug$X_TIME_TO_EVENT, os_stage_drug$X_EVENT)
                                   ~ os_stage_drug$drug_name + os_stage_drug$tumor_stage.diagnoses)
os_stage_drug_survival

###Race/OS
os_race1 <- merge(race, os,
                  by = "sample")
os_race <- na.omit(os_race1)
dim(os_race)

os_race_survival <- survdiff(Surv(os_race$X_TIME_TO_EVENT, os_race$X_EVENT)
                             ~ os_race$race.demographic)
os_race_survival

###Gender/OS
os_gender1 <- merge(gender, os,
                    by = "sample")
os_gender <- na.omit(os_gender1)
dim(os_gender)

os_gender_survival <- survdiff(Surv(os_gender$X_TIME_TO_EVENT, os_gender$X_EVENT)
                               ~ os_gender$gender.demographic)
os_gender_survival

###Cigarettes per Day/OS
os_cpd1 <- merge(cpd, os,
                 by = "sample")
os_cpd <- na.omit(os_cpd1)
dim(os_cpd)

os_cpd_survival <- survdiff(Surv(os_cpd$X_TIME_TO_EVENT, os_cpd$X_EVENT)
                            ~ os_cpd$cigarettes_per_day.exposures)
os_cpd_survival

###Days to Drug Therapy/OS
os_ddt1 <- merge(ddt, os,
                 by = "sample")
os_ddt <- na.omit(os_ddt1)
dim(os_ddt)

os_ddt_survival <- survdiff(Surv(os_ddt$X_TIME_TO_EVENT, os_ddt$X_EVENT)
                            ~ os_ddt$days_to_drug_therapy_start)
os_ddt_survival

###Years Smoked/OS
os_ys1 <- merge(ys, os,
                by = "sample")
os_ys <- na.omit(os_ys1)
dim(os_ys)

os_ys_survival <- survdiff(Surv(os_ys$X_TIME_TO_EVENT, os_ys$X_EVENT)
                           ~ os_ys$years_smoked.exposures)
os_ys_survival

###Race + Stage/OS
os_race_stage1 <- merge(os_stage, race, by = "sample")
os_race_stage <- na.omit(os_race_stage1)
dim(os_race_stage)

os_race_stage_survival <- survdiff(Surv(os_race_stage$X_TIME_TO_EVENT, os_race_stage$X_EVENT)
                                   ~ os_race_stage$tumor_stage.diagnoses + os_race_stage$race.demographic)
os_race_stage_survival

###Race + Drug/OS
os_race_drug1 <- merge(os_race, drug, by = "sample")
os_race_drug <- na.omit(os_race_drug1)
dim(os_race_drug)

os_race_drug_survival <- survdiff(Surv(os_race_drug$X_TIME_TO_EVENT, os_race_drug$X_EVENT)
                                  ~ os_race_drug$drug_name + os_race_drug$race.demographic)
os_race_drug_survival

###CPD + YS/OS
os_cpd_ys1 <- merge(os_cpd, ys, by="sample")
os_cpd_ys <- na.omit(os_cpd_ys1)
dim(os_cpd_ys)

os_cpd_ys_survival <- survdiff(Surv(os_cpd_ys$X_TIME_TO_EVENT, os_cpd_ys$X_EVENT)
                               ~ os_cpd_ys$cigarettes_per_day.exposures + os_cpd_ys$years_smoked.exposures)
os_cpd_ys_survival

###Race + DDT/OS
os_race_ddt1 <- merge(os_race, ddt, by = "sample")
os_race_ddt <- na.omit(os_race_ddt1)
dim(os_race_ddt)

os_race_ddt_survival <- survdiff(Surv(os_race_ddt$X_TIME_TO_EVENT, os_race_ddt$X_EVENT)
                                 ~ os_race_ddt$days_to_drug_therapy_start + os_race_ddt$race.demographic)
os_race_ddt_survival

###Race + CPD/OS
os_race_cpd1 <- merge(os_race, cpd, by = "sample")
os_race_cpd <- na.omit(os_race_cpd1)
dim(os_race_cpd)

os_race_cpd_survival <- survdiff(Surv(os_race_cpd$X_TIME_TO_EVENT, os_race_cpd$X_EVENT)
                                 ~ os_race_cpd$cigarettes_per_day.exposures + os_race_cpd$race.demographic)
os_race_cpd_survival

###Race + YS/OS
os_race_ys1 <- merge(os_race, ys, by = "sample")
os_race_ys <- na.omit(os_race_ys1)
dim(os_race_ys)

os_race_ys_survival <- survdiff(Surv(os_race_ys$X_TIME_TO_EVENT, os_race_ys$X_EVENT)
                                ~ os_race_ys$years_smoked.exposures + os_race_ys$race.demographic)
os_race_ys_survival

##Hazard Ratios for Race (OS and NT)

## New Tumor Event Survival Analysis

#New Tumor Survival Data
nt1 <- all.dat[ , c(1, 28, 55)]
nt <-na.omit(nt1)

nt$new_tumor_event_after_initial_treatment <- gsub("NO", "0", nt$new_tumor_event_after_initial_treatment)
nt$new_tumor_event_after_initial_treatment <- gsub("YES", "1", nt$new_tumor_event_after_initial_treatment)
nt$new_tumor_event_after_initial_treatment <- as.numeric(as.character(nt$new_tumor_event_after_initial_treatment))

###Drug/nt
nt_drug1 <- merge(drug, nt,
                  by = "sample")
nt_drug <- na.omit(nt_drug1)
dim(nt_drug)

nt_drug_survival <- survdiff(Surv(nt_drug$days_to_new_tumor_event_after_initial_treatment, nt_drug$new_tumor_event_after_initial_treatment)
                             ~ nt_drug$drug_name)
nt_drug_survival

###Stage/nt
nt_stage1 <- merge(stage, nt,
                   by = "sample")
nt_stage <- na.omit(nt_stage1)
dim(nt_stage)

nt_stage_survival <- survdiff(Surv(nt_stage$days_to_new_tumor_event_after_initial_treatment, nt_stage$new_tumor_event_after_initial_treatment)
                              ~ nt_stage$tumor_stage.diagnoses)
nt_stage_survival

###Stage+Drug/nt
nt_stage_drug1 <- merge(nt_stage, drug,
                        by = "sample")
nt_stage_drug <- na.omit(nt_stage_drug1)
dim(nt_stage_drug)

nt_stage_drug_survival <- survdiff(Surv(nt_stage_drug$days_to_new_tumor_event_after_initial_treatment, nt_stage_drug$new_tumor_event_after_initial_treatment)
                                   ~ nt_stage_drug$drug_name + nt_stage_drug$tumor_stage.diagnoses)
nt_stage_drug_survival

###Race/nt
nt_race1 <- merge(race, nt,
                  by = "sample")
nt_race <- na.omit(nt_race1)
dim(nt_race)

nt_race_survival <- survdiff(Surv(nt_race$days_to_new_tumor_event_after_initial_treatment, nt_race$new_tumor_event_after_initial_treatment)
                             ~ nt_race$race.demographic)
nt_race_survival

###Gender/nt
nt_gender1 <- merge(gender, nt,
                    by = "sample")
nt_gender <- na.omit(nt_gender1)
dim(nt_gender)

nt_gender_survival <- survdiff(Surv(nt_gender$days_to_new_tumor_event_after_initial_treatment, nt_gender$new_tumor_event_after_initial_treatment)
                               ~ nt_gender$gender.demographic)
nt_gender_survival

###Days to Drug Therapy/nt
nt_ddt1 <- merge(ddt, nt,
                 by = "sample")
nt_ddt <- na.omit(nt_ddt1)
dim(nt_ddt)

nt_ddt_survival <- survdiff(Surv(nt_ddt$days_to_new_tumor_event_after_initial_treatment, nt_ddt$new_tumor_event_after_initial_treatment)
                            ~ nt_ddt$days_to_drug_therapy_start)
nt_ddt_survival

###Years Smoked/nt
nt_ys1 <- merge(ys, nt, by = "sample")
nt_ys <- na.omit(nt_ys1)
dim(nt_ys)

nt_ys_survival <- survdiff(Surv(nt_ys$days_to_new_tumor_event_after_initial_treatment, nt_ys$new_tumor_event_after_initial_treatment)
                           ~ nt_ys$years_smoked.exposures)
nt_ys_survival

###Cigarettes per Day/nt
nt_cpd1 <- merge(cpd, nt, by = "sample")
nt_cpd <- na.omit(nt_cpd1)
dim(nt_cpd)

nt_cpd_survival <- survdiff(Surv(nt_cpd$days_to_new_tumor_event_after_initial_treatment, nt_cpd$new_tumor_event_after_initial_treatment)
                            ~ nt_cpd$cigarettes_per_day.exposures)
nt_cpd_survival

###Race + Stage/nt
nt_race_stage1 <- merge (nt_race, stage, by="sample")
nt_race_stage <- na.omit(nt_race_stage1)
dim(nt_race_stage)

nt_race_stage_survival <- survdiff(Surv(nt_race_stage$days_to_new_tumor_event_after_initial_treatment, nt_race_stage$new_tumor_event_after_initial_treatment)
                                   ~ nt_race_stage$race.demographic + nt_race_stage$tumor_stage.diagnoses)
nt_race_stage_survival

###Race + Drug/nt
nt_race_drug1 <- merge (nt_race, drug, by="sample")
nt_race_drug <- na.omit(nt_race_drug1)
dim(nt_race_drug)

nt_race_drug_survival <- survdiff(Surv(nt_race_drug$days_to_new_tumor_event_after_initial_treatment, nt_race_drug$new_tumor_event_after_initial_treatment)
                                  ~ nt_race_drug$race.demographic + nt_race_drug$drug_name)
nt_race_drug_survival

###Race + DDT/nt
nt_race_ddt1 <- merge (nt_race, ddt, by="sample")
nt_race_ddt <- na.omit(nt_race_ddt1)
dim(nt_race_ddt)

nt_race_ddt_survival <- survdiff(Surv(nt_race_ddt$days_to_new_tumor_event_after_initial_treatment, nt_race_ddt$new_tumor_event_after_initial_treatment)
                                 ~ nt_race_ddt$race.demographic + nt_race_ddt$days_to_drug_therapy_start)
nt_race_ddt_survival

###Race + CPD/nt
nt_race_cpd1 <- merge (nt_race, cpd, by="sample")
nt_race_cpd <- na.omit(nt_race_cpd1)
dim(nt_race_cpd)

nt_race_cpd_survival <- survdiff(Surv(nt_race_cpd$days_to_new_tumor_event_after_initial_treatment, nt_race_cpd$new_tumor_event_after_initial_treatment)
                                 ~ nt_race_cpd$race.demographic + nt_race_cpd$cigarettes_per_day.exposures)
nt_race_cpd_survival

###Race + YS/nt
nt_race_ys1 <- merge (nt_race, ys, by="sample")
nt_race_ys <- na.omit(nt_race_ys1)
dim(nt_race_ys)

nt_race_ys_survival <- survdiff(Surv(nt_race_ys$days_to_new_tumor_event_after_initial_treatment, nt_race_ys$new_tumor_event_after_initial_treatment)
                                ~ nt_race_ys$race.demographic + nt_race_ys$years_smoked.exposures)
nt_race_ys_survival

###CPD + YS/nt
nt_cpd_ys1 <- merge (nt_cpd, ys, by="sample")
nt_cpd_ys <- na.omit(nt_cpd_ys1)
dim(nt_cpd_ys)

nt_cpd_ys_survival <- survdiff(Surv(nt_cpd_ys$days_to_new_tumor_event_after_initial_treatment, nt_cpd_ys$new_tumor_event_after_initial_treatment)
                               ~ nt_cpd_ys$cigarettes_per_day.exposures + nt_cpd_ys$years_smoked.exposures)
nt_cpd_ys_survival

#######MORE THIGNS?!?!?