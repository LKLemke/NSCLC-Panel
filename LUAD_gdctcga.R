## Uploading Data
# Clinical Data
read.csv(file = "xena_gdctcga_adenocarcinoma_phenotype.csv", header = TRUE)
phen.dat <- read.csv(file = "xena_gdctcga_adenocarcinoma_phenotype.csv", header = TRUE)

phen.dat[phen.dat==""] <- NA

# Survival Data
read.csv(file = "xena_gdctcga_adenocarcinoma_survival.csv", header = TRUE)
os1 <- read.csv(file = "xena_gdctcga_adenocarcinoma_survival.csv", header = TRUE)
os <- os1[ , c(1, 2, 4)]
colnames(os)[which(colnames(os) == 'ï..sample')] <- 'sample'


# Merge Phenotype and Survivial Data
names(phen.dat)[1]<-"sample"
all.dat <- merge(phen.dat, os,
                 by = "sample")

# HTSeq-FPKM Data
##exp.dat <- read.delim("tcga-LUAD.htseq_fpkm.tsv", header = TRUE)
##exp1 <- t(exp.dat)
##colnames(exp1)[which(colnames(exp1) == 'Ensembl_ID')] <- 'sample'
##exp1[exp1==""] <- NA
##exp <- na.omit(exp1)


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
race1[race1==""] <- NA
race <- na.omit(race1)

counts_race <- table(race)
barplot(counts_race,
        main="Race", ylab="Frequency")

# Gender
gender1 <- all.dat[ , c(1,106)]
gender1$gender.demographic <- gsub("not reported", "", gender1$gender.demographic)
gender1[gender1==""] <- NA
gender <- na.omit(gender1)

counts_gender <- table(gender)
barplot(counts_gender,
        main="Gender", ylab="Frequency")

# Age at Pathologic Diagnosis
age1 <- all.dat[ , c(1,6)]
age1[age1==""] <- NA
age <- na.omit(age1)

range(age$age_at_initial_pathologic_diagnosis)
mean(age$age_at_initial_pathologic_diagnosis)
median(age$age_at_initial_pathologic_diagnosis)

# Days to Start Drug Therapy
ddt1 <- all.dat[ , c(1,25)]
ddt1[ddt1==""] <- NA
ddt <- na.omit(ddt1)

range(ddt$days_to_drug_therapy_start)

# Cigarettes/Day
cpd1 <- all.dat[ , c(1,126)]
cpd1[cpd1==""] <- NA
cpd <- na.omit(cpd1)

range(cpd$cigarettes_per_day.exposures)
mean(cpd$cigarettes_per_day.exposures)
median(cpd$cigarettes_per_day.exposures)
boxplot(cpd$cigarettes_per_day.exposures, main="Cigarettes Per Day")

# Years Smoked
ys1 <- all.dat[ , c(1,127)]
ys1[ys1==""] <- NA
ys <- na.omit(ys1)

range(ys$years_smoked.exposures)
mean(ys$years_smoked.exposures)
median(ys$years_smoked.exposure)
boxplot(ys$years_smoked.exposure, main="Years Smoked")

# Pack Years Smoked
pys1 <- all.dat[ , c(1,57)]
pys1[pys1==""] <- NA
pys <- na.omit(pys1)

range(pys$number_pack_years_smoked)
mean(pys$number_pack_years_smoked)
median(pys$number_pack_years_smoked)
boxplot(pys$number_pack_years_smoked, main="Pack Years Smoked")

#Prognostic Biomarkers
##EGFR
##ERBB2
##RET
##MET (amplification)
##FGFR1 (amplification)
##KRAS
##BRAF
##PIK3CA
##PTEN
##H3F3A

# Mutation Data
mutdat1 <- read.delim("tcga-luad.mutations.tsv.txt", header = TRUE, stringsAsFactors = FALSE)
colnames(mutdat1)[which(colnames(mutdat1) == 'SAMPLE_ID')] <- 'sample'
mutdat <- mutdat1[ ,-1]
length(unique(mutdat$EGFR))

mut_columns <- as.matrix(mutdat[,-1])
mut_columns[! is.na(mut_columns)] <- 1
mut_columns[is.na(mut_columns)] <- 0
colnames(mut_columns) <- paste0(colnames(mut_columns), "_mutated")
mutdat <- cbind(mutdat, mut_columns, stringsAsFactors = F)

os_sample <- os
os_sample$sample <- substr(os_sample$sample, start = 1, stop = 15)
table(duplicated(os_sample$sample))
duplicates <- unique(os_sample$sample[duplicated(os_sample$sample)])
dup.df <- os_sample[os_sample$sample %in% duplicates,]
os_sample <- os_sample[! duplicated(os_sample$sample),]
table(duplicated(os_sample$sample))

table(mutdat$sample %in% os_sample$sample)
os_mut <- merge(mutdat, os_sample, by = "sample", all.x = F, all.y = F, stringsAsFactors = F)


# CNV Data (MET, FGFR1)
cnvdat1 <- read.delim("tcga-luad.cnv.tsv.txt", header = TRUE)
colnames(cnvdat1)[which(colnames(cnvdat1) == 'SAMPLE_ID')] <- 'sample'
cnvdat <- cnvdat1[ ,-1]

os_cnv <- merge(cnvdat, os_sample, by = "sample")

#SURVIVAL ANALYSIS
library(survival)

Surv.luad <- Surv(event=os_mut$X_EVENT, time=os_mut$X_TIME_TO_EVENT)

library(glmnet)
luadCOX.cv <- cv.glmnet(t(luad.clara$medoids), Surv.luad, family="cox", alpha=1, nfolds=10)
luadCOX <- glmnet(t(luad.clara$medoids), Surv.luad, family="cox", alpha=1, lambda=luadCOX.cv$lambda.min)

library(randomForestSRC)
temp <- cbind.data.frame("event"=luad.surv.2$X_EVENT, "time"=luad.surv.2$X_TIME_TO_EVENT, t(luad.clara$medoids))
luadRF <- rfsrc(Surv(time, event) ~ ., data=temp, ntree = 500, block.size = 10, importance=TRUE)

par(cex=0.25, mar=c(1,14,2,2))
barplot(sort(luadRF$importance), horiz=TRUE, las=1)


##OVERALL SURVIVAL

###EGFR/os
os_mut_survival_EGFR <- survdiff(Surv(os_mut$X_TIME_TO_EVENT, os_mut$X_EVENT) ~ os_mut$EGFR_mutated)
os_mut_survival_EGFR

###ERBB2/os
os_mut_survival_ERBB2 <- survdiff(Surv(os_mut$X_TIME_TO_EVENT, os_mut$X_EVENT) ~ os_mut$ERBB2_mutated)
os_mut_survival_ERBB2

###RET/os
os_mut_survival_RET <- survdiff(Surv(os_mut$X_TIME_TO_EVENT, os_mut$X_EVENT) ~ os_mut$RET_mutated)
os_mut_survival_RET

###MET/os
os_cnv_survival_MET <- survdiff(Surv(os_cnv$X_TIME_TO_EVENT, os_cnv$X_EVENT) ~ os_cnv$MET)
os_cnv_survival_MET

###FGFR1/os
os_cnv_survival_FGFR1 <- survdiff(Surv(os_cnv$X_TIME_TO_EVENT, os_cnv$X_EVENT) ~ os_cnv$FGFR1)
os_cnv_survival_FGFR1

###KRAS/os
os_mut_survival_KRAS <- survdiff(Surv(os_mut$X_TIME_TO_EVENT, os_mut$X_EVENT) ~ os_mut$KRAS_mutated)
os_mut_survival_KRAS

###BRAF/os
os_mut_survival_BRAF <- survdiff(Surv(os_mut$X_TIME_TO_EVENT, os_mut$X_EVENT) ~ os_mut$BRAF_mutated)
os_mut_survival_BRAF

###PIK3CA/os
os_mut_survival_PIK3CA <- survdiff(Surv(os_mut$X_TIME_TO_EVENT, os_mut$X_EVENT) ~ os_mut$PIK3CA_mutated)
os_mut_survival_PIK3CA

###PTEN/os
os_mut_survival_PTEN <- survdiff(Surv(os_mut$X_TIME_TO_EVENT, os_mut$X_EVENT) ~ os_mut$PTEN_mutated)
os_mut_survival_PTEN

###H3F3A/os
      #no variants :(


###Drug/OS
os_drug1 <- merge(drug, os,
                  by = "sample")
os_drug <- na.omit(os_drug1)
dim(os_drug)

os_drug_survival <- survdiff(Surv(os_drug$X_TIME_TO_EVENT, os_drug$X_EVENT) ~ os_drug$drug_name)
os_drug_survival


###Stage/OS
os_stage1 <- merge(stage, os,
                   by = "sample")
os_stage <- na.omit(os_stage1)

os_stage_survival <- survdiff(Surv(os_stage$X_TIME_TO_EVENT, os_stage$X_EVENT)
                              ~ os_stage$tumor_stage.diagnoses)
os_stage_survival


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

###Pack Years Smoked/OS
os_pys1 <- merge(pys, os,
                by = "sample")
os_pys <- na.omit(os_pys1)
dim(os_pys)

os_pys_survival <- survdiff(Surv(os_pys$X_TIME_TO_EVENT, os_pys$X_EVENT)
                           ~ os_pys$number_pack_years_smoked)
os_pys_survival

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

###Biomarker/OS
os_mutdat1 <- merge(os, mutdat, by = "sample")


##Hazard Ratios for Race (OS and NT)????


## New Tumor Event Survival Analysis

#New Tumor Survival Data
nt1 <- all.dat[ , c(1, 28, 55)]
nt <-na.omit(nt1)

nt$new_tumor_event_after_initial_treatment <- gsub("NO", "0", nt$new_tumor_event_after_initial_treatment)
nt$new_tumor_event_after_initial_treatment <- gsub("YES", "1", nt$new_tumor_event_after_initial_treatment)
nt$new_tumor_event_after_initial_treatment <- as.numeric(as.character(nt$new_tumor_event_after_initial_treatment))

nt_sample <- nt
nt_sample$sample <- substr(nt_sample$sample, start = 1, stop = 15)
table(duplicated(nt_sample$sample))
duplicates_nt <- unique(nt_sample$sample[duplicated(nt_sample$sample)])
dup.df_nt <- nt_sample[nt_sample$sample %in% duplicates_nt,]
nt_sample <- nt_sample[! duplicated(nt_sample$sample),]
table(duplicated(nt_sample$sample))

nt_mut <- merge(mutdat, nt_sample, by = 'sample')
nt_cnv <- merge(cnvdat, nt_sample, by = "sample")


###EGFR/nt
nt_mut_survival_EGFR <- survdiff(Surv(nt_mut$days_to_new_tumor_event_after_initial_treatment,
                                      nt_mut$new_tumor_event_after_initial_treatment)~ nt_mut$EGFR_mutated)
nt_mut_survival_EGFR

###ERBB2/nt
nt_mut_survival_ERBB2 <- survdiff(Surv(nt_mut$days_to_new_tumor_event_after_initial_treatment,
                                       nt_mut$new_tumor_event_after_initial_treatment) ~ nt_mut$ERBB2_mutated)
nt_mut_survival_ERBB2

###RET/nt
nt_mut_survival_RET <- survdiff(Surv(nt_mut$days_to_new_tumor_event_after_initial_treatment,
                                     nt_mut$new_tumor_event_after_initial_treatment) ~ nt_mut$RET_mutated)
nt_mut_survival_RET

###MET/nt
nt_cnv_survival_MET <- survdiff(Surv(nt_cnv$days_to_new_tumor_event_after_initial_treatment,
                                     nt_cnv$new_tumor_event_after_initial_treatment) ~ nt_cnv$MET)
nt_cnv_survival_MET

###FGFR1/nt
nt_cnv_survival_FGFR1 <- survdiff(Surv(nt_cnv$days_to_new_tumor_event_after_initial_treatment, 
                                       nt_cnv$new_tumor_event_after_initial_treatment) ~ nt_cnv$FGFR1)
nt_cnv_survival_FGFR1

###KRAS/nt
nt_mut_survival_KRAS <- survdiff(Surv(nt_mut$days_to_new_tumor_event_after_initial_treatment, 
                                      nt_mut$new_tumor_event_after_initial_treatment) ~ nt_mut$KRAS_mutated)
nt_mut_survival_KRAS

###BRAF/nt
nt_mut_survival_BRAF <- survdiff(Surv(nt_mut$days_to_new_tumor_event_after_initial_treatment, 
                                      nt_mut$new_tumor_event_after_initial_treatment) ~ nt_mut$BRAF_mutated)
nt_mut_survival_BRAF

###PIK3CA/os
nt_mut_survival_PIK3CA <- survdiff(Surv(nt_mut$days_to_new_tumor_event_after_initial_treatment, 
                                        nt_mut$new_tumor_event_after_initial_treatment) ~ nt_mut$PIK3CA_mutated)
nt_mut_survival_PIK3CA

###PTEN/os
nt_mut_survival_PTEN <- survdiff(Surv(nt_mut$days_to_new_tumor_event_after_initial_treatment, 
                                      nt_mut$new_tumor_event_after_initial_treatment) ~ nt_mut$PTEN_mutated)
nt_mut_survival_PTEN

###H3F3A/nt
#no variants :(

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

