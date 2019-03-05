# Uploading Data
## Phenotype Data
setwd("~/Desktop/nsclc/")
read.csv(file = "xena_gdctcga_pancancer_phenotype.csv", header = TRUE)
panphen.dat <- read.csv(file = "xena_gdctcga_pancancer_phenotype.csv", header = TRUE)

panphen.dat[panphen.dat==""] <- NA

## Survival Data
read.csv(file = "xena_gdctcga_pancancer_survival.csv", header = TRUE)
panos1 <- read.csv(file = "xena_gdctcga_adenocarcinoma_survival.csv", header = TRUE)
panos <- panos1[ , c(1, 3, 4)]

## Merge Phenotype and Survivial Data
names(panphen.dat)[1]<-"sample"
panall.dat <- merge(panphen.dat, panos,
                    by = "sample")
## Stages
panstage1 <- panall.dat[ , c(1,785)]
panstage1$tumor_stage.diagnoses <- gsub("not reported", "", panstage1$tumor_stage.diagnoses)
panstage1[panstage1==""] <- NA
panstage <- na.omit(panstage1)

### Compiled Stages (stage i = stage i, stage ia, stage ib, etc)

panstage$tumor_stage.diagnoses <- gsub("stage ia", "stage i", panstage$tumor_stage.diagnoses)
panstage$tumor_stage.diagnoses <- gsub("stage ib", "stage i", panstage$tumor_stage.diagnoses)

panstage$tumor_stage.diagnoses <- gsub("stage iia", "stage ii", panstage$tumor_stage.diagnoses)
panstage$tumor_stage.diagnoses <- gsub("stage iib", "stage ii", panstage$tumor_stage.diagnoses)

panstage$tumor_stage.diagnoses <- gsub("stage iiia", "stage iii", panstage$tumor_stage.diagnoses)
panstage$tumor_stage.diagnoses <- gsub("stage iiib", "stage iii", panstage$tumor_stage.diagnoses)

counts <- table(panstage)

barplot(counts, 
        main="Staging Data",
        xlab="Tumor Stage", 
        ylab="Frequency")
##Drug Therapy Data
pandrug1 <- panall.dat[ , c(1, 211)]
pandrug <- na.omit(pandrug1)
pandrugfreq <- table(pandrug$drug_name)
View(pandrugfreq)

drugs <- na.omit(panall.dat$drug_name)
drugs1 <-table(drugs)
View(drugs1)