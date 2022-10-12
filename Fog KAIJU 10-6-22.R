#Impact of HIV infection or exposure on the nasopharyngeal microbiome of children in Botswana
#Redoing with KAIJU output
#Compiled February 17, 2022 by Sweta Patel
#Redone Feb 22, 2022 by Sweta Patel (using MK Kaiju pipeline ID'ing with progenome database)
#Redone April 19, 2022 by Sweta Patel (using decontaminated kaiju data)
#Redone May 6, 2022 by Sweta Patel incorporating 5/5 lab meeting feedback
#Updated May 26, 2022 by Sweta Patel to clean up code 
#Updated June 9, 2022 by Sweta Patel to incorporate age into all compositional analyses
#Updated Aug 3, 2022 by Sweta Patel to update figures and analyses after MK manuscript revision #1
#Updated Aug 28, 2022 by Sweta Patel to update figures and analysis after MK manuscript revision #2
#Updated Oct 6, 2022 by Sweta Patel to incorporate co-author feedback

#***GENERAL ORDER OF STEPS***
#1. Clean up taxtable (+ ID all corynebacterium spp and other common sp) and agglomerate at species level
#2. Look at observed reads and see where they plateau +/-construct rarefaction curves using the raw data (i.e. before filtering or normalization)
#3. Use the plateau/rarefaction curves to inform your cutoff for sample pruning
#4. Prune
#5. Basic demographics
#6. Alpha diversity analyses
#7. Transform data (CLR) and do Beta-diversity analyses (PCoA plots, PERMANOVA)
  #UPDATE 5/5: the order covariables are listed in PERMANOVA affects the results
#8. Generate relative abundance plots using NON-transformed data
#9. UPDATE 5/21: Filter on mean relative abundance of 50 and remove NOS species
#10. Differential abundance analyses using Maaslin and FILTERED data
    #set Q value to 0.1-0.2
#11. Logistic regression analyses (Aim 2)
                                                                           
#NOTE: Kaiju results are in counts, not rel abundances like metaphlan data

#***LOADING PACKAGES***
library(plyr)
library(dplyr)
library(ggplot2)
library(phyloseq)
library(tidyr)
library(vegan)
library(metagenomeSeq)
library(httr)
library(gridExtra)
library(data.table)
library(RColorBrewer)
library(scales)
library(gsubfn)
library(readxl)
library(microbiome)
library(MASS)
library(Maaslin2)
library(cowplot)

set.seed(1234)

setwd("OneDrive - Duke University/Fogarty coding/Sequencing/")

###############################################################################################
#INITIAL STEPS: SKIP TO LINE 181 TO IMPORT PSEQ WITH NEW VARIABLES AND LINE 396 FOR PRUNED PSEQ
###############################################################################################
#***LOADING DATA***
#importing RDS into R
meta <- readRDS("phy.kaiju.decontam.04192022.rds")

#creating metadata from our phyloseq object
metadata <- data.frame(sample_data(meta))

#Number of reads per sample #NOTE: this is from 02-17-22
reads <- read.csv("HIV_NGS_kneaddata_counts.csv")
View(reads)
reads$Sample <- gsubfn("_R1_kneaddata", "", reads$Sample)
reads$Sample <- gsubfn("-", ".", reads$Sample)
reads <- data.frame(reads$Sample, reads$final.pair1)

#rename Sample variable to sample_id to facilitate merging with metadata file
reads <- rename(reads, sample_id = reads.Sample)
reads <- rename(reads, paired.reads = reads.final.pair1)

##############
#NEW VARIABLES
##############
#sample_id
#rownames as new variable: sample_id
names(metadata)
metadata$sample_id <- rownames(metadata)

#subject
#grepl to create mom, child, sib as subject
#source (grepl): https://community.rstudio.com/t/creating-a-string-variable-based-on-another-variable/27270/4
metadata$subject <- "99"
metadata <- mutate(metadata, subject = case_when(grepl("MOM", sample_id) ~ "mom", grepl("CHI", sample_id) ~ "child",
                                                 grepl("SIB", sample_id) ~ "sib"))
table(metadata$subject)

#child HIV status
metadata$hiv <- NA
metadata$hiv [metadata$child_pcr_1 == 1 | metadata$child_pcr_2 == 1] <- "Infected"
metadata$hiv [metadata$child_pcr_1 == 2 & metadata$mat_hiv == 1] <- "Exposed Uninfected"
metadata$hiv [metadata$child_pcr_1 == 4] <- "Unexposed"
table(metadata$hiv)

#child HIV status, collapsed into HIV/HEU vs HUU
metadata$hivcat <- NA
metadata$hivcat [metadata$child_pcr_1 == 1 | metadata$child_pcr_2 == 1] <- "EIEU"
metadata$hivcat [metadata$child_pcr_1 == 2 & metadata$mat_hiv == 1] <- "EIEU"
metadata$hivcat [metadata$child_pcr_1 == 4] <- "HUU"
table(metadata$hivcat)

#season
  #For consistency with MK: Rainy season = Nov to March. Dry season = April to Oct
metadata$season <- "99"
metadata$season [metadata$month == 1 | metadata$month == 2 | metadata$month == 3 | metadata$month == 11 | metadata$month == 12] <- "Rainy"
metadata$season [metadata$month == 4 | metadata$month == 5 | metadata$month == 6 | metadata$month == 7 | metadata$month == 8 | metadata$month == 9 | metadata$month == 10] <- "Dry"
table(metadata$season)

#people in household
metadata$hhsize <- metadata$num_adol + metadata$num_adults + metadata$num_child 

#maternal education
table(metadata$mat_educ) #None (1), Primary (14), Secondary (109), Tertiary (17), Unknown (1)
metadata$mat_educ2 <- NA
metadata$mat_educ2 [metadata$mat_educ == 1 | metadata$mat_educ == 2] <- "None or Primary"
metadata$mat_educ2 [metadata$mat_educ == 3 | metadata$mat_educ == 4] <- "Secondary"
metadata$mat_educ2 [metadata$mat_educ == 5] <- "Tertiary"
table(metadata$mat_educ2)

#Adding binary variable for recent/current URI
metadata$uri_cur2 <- NA
metadata$uri_cur2 [metadata$uri_current == 1] <- "Yes"
metadata$uri_cur2 [metadata$uri_current == 2] <- "No"
table(metadata$uri_cur2)
table(metadata$uri_current)

metadata$uri_rec2 <- NA
metadata$uri_rec2 [metadata$uri_recent == 1] <- "Yes"
metadata$uri_rec2 [metadata$uri_recent == 2] <- "No"
table(metadata$uri_rec2)

#Antibiotic exposure
table(metadata$meds_name)  #referring to "REDCap recoding notes - 9-8-21.docx" will create new abx variable
  #Only include systemically received abx for now (exclude otic drops and topical abx)
  #ATT, penicillin, amoxicillin, cloxacillin, metronidazole, TMP-SMX, nalidixic acid, erythromycin, ceftriaxone, cefotaxime , amox-clav
#d/w MK 12-9-21: remove all non-abx in excel and then re-import into REDCap --> Done
  #will merge the clean meds_name data with our current metadata dataset
abx <- read.csv("systemicabx_2021-12-09_0921.csv")
  #need to add NAs? 
#drop current med_names variable from metadata and then merge with abx df by study_id
  #standardize study_id
toString(abx$study_id)
abx$study_id <- gsub("-", ".", abx$study_id)
metadata$meds_name <- NULL
metadata <- as.data.frame(metadata)
metadata <- merge(metadata, abx, by ="study_id", sort = TRUE)
table(metadata$meds_name)
metadata$abx <- "Yes"
metadata$abx [metadata$meds_name == ""] <- "No"
table(metadata$abx) #visually compared to meds_name and confirmed it looks good!

#will drop arv variables from the meta dataset and then merge with arv dataset by studyID
arvs <- read.csv("TheImpactOfHIVInfect-ARVs_DATA_2021-10-07_1127.csv")
toString(arvs$study_id)
arvs$study_id <- gsub("-", ".", arvs$study_id)
metadata$mat_arv_meds <- NULL
metadata$child_arv_meds <-NULL
metadata <- merge(metadata, arvs, by = "study_id", sort = TRUE)
table(metadata$mat_arv_meds)
table(metadata$child_arv_meds)
#looks like it matches up! 

#number of reads
metadata <- merge(metadata, reads, by = "sample_id", sort = TRUE)

#import the new subject var back into pseq
#source (importing metadata into pseq): https://github.com/joey711/phyloseq/issues/1106
  #NOTE: new problems doing this, likely related to changes in row names when we merged abx and arv datasets
row.names(metadata) <- metadata$sample_id
metadata <- sample_data(metadata)
meta <- merge_phyloseq(meta, metadata)
#now see if it worked:
bloop <- data.frame(sample_data(meta))
#success!
remove(bloop)

#clean up other extra datasets that we don't need anymore: abx, arvs, reads
remove(abx, arvs, reads)

#save updated pseq object with the new variables so we don't have to re-do them every time
saveRDS(meta, file = "phy.kaiju041922.rds", ascii = FALSE, version = NULL,
        compress = FALSE, refhook = NULL)

meta <- readRDS("phy.kaiju041922.rds")
metadata <- data.frame(sample_data(meta))

###########################################
#NEW PHYLOSEQ OBJECT WITH DUPLICATE REMOVED
###########################################
otu <- as.data.frame(otu_table(meta, taxa_are_rows=TRUE))
#to ID any study_ids that are missing
csum <- colSums(otu)
# check if any are missing
any(is.na(csum))
#NONE missing in kaiju set
remove(csum)

#new pseq object excluding specific samples (source: https://github.com/joey711/phyloseq/issues/618)
#NOTE: will exclude B24 since that was the second enrollment for one of the pairs

meta_new = subset_samples(meta, sample_id != "B.24.CHI" & sample_id != "B.24.MOM")
remove(meta)
metanew <- data.frame(sample_data(meta_new))

#######################
#MOST ABUNDANT SPECIES?
#######################
#Source: https://github.com/joey711/phyloseq/issues/1487

# Merge species
gp = tax_glom(meta_new, taxrank = "Species")

# Calculate taxa sum of the selected samples
top100 = head(sort(rowSums(otu_table(gp)), decreasing = TRUE), 100)

# Combine count and taxonomyTable
top100 = cbind(as.data.frame(tax_table(gp)[names(top100),]), Count = top100)

table(top100$Genus)

#Top 10 species in our data (overall, NOT by subject):
  #1. D. pigrum
  #2. Corynebacterium propinquum
  #3. Coryne pseudodiphtheriticum
  #4. Moraxella nonliquefaciens
  #5. M. cattarhalis
  #6. Coryne accolens
  #7. Staph epi
  #8. Strep pneumo
  #9. Coryne 1859 (unidentified species)
  #10. M. lincolnii

###################
#CLEANING TAXTABLE
###################
taxtable <- as.data.frame(as(tax_table(meta_new),"matrix"),stringsAsFactors=FALSE)
#We know from MK that we need to edit our taxtable and import back into phyloseq obj
table(taxtable$Kingdom)  #All from kingdom Bacteria
table(taxtable$Phylum)    #Phyla look ok
table(taxtable$Class) 
table(taxtable$Order)   #Looks good, no changes
table(taxtable$Family) 
table(taxtable$Genus)

#10 species identified as "symbiont of" and 4 identified as "endosymbiont of"
  #Per 3/3/22 call with MK, will recode to NA
#1829 species identified with strain numbers or other identifier, not species name
  #ex "Haemophilus sp. C1" or "Corynebacterium sp. KPL1995"
#For all Coryne species and the NON-Coryne species in the top 100 most common, took the following steps:
 #1. Obtained GenBank accession number for each strain using https://www.ncbi.nlm.nih.gov/nucleotide/
      #(type in strain name ex KPL1995, and accession number will display on results page: Accession: KI515715.1)
 #2. Used TYGS website to identify strains using GenBank accession number: https://tygs.dsmz.de/user_requests/new
 #3. Saved all TYGS results as PDF files
#For the other NON-CORYNE species that are not identified by name and were not in the 100 most common species, will recode to NA

class(taxtable$Species) #character variable
taxtable$Species <- gsub("symbiont", NA, taxtable$Species)
#We have 16 Corynebacterium species to identify
  #No match in NCBI Assembly or TYGS found for: KPL1859, NML98-0116, HMSC29G08, HMSC074E01, HMSC055A01, 
    #HMSC078H07, HMSC04H06, HMSC034A01
taxtable$Species <- gsub("Corynebacterium sp. KPL1995", "Corynebacterium pseudodiphtheriticum", taxtable$Species)
taxtable$Species <- gsub("Corynebacterium sp. KPL1986", "Corynebacterium accolens", taxtable$Species)
taxtable$Species <- gsub("Corynebacterium sp. ATCC 6931", "Corynebacterium amycolatum", taxtable$Species)
taxtable$Species <- gsub("Corynebacterium sp. HMSC28B08", "Corynebacterium auriscanis", taxtable$Species)
taxtable$Species <- gsub("Corynebacterium sp. HMSC11E11", "Corynebacterium freneyi", taxtable$Species)

#Non-Corynebacterium species
taxtable$Species <- gsub("Haemophilus sp. C1", "Haemophilus haemolyticus", taxtable$Species)
taxtable$Species <- gsub("Streptomyces sp. e14", "Streptomyces carpinensis", taxtable$Species)

#NML98-0116 --> inconclusive by NCBI Assembly and TYGS
#HMSC28B08 --> Corynebacterium auriscanis by NCBI assembly; TYGS also ID'd auriscanis but with c/f potentially unreliable ID
#ATCC 6931 --> Corynebacterium amycolatum by NCBI assembly and TYGS
#HMSC29G08 --> inconclusive by NCBI Assembly and TYGS
#HMSC11E11 --> Corynebacterium freneyi by NCBI assembly; TYGS also ID'd freneyi but with c/f potentially unreliable ID
#HMSC074E01 --> inconclusive by NCBI Assembly and TYGS
#HMSC055A01 --> inconclusive by NCBI Assembly and TYGS
#HMSC078H07 --> inconclusive by NCBI Assembly and TYGS
#HMSC04H06 --> inconclusive by NCBI Assembly and TYGS
#HMSC034A01 --> inconclusive by NCBI Assembly and TYGS
#NML130628 --> inconclusive by NCBI Assembly and TYGS
#CNJ-954 --> inconclusive by NCBI Assembly and TYGS
#NML140438 --> inconclusive by NCBI Assembly and TYGS

#there are 8 non-Corynebacterium species that are not identified in our top 100
#Propionibacterium sp. KPL1844 --> anomalous assembly by NCBI; Cutibacterium granulosum by TYGS
#Haemophilus sp. C1 --> Haemophilus haemolyticus by assembly but inconclusive by TYGS
#Streptococcus sp. SK643 --> inconclusive by assembly and TYGS
#Streptococcus sp. oral taxon 431 --> inconclusive by assembly and TYGS
#Rhodococcus sp. RD6.2 --> inconclusive by assembly and TYGS
#Streptomyces sp. e14 --> Streptomyces carpinensis by NCBI but Streptomyces sennicomposti by TYGS
#Microbacterium sp. 70-16 --> inconclusive by assembly and TYGS
#Streptococcus sp. DD10 --> inconclusive by assembly and TYGS

# taxtable$Species <- gsub("sp.", NA, taxtable$Species)
#using "sp." removes all species names that contain "sp" like H. sputorum and Paracoccus sphaerophysae
taxtable$Species <- gsub("sp. ", NA, taxtable$Species)  #adding the space after fixes it!
taxtable$Species <- gsub("uncultured ", NA, taxtable$Species)
taxtable$Species <- gsub("Bradyrhizobium sp.", NA, taxtable$Species)

#FOR MK: what do i do with "candidatus" species? Per 4/10 email from MK, can leave or delete. Will leave for now.

#To get rid of all NAs in the taxtable (code from MK):
# Replace NA with data from higher taxonomic rank (identifies taxonomic rank with letter)
rplc <- which(taxtable$Phylum=="" | is.na(taxtable$Phylum))
taxtable$Phylum[rplc] <- paste(taxtable$Kingdom[rplc],";p",sep="")
rplc <- which(taxtable$Class=="" | is.na(taxtable$Class))
taxtable$Class[rplc] <- paste(taxtable$Phylum[rplc],";c",sep="")
rplc <- which(taxtable$Order=="" | is.na(taxtable$Order))
taxtable$Order[rplc] <- paste(taxtable$Class[rplc],";o",sep="")
rplc <- which(taxtable$Family=="" | is.na(taxtable$Family))
taxtable$Family[rplc] <- paste(taxtable$Order[rplc],";f",sep="")
rplc <- which(taxtable$Genus=="NA" | is.na(taxtable$Genus) | taxtable$Genus=="NA")
taxtable$Genus[rplc] <- paste(taxtable$Family[rplc],";g",sep="")
# rplc <- which(taxtable$Genus=="NA;c;o;f;g")
# taxtable$Genus[rplc]<-"Unassigned"
rplc <- which(taxtable$Species=="NA" | is.na(taxtable$Species) | taxtable$Species=="NA")
taxtable$Species[rplc] <- paste(taxtable$Genus[rplc],";s",sep="")
write.csv(taxtable, "taxtable_050922.csv", row.names=F)
remove(rplc)

#now import updated taxtable into pseq object (convert to matrix first) and take a look
TT <- as.matrix(taxtable)
tax_table(meta_new) <- TT
#can we test it to see if it worked? 
test <- as.data.frame(as(tax_table(meta_new),"matrix"),stringsAsFactors=FALSE) #it worked
remove(TT, test, taxtable)

#agglomerate at species level BEFORE looking at read cutoffs
#source: https://rdrr.io/bioc/phyloseq/man/tax_glom.html
#How many taxa before/after agglomeration?
ntaxa(meta_new) #7998
meta_new <- tax_glom(meta_new, taxrank="Species")
ntaxa(meta_new) #6709 taxa after agglom

########
#PRUNING
########
#Create dataframe with diversity indices for each specimen (so we can plot #observed spp by #reads)
diversity_mn <- estimate_richness(meta_new, measures = c("Shannon", "Simpson", "Chao1", "Observed"))

setDT(diversity_mn, keep.rownames = TRUE)[]
colnames(diversity_mn) <- c("sample_id", "Observed", "Chao1", "se.chao1", "Shannon", "Simpson")
rownames(diversity_mn) <- diversity_mn$sample_id
diversity_mn <- merge(diversity_mn, metanew, by="sample_id")
summary(diversity_mn$Observed)
sd(diversity_mn$Observed)
#Median (IQR) # spec: 82 (33, 232)
#Mean (SD): 191.7 (289.1)

#Plot # species by # reads
png(file="Num_spec by num_reads decontam 041922.png",
    width = 7, height = 4, units = 'in', res = 600)
ggplot(diversity_mn, aes(x=paired.reads, y=Observed)) + geom_point() +
  labs(title="Number observed species by # reads",x="Number of reads", y = "# Observed species")+
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, hjust = 0.5),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12))
dev.off()

#To zoom:
p <- ggplot(diversity_mn, aes(x=paired.reads, y=Observed)) + geom_point() +
  labs(title="#Observed species by # reads",x="Number of reads", y = "# Observed species")+
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, hjust = 0.5),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12))
p + coord_cartesian(xlim = c(0, 500000))
p + coord_cartesian(xlim = c(0, 50000), ylim = c(0,100))
png(file="Num_spec by num_reads 20K 041922.png",
    width = 7, height = 4, units = 'in', res = 600)
p + coord_cartesian(xlim = c(0, 20000), ylim = c(0,100))
dev.off()
png(file="Num_spec by num_reads 10K 041922.png",
    width = 7, height = 4, units = 'in', res = 600)
p + coord_cartesian(xlim = c(0, 10000), ylim = c(0,100))
dev.off()

remove(p)

#I don't really see a plateau; discussed with MK and based on 10K and 20K magnified figures, will prune at 2.5K (see 4/10 email)
#"It's a balance between losing samples and ensuring depth isn't influencing or driving your findings.
#"I think with the plots that you showed, you could prune at either 2500 or 5000 reads."

#Pruning
pruned = prune_samples(sample_sums(meta_new)>=2500, meta_new)
#how many ASVs do we have in our pruned sample? 
ntaxa(pruned)  #6709; unchanged from agglom, but we are down to 272 participants
#metadata for 272 participants
meta_prune <- data.frame(sample_data(pruned))

#save pruned pseq object
saveRDS(pruned, file = "phy.pruned05092022.rds", ascii = FALSE, version = NULL,
        compress = FALSE, refhook = NULL)

pruned <- readRDS("phy.pruned05092022.rds")
#confirm the taxtable, etc looks good:
test <- as.data.frame(as(tax_table(pruned),"matrix"),stringsAsFactors=FALSE)  #success!
meta_prune <- data.frame(sample_data(pruned))
meta_prune <- within(meta_prune, rm(sib_weight, sib_sex, sib_height, sib_muac, sib_race, sib_bfeed_1,
                                  sib_bfeed_2, sib_uri_recent, sib_uri_current, sib_pcr_1, sib_pcr_2,
                                  sib_bcg, sib_hepb, sib_clinic, sib_clinic_dx___5, sib_clinic_dx_oth,
                                  sib_hosp, sib_meds, sib_meds_name, sibmo, sibyr, sib_vl1_weeks, sib_vl_weeks2,
                                  sib_med_days, sib_dpt, sib_pcv, sib_rota, sib_polio, sib_measles, sib_age))
meta_prune <- within(meta_prune, rm(samp_storage, enr_clinic, enr_clinic_other, mat_np))

remove(meta_new, test)

#############################################
#DEMOGRAPHICS OF PARTICIPANTS WITH SEQUENCING
#############################################
#Thinking about Table 1: all HEU and HUU kids were age and sex matched
#start with general demographics first using pruned group (n=272)
child <- subset(meta_prune, subject == "child") #n=143
table(child$hiv)  #50 HUU kids, 49 HEU kids, 44 HIV+ kids

#Age
hist(child$age)
summary(child$age)
tapply(child$age, child$hiv, summary)
kruskal.test(child$age, child$hiv)

#Sex
table(child$sex)  #74 girls
table(child$sex, child$hiv)
prop.table(table(child$sex, child$hiv), 2)
summary(table(child$sex, child$hiv))

#Maternal age
hist(child$mat_age2)  #looks quite normal, but will still present median and IQR for consistency
summary(child$mat_age2)
tapply(child$mat_age2, child$hiv, summary)
kruskal.test(child$mat_age2, child$hiv)

#Season
table(child$season) 
table(child$season, child$hiv)
prop.table(table(child$season, child$hiv), 2)
summary(table(child$season, child$hiv))

#Maternal education
table(child$mat_educ)
table(child$mat_educ2)
table(child$mat_educ2, child$hiv)
prop.table(table(child$mat_educ2, child$hiv), 2)
summary(table(child$mat_educ2, child$hiv))
fisher.test(child$mat_educ2, child$hiv)

#Electricity
table(child$elec) #106 children have electricity
table(child$elec, child$hiv)
prop.table(table(child$elec, child$hiv), 2)
summary(table(child$elec, child$hiv))

#Wood
table(child$wood) #90 children in households using wood
table(child$wood, child$hiv)
prop.table(table(child$wood, child$hiv), 2)
summary(table(child$wood, child$hiv))

#Num household members
hist(child$hhsize)
summary(child$hhsize)
tapply(child$hhsize, child$hiv, summary)
kruskal.test(child$hhsize, child$hiv)

#Num children <5
hist(child$num_child)
summary(child$num_child)
tapply(child$num_child, child$hiv, summary)
kruskal.test(child$num_child, child$hiv)

#antibiotics
table(child$abx) #24 children received abx
table(child$abx, child$hiv)
prop.table(table(child$abx, child$hiv), 2)
summary(table(child$abx, child$hiv))

#Clinic in past 3 months
table(child$clinic)
table(child$clinic, child$hiv)
prop.table(table(child$clinic, child$hiv), 2)
summary(table(child$clinic, child$hiv))

#URI in past 1 month
table(child$uri_recent)
table(child$uri_recent, child$hiv)
prop.table(table(child$uri_recent, child$hiv), 2)
summary(table(child$uri_recent, child$hiv))

#URI currently (don't need both this and uri recent in table 1; need to choose)
table(child$uri_current)
table(child$uri_current, child$hiv)
prop.table(table(child$uri_current, child$hiv), 2)
summary(table(child$uri_current, child$hiv))

#Hospital in past 3 months -->> 4 children hospitalized, all 4 HIV+
  #*MENTION IN MANUSCRIPT
table(child$hosp)
table(child$hosp, child$hiv)
table(child$hosp_dx)  #1 with lower respiratory infx, 1 with gastroenteritis, 2 with other
table(child$hosp_dx_oth)  #malnutrition, sepsis + chronic suppurative OM

#PCV doses
table(child$pcv)
table(child$pcv, child$hiv)
prop.table(table(child$pcv, child$hiv), 2)
summary(table(child$pcv, child$hiv))
fisher.test(child$pcv, child$hiv)

#HiB
table(child$dpt)
table(child$dpt, child$hiv)
prop.table(table(child$dpt, child$hiv), 2)
summary(table(child$dpt, child$hiv))
fisher.test(child$dpt, child$hiv)

#breakdown of subjects (not for table 1, but for sequencing figure)
table(meta_prune$subject)

#######################
#HIV SPECIFIC VARIABLES
#######################
chi_hiv <- subset(child, hiv == "Infected")
#restrict to just HIV variables for easier cleaning and analysis
chi_hiv <- chi_hiv[c("child_pcr_1", "child_vl_1", "child_pcr_2", "child_vl_2", "chi_cd4_num_1", 
                     "chi_cd4_perc_1", "chi_cd4_num_2", "chi_cd4_perc_2", "chi_cd4_weeks1", 
                     "chi_cd4_weeks2", "child_vl1_weeks", "child_vl2_weeks",
                     "child_arv", "child_arv_meds", "child_tmpsmx", "abx")]

#How many kids receiving TMP-SMX? NOTE: we did not include these kids in the abx category. Need to re-do analyses with this included?
table(chi_hiv$child_tmpsmx) #15 kids on TMP-SMX
#manual inspection: 10 kids on TMP-SMX without other abx exposure

#I want the most recent viral load values for my analysis. Will try creating a new variable to ID which value is more recent first:
#missing 3 vl1 and 14vl2's
chi_hiv$newest_vL <- NA
chi_hiv$newest_vL [chi_hiv$child_vl1_weeks < chi_hiv$child_vl2_weeks] <- 1
chi_hiv$newest_vL [chi_hiv$child_vl2_weeks < chi_hiv$child_vl1_weeks] <- 2
chi_hiv$newest_vL [chi_hiv$child_vl2_weeks < 0] <- 1  #we want to include the viral loads collected BEFORE enrollment where possible
chi_hiv$newest_vL [is.na(chi_hiv$child_vl_2)] <- 1    #if we are missing vl_2, we will use vl_1
chi_hiv$newest_vL [is.na(chi_hiv$child_vl_1) & is.na(chi_hiv$child_vl_2)] <- NA
table(chi_hiv$newest_vL)

#Next, can we pull the viral load that was most recent using our newest_vL variable?
  #source: https://stackoverflow.com/questions/67944727/creating-new-variable-by-selecting-column-based-on-value-of-another-column
  #Need to do this with rows that have values (aka need to exclude the NA viral load children for this to work)
#For future phyloseq analyses: will need to remove B-04, B-07, B-52
chi_viral <- subset(chi_hiv, !is.na(newest_vL))
chi_viral$Var1 <- chi_viral$child_vl_1
chi_viral$Var2 <- chi_viral$child_vl_2

chi_viral <- chi_viral %>% rowwise() %>%
  mutate(newest_vlnum = get(paste0('Var', newest_vL)))
#now we have a dataframe of the 41 children with viral load measurements and 1 column containing the most recent viral loads for each
summary(chi_viral$newest_vlnum)
table(chi_viral$newest_vlnum) #27 with viral suppression, 14 without viral suppression

#make new variable for viral suppression or not, with vL < 400 = suppressed
chi_viral$vl_suppr <- NA
chi_viral$vl_suppr [chi_viral$newest_vlnum < 400] <- "1"
chi_viral$vl_suppr [chi_viral$newest_vlnum >= 400] <- "0"
table(chi_viral$vl_suppr)

#REPEAT the above steps to get date associated with most recent viral load, most recent CD4, and date a/w most recent CD4
#Viral load collection date
chi_viral$Var1 <- chi_viral$child_vl1_weeks
chi_viral$Var2 <- chi_viral$child_vl2_weeks

chi_viral <- chi_viral %>% rowwise() %>%
  mutate(newest_vldate = get(paste0('Var', newest_vL)))   #now have info for when vL were checked in relation to enrollment

summary(chi_viral$newest_vldate)

chi_viral <- within(chi_viral, rm(Var1, Var2))

#CD4
  #missing CD4 data for B-45 and B-52
  #will focus on percentages since this is what is presented for kids
chi_hiv$newest_cd4 <- NA
chi_hiv$newest_cd4 [chi_hiv$chi_cd4_weeks1 < chi_hiv$chi_cd4_weeks2] <- 1
chi_hiv$newest_cd4 [chi_hiv$chi_cd4_weeks2 < chi_hiv$chi_cd4_weeks1] <- 2
chi_hiv$newest_cd4 [chi_hiv$chi_cd4_weeks2 < 0] <- 1  #we want to include the viral loads collected BEFORE enrollment where possible
chi_hiv$newest_cd4 [is.na(chi_hiv$chi_cd4_perc_2)] <- 1    #if we are missing cd4_2, we will use cd4_1
chi_hiv$newest_cd4 [is.na(chi_hiv$chi_cd4_perc_1)] <- 2    #one child is missing cd4_1 but has cd4_2
chi_hiv$newest_cd4 [is.na(chi_hiv$chi_cd4_perc_1) & is.na(chi_hiv$chi_cd4_perc_2)] <- NA
table(chi_hiv$newest_cd4)

#create new dataframe without B-45 and B-52 (aka new df with children missing CD4 data removed)
chi_cd4 <- subset(chi_hiv, !is.na(newest_cd4))
chi_cd4$Var1 <- chi_cd4$chi_cd4_perc_1
chi_cd4$Var2 <- chi_cd4$chi_cd4_perc_2

chi_cd4 <- chi_cd4 %>% rowwise() %>%
  mutate(newest_cd4perc = get(paste0('Var', newest_cd4)))

summary(chi_cd4$newest_cd4perc)
table(chi_cd4$newest_cd4perc)

#make new variable for immune competent CD4 or not, with CD4 percentage *** = immune competent
  #See sub-analyses file for code

#arvs
table(chi_hiv$child_arv)
table(chi_hiv$child_arv_meds)

################################################
#NUMBER OF READS BY SUBJECT CLASS AND HIV STATUS
################################################
###PRUNED and agglomerated###
summary(meta_prune$paired.reads)   #Median (IQR) of 73,166 (29,605; 217,115) reads per sample
#mean 183,726 reads per sample
hist(meta_prune$paired.reads)

#Did maternal samples contain more reads than child or sib reads? --> the opposite! Why?
#children = petri dishes
tapply(meta_prune$paired.reads, meta_prune$subject, summary)
kruskal.test(meta_prune$paired.reads, meta_prune$subject)

#number of reads
hist(child$paired.reads)
summary(child$paired.reads)
tapply(child$paired.reads, child$hiv, summary)
kruskal.test(paired.reads ~ hiv, data = child) 

#Sum of sequences prior to filtering out singletons
sum(sample_sums(pruned))     #35,101,748 sequences obtained 
mean(sample_sums(pruned))    #Mean of 129,050.5 sequences per sample
median(sample_sums(pruned))  #Median of 45,235 sequences per sample
summary(sample_sums(pruned)) #IQR of sequencing reads per sample: 15730, 146428
nsamples(pruned)        #272 NP samples 
ntaxa(pruned) #6709 taxa

################
#ALPHA DIVERSITY
################
#Create dataframe with diversity indices for each specimen using PRUNED pseq object
diversity <- estimate_richness(pruned, measures = c("Shannon", "Simpson", "Chao1", "Observed"))

setDT(diversity, keep.rownames = TRUE)[]
colnames(diversity) <- c("sample_id", "Observed", "Chao1", "se.chao1", "Shannon", "Simpson")
rownames(diversity) <- diversity$sample_id
diversity <- merge(diversity, meta_prune, by="sample_id")

#
#OBSERVED SPECIES
#
#How many species on average per sample? Did that differ by subject or by HIV status (among children only)?
hist(diversity$Observed)
summary(diversity$Observed) 
sd(diversity$Observed)
#Median (IQR) # spec: 106 (53, 260)
#Mean (SD): 216 (300)

#what does it look like by subject? Moms have significantly fewer species per sample
tapply(diversity$Observed, diversity$subject, summary)   
kruskal.test(Observed ~ subject, data = diversity)

#what about by HIV status of the child? No signif diff
div_child <- subset(diversity, subject == "child")
tapply(div_child$Observed, div_child$hiv, summary)
kruskal.test(Observed ~ hiv, data = div_child) 

#
#SHANNON
#

#Is our SDI data normally distributed? 
shapiro.test(diversity$Shannon)  #technically no, because p<0.05; looks pretty good on histogram but QQ is wonky
histogram(diversity$Shannon)
qqnorm(diversity$Shannon, pch = 1, frame = FALSE)
qqline(diversity$Shannon, col = "steelblue", lwd = 2)

#Shannon diversity median and IQR for all samples
summary(diversity$Shannon)
sd(diversity$Shannon)
#Median (IQR) SDI: 1.85 (1.51, 2.26)
#Mean (SD) 1.93 (0.66)  

#UNPRUNED VALUES
summary(diversity_mn$Shannon)
sd(diversity_mn$Shannon)
  #Median (IQR) SDI: 1.74 (1.42, 2.17)
  #Mean (SD) 1.83 (0.67)

#what does it look like by subject? 
tapply(diversity$Shannon, diversity$subject, summary)
kruskal.test(Shannon ~ subject, data = diversity) #child signif lower (p=0.01)
#Since our data is normal-ish in distribution, what does the ANOVA look like?
subj.aov <- aov(Shannon ~ subject, data = diversity)
summary(subj.aov)   #P value 0.007; SIGNIFICANT
table(diversity$subject)

#what about by HIV status of the child?
tapply(div_child$Shannon, div_child$hiv, summary)
kruskal.test(Shannon ~ hiv, data = div_child)   #no significant difference by HIV status using KW test (P=0.71)
  #ANOVA
hiv.aov <- aov(Shannon ~ hiv, data = div_child)
summary(hiv.aov)  #P=0.37, not significant

#Age of the child
div_age <- lm(Shannon ~ age, data = div_child)
summary(div_age)  #p=0.98

#Wood smoke exposure
tapply(div_child$Shannon, div_child$wood, summary)
wilcox.test(Shannon ~ wood, data = div_child)   #p=0.30
t.test(Shannon ~ wood, data = div_child)        #P=0.73

#Season
tapply(div_child$Shannon, div_child$season, summary)
wilcox.test(Shannon ~ season, data = div_child)   #p=0.69
t.test(Shannon ~ season, data = div_child)  #P=0.70

#Receipt of abx in prior 3 months 
tapply(div_child$Shannon, div_child$abx, summary)
wilcox.test(Shannon ~ abx, data = div_child)   #p=0.86
t.test(Shannon ~ abx, data = div_child)   #P=0.71

#Current URI
tapply(div_child$Shannon, div_child$uri_cur2, summary)
wilcox.test(Shannon ~ uri_cur2, data = div_child) #p=0.65
t.test(Shannon ~ uri_cur2, data = div_child) #p=0.37

#Recent URI
tapply(div_child$Shannon, div_child$uri_rec2, summary)
wilcox.test(Shannon ~ uri_rec2, data = div_child) #p=0.30
t.test(Shannon ~ uri_rec2, data = div_child) #p=0.31

#Household members
summary(div_hh <- lm(Shannon ~ hhsize, data = div_child)) #p=0.47

#MULTIVARIABLE MODEL
#Need to set reference for HIV variable
class(div_child$hiv) #character; need to make factor to relevel
div_child$hiv <- as.factor(div_child$hiv)
div_child$hiv <- relevel(div_child$hiv, ref = "Unexposed")
summary(div_mv <- lm(Shannon ~ age + hiv + wood + season + abx + uri_rec2, data = div_child))
summary(div_mv <- lm(Shannon ~ age + hiv + wood + season + abx + uri_rec2 + hhsize, data = div_child))

#
#CHAO1
#

#Is our Chao1 data normally distributed? 
shapiro.test(diversity$Chao1)  #no, because p<0.05; also looks skewed on histogram and QQ is wonky
histogram(diversity$Chao1)
qqnorm(diversity$Chao1, pch = 1, frame = FALSE)
qqline(diversity$Chao1, col = "steelblue", lwd = 2)

#What does it look like if we log transform it? Will this allow for parametric testing?
diversity$Chao1_log <- log(diversity$Chao1)
shapiro.test(diversity$Chao1_log) 
histogram(diversity$Chao1_log)
qqnorm(diversity$Chao1_log, pch = 1, frame = FALSE)
qqline(diversity$Chao1_log, col = "steelblue", lwd = 2)
  #Looks very normal! Can use t.test and linear regression on transformed data
div_child <- subset(diversity, subject == "child")

#Chao1 median and IQR for all samples
summary(diversity$Chao1)
sd(diversity$Chao1)
#Median (IQR) Chao1: 106 (53, 260)
#Mean (SD) 216 (300.4)

#what does it look like by subject? 
tapply(diversity$Chao1, diversity$subject, summary)
tapply(diversity$Chao1_log, diversity$subject, summary) #log mean of moms is lower
summary(aov(Chao1_log ~ subject, data = diversity)) #significant difference (p=0.03)

#what about by HIV status of the child?
tapply(div_child$Chao1, div_child$hiv, summary)
tapply(div_child$Chao1_log, div_child$hiv, summary)
summary(aov(Chao1_log ~ hiv, data = div_child)) #no significant difference by HIV status using log values and ANOVA (p=0.2)

#Age of the child
summary(div_age <- lm(Chao1_log ~ age, data = div_child)) #p=0.38

#Wood smoke exposure
tapply(div_child$Chao1, div_child$wood, summary)
t.test(Chao1_log ~ wood, data = div_child)    #p=0.93

#Season
tapply(div_child$Chao1, div_child$season, summary)
t.test(Chao1_log ~ season, data = div_child)    #p=0.43

#Receipt of abx in prior 3 months 
tapply(div_child$Chao1, div_child$abx, summary)
t.test(Chao1_log ~ abx, data = div_child)    #p=0.92

#Current URI  #p=0.06
tapply(div_child$Chao1_log, div_child$uri_cur2, summary)
t.test(Chao1_log ~ uri_cur2, data = div_child)

#Recent URI 
tapply(div_child$Chao1_log, div_child$uri_rec2, summary)
t.test(Chao1_log ~ uri_rec2, data = div_child)  #p=0.69

#Household members
summary(div_hh <- lm(Chao1_log ~ hhsize, data = div_child)) #p=0.42

#MULTIVARIABLE MODEL
#log with HUU as reference value
class(div_child$hiv) #character; need to make factor to relevel
div_child$hiv <- as.factor(div_child$hiv)
div_child$hiv <- relevel(div_child$hiv, ref = "Unexposed")
summary(div_mv <- lm(Chao1_log ~ age + hiv + wood + season + abx + uri_rec2, data = div_child))
summary(div_mv <- lm(Chao1_log ~ age + hiv + wood + season + abx + uri_cur2, data = div_child))
summary(div_mv <- lm(Chao1_log ~ age + hiv + wood + season + abx + uri_rec2 + hhsize, data = div_child))

#Good boxplot resource: http://www.sthda.com/english/wiki/ggplot2-box-plot-quick-start-guide-r-software-and-data-visualization
#grayscale boxplots
  #SDI by subject
png(file="/Users/swetapatel/OneDrive - Duke University/Fogarty coding/Shannon by subject KAIJU_pruned.png",
    width = 4.5, height = 5, units = 'in', res = 600)
ggplot(diversity, aes(x=subject, y=Shannon)) + 
  geom_boxplot(fill="gray") + labs(title="Shannon diversity index by Subject",x="Subject", y = "Shannon diversity index")+
  theme_classic() +
  theme(
    plot.title = element_text(size = 16, hjust = 0.5),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.position = "right",
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14))
dev.off()

  #SDI by HIV 
# div_child$hiv <- as.factor(div_child$hiv)
div_child$hiv <- reorder(div_child$hiv, new.order=c("Infected", "Exposed Uninfected", "Unexposed"))
png(file="/Users/swetapatel/OneDrive - Duke University/Fogarty coding/Shannon by HIV KAIJU_pruned.png",
    width = 4.5, height = 5, units = 'in', res = 600)
ggplot(div_child, aes(x=hiv, y=Shannon)) + 
  geom_boxplot(fill="gray") + labs(title="Shannon diversity index by HIV status",x="HIV status", y = "Shannon diversity index")+
  theme_classic() +
  theme(
    plot.title = element_text(size = 16, hjust = 0.5),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.position = "right",
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14))
dev.off()

  #Chao1 by subject
png(file="/Users/swetapatel/OneDrive - Duke University/Fogarty coding/Chao1 by subject KAIJU_pruned.png",
    width = 4.5, height = 5, units = 'in', res = 600)
ggplot(diversity, aes(x=subject, y=Chao1)) + 
  geom_boxplot(fill="gray") + labs(title="Chao1 richness by Subject",x="Subject", y = "Chao1 richness")+
  theme_classic() +
  theme(
    plot.title = element_text(size = 16, hjust = 0.5),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.position = "right",
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14))
dev.off()

  #Chao1 by HIV
png(file="/Users/swetapatel/OneDrive - Duke University/Fogarty coding/Chao1 by HIV KAIJU_pruned.png",
    width = 5, height = 5, units = 'in', res = 600)
ggplot(div_child, aes(x=hiv, y=Chao1)) + 
  geom_boxplot(fill="gray") + labs(title="Chao1 richness by HIV status",x="HIV status", y = "Chao1 richness")+
  theme_classic() +
  theme(
    plot.title = element_text(size = 16, hjust = 0.5),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.position = "right",
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14))
dev.off()

#Alternative visualization of SDI
plot_richness(meta_new, x="subject", measures=c("Shannon"))

###UPDATE 7-4-22: redoing boxplots for manuscript per MK suggestion
###UPDATE 8-17-22: redoing boxplots as 1 supplemental figure using cowplot
#SDI by HIV 
#make new variable with HIV status abbreviated
div_child$hiv2 <- NA
div_child$hiv2 [div_child$hiv == "Infected"] <- "Children with HIV"
div_child$hiv2 [div_child$hiv == "Exposed Uninfected"] <- "HEU children"
div_child$hiv2 [div_child$hiv == "Unexposed"] <- "HUU children"
table(div_child$hiv2)
#matches up with table(relative_df$hiv)
div_child$hiv2 <- as.factor(div_child$hiv2)
div_child$hiv2 <- reorder(div_child$hiv2, new.order=c("Children with HIV", "HEU children", "HUU children"))

# png(file="/Users/swetapatel/OneDrive - Duke University/Fogarty coding/Shannon_HIV_070422.png",
    # width = 4.5, height = 5, units = 'in', res = 600)
s1a <- ggplot(div_child, aes(x=hiv2, y=Shannon, fill=hiv2)) + 
  geom_boxplot() + labs(x="HIV status", y = "Shannon diversity index")+
  scale_fill_manual(values=c('#800000FF','#155F83FF', '#FFA319FF')) +
  theme_classic() +
  theme(
    legend.text = element_text(size = 12),
    legend.position = "none",
    legend.title = element_blank(),
    axis.text = element_text(size = 10),
    axis.title.y = element_text(size = 12),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size=10, angle = 45, vjust = 1, hjust = 1, colour = "black")) 
# dev.off()

#what if we try with no x axis labels and legend only? 
s1a2 <- ggplot(div_child, aes(x=hiv2, y=Shannon, fill=hiv2)) + 
  geom_boxplot() + labs(x="HIV status", y = "Shannon diversity index")+
  scale_fill_manual(values=c('#800000FF','#155F83FF', '#FFA319FF')) +
  theme_classic() +
  theme(
    legend.text = element_text(size = 10),
    legend.position = "none",
    legend.title = element_blank(),
    axis.text = element_text(size = 10),
    axis.title.y = element_text(size = 12),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank())

#Chao1 by HIV
#8-17-22: with legend
# png(file="/Users/swetapatel/OneDrive - Duke University/Fogarty coding/Chao1_HIV_070422.png",
#     width = 5, height = 5, units = 'in', res = 600)
s1b <- ggplot(div_child, aes(x=hiv2, y=Chao1, fill=hiv2)) + 
  geom_boxplot() + labs(x="HIV status", y = "Chao1 richness")+
  scale_fill_manual(values=c('#800000FF','#155F83FF', '#FFA319FF')) +
  theme_classic() +
  theme(
    legend.text = element_text(size = 10),
    legend.position = "none",
    legend.title = element_blank(),
    axis.text = element_text(size = 10),
    axis.title.y = element_text(size = 12),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size=10, angle = 45, vjust = 1, hjust = 1, colour = "black"))
# dev.off()

#what if we try with no x axis labels and legend only? 
s1b2 <- ggplot(div_child, aes(x=hiv2, y=Chao1, fill=hiv2)) + 
  geom_boxplot() + labs(x="HIV status", y = "Chao1 richness")+
  scale_fill_manual(values=c('#800000FF','#155F83FF', '#FFA319FF')) +
  theme_classic() +
  theme(
    legend.text = element_text(size = 10),
    legend.position = "none",
    legend.title = element_blank(),
    axis.text = element_text(size = 10),
    axis.title.y = element_text(size = 12),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank())

#try getting just legend to make it easier to size all parts of the figure
legend <- get_legend(s1b2)
  # s1b + theme(legend.box.margin = (0, 0, 0, 12))
  # ) 


#combine the versions using cowplot:
#with x axis
prow <- plot_grid(s1a, s1b, labels = "AUTO", nrow = 1, ncol = 2, align = "h", axis = "b")
FigS1_v1 <- plot_grid(prow, legend, rel_widths = c(3, 1))
png(file="/Users/swetapatel/OneDrive - Duke University/Fogarty coding/FigS1_v1.png",
    width = 7, height = 4, units = 'in', res = 600)
FigS1_v1
dev.off()

png(file="/Users/swetapatel/OneDrive - Duke University/Fogarty coding/FigS1_v1.png",
    width = 7, height = 4, units = 'in', res = 600)
FigS1_v1
dev.off()

#without x axis
prow2 <- plot_grid(s1a2, s1b2, labels = "AUTO", nrow = 1, ncol = 2, align = "h", axis = "b", rel_widths = c(1,1))
FigS1_v2 <- plot_grid(prow2, legend, rel_widths = c(3, 1))

png(file="/Users/swetapatel/OneDrive - Duke University/Fogarty coding/FigS1_v2.png",
    width = 7, height = 4, units = 'in', res = 600)
FigS1_v2
dev.off()

#clean up environment
remove(div_age, div_mv, hiv.aov, m1, subj.aov)

###################
#CLR TRANSFORMATION
###################
#Need to generate centered log-ratio (CLR)-transformed sample counts
  #Sources: https://rdrr.io/github/microbiome/microbiome/man/transform.html
    #https://microbiome.github.io/tutorials/Preprocessing.html

pruned_clr <- microbiome::transform(pruned, "clr")

#####################
#BETA DIVERSITY PLOTS
#####################
#What do the non-transformed plots look like? 
#Subject
meta_pruned_nmds <- ordinate(pruned, method = "NMDS", distance = "bray")
lines <- c("solid", "solid", "solid")
bray_subj <- plot_ordination(pruned, meta_pruned_nmds, color = "subject", title = "Bray-Curtis NMDS Plot") +
  geom_point(size = 2) + theme_classic() + theme(plot.title = element_text(size = 14, hjust = 0.5)) +
  stat_ellipse(aes(linetype=subject), geom="polygon", alpha=0, type="t", level=0.8, size=0.7) +
  scale_linetype_manual(values=lines)
bray_subj

#HIV status (kids only)
chi_prune = subset_samples(pruned, subject == "child")
meta_chi_nmds <- ordinate(chi_prune, method = "NMDS", distance = "bray")
lines <- c("solid", "solid", "solid")
plot_ordination(chi_prune, meta_chi_nmds, color = "hiv", title = "Bray-Curtis NMDS Plot") +
  geom_point(size = 2) + theme_classic() + theme(plot.title = element_text(size = 14, hjust = 0.5))+
  stat_ellipse(aes(linetype=subject), geom="polygon", alpha=0, type="t", level=0.8, size=0.7) +
  scale_linetype_manual(values=lines)

#Transformed plots --> because we have negative values in our transformed data, we CAN'T use bray-curtis or jaccard
  #Also can't use NMDS plots
  #Need to generate a PCoA plot using euclidean distances 
#SUBJECT
clr_pcoa <- ordinate(pruned_clr, method = "PCoA", distance = "euclidean")
lines <- c("solid", "solid", "solid")
clr_subj <- plot_ordination(pruned_clr, clr_pcoa, color = "subject", title = "Euclidean PCoA Plot by Subject") +
  geom_point(size = 2) + theme_classic() + theme(plot.title = element_text(size = 14, hjust = 0.5)) +
  stat_ellipse(aes(linetype=subject), geom="polygon", alpha=0, type="t", level=0.8, size=0.7) +
  xlab("PC1 (13.4%)") + ylab("PC2 (10.6%)") +
  scale_linetype_manual(values=lines)
png(file="clr_subj.png",
    width = 5, height = 4, units = 'in', res=600)
clr_subj
dev.off()

#CHILD
#HIV status: need to relevel so legend is HUU, HEU, HEI
chi_clr = subset_samples(pruned_clr, subject == "child")
clr_pcoa <- ordinate(chi_clr, method = "PCoA", distance = "euclidean")
lines <- c("solid", "solid", "solid")
clr_hiv <- plot_ordination(chi_clr, clr_pcoa, color = "hiv", title = "Euclidean PCoA Plot by HIV Status") +
  geom_point(size = 1) + theme_classic() + theme(
    plot.title = element_text(size = 14, hjust = 0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 10),
    legend.position = "right",
    legend.background = element_blank(),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12)) +
  scale_color_manual(values=c('#800000FF','#155F83FF', '#FFA319FF'))+
  stat_ellipse(aes(linetype=hiv), geom="polygon", alpha=0, type="t", level=0.8, size=0.4) +
  xlab("PC1 (14.3%)") + ylab("PC2 (8.9%)") +
  scale_linetype_manual(values=lines)
png(file="clr_hiv.png",
    width = 5, height = 4, units = 'in', res=600)
clr_hiv
dev.off()

#trying out different options: dashed ellipses 
  #source: https://stackoverflow.com/questions/50293125/r-change-ggplot-legend-names-with-scale-linetype-manual
library(ggplot2)
png(file="clr_hiv.png",
    width = 5, height = 4, units = 'in', res=600)
plot_ordination(chi_clr, clr_pcoa, color = "hiv", title = "Euclidean PCoA Plot by HIV Status") +
  geom_point(size = 0.005) + theme_classic() + theme(
    plot.title = element_text(size = 14, hjust = 0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 10),
    legend.position = "right",
    legend.background = element_blank(),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12)) +
  scale_color_manual(values=c('#800000FF','#155F83FF', '#FFA319FF'))+
  stat_ellipse(aes(linetype=hiv), geom="polygon", alpha=0, type="t", level=0.8, size=0.5) +
  xlab("PC1 (14.3%)") + ylab("PC2 (8.9%)") +
  scale_linetype_manual(values= c("dashed", "dashed", "dashed"))
dev.off()

#differing shapes: makes lines go away
plot_ordination(chi_clr, clr_pcoa, shape = "hiv", title = "Euclidean PCoA Plot by HIV Status") +
  geom_point(aes(color=hiv), size = 1) + theme_classic() + theme(
    plot.title = element_text(size = 14, hjust = 0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 10),
    legend.position = "right",
    legend.background = element_blank(),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12)) +
  scale_color_manual(values=c('#800000FF','#155F83FF', '#FFA319FF'))+
  stat_ellipse(aes(linetype=hiv), geom="polygon", alpha=0, type="t", level=0.8, size=0.5) +
  xlab("PC1 (14.3%)") + ylab("PC2 (8.9%)") +
  scale_linetype_manual(values=lines)

#Plotting out the significant categorical variables from PERMANOVA: HIV plotted above, now will plot season and URI (current and recent)
lines <- c("solid", "solid")
clr_seas <- plot_ordination(chi_clr, clr_pcoa, color = "season", title = "Euclidean PCoA Plot by Season") +
  geom_point(size = 2) + theme_classic() + theme(plot.title = element_text(size = 14, hjust = 0.5)) +
  stat_ellipse(aes(linetype=season), geom="polygon", alpha=0, type="t", level=0.8, size=0.7) +
  xlab("PC1 (14.3%)") + ylab("PC2 (8.9%)") +
  scale_linetype_manual(values=lines)
png(file="clr_seas.png",
    width = 5, height = 4, units = 'in', res=600)
clr_seas
dev.off()

#Plot with dashed lines and new colors
lines <- c("dashed", "dashed")
clr_seas <- plot_ordination(chi_clr, clr_pcoa, color = "season", title = "Euclidean PCoA Plot by Season") +
  geom_point(size = 2) + theme_classic() + theme(
    plot.title = element_text(size = 14, hjust = 0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 10),
    legend.position = "right",
    legend.background = element_blank(),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12)) +
  scale_color_manual(values=c('#800000FF','#155F83FF'))+
  stat_ellipse(aes(linetype=season), geom="polygon", alpha=0, type="t", level=0.8, size=0.7) +
  xlab("PC1 (14.3%)") + ylab("PC2 (8.9%)") +
  scale_linetype_manual(values=lines)
png(file="clr_seas.png",
    width = 5, height = 4, units = 'in', res=600)
clr_seas
dev.off()

clr_uricur <- plot_ordination(chi_clr, clr_pcoa, color = "uri_cur2", title = "Euclidean PCoA Plot by Current URI Status") +
  geom_point(size = 2) + theme_classic() + theme(plot.title = element_text(size = 14, hjust = 0.5)) +
  stat_ellipse(aes(linetype=uri_cur2), geom="polygon", alpha=0, type="t", level=0.8, size=0.7) +
  xlab("PC1 (14.3%)") + ylab("PC2 (8.9%)") +
  scale_linetype_manual(values=lines)
png(file="clr_uricur.png",
    width = 5, height = 4, units = 'in', res=600)
clr_uricur
dev.off()

#Plot with dashed lines and new colors
clr_uricur <- plot_ordination(chi_clr, clr_pcoa, color = "uri_cur2", title = "Euclidean PCoA Plot by Current URI Status") +
  geom_point(size = 2) + theme_classic() + theme(
    plot.title = element_text(size = 14, hjust = 0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 10),
    legend.position = "right",
    legend.background = element_blank(),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12)) +
  scale_color_manual(values=c('#800000FF','#155F83FF'))+
  stat_ellipse(aes(linetype=uri_cur2), geom="polygon", alpha=0, type="t", level=0.8, size=0.7) +
  xlab("PC1 (14.3%)") + ylab("PC2 (8.9%)") +
  scale_linetype_manual(values=lines)
png(file="clr_uricur.png",
    width = 5, height = 4, units = 'in', res=600)
clr_uricur
dev.off()

clr_urirec <- plot_ordination(chi_clr, clr_pcoa, color = "uri_rec2", title = "Euclidean PCoA Plot by Recent URI Status") +
  geom_point(size = 2) + theme_classic() + theme(plot.title = element_text(size = 14, hjust = 0.5)) +
  stat_ellipse(aes(linetype=uri_rec2), geom="polygon", alpha=0, type="t", level=0.8, size=0.7) +
  xlab("PC1 (14.3%)") + ylab("PC2 (8.9%)") +
  scale_linetype_manual(values=lines)
png(file="clr_urirec.png",
    width = 5, height = 4, units = 'in', res=600)
clr_urirec
dev.off()

#Plot with dashed lines and new colors
clr_urirec <- plot_ordination(chi_clr, clr_pcoa, color = "uri_rec2", title = "Euclidean PCoA Plot by Recent URI Status") +
  geom_point(size = 2) + theme_classic() + theme(
    plot.title = element_text(size = 14, hjust = 0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 10),
    legend.position = "right",
    legend.background = element_blank(),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12)) +
  scale_color_manual(values=c('#800000FF','#155F83FF'))+
  stat_ellipse(aes(linetype=uri_rec2), geom="polygon", alpha=0, type="t", level=0.8, size=0.7) +
  xlab("PC1 (14.3%)") + ylab("PC2 (8.9%)") +
  scale_linetype_manual(values=lines)
png(file="clr_urirec.png",
    width = 5, height = 4, units = 'in', res=600)
clr_urirec
dev.off()

#What is a good way to display clustering by age? 
bray_age <- plot_ordination(meta_chi, ord, color="age") +
  geom_point(size=2) + theme_classic() +
  stat_ellipse(aes(linetype=age), geom="polygon", alpha=0, type="t", level=0.8, size=0.7) +
  xlab("PC1 (19.9%)") + ylab("PC2 (10.7%)") +
  scale_linetype_manual(values=lines)
bray_age
png(file="R Plots/Figure_2a.png",
    width = 10, height = 8.3, units = 'in', res=600)
bray_age
dev.off()

##UPDATE AUG 3: REVISING FIGURES, COMBINING USING COWPLOT FOR:
  #HIV PCOA, HIV rel abundance, CD4 PCOA, CD4 rel abundance
  #make each figure first and label as p1 - p4, then use cowplot to organize into grids
#For HIV PCoA, need to: remove title, change legend headings, make plot square, change order of HIV cat
#Try creating new HIV var as factor, then relevel
table(child$hiv)
child$hiv2 <- NA
child$hiv2 [child$hiv == "Infected"] <- "Children with HIV"
child$hiv2 [child$hiv == "Exposed Uninfected"] <- "HEU children"
child$hiv2 [child$hiv == "Unexposed"] <- "HUU children"
table(child$hiv2)
class(child$hiv2)
child$hiv2 <- as.factor(child$hiv2)
child$hiv2 <- reorder(child$hiv2, new.order=c("Children with HIV", "HEU children", "HUU children"))

#Import back into pseq (source: https://jacobrprice.github.io/2017/08/26/phyloseq-to-vegan-and-back.html)
sample_data(chi_clr) <- as.data.frame(child)
#now see if it worked:
bloop <- data.frame(sample_data(chi_clr))
#success!
remove(bloop)

# png(file="clr_hiv.png",
#     width = 5, height = 4, units = 'in', res=600)
p1 <- plot_ordination(chi_clr, clr_pcoa, color = "hiv2") +
  geom_point(size = 0.005) + theme_classic() + theme(
    plot.title = element_text(size = 14, hjust = 0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 10),
    legend.position = "right",
    legend.background = element_blank(),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12)) +
  scale_color_manual(values=c('#800000FF','#155F83FF', '#FFA319FF'))+
  stat_ellipse(aes(linetype=hiv2), geom="polygon", alpha=0, type="t", level=0.8, size=0.5) +
  xlab("PC1 (14.3%)") + ylab("PC2 (8.9%)") +
  scale_linetype_manual(values= c("dashed", "dashed", "dashed"))
# dev.off()

#Clean up environment
remove(clr_hiv, clr_pcoa, clr_seas, clr_subj, clr_uricur, clr_urirec) 

##########
#PERMANOVA
##########
otu <- as.data.frame(otu_table(pruned_clr, taxa_are_rows=TRUE))
#Need species in columns and observations in rows: transpose using t() function and save as a dataframe
#Source: https://github.com/edamame-course/2015-tutorials/blob/master/final/Ashley_Intro_to_R.md
otu_trans <- data.frame(t(otu))

#Calculate distance and save as a matrix
BC.dist=vegdist(otu_trans, method="euclidean")
#Run PERMANOVA on distances, unadjusted
adonis(BC.dist ~ hiv, data = meta_prune, permutations = 1000)
  #Significant difference in overall composition by hiv status (P = 0.004) among all participants; R^2 = 0.013

#What about moms vs kids vs sibs?
adonis(BC.dist ~ subject, data = meta_prune, permutations = 1000)
  #Significant difference in overall composition by subject classification (P = 0.001, R^2 = 0.081)

#Restricting to just the kids: must redo all steps above for calculating distances
otu <- as.data.frame(otu_table(chi_clr, taxa_are_rows=TRUE))
#Need species in columns and observations in rows: transpose using t() function and save as a dataframe
#Source: https://github.com/edamame-course/2015-tutorials/blob/master/final/Ashley_Intro_to_R.md
otu_trans <- data.frame(t(otu))
child <- data.frame(sample_data(chi_clr))
BC.dist=vegdist(otu_trans, method="euclidean")
#HIV
adonis(BC.dist ~ hiv, data = child, permutations = 1000)
  #Significant difference by HIV status among children (P = 0.043, R^2 = 0.019)
  #HIV + age (per 6/9 MK discussion to account for age in all PERMANOVAs)
adonis(BC.dist ~ hiv + age, data = child, permutations = 1000)
  #Significant difference by HIV status among children (P = 0.042, R^2 = 0.019)
    #Significant difference by age as well (P = 0.03, R^2 = 0.012)

#Age
adonis(BC.dist ~ age, data = child, permutations = 1000)
  #Significant difference in overall composition by age (P = 0.03, R2 = 0.012)

#Season
  #Significant difference in overall composition by season (P = 0.005, R2 = 0.016)
adonis(BC.dist ~ season, data = child, permutations = 1000)
  #Season + age (per 6/9 MK discussion to account for age in all PERMANOVAs)
adonis(BC.dist ~ season + age, data = child, permutations = 1000)
#Significant difference by season (P = 0.002, R^2 = 0.016)
#Significant difference by age as well (P = 0.02, R^2 = 0.011)

#Wood
adonis(BC.dist ~ wood, data = child, permutations = 1000)
  #No significant difference by wood smoke exposure among children (P = 0.10)
#Wood + age (per 6/9 MK discussion to account for age in all PERMANOVAs)
adonis(BC.dist ~ wood + age, data = child, permutations = 1000)
#No significant difference by wood smoke exposure (P = 0.11)
#Significant difference by age  (P = 0.02, R^2 = 0.012)

#Abx
adonis(BC.dist ~ abx, data = child, permutations = 1000)
  #Trend but no significant difference in overall composition by abx exposure (P = 0.08)
#Abx + age (per 6/9 MK discussion to account for age in all PERMANOVAs)
adonis(BC.dist ~ abx + age, data = child, permutations = 1000)
#Trend but no significant difference in overall composition by abx exposure (P = 0.08)
#Significant difference by age as well (P = 0.03, R^2 = 0.012)


#Recent/current URI?
adonis(BC.dist ~ uri_cur2, data = child, permutations = 1000)
#Significant difference in composition among children with current URI (P = 0.001, R2=0.023)
  #Current URI + age (per 6/9 MK discussion to account for age in all PERMANOVAs)
adonis(BC.dist ~ uri_cur2 + age, data = child, permutations = 1000)
#Significant difference in composition among children with current URI (P = 0.001, R2=0.023)
#Significant difference by age as well (P = 0.03, R^2 = 0.011)

adonis(BC.dist ~ uri_rec2, data = child, permutations = 1000)
#Significant difference in composition among children with recent URI (P = 0.006, R2=0.015)
  #Recent URI + age (per 6/9 MK discussion to account for age in all PERMANOVAs)
adonis(BC.dist ~ uri_rec2 + age, data = child, permutations = 1000)
#Significant difference in composition among children with current URI (P = 0.003, R2=0.015)
#Significant difference by age as well (P = 0.03, R^2 = 0.011)

#Multivariable
  #NOTE: order of variables affects results: https://stats.stackexchange.com/questions/188519/adonis-in-vegan-order-of-variables-or-use-of-strata
    #https://stats.stackexchange.com/questions/312302/adonis-in-vegan-order-of-variables-non-nested-with-one-degree-of-freedom-for
#one approach: include only significant variables in MV analysis and vary the order to see what remains significant
  #alternative: strata (ex difference by HIV status among children enrolled within rainy season)
    #can't do an age strata right? unless created ordinal or binary variable
#start with model that lists covariables by effect size from univariable analysis (signif variables only)
  #HIV R2 = 0.019, season R2 = 0.016, recent URI R2=0.015, age: R2 = 0.012
adonis(BC.dist ~ hiv + season + uri_rec2 + age, data = child, permutations = 1000)
  #In a multivariable model, HIV status (P=0.04, R2=0.019), season (P=0.003, R2=0.015), recent URI (P=0.005, R2=0.014), and age (P=0.04, R2=0.011) were a/w differences in overall microbiome composition

#How does this compare to listing the variables in the opposite order (smallest -> largest R2)
adonis(BC.dist ~ age + uri_rec2 + season + hiv, data = child, permutations = 1000)  #interestingly not that different
  #Age (P=0.03, R2=0.012), recent URI (P=0.006, R2=0.014), season (P=0.003, R2=0.016), and HIV status (P=0.046, R2=0.019) signif associated with composition

#*stick with results of the first model (largest -> smallest R2)*

#We have a number of factors significantly a/w composition --> reliable?
  #Source: https://microbiome.github.io/tutorials/PERMANOVA.html
#QUESTION: is this still done with transformed data? If so, age and HIV have significant P-values
#Assuming that the non-significant P value for the F statistic implies that variance is similar (enough) across all groups
anova(betadisper(BC.dist, child$age))
anova(betadisper(BC.dist, child$season))
anova(betadisper(BC.dist, child$hiv))
anova(betadisper(BC.dist, child$uri_rec2))

#Clean up environment
remove(otu, otu_trans, pruned_clr, chi_clr)

##############################################################################
#CREATING DATAFRAMES FOR VISUALIZING SAMPLE COMPOSITION (NON-TRANSFORMED DATA)
##############################################################################
#Goal: generate stacked barplot of most common taxa by HIV, etc --> using NON TRANSFORMED pruned data
#(AGGREGATED AT THE SPECIES LEVEL)
# Create dataframe with overall relative abundances of OTUs (needs to be done before agglomeration --> does it? we agglom'd by species a while back)
#Using chi_pruned dataframe
chi_prune = subset_samples(pruned, subject == "child")
phy.otus <- transform_sample_counts(chi_prune, function(Abundance) Abundance/sum(Abundance))
otu_df <- psmelt(phy.otus)
otu_df$OTU <- as.character(otu_df$OTU)
otu_abundances <- aggregate(otu_df$Abundance, by=list(OTU=otu_df$OTU), FUN=sum)
otu_abundances$x <- (otu_abundances$x)/(nsamples(chi_prune))
names(otu_abundances)[names(otu_abundances) == 'x'] <- 'otu_Ab'
sum(otu_abundances$otu_Ab)  #sums to 1
nrow(otu_abundances)        #6709 species
remove(phy.otus, otu_df)

# Agglomerate taxa into species --> already done, so no need to do again
# Transform to relative abundances
chi.rel <- transform_sample_counts(chi_prune, function(Abundance) Abundance/sum(Abundance))
sample_sums(chi.rel)  # This is a sanity check to make sure that relative abundance was calculated for each sample prior to pooling data
#All samples add to 1
melted_df <- psmelt(chi.rel)
#clean up melted_df
melted_df <- within(melted_df, rm(sib_weight, sib_sex, sib_height, sib_muac, sib_race, sib_bfeed_1,
                                    sib_bfeed_2, sib_uri_recent, sib_uri_current, sib_pcr_1, sib_pcr_2,
                                    sib_bcg, sib_hepb, sib_clinic, sib_clinic_dx___5, sib_clinic_dx_oth,
                                    sib_hosp, sib_meds, sib_meds_name, sibmo, sibyr, sib_vl1_weeks, sib_vl_weeks2,
                                    sib_med_days, sib_dpt, sib_pcv, sib_rota, sib_polio, sib_measles, sib_age))
melted_df <- within(melted_df, rm(samp_storage, enr_clinic, enr_clinic_other, mat_np))

# Create dataframes with overall relative abundances of phyla and genera
melted_df$Phylum <- as.character(melted_df$Phylum)
phyla_abundances <- aggregate(melted_df$Abundance, by=list(Phylum=melted_df$Phylum), FUN=sum)
phyla_abundances$x <- (phyla_abundances$x)/(nsamples(chi.rel))
names(phyla_abundances)[names(phyla_abundances) == 'x'] <- 'phyla_Ab'
sum(phyla_abundances$phyla_Ab)    # Should sum to 1 (sum of relative abundances of phyla)
nrow(phyla_abundances)            # Corresponds to # of unique phyla: 116
melted_df$Genus <- as.character(melted_df$Genus)
genus_abundances <- aggregate(melted_df$Abundance, by=list(Phylum=melted_df$Phylum, Genus=melted_df$Genus,
                                                           OTU=melted_df$OTU), FUN=mean)
genus_abundances <- aggregate(genus_abundances$x, by=list(Phylum=genus_abundances$Phylum,
                                                          Genus=genus_abundances$Genus), FUN=sum)
names(genus_abundances)[names(genus_abundances) == 'x'] <- 'genus_Ab'
sum(genus_abundances$genus_Ab)    # Should sum to 1 (sum of relative abundances of genera)
nrow(genus_abundances)            # Corresponds to # of unique genera: 2146
abundances <- merge(genus_abundances, phyla_abundances, by="Phylum")

melted_df$Species <- as.character(melted_df$Species)
species_abundances <- aggregate(melted_df$Abundance, by=list(Phylum=melted_df$Phylum, Genus=melted_df$Genus,Species=melted_df$Species,
                                                           OTU=melted_df$OTU), FUN=mean)
species_abundances <- aggregate(species_abundances$x, by=list(Phylum=species_abundances$Phylum,
                                                              Genus=species_abundances$Genus, Species=species_abundances$Species), FUN=sum)
names(species_abundances)[names(species_abundances) == 'x'] <- 'species_Ab'
sum(species_abundances$species_Ab)    # Should sum to 1 (sum of relative abundances of genera)
nrow(species_abundances)            # Corresponds to # of unique species: 6708 ***WHY is this different from 6709?###
abundances <- merge(abundances, species_abundances, by=c("Phylum", "Genus"))

phyla_abundances <- arrange(phyla_abundances, desc(phyla_Ab)) 
# Create a df with the top 5 phyla
TOPPhyla <- unique(phyla_abundances$Phylum[1:5])
phylum_df <- phyla_abundances[phyla_abundances$Phylum %in% TOPPhyla,]
phylum_df$Phylum <- factor(phylum_df$Phylum, levels = phylum_df$Phylum[order(-phylum_df$phyla_Ab)])

genus_abundances <- arrange(genus_abundances, desc(genus_Ab))  
TOPGenera <- unique(genus_abundances$Genus[1:10])
genus_df <- genus_abundances[genus_abundances$Genus %in% TOPGenera,]
genus_df$Genus <- factor(genus_df$Genus, levels = genus_df$Genus[order(-genus_df$genus_Ab)])

#create df with top 10 species
species_abundances <- arrange(species_abundances, desc(species_Ab))  
TOPSpecies <- unique(species_abundances$Species[1:10])
species_df <- species_abundances[species_abundances$Species %in% TOPSpecies,]
species_df$Species <- factor(species_df$Species, levels = species_df$Species[order(-species_df$species_Ab)])

# Rename species other than top10 as "Other" in creating dataframe relative_df
# Rename genera other than top10 as "Other" in creating dataframe relative_df 
relative_df <- merge(melted_df, abundances, by=c("Genus", "Species"))
#keep getting vector memory exhausted error; need to reboot R so will save melted_df and abundances as csv files
  #source for fix: https://stackoverflow.com/questions/51295402/r-on-macos-error-vector-memory-exhausted-limit-reached
# write.csv(melted_df, 'melted_df.csv')
# write.csv(abundances, 'abundances.csv')
# abundances <- read.csv("abundances.csv")
# melted_df <- read.csv("melted_df.csv")
#note: using the csv files led to very weird errors (sample sums 0.99 instead of 1, 1 phyla missing), 
  #so fixed memory error and reran the code instead of using the abundances and melted_df csv files

relative_df$Species[!(relative_df$Species %in% TOPSpecies)] <- "Other"
sum(relative_df$Abundance)#sums to 143.0008; n=143 so tracks?
table(relative_df$Species) #We have 11 levels for the Species variable, which are consistent with our top 10 + other

#Clean up relative_df and save as a csv file so we can come back to it later
write.csv(relative_df, 'rel_df.csv')

#########################
#RELATIVE ABUNDANCE PLOTS
#########################
#To generate color palette of 11 colors used previously: 
library(RColorBrewer)
library(scales)
chisp2 <- c("#293352", "#800000FF", "#767676FF", "#FFA319FF", "#8A9045FF", "#155F83FF", "#C16622FF", "#4D004B", "#8F3931FF", "#58593FFF", "#350E20FF")
show_col(chisp2)

#Creating relative abundance figure: HIV
#CREATING RELATIVE ABUNDANCE PLOT
  #INCORPORATING CHANGES FROM AUG 3: change Y axis to mean rel abundance, update categories for consistency with PCoA plot
#source for rotating x labels: https://stackoverflow.com/questions/1330989/rotating-and-spacing-axis-labels-in-ggplot2
relative_df$subject_abundances <- (relative_df$Abundance)/(length(unique(relative_df$Sample))) 
#abbreviate HIV status
relative_df$hiv2 <- NA
relative_df$hiv2 [relative_df$hiv == "Infected"] <- "Children with HIV"
relative_df$hiv2 [relative_df$hiv == "Exposed Uninfected"] <- "HEU children"
relative_df$hiv2 [relative_df$hiv == "Unexposed"] <- "HUU children"
table(relative_df$hiv2)
#matches up with table(relative_df$hiv)
relative_df$hiv2 <- as.factor(relative_df$hiv2)
relative_df$hiv2 <- reorder(relative_df$hiv2, new.order=c("Children with HIV", "HEU children", "HUU children"))
relative_df$Species <- factor(relative_df$Species)
relative_df$Species <- reorder(relative_df$Species, new.order=c("Corynebacterium accolens", "Corynebacterium propinquum", "Corynebacterium pseudodiphtheriticum", "Dolosigranulum pigrum",  
                                                                "Haemophilus influenzae", "Moraxella catarrhalis", "Moraxella lincolnii", "Moraxella nonliquefaciens", 
                                                                "Staphylococcus aureus", "Streptococcus pneumoniae", "Other"))
subject_plot2 <- ggplot(relative_df[order(relative_df$Species, decreasing = TRUE),], aes(x=hiv2, y=subject_abundances, fill=Species)) + 
  geom_bar(stat="identity", position="fill") + theme_classic() + theme(legend.position="right", 
                                                                       legend.text=element_text(size=10, face="italic"),
                                                                       axis.title.y = element_text(size=12, margin=margin(0,20,0,0)), axis.text.y = element_text(size=10), 
                                                                       axis.title.x = element_blank(), axis.text.x = element_text(size=10, angle = 45, vjust = 1, hjust = 1, colour = "black")) + 
  scale_fill_manual(values = chisp2) +
  # scale_x_discrete(labels=c("No", "Yes")) +
  xlab("HIV Status") + ylab("Mean relative abundance") 

#element_text(size=12, margin=margin(10,0,0,0))

p2 <- plot(subject_plot2)
p2
# png(file="genplot_kaiju_05122022.png",
#     width = 9, height = 5, units = 'in', res = 600)
# plot(subject_plot2)
# dev.off()

###UPDATE AUG 3: try cowplot with just p1 and p2
plot_grid(p1, p2, labels = c('A', 'B'), label_size = 12)
#for p3 and p4: ran code from file "Fog Kaiju sub-analyses 5-26-22.R"

################
#PLOTS BY SAMPLE [SKIPPED 5/4/22]
################
#Stacked barplots using chi_pruned dataframe
  #Source (again): https://mvuko.github.io/meta_phyloseq/#321_ordination
  #Source merges counts together (since each ID has multiple entries) --> skip this step as each study ID has only one sample_id associated
#check taxonomy level names
rank_names(meta_chi_rel)
meta_chi_rel  #137 taxa; as a test, will agglom by genus and see what happens
meta_chi_rel_genus<- tax_glom(meta_chi_rel, taxrank = "Genus")
#check to how many genera were these assigned to
meta_chi_rel_genus  #56 taxa
#check sample sums to make sure nothing was deleted due to NAs (should sum up 100%)
sample_sums(meta_chi_rel_genus)   #all sums to 100!
plot_bar(meta_chi_rel_genus, fill = "Genus")  #too many genera
#Since 56 genera are too many to plot in a stacked barplot, we will filter the low abundant ones and put them into one category. 
  #To do this, we will use the power of tidyverse again. 
  #First, we will create a normal data frame out of the phyloseq object and then add another column where all taxa with abundance lower than 1% will be renamed to “< 1%”.

#transform phyloseq object to a data frame
chigenus<- psmelt(meta_chi_rel_genus)
str(chigenus)

#make the genera characters, not factors
chigenus$Genus<- as.character(chigenus$Genus)

#add new column with renamed low abundant taxa
chigenus<- chigenus %>% 
  mutate(Genus2 = replace(Genus, Abundance < 1, "< 1%"))

#check all phyla names
unique(chigenus$Genus2)

#reorder the phyla so that they are stacked according to abundance
chigenus$Genus2<- reorder(chigenus$Genus2, chigenus$Abundance)

ggplot(chigenus, aes(Sample, Abundance, fill=Genus2)) +
  geom_bar(stat = "identity") +
  labs(y= "Relative abundance [%]",
       fill= "Genera") +
  theme_bw()

#that's a lot of genera; what if we clump all genera that are <5% together? 
chigenus<- chigenus %>% 
  mutate(Genus3 = replace(Genus, Abundance < 5, "< 5%"))
unique(chigenus$Genus3)   #still 18 different genera; what about <10%?

chigenus<- chigenus %>% 
  mutate(Genus4 = replace(Genus, Abundance < 10, "< 10%"))
unique(chigenus$Genus4) #16 genera

chigenus<- chigenus %>% 
  mutate(Genus5 = replace(Genus, Abundance < 20, "< 20%"))
unique(chigenus$Genus5) #15 genera; will stick with genus 4 and see what that looks like

#reorder the phyla so that they are stacked according to abundance
chigenus$Genus4<- reorder(chigenus$Genus4, chigenus$Abundance)
  #get unique set of colors: https://medialab.github.io/iwanthue/
ggplot(chigenus, aes(Sample, Abundance, fill=Genus4)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#a256c1",
                               "#5db645",
                               "#d14799",
                               "#5bbe7d",
                               "#646bd0",
                               "#b4b536",
                               "#638dca",
                               "#dd9a4c",
                               "#c786c5",
                               "#427a3f",
                               "#dc6072",
                               "#4bb9ac",
                               "#c95535",
                               "#8f9b4b",
                               "#9c4864",
                               "#8f6c30")) +
  labs(y= "Relative abundance [%]",
       fill= "Genera") +
  theme_bw()

#since we have metagenomic data, what if we look at species and clump everything <20% together?
meta_chi_rel_species<- tax_glom(meta_chi_rel, taxrank = "Species")
#check to how many genera were these assigned to
meta_chi_rel_species  #137 taxa
#check sample sums to make sure nothing was deleted due to NAs (should sum up 100%)
sample_sums(meta_chi_rel_species)   #all sums to 100!
#Since 137 species are too many to plot in a stacked barplot, we will filter the low abundant ones and put them into one category. 
#To do this, we will use the power of tidyverse again. 
#First, we will create a normal data frame out of the phyloseq object and then add another column where all taxa with abundance lower than 1% will be renamed to “< 1%”.

#transform phyloseq object to a data frame
chisp<- psmelt(meta_chi_rel_species)

#make the species characters, not factors
chisp$Species<- as.character(chisp$Species)

#add new column with renamed low abundant taxa
chisp<- chisp %>% 
  mutate(Species2 = replace(Species, Abundance < 50, "< 50%"))
unique(chisp$Species2)  #13 categories
#reorder the species so that they are stacked according to abundance
chisp$Species2<- reorder(chisp$Species2, chisp$Abundance)

species <- ggplot(chisp, aes(Sample, Abundance, fill=Species2)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#c688ca",
                               "#68b647",
                               "#be51b8",
                               "#beae47",
                               "#6f66c9",
                               "#657a35",
                               "#d9457a",
                               "#56b189",
                               "#cc4d34",
                               "#5b97cf",
                               "#c0803e",
                               "#9d486d",
                               "#d47573")) +
  labs(y= "Relative abundance [%]",
       fill= "Species") +
  theme(axis.text.x = element_text(angle = 90))
  theme_bw()
  
species

png(file="Species1.png",
    width = 15, height = 5, units = 'in', res=600)
species
dev.off()


################
#MAASLIN2 MODELS
################
#Maaslin2 tutorial: https://github.com/biobakery/biobakery/wiki/maaslin2
# if(!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("Maaslin2")

###NOTE: I ran several different models with varying parameters that can be found in the prior version of this code (Fog KAIJU 4-10-22 - decontam.R)
  #deleted here to save space
#what if we filter WAY DOWN to just the top 50 species and re-run our model?
chi.filter3 <- filter_taxa(chi_prune, function(Abundance) mean(Abundance)>=48.2, TRUE)
ntaxa(chi.filter3) #50

#UPDATE per d/w MK 5/21: filter based on round abundance of 50 and re-run the model
chi.filter3 <- filter_taxa(chi_prune, function(Abundance) mean(Abundance)>=50, TRUE)
ntaxa(chi.filter3) #43

otu <- as.data.frame(otu_table(chi.filter3, taxa_are_rows=TRUE))
#need to remove the unidentified species from here: row 11, 25, 37, 38, 40
  #source: https://stackoverflow.com/questions/12328056/how-do-i-delete-rows-in-a-data-frame
#corresponding respectively to: Corynebacterium sp. KPL1859, Streptococcus sp. SK643, Alkalibacterium sp. 20, uncultured Clostridium sp., Nocardioides sp. Root122
otu <- otu[-c(11, 25, 37, 38, 40),]
fit_data2 = Maaslin2(
  input_data = otu, 
  input_metadata = child, 
  min_prevalence = 0.2,
  max_significance = 0.2,
  output = "maaslin_mv_kaiju_filt50_v4", 
  fixed_effects = c("hiv", "season", "age", "abx", "wood", "uri_rec2"),
  reference = c("hiv,Unexposed"))
#NOTE: maaslin_mv_kaiju_filt50_v2 = results based on chi.filter with abundance >47
#NOTE: maaslin_mv_kaiju_filt50_v3 = results based on chi.filter with abundance >48.2 to give us the top 50 species AND with NOS species removed from OTU table (so n=45 identified species)
#NOTE: maaslin_mv_kaiju_filt50_v4 = results based on chi.filter with abundance >50 AND with NOS species removed from OTU table (so n=38 identified species)

#UPDATE AUG 3: use Maaslin2 results to create barplot
maaslin <- read.csv("maaslin_mv_kaiju_filt50_v4/significant_results.tsv", sep = "\t")
#clean up data: remove "." from species name and make names factors so we can reorder them
class(maaslin$feature)
maaslin$feature <- gsub(".", " ", maaslin$feature, fixed = TRUE)
maaslin$feature <- gsub("sp  KPL1995", "pseudodiphtheriticum", maaslin$feature, fixed = TRUE)
maaslin$feature <- as.factor(maaslin$feature)
#rename feature to species
maaslin <- rename(maaslin, Species = 'feature')
#create new variable for positive or neg coefficient that we can use to add color to our figures
maaslin$posneg <- NA
maaslin$posneg [maaslin$coef >= 0] <- "Pos"
maaslin$posneg [maaslin$coef < 0] <- "Neg"
table(maaslin$posneg)

#will create 4 horizontal barplots and use cowplot to align them into 1 figure
#MAKE NEW VAR for coef + or neg and then use that for your barplot fill
#season
m_seas <- subset(maaslin, metadata == "season")
#relevel so all negative values grouped together
m_seas$Species <- reorder(m_seas$Species, new.order=c("Micrococcus luteus", "Kocuria polaris", "Moraxella bovis", "Moraxella caprae",
                                                      "Moraxella catarrhalis", "Moraxella caviae",  "Moraxella cuniculi",
                                                      "Moraxella lacunata", "Moraxella lincolnii", "Moraxella macacae", 
                                                      "Moraxella nonliquefaciens",  "Moraxella oblonga", "Moraxella ovis"))

m_s <- ggplot(m_seas, aes(x=coef, y=Species, fill=posneg)) +
  geom_bar(stat="identity") +
  ggtitle("Rainy season") +
  scale_fill_manual(values=c('#800000FF', '#155F83FF')) +
  scale_x_continuous(breaks = seq(-6, 4, 1), limits = c(-5, 3)) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 12, hjust = 0.5, face = "bold"),
    # legend.text = element_text(size = 12),
    legend.position = "none",
    # legend.title = element_blank(),
    axis.text = element_text(size = 10),
    axis.text.y = element_text(face = "italic"),
    axis.title = element_blank()) +
  geom_errorbar(aes(xmin=coef-stderr, xmax=coef+stderr), width=.2) 

#age
m_age <- subset(maaslin, metadata == "age")

m_a <- ggplot(m_age, aes(x=coef, y=Species, fill=posneg)) +
  geom_bar(stat="identity") +
  ggtitle("Age") +
  scale_fill_manual(values=c('#800000FF', '#155F83FF')) +
  scale_x_continuous(breaks = seq(-6, 4, 1), limits = c(-5, 3)) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 12, hjust = 0.5, face = "bold"),
    # legend.text = element_text(size = 12),
    legend.position = "none",
    # legend.title = element_blank(),
    axis.text = element_text(size = 10),
    axis.text.y = element_text(face = "italic"),
    axis.title = element_blank()) +
  geom_errorbar(aes(xmin=coef-stderr, xmax=coef+stderr), width=.2) 

#recent URI
m_uri <- subset(maaslin, metadata == "uri_rec2")

m_u <- ggplot(m_uri, aes(x=coef, y=Species, fill=posneg)) +
  geom_bar(stat="identity") +
  ggtitle("Recent URI symptoms") +
  scale_fill_manual(values=c('#155F83FF')) +
  scale_x_continuous(breaks = seq(-6, 4, 1), limits = c(-5, 3)) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 12, hjust = 0.5, face = "bold"),
    # legend.text = element_text(size = 12),
    legend.position = "none",
    # legend.title = element_blank(),
    axis.text = element_text(size = 10),
    axis.text.y = element_text(face = "italic"),
    axis.title = element_blank()) +
  geom_errorbar(aes(xmin=coef-stderr, xmax=coef+stderr), width=.2) 

#hiv
m_hiv <- subset(maaslin, metadata == "hiv")

#Source to split C. pseudo into 2 lines on the Y axis for space: https://stackoverflow.com/questions/20123147/add-line-break-to-axis-labels-and-ticks-in-ggplot
m_h <- ggplot(m_hiv, aes(x=coef, y=Species, fill=posneg)) +
  geom_bar(stat="identity") +
  ggtitle("HIV infection") +
  scale_fill_manual(values=c('#800000FF')) +
  scale_x_continuous(breaks = seq(-6, 4, 1), limits = c(-5, 3)) +
  scale_y_discrete(labels=function(x){sub("\\s", "\n", x)}) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 12, hjust = 0.5, face = "bold"),
    # legend.text = element_text(size = 12),
    legend.position = "none",
    # legend.title = element_blank(),
    axis.text = element_text(size = 10),
    axis.text.y = element_text(face = "italic"),
    axis.title = element_blank()) +
  geom_errorbar(aes(xmin=coef-stderr, xmax=coef+stderr), width=.2) 

#combine into 1 cowplot
right <- plot_grid(m_h, m_u, m_a,
          labels = c('B', 'C', 'D'),
          nrow = 3, 
          label_size = 12,
          rel_heights = c(1, 1, 1.75),
          align = "v")

fig3 <- plot_grid(m_s, right,
                  labels = c('A', '', '', ''),
                  nrow = 1,
                  label_size = 14, 
                  rel_heights = c(1, 0.8))
                  # align = "v", axis = "b")

png(file="fig3.png",
    width = 12, height = 4, units = 'in', res = 600)
fig3
dev.off()

#what if you make it a standard 4 part figure?
fig_3 <- plot_grid(m_h, m_a, m_u, m_s, labels = "AUTO", nrow = 2, ncol = 2, align = "hv", axis = "b")
png(file="fig3.png",
    width = 12, height = 9, units = 'in', res = 600)
fig_3
dev.off()


#Univariable
fit_data2 = Maaslin2(
  input_data = otu, 
  input_metadata = child, 
  min_prevalence = 0.2,
  max_significance = 0.2,
  output = "maaslin_hiv_kaiju_filt50", 
  fixed_effects = c("hiv"),
  reference = c("hiv,Unexposed"))
#No significant assoc

#HEI/HEU vs HUU
fit_data2 = Maaslin2(
  input_data = otu, 
  input_metadata = child,
  min_abundance = 0.0005,
  min_prevalence = 0.2,
  max_significance = 0.2,
  output = "maaslin_mv_kaiju_hivcat", 
  fixed_effects = c("hivcat", "season", "age", "abx", "wood", "uri_rec2"),
  reference = c("hivcat,HUU"))
#Signif association with recent URI and season

############################################
#INDIVIDUAL SPECIES ABUNDANCES BY HIV STATUS
############################################
#We know what our top 10 most common species are; is there an association with HIV status in univariable analysis?
  #Note: using data from chi_pruned (NOT transformed); need to do with transformed?
#Corynebacterium accolens
cor_acc <- filter(melted_df, Species == "Corynebacterium accolens")
summary(cor_acc$Abundance)
hist(cor_acc$Abundance)
tapply(cor_acc$Abundance, cor_acc$hiv, summary) #HIV+ kids have lowest median rel abundance
kruskal.test(cor_acc$Abundance, cor_acc$hiv)  #P=0.02

#Corynebacterium propinquum
cor_p <- filter(melted_df, Species == "Corynebacterium propinquum")
summary(cor_p$Abundance)
hist(cor_p$Abundance)
tapply(cor_p$Abundance, cor_p$hiv, summary)
kruskal.test(cor_p$Abundance, cor_p$hiv)  #trend in abundance by HIV status (P=0.06)

#Corynebacterium pseudodiphtheriticum
cor_ps <- filter(melted_df, Species == "Corynebacterium pseudodiphtheriticum")
summary(cor_ps$Abundance)
hist(cor_ps$Abundance)
tapply(cor_ps$Abundance, cor_ps$hiv, summary) #HIV+ kids have lowest median rel abundance
kruskal.test(cor_ps$Abundance, cor_ps$hiv)  #P=0.02

#Dolosigranulum pigrum
dpig <- filter(melted_df, Species == "Dolosigranulum pigrum")
hist(dpig$Abundance)
summary(dpig$Abundance)
tapply(dpig$Abundance, dpig$hiv, summary)
kruskal.test(dpig$Abundance, dpig$hiv)  #P=0.87

#Haemophilus influenzae
hflu <- filter(melted_df, Species == "Haemophilus influenzae")
hist(hflu$Abundance)  #lots of zeros
summary(hflu$Abundance)
tapply(hflu$Abundance, hflu$hiv, summary) #limited by zeros (all medians = 0)
kruskal.test(hflu$Abundance, hflu$hiv)  #P=0.50

#Moraxella cattarhalis
mcat <- filter(melted_df, Species == "Moraxella catarrhalis")
hist(mcat$Abundance)  #lots of zeros
summary(mcat$Abundance)
tapply(mcat$Abundance, mcat$hiv, summary) 
kruskal.test(mcat$Abundance, mcat$hiv)  #P=0.44

#Moraxella lincolnii
mlin <- filter(melted_df, Species == "Moraxella lincolnii")
hist(mlin$Abundance)  #lots of zeros
summary(mlin$Abundance)
tapply(mlin$Abundance, mlin$hiv, summary) 
kruskal.test(mlin$Abundance, mlin$hiv)  #P=0.16

#Moraxella nonliquefaciens
mnon <- filter(melted_df, Species == "Moraxella nonliquefaciens")
hist(mnon$Abundance)  #lots of zeros but more evenly distrib than other Moraxellas
summary(mnon$Abundance)
tapply(mnon$Abundance, mnon$hiv, summary)
kruskal.test(mnon$Abundance, mnon$hiv)  #P=0.43

#Staphylococcus aureus
sa <- filter(melted_df, Species == "Staphylococcus aureus")
hist(sa$Abundance)  #lots of zeros
summary(sa$Abundance)
tapply(sa$Abundance, sa$hiv, summary) #limited by zeros (all medians = 0)
kruskal.test(sa$Abundance, sa$hiv)  #P=0.049

#Streptococcus pneumoniae
spneum <- filter(melted_df, Species == "Streptococcus pneumoniae")
hist(spneum$Abundance)  #lots of zeros
summary(spneum$Abundance)
tapply(spneum$Abundance, spneum$hiv, summary)
kruskal.test(spneum$Abundance, spneum$hiv)  #P=0.55

#clean up environment (but KEEP pathobiont df for aim 2 analysis)
remove(cor_acc, cor_p, cor_ps, dpig, mlin, mnon)

########################################
#AIM 2. LOGISTIC REGRESSION: PATHOBIONTS
########################################
#We are interested in understanding factors associated with the odds of colonization by our 4 pathobionts of interest:
  #Strep pneumo, S aureus, M cat, H flu
  #Logistic regression has a binary outcome --> create new binary vars in our org-specific datasets to use as our outcome

#BINARY VARIABLES AND INITIAL CHI SQUARED TESTS
spneum$present <- NA
spneum$present [spneum$Abundance == 0] <- 0
spneum$present [spneum$Abundance > 0] <- 1
table(spneum$present)
table(spneum$present, spneum$hiv)
prop.table(table(spneum$present, spneum$hiv), 2)
summary(table(spneum$present, spneum$hiv))  #Trend towards higher %HEU kids with S.pneumo col (96% vs 82% in HIV+ and 86% in HUU) but not signif; p=0.09

sa$present <- NA
sa$present [sa$Abundance == 0] <- 0
sa$present [sa$Abundance > 0] <- 1
table(sa$present)
table(sa$present, sa$hiv)
prop.table(table(sa$present, sa$hiv), 2)
summary(table(sa$present, sa$hiv))

mcat$present <- NA
mcat$present [mcat$Abundance == 0] <- 0
mcat$present [mcat$Abundance > 0] <- 1
table(mcat$present)
table(mcat$present, mcat$hiv)
prop.table(table(mcat$present, mcat$hiv), 2)
summary(table(mcat$present, mcat$hiv))

hflu$present <- NA
hflu$present [hflu$Abundance == 0] <- 0
hflu$present [hflu$Abundance > 0] <- 1
table(hflu$present)
table(hflu$present, hflu$hiv)
prop.table(table(hflu$present, hflu$hiv), 2)
summary(table(hflu$present, hflu$hiv))

#LOGIT MODELS
  #Reference: https://stats.oarc.ucla.edu/r/dae/logit-regression/

#STREP PNEUMO
#S.pneumo univariable: HIV
spneum$hiv <- as.factor(spneum$hiv)
spneum$hiv <- relevel(spneum$hiv, ref = "Unexposed")
summary(mylogit <- glm(present ~ hiv, data = spneum, family = "binomial"))
summary(mylogit)  #No association between HEU (P=0.11) or HEI (P=0.58) status and S. pneumo col
#AIC 104.93

#season: no association
summary(mylogit <- glm(present ~ season, data = spneum, family = "binomial"))
#Age: no association
summary(mylogit <- glm(present ~ age, data = spneum, family = "binomial"))
#abx: no association
summary(mylogit <- glm(present ~ abx, data = spneum, family = "binomial"))
#wood: no association
summary(mylogit <- glm(present ~ wood, data = spneum, family = "binomial"))
#recent URI: no association
summary(mylogit <- glm(present ~ uri_rec2, data = spneum, family = "binomial"))
#current URI (out of curiosity): no association
summary(mylogit <- glm(present ~ uri_cur2, data = spneum, family = "binomial"))

#S. pneumo multivariable: no associations with any covariables and odds of S.pneumo colonization
summary(mylogit <- glm(present ~ hiv + season + age + abx + wood + uri_rec2, data = spneum, family = "binomial"))

#STAPH AUREUS
#S. aureus univariable
sa$hiv <- as.factor(sa$hiv)
sa$hiv <- relevel(sa$hiv, ref = "Unexposed")
summary(mylogit <- glm(present ~ hiv, data = sa, family = "binomial"))
#Significant increase in odds of S. aureus colonization among HEU children (P=0.009)
exp(mylogit$coefficients) #OR (95% CI) for HEU kids: 3.26 (1.36, 8.23)
exp(confint(mylogit))

#season: no association
summary(mylogit <- glm(present ~ season, data = sa, family = "binomial"))
#Age: trend towards reduction in odds with age but not significant (p=0.09)
summary(mylogit <- glm(present ~ age, data = sa, family = "binomial"))
#abx: no association
summary(mylogit <- glm(present ~ abx, data = sa, family = "binomial"))
#wood: no association
summary(mylogit <- glm(present ~ wood, data = sa, family = "binomial"))
#recent URI: no association
summary(mylogit <- glm(present ~ uri_rec2, data = sa, family = "binomial"))
#current URI (out of curiosity): no association
summary(mylogit <- glm(present ~ uri_cur2, data = sa, family = "binomial"))

#S. aureus multivariable
summary(mylogit <- glm(present ~ hiv + season + age + abx + wood + uri_rec2, data = sa, family = "binomial"))
#Significant increase in odds of S.aureus colonization among HEU children
exp(mylogit$coefficients) #OR (95% CI) for HEU kids: 3.13 (1.28, 8.05)
exp(confint(mylogit))

#MORAXELLA CATARRHALIS
#HIV: no association
mcat$hiv <- as.factor(mcat$hiv)
mcat$hiv <- relevel(mcat$hiv, ref = "Unexposed")
summary(mylogit <- glm(present ~ hiv, data = mcat, family = "binomial"))
#season: significantly reduced odds of M.cat colonization in rainy season
summary(mylogit <- glm(present ~ season, data = mcat, family = "binomial"))
exp(mylogit$coefficients) #OR (95% CI) for rainy season: 0.27 (0.09, 0.66)
exp(confint(mylogit))

#Age: no association
summary(mylogit <- glm(present ~ age, data = mcat, family = "binomial"))
#abx: no association
summary(mylogit <- glm(present ~ abx, data = mcat, family = "binomial"))
#wood: no association
summary(mylogit <- glm(present ~ wood, data = mcat, family = "binomial"))
#recent URI: no association
summary(mylogit <- glm(present ~ uri_rec2, data = mcat, family = "binomial"))
#current URI (out of curiosity): no association
summary(mylogit <- glm(present ~ uri_cur2, data = mcat, family = "binomial"))

#M. catarrhalis multivariable
summary(mylogit <- glm(present ~ hiv + season + age + abx + wood + uri_rec2, data = mcat, family = "binomial"))
#Significantly reduced odds of M.cat colonization in rainy season
exp(mylogit$coefficients) #OR (95% CI) for rainy season: 0.24 (0.08, 0.62)
exp(confint(mylogit))

#HAEMOPHILUS INFLUENZAE
#HIV: no association
hflu$hiv <- as.factor(hflu$hiv)
hflu$hiv <- relevel(hflu$hiv, ref = "Unexposed")
summary(mylogit <- glm(present ~ hiv, data = hflu, family = "binomial"))

#season: no association
summary(mylogit <- glm(present ~ season, data = hflu, family = "binomial"))
#Age: no association
summary(mylogit <- glm(present ~ age, data = hflu, family = "binomial"))
#abx: no association
summary(mylogit <- glm(present ~ abx, data = hflu, family = "binomial"))
#wood: no association
summary(mylogit <- glm(present ~ wood, data = hflu, family = "binomial"))
#recent URI: no association
summary(mylogit <- glm(present ~ uri_rec2, data = hflu, family = "binomial"))
#current URI (out of curiosity): no association
summary(mylogit <- glm(present ~ uri_cur2, data = hflu, family = "binomial"))

#H. flu multivariable
summary(mylogit <- glm(present ~ hiv + season + age + abx + wood + uri_rec2, data = hflu, family = "binomial"))
#No association between HIV status, season, age, abx exposure, wood smoke, or recent URI and odds of H.flu col

###############################################
#ADDITIONAL ANALYSES: INTERSPECIES CORRELATIONS
###############################################
#Need to create new dataframe that contains the relative abundances from our species-specific dataframes above
  #Then do correlations between the abundances
#Rename abundance variables in each dataframe
cor_ps <- rename(cor_ps, ab.cor_ps = Abundance)
cor_p <- rename(cor_p, ab.cor_p = Abundance)
dpig <- rename(dpig, ab.dpig = Abundance)
hflu <- rename(hflu, ab.hflu = Abundance)
mcat <- rename(mcat, ab.mcat = Abundance)
spneum <- rename(spneum, ab.spneum = Abundance)
sa <- rename(sa, ab.sa = Abundance)

#Create simple dataframe with abundance variables, sample_id, hiv, age, and merge on sample_id
corr <- cor_ps[c("sample_id", "age", "hiv", "pcv", "dpt", "abx", "uri_rec2", "wood", "season", "ab.cor_ps")]
test <- cor_p[c("sample_id", "ab.cor_p")]
corr <- merge(corr, test, by ="sample_id", sort = TRUE)
test <- hflu[c("sample_id", "ab.hflu")]
corr <- merge(corr, test, by ="sample_id", sort = TRUE)
test <- mcat[c("sample_id", "ab.mcat")]
corr <- merge(corr, test, by ="sample_id", sort = TRUE)
test <- spneum[c("sample_id", "ab.spneum")]
corr <- merge(corr, test, by ="sample_id", sort = TRUE)
test <- sa[c("sample_id", "ab.sa")]
corr <- merge(corr, test, by ="sample_id", sort = TRUE)
test <- dpig[c("sample_id", "ab.dpig")]
corr <- merge(corr, test, by ="sample_id", sort = TRUE)

#Per 6/9 MK meeting: create binary variables for presence/absence spneumo and s.aureus for logistic regression models
corr$sp_present <- NA
corr$sp_present [corr$ab.spneum == 0] <- 0
corr$sp_present [corr$ab.spneum > 0] <- 1
table(corr$sp_present)
#matches our Aim 2 numbers

corr$sa_present <- NA
corr$sa_present [corr$ab.sa == 0] <- 0
corr$sa_present [corr$ab.sa > 0] <- 1
table(corr$sa_present)
#matches our Aim 2 numbers

#Per 7/14 MK meeting, will log transform SP and SA rel abundances for linear regression model
  #does this mean we also need to transform our cor_ps and dpig abundances? 
#UPDATE 8-30: skipping per MK edits of manuscript (correlations only)
# corr$sp_log <- log(corr$ab.spneum)
# hist(corr$sp_log)
# 
# corr$sa_log <- log(corr$ab.sa)
# hist(corr$sa_log)

#save as csv in case we want to come back to it
# write.csv(corr, "corrtable_all.csv", row.names=F)
write.csv(corr, "corrtable_updatedAug22.csv", row.names = F)
corr <- read.csv("corrtable_updatedAug22.csv")

#Correlations
  #source: https://stats.stackexchange.com/questions/8071/how-to-choose-between-pearson-and-spearman-correlation
#based on the above, will run both tests to look at how they compare (can help us see if the relationship is monotonic and/or linear)

###CORYNE PSEUDO VS STREP PNEUMO
cor.test(corr$ab.cor_ps, corr$ab.spneum, method = "pearson")  #P=0.026, cor = -0.186
cor.test(corr$ab.cor_ps, corr$ab.spneum, method = "spearman")  #P=0.009, estimate = -0.217 BUT warning message
  #warning:Cannot compute exact p-value with ties
#need to plot this out and see what it looks like (overall and then by HIV status)
  #source: https://bookdown.org/dli/rguide/scatterplots-and-best-fit-lines-two-sets.html#two-scatterplots-using-ggplot2
#what to do about the ties warning message? 
  #Source: https://stackoverflow.com/questions/10711395/spearman-correlation-and-ties
#try with (exact = FALSE) and see what happens 
cor.test(corr$ab.cor_ps, corr$ab.spneum, method = "spearman", exact = FALSE)
  #get the same result; per link above R calculates the Spearman coefficient corrected for ties

#Correlations by HIV status:
corr_hei <- subset(corr, hiv == "Infected")
corr_heu <- subset(corr, hiv == "Exposed Uninfected")
corr_huu <- subset(corr, hiv == "Unexposed")

cor.test(corr_hei$ab.cor_ps, corr_hei$ab.spneum, method = "pearson")  #P=0.34, cor = -0.148
cor.test(corr_hei$ab.cor_ps, corr_hei$ab.spneum, method = "spearman") #P=0.34, cor = -0.145 with warning re: ties

cor.test(corr_heu$ab.cor_ps, corr_heu$ab.spneum, method = "pearson")  #P=0.15, cor = -0.209
cor.test(corr_heu$ab.cor_ps, corr_heu$ab.spneum, method = "spearman") #P=0.002, cor = -0.426 with warning

cor.test(corr_huu$ab.cor_ps, corr_huu$ab.spneum, method = "pearson")  #P=0.16, cor = -0.201
cor.test(corr_huu$ab.cor_ps, corr_huu$ab.spneum, method = "spearman") #P=0.39, cor = -0.124 with warning re: ties

#Plot out
#UPDATED 8-19-22: plot correlations in same figure using cowplot
#abbreviate HIV status and relevel
corr$hiv2 <- NA
corr$hiv2 [corr$hiv == "Infected"] <- "HEI"
corr$hiv2 [corr$hiv == "Exposed Uninfected"] <- "HEU"
corr$hiv2 [corr$hiv == "Unexposed"] <- "HUU"
table(corr$hiv2)
#matches up with known numbers
corr$hiv2 <- as.factor(corr$hiv2)
corr$hiv2 <- reorder(corr$hiv2, new.order=c("HUU", "HEU", "HEI"))

#C. pseudodiphtheriticum vs S. pneumo plot: all participants together
#Source for adding text: http://www.sthda.com/english/wiki/ggplot2-texts-add-text-annotations-to-a-graph-in-r-software
library(ggtext)
# png(file="/Users/swetapatel/OneDrive - Duke University/Fogarty coding/Coryne_SP_all.png",
#     width = 5, height = 5, units = 'in', res = 600)
csp <- ggplot(data = corr, aes(x = ab.cor_ps, y = ab.spneum)) + 
  geom_point() + 
  geom_smooth(method = "lm", color="red", size=0.5, se = FALSE) +
  labs(title="Interspecies correlations",x="*C. pseudodiphtheriticum* relative abundance", y = "*S. pneumoniae* relative abundance")+
  geom_text(x=0.4, y=0.55, label="r = - 0.22, p = 0.009") +
  theme_bw() +
  theme(
    # plot.title = element_text(size = 16, hjust = 0.5),
    plot.title = element_blank(),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.position = "right",
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    axis.title.x = ggtext::element_markdown(),
    axis.title.y = ggtext::element_markdown())
# dev.off()

csp

#C. pseudodiphtheriticum vs S. pneumo plot: all participants by HIV status
ggplot(data = corr, aes(x = ab.cor_ps, y = ab.spneum, color = hiv2)) + 
  geom_point(aes(shape=hiv2)) + 
  scale_color_manual(values=c('#800000FF','#155F83FF', '#FFA319FF'))+
  labs(title="Interspecies Correlations",x="*C. pseudodiphtheriticum*", y = "*S. pneumoniae*")+
  theme_bw() +
  theme(
    plot.title = element_text(size = 16, hjust = 0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    legend.position = c(0.83,0.86),
    legend.background = element_blank(),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    axis.title.x = ggtext::element_markdown(),
    axis.title.y = ggtext::element_markdown())

#What about 3 plots, one for each HIV category? (Facet wrap)
png(file="/Users/swetapatel/OneDrive - Duke University/Fogarty coding/Coryne_SP_facet.png",
    width = 5, height = 5, units = 'in', res = 600)
corr$hiv <- as.factor(corr$hiv)
corr$hiv <- reorder(corr$hiv, new.order=c("Unexposed", "Exposed Uninfected", "Infected"))
ggplot(data = corr, aes(x = ab.cor_ps, y = ab.spneum)) + 
  geom_point() + 
  geom_smooth(method = "lm", color="red", size=0.5, se = FALSE) +
  facet_grid(~hiv) +
  labs(title="Interspecies correlations",x="*C. pseudodiphtheriticum*", y = "*S. pneumoniae*")+
  theme_bw() +
  theme(
    plot.title = element_text(size = 16, hjust = 0.5),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.position = "right",
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    axis.title.x = ggtext::element_markdown(),
    axis.title.y = ggtext::element_markdown())
dev.off()
#manually annotate for now after export 

#Logistic regression models with presence/absence S. pneumo as our outcome variable
  #Will look at Coryne pseudo abundance, age, recent URI, abx, +/- HIV (a priori expect these things to affect S.pneumo col)
#UPDATE 8-19-22: analyses for manuscript need to be uniform; adjust for all vars (age, HIV, URI, abx, wood, season)
  #NOTE: as part of Aim 2 above (lines 1698-1746), we've already looked at these variables EXCEPT Coryne; none associated
summary(mylogit <- glm(sp_present ~ ab.cor_ps, data = corr, family = "binomial")) #No association; P=0.871
  #age-adjusted with Coryne
summary(mylogit <- glm(sp_present ~ ab.cor_ps + age, data = corr, family = "binomial")) #No assoc; p=0.85 for Coryne

#S. pneumo multivariable: no associations with any covariables and odds of S.pneumo colonization
summary(mylogit <- glm(sp_present ~ ab.cor_ps + age + abx + uri_rec2, data = corr, family = "binomial"))
  #Multivariable including HIV: none associated
summary(mylogit <- glm(sp_present ~ ab.cor_ps + age + abx + uri_rec2 + hiv, data = corr, family = "binomial"))
  #Multivariable including all 6 vars: none associated
#relevel HIV variable first
class(corr$hiv) #character; need to make factor to relevel
corr$hiv <- as.factor(corr$hiv)
corr$hiv <- relevel(corr$hiv, ref = "Unexposed")
summary(mylogit <- glm(sp_present ~ ab.cor_ps + age + abx + uri_rec2 + hiv + wood + season, data = corr, family = "binomial"))
exp(mylogit$coefficients) #OR (95% CI) for C. pseudo rel abundance: 1.44 (0.05, 86.96)
exp(confint(mylogit))

#linear regression model with log txformed SP abundance as our outcome var
summary(my_mv <- lm(sp_log ~ age + hiv + wood + season + abx + uri_rec2, data = corr))
  #error message: can't run model because we have values of sp_log that are "-Inf"
  #Skip for now, return pending MK discussion

###CORYNE PSEUDO VS STAPH AUREUS
cor.test(corr$ab.cor_ps, corr$ab.sa, method = "pearson")  #P=0.08, cor = -0.148
cor.test(corr$ab.cor_ps, corr$ab.sa, method = "spearman")  #P=0.0001, estimate = -0.317 BUT warning message
#warning:Cannot compute exact p-value with ties
cor.test(corr$ab.cor_ps, corr$ab.sa, method = "spearman", exact = FALSE) #same as above but without the warning

#Correlations by HIV status:
cor.test(corr_hei$ab.cor_ps, corr_hei$ab.sa, method = "pearson")  #P=0.37, cor = -0.139
cor.test(corr_hei$ab.cor_ps, corr_hei$ab.sa, method = "spearman") #P=0.005, cor = -0.413 with warning re: ties

cor.test(corr_heu$ab.cor_ps, corr_heu$ab.sa, method = "pearson")  #P=0.21, cor = -0.184
cor.test(corr_heu$ab.cor_ps, corr_heu$ab.sa, method = "spearman") #P=0.04, cor = -0.302 with warning

cor.test(corr_huu$ab.cor_ps, corr_huu$ab.sa, method = "pearson")  #P=0.26, cor = -0.161
cor.test(corr_huu$ab.cor_ps, corr_huu$ab.sa, method = "spearman") #P=0.03, cor = -0.308 with warning re: ties

#Plot out
#C. pseudodiphtheriticum vs S. aureus plot: all participants together
#Source for adding text: http://www.sthda.com/english/wiki/ggplot2-texts-add-text-annotations-to-a-graph-in-r-software
# png(file="/Users/swetapatel/OneDrive - Duke University/Fogarty coding/Coryne_SA_all.png",
#     width = 5, height = 5, units = 'in', res = 600)
csa <- ggplot(data = corr, aes(x = ab.cor_ps, y = ab.sa)) + 
  geom_point() + 
  geom_smooth(method = "lm", color="red", size=0.5, se = FALSE) +
  labs(title="Interspecies correlations",x="*C. pseudodiphtheriticum* relative abundance", y = "*S. aureus* relative abundance")+
  geom_text(x=0.4, y=0.55, label="r = - 0.32, p = 0.0001") +
  scale_y_continuous(breaks = seq(0, 0.8, 0.2), limits = c(0, 0.8)) +
  theme_bw() +
  theme(
    # plot.title = element_text(size = 16, hjust = 0.5),
    plot.title = element_blank(),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.position = "right",
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    axis.title.x = ggtext::element_markdown(),
    axis.title.y = ggtext::element_markdown())
# dev.off()

#C. pseudodiphtheriticum vs S. aureus plot: all participants by HIV status
ggplot(data = corr, aes(x = ab.cor_ps, y = ab.sa, color = hiv2)) + 
  geom_point(aes(shape=hiv2)) + 
  scale_color_manual(values=c('#800000FF','#155F83FF', '#FFA319FF'))+
  labs(title="Interspecies Correlations",x="*C. pseudodiphtheriticum*", y = "*S. aureus*")+
  theme_bw() +
  theme(
    plot.title = element_text(size = 16, hjust = 0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    legend.position = c(0.83,0.86),
    legend.background = element_blank(),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    axis.title.x = ggtext::element_markdown(),
    axis.title.y = ggtext::element_markdown())

#3 plots, one for each HIV category (Facet wrap)
png(file="/Users/swetapatel/OneDrive - Duke University/Fogarty coding/Coryne_SA_facet.png",
    width = 5, height = 5, units = 'in', res = 600)
corr$hiv <- as.factor(corr$hiv)
corr$hiv <- reorder(corr$hiv, new.order=c("Unexposed", "Exposed Uninfected", "Infected"))
ggplot(data = corr, aes(x = ab.cor_ps, y = ab.sa)) + 
  geom_point() + 
  geom_smooth(method = "lm", color="red", size=0.5, se = FALSE) +
  facet_grid(~hiv) +
  labs(title="Interspecies correlations",x="*C. pseudodiphtheriticum*", y = "*S. aureus*")+
  theme_bw() +
  theme(
    plot.title = element_text(size = 16, hjust = 0.5),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.position = "right",
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    axis.title.x = ggtext::element_markdown(),
    axis.title.y = ggtext::element_markdown())
dev.off()
#manually annotate for now after export 

#Logistic regression models with presence/absence S. aureus as our outcome variable
#Will look at Coryne pseudo abundance, age, recent URI, abx, +/- HIV (a priori expect these things to affect S.pneumo col)
#NOTE: as part of Aim 2 above (lines 1698-1746), we've already looked at these variables EXCEPT Coryne
  #Only HIV status associated
summary(mylogit <- glm(sa_present ~ ab.cor_ps, data = corr, family = "binomial")) #Associated!! P = 0.03; AIC 176.33
#age-adjusted with Coryne
summary(mylogit <- glm(sa_present ~ ab.cor_ps + age, data = corr, family = "binomial")) #Still assoc; P = 0.03
exp(mylogit$coefficients) #OR (95% CI) for C. pseudo rel abundance: 0.04 (0.001, 0.60)
exp(confint(mylogit))

#S. aureus multivariable: 
summary(mylogit <- glm(sa_present ~ ab.cor_ps + age + abx + uri_rec2, data = corr, family = "binomial"))
exp(mylogit$coefficients) #OR (95% CI) for C. pseudo rel abundance: 0.04 (0.001, 0.69)
exp(confint(mylogit))
  #multivariable including HIV (which was prev assoc)
summary(mylogit <- glm(sa_present ~ ab.cor_ps + age + abx + uri_rec2 + hiv, data = corr, family = "binomial"))
exp(mylogit$coefficients) #OR (95% CI) for C. pseudo rel abundance: 0.03 (0.001, 0.68)
exp(confint(mylogit))
#OR (95% CI) for HEU = 3.13 (1.27, 8.17); not far off from our Aim 2 results
  #multivariable including all 6 covars
summary(mylogit <- glm(sa_present ~ ab.cor_ps + age + abx + uri_rec2 + hiv + wood + season, data = corr, family = "binomial"))
exp(mylogit$coefficients) #OR (95% CI) for C. pseudo rel abundance: 0.04 (0.001, 0.70)
exp(confint(mylogit))
  #OR (95% CI) for HEU: 3.10 (1.25, 8.12)

###CORYNE PSEUDO VS MORAXELLA CATARRHALIS
cor.test(corr$ab.cor_ps, corr$ab.mcat, method = "pearson")  #P=0.03, cor = -0.186
cor.test(corr$ab.cor_ps, corr$ab.mcat, method = "spearman", exact = FALSE)  #P=0.48, estimate = -0.059 BUT warning message
#warning:Cannot compute exact p-value with ties

#Correlations by HIV status:
cor.test(corr_hei$ab.cor_ps, corr_hei$ab.mcat, method = "pearson")  #P=0.26, cor = -0.172
cor.test(corr_hei$ab.cor_ps, corr_hei$ab.mcat, method = "spearman") #P=0.33, cor = -0.151 with warning re: ties

cor.test(corr_heu$ab.cor_ps, corr_heu$ab.mcat, method = "pearson")  #P=0.17, cor = -0.201
cor.test(corr_heu$ab.cor_ps, corr_heu$ab.mcat, method = "spearman") #P=0.47, cor = 0.105 with warning

cor.test(corr_huu$ab.cor_ps, corr_huu$ab.mcat, method = "pearson")  #P=0.13, cor = -0.217
cor.test(corr_huu$ab.cor_ps, corr_huu$ab.mcat, method = "spearman") #P=0.36, cor = -0.132 with warning re: ties

#Plot out
#C. pseudodiphtheriticum vs M. catarrhalis plot: all participants together
#Source for adding text: http://www.sthda.com/english/wiki/ggplot2-texts-add-text-annotations-to-a-graph-in-r-software
png(file="/Users/swetapatel/OneDrive - Duke University/Fogarty coding/Coryne_MC_all.png",
    width = 5, height = 5, units = 'in', res = 600)
ggplot(data = corr, aes(x = ab.cor_ps, y = ab.mcat)) + 
  geom_point() + 
  geom_smooth(method = "lm", color="red", size=0.5, se = FALSE) +
  labs(title="Interspecies correlations",x="*C. pseudodiphtheriticum*", y = "*M. cattarhalis*")+
  geom_text(x=0.4, y=0.55, label="r = - 0.059, P = 0.48") +
  theme_bw() +
  theme(
    plot.title = element_text(size = 16, hjust = 0.5),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.position = "right",
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    axis.title.x = ggtext::element_markdown(),
    axis.title.y = ggtext::element_markdown())
dev.off()

#C. pseudodiphtheriticum vs M. cattarhalis plot: all participants by HIV status
ggplot(data = corr, aes(x = ab.cor_ps, y = ab.mcat, color = hiv2)) + 
  geom_point(aes(shape=hiv2)) + 
  scale_color_manual(values=c('#800000FF','#155F83FF', '#FFA319FF'))+
  labs(title="Interspecies Correlations",x="*C. pseudodiphtheriticum*", y = "*M. cattarhalis*")+
  theme_bw() +
  theme(
    plot.title = element_text(size = 16, hjust = 0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    legend.position = c(0.83,0.86),
    legend.background = element_blank(),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    axis.title.x = ggtext::element_markdown(),
    axis.title.y = ggtext::element_markdown())

#3 plots, one for each HIV category (Facet wrap)
png(file="/Users/swetapatel/OneDrive - Duke University/Fogarty coding/Coryne_MC_facet.png",
    width = 5, height = 5, units = 'in', res = 600)
corr$hiv <- as.factor(corr$hiv)
corr$hiv <- reorder(corr$hiv, new.order=c("Unexposed", "Exposed Uninfected", "Infected"))
ggplot(data = corr, aes(x = ab.cor_ps, y = ab.mcat)) + 
  geom_point() + 
  geom_smooth(method = "lm", color="red", size=0.5, se = FALSE) +
  facet_grid(~hiv) +
  labs(title="Interspecies correlations",x="*C. pseudodiphtheriticum*", y = "*M. cattarhalis*")+
  theme_bw() +
  theme(
    plot.title = element_text(size = 16, hjust = 0.5),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.position = "right",
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    axis.title.x = ggtext::element_markdown(),
    axis.title.y = ggtext::element_markdown())
dev.off()
#manually annotate after export 

###CORYNE PSEUDO VS HAEMOPHILUS INFLUENZAE
cor.test(corr$ab.cor_ps, corr$ab.hflu, method = "pearson")  #P=0.04, cor = -0.169
cor.test(corr$ab.cor_ps, corr$ab.hflu, method = "spearman", exact = FALSE)  #P=0.301, estimate = -0.087

#Correlations by HIV status:
cor.test(corr_hei$ab.cor_ps, corr_hei$ab.hflu, method = "pearson")  #P=0.21, cor = -0.192
cor.test(corr_hei$ab.cor_ps, corr_hei$ab.hflu, method = "spearman", exact = FALSE) #P=0.09, cor = -0.260

cor.test(corr_heu$ab.cor_ps, corr_heu$ab.hflu, method = "pearson")  #P=0.22, cor = -0.179
cor.test(corr_heu$ab.cor_ps, corr_heu$ab.hflu, method = "spearman", exact = FALSE) #P=0.36, cor = -0.133

cor.test(corr_huu$ab.cor_ps, corr_huu$ab.hflu, method = "pearson")  #P=0.37, cor = -0.131
cor.test(corr_huu$ab.cor_ps, corr_huu$ab.hflu, method = "spearman", exact = FALSE) #P=0.47, cor = 0.104

#Plot out
#C. pseudodiphtheriticum vs H. influenzae plot: all participants together
#Source for adding text: http://www.sthda.com/english/wiki/ggplot2-texts-add-text-annotations-to-a-graph-in-r-software
png(file="/Users/swetapatel/OneDrive - Duke University/Fogarty coding/Coryne_HI_all.png",
    width = 5, height = 5, units = 'in', res = 600)
ggplot(data = corr, aes(x = ab.cor_ps, y = ab.hflu)) + 
  geom_point() + 
  geom_smooth(method = "lm", color="red", size=0.5, se = FALSE) +
  labs(title="Interspecies correlations",x="*C. pseudodiphtheriticum*", y = "*H. influenzae*")+
  geom_text(x=0.4, y=0.55, label="r = - 0.087, P = 0.301") +
  theme_bw() +
  theme(
    plot.title = element_text(size = 16, hjust = 0.5),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.position = "right",
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    axis.title.x = ggtext::element_markdown(),
    axis.title.y = ggtext::element_markdown())
dev.off()

#C. pseudodiphtheriticum vs M. cattarhalis plot: all participants by HIV status
ggplot(data = corr, aes(x = ab.cor_ps, y = ab.hflu, color = hiv2)) + 
  geom_point(aes(shape=hiv2)) + 
  scale_color_manual(values=c('#800000FF','#155F83FF', '#FFA319FF'))+
  labs(title="Interspecies Correlations",x="*C. pseudodiphtheriticum*", y = "*H. influenzae*")+
  theme_bw() +
  theme(
    plot.title = element_text(size = 16, hjust = 0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    legend.position = c(0.83,0.86),
    legend.background = element_blank(),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    axis.title.x = ggtext::element_markdown(),
    axis.title.y = ggtext::element_markdown())

#3 plots, one for each HIV category (Facet wrap)
png(file="/Users/swetapatel/OneDrive - Duke University/Fogarty coding/Coryne_HI_facet.png",
    width = 5, height = 5, units = 'in', res = 600)
corr$hiv <- as.factor(corr$hiv)
corr$hiv <- reorder(corr$hiv, new.order=c("Unexposed", "Exposed Uninfected", "Infected"))
ggplot(data = corr, aes(x = ab.cor_ps, y = ab.hflu)) + 
  geom_point() + 
  geom_smooth(method = "lm", color="red", size=0.5, se = FALSE) +
  facet_grid(~hiv) +
  labs(title="Interspecies correlations",x="*C. pseudodiphtheriticum*", y = "*H. influenzae*")+
  theme_bw() +
  theme(
    plot.title = element_text(size = 16, hjust = 0.5),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.position = "right",
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    axis.title.x = ggtext::element_markdown(),
    axis.title.y = ggtext::element_markdown())
dev.off()
#manually annotate after export 

#S.pneumo vs S.aureus (known antagonists) and see what we get
cor.test(corr$ab.sa, corr$ab.spneum, method = "pearson")  #P=0.396, estimate = -0.072
cor.test(corr$ab.sa, corr$ab.spneum, method = "spearman")  #P=0.0005, estimate = 0.286 BUT warning message about ties

###D. PIGRUM tests (not stratified by HIV); added 8-19-22
  #run all 4 spearman corrs and then plot out whatever is significant
cor.test(corr$ab.dpig, corr$ab.spneum, method = "spearman", exact = FALSE) #P<0.0001, estimate = -0.35
cor.test(corr$ab.dpig, corr$ab.sa, method = "spearman", exact = FALSE)    #P=0.0003, estimate = -0.30
cor.test(corr$ab.dpig, corr$ab.mcat, method = "spearman", exact = FALSE)  #P=0.10, estimate = -0.14
cor.test(corr$ab.dpig, corr$ab.hflu, method = "spearman", exact = FALSE)  #P=0.16, estimate = -0.12

#S pneumo and S aureus significant, so will plot these out as fig 4C and 4D for manuscript
#S. pneumoniae
dsp <- ggplot(data = corr, aes(x = ab.dpig, y = ab.spneum)) + 
  geom_point() + 
  geom_smooth(method = "lm", color="red", size=0.5, se = FALSE) +
  labs(title="Interspecies correlations",x="*D. pigrum* relative abundance", y = "*S. pneumoniae* relative abundance")+
  geom_text(x=0.4, y=0.55, label="r = - 0.35, p < 0.0001") +
  theme_bw() +
  theme(
    plot.title = element_blank(),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.position = "right",
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    axis.title.x = ggtext::element_markdown(),
    axis.title.y = ggtext::element_markdown())

#S. aureus
dsa <- ggplot(data = corr, aes(x = ab.dpig, y = ab.sa)) + 
  geom_point() + 
  geom_smooth(method = "lm", color="red", size=0.5, se = FALSE) +
  labs(title="Interspecies correlations",x="*D. pigrum* relative abundance", y = "*S. aureus* relative abundance")+
  geom_text(x=0.4, y=0.55, label="r = - 0.30, p = 0.0003") +
  scale_y_continuous(breaks = seq(0, 0.8, 0.2), limits = c(0, 0.8)) +
  theme_bw() +
  theme(
    plot.title = element_blank(),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.position = "right",
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    axis.title.x = ggtext::element_markdown(),
    axis.title.y = ggtext::element_markdown())

###C. PROPINQUUM tests (not stratified by HIV); added 8-31-22
#run all 4 spearman corrs and then plot out whatever is significant
cor.test(corr$ab.cor_p, corr$ab.spneum, method = "spearman", exact = FALSE) #P=0.0005, estimate = -0.29
cor.test(corr$ab.cor_p, corr$ab.sa, method = "spearman", exact = FALSE)    #P=0.0002, estimate = -0.31
cor.test(corr$ab.cor_p, corr$ab.mcat, method = "spearman", exact = FALSE)  #P=0.002, estimate = -0.25
cor.test(corr$ab.cor_p, corr$ab.hflu, method = "spearman", exact = FALSE)  #P=0.25, estimate = -0.10

#S pneumo, S aureus, AND MCAT are significant --> plot all SA and SP out as figures 4E and 4F for manuscript?
#S. pneumoniae
cpsp <- ggplot(data = corr, aes(x = ab.cor_p, y = ab.spneum)) + 
  geom_point() + 
  geom_smooth(method = "lm", color="red", size=0.5, se = FALSE) +
  labs(title="Interspecies correlations",x="*C. propinquum* relative abundance", y = "*S. pneumoniae* relative abundance")+
  geom_text(x=0.4, y=0.55, label="r = - 0.29, p = 0.0005") +
  theme_bw() +
  theme(
    plot.title = element_blank(),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.position = "right",
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    axis.title.x = ggtext::element_markdown(),
    axis.title.y = ggtext::element_markdown())

#S. aureus
cpsa <- ggplot(data = corr, aes(x = ab.cor_p, y = ab.sa)) + 
  geom_point() + 
  geom_smooth(method = "lm", color="red", size=0.5, se = FALSE) +
  labs(title="Interspecies correlations",x="*C. propinquum* relative abundance", y = "*S. aureus* relative abundance")+
  geom_text(x=0.4, y=0.55, label="r = - 0.31, p = 0.0002") +
  scale_y_continuous(breaks = seq(0, 0.8, 0.2), limits = c(0, 0.8)) +
  theme_bw() +
  theme(
    plot.title = element_blank(),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.position = "right",
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    axis.title.x = ggtext::element_markdown(),
    axis.title.y = ggtext::element_markdown())

#now put all 4 figs together to make Fig 4 using cowplot
fig4 <- plot_grid(csp, csa, dsp, dsa, labels = "AUTO", nrow = 2, ncol = 2, align = "h", axis = "b")
#UPDATED WITH C. PROPINQUUM (6 part figure)
fig4 <- plot_grid(cpsp, cpsa, csp, csa, dsp, dsa, labels = "AUTO", nrow = 3, ncol = 2, align = "h", axis = "b")

png(file="fig4.png",
    width = 12, height = 12, units = 'in', res = 600)
fig4
dev.off()

#D. pig logistic regression models: S. pneumo
  #univariable
summary(mylogit <- glm(sp_present ~ ab.dpig, data = corr, family = "binomial")) #No assoc (p=0.864)
  #multivariable including all 6 covars
summary(mylogit <- glm(sp_present ~ ab.dpig + age + abx + uri_rec2 + hiv + wood + season, data = corr, family = "binomial"))
exp(mylogit$coefficients) #OR (95% CI) for D.pig rel abundance: 1.12 (0.04, 36.27)
exp(confint(mylogit))

#D. pig logistic regression models: S. aureus
summary(mylogit <- glm(sa_present ~ ab.dpig, data = corr, family = "binomial")) #Assoc! P=0.003
#multivariable including all 6 covars
summary(mylogit <- glm(sa_present ~ ab.dpig + age + abx + uri_rec2 + hiv + wood + season, data = corr, family = "binomial"))
exp(mylogit$coefficients) #OR (95% CI) for D.pig rel abundance: 0.019 (0.001, 0.24)
exp(confint(mylogit))




##########################################################################
#ADDITIONAL ANALYSES: FACTORS A/W CORYNE PSEUDODIPTHERITICUM REL ABUNDANCE
##########################################################################
hist(cor_ps$ab.cor_ps)
#sex
tapply(cor_ps$ab.cor_ps, cor_ps$sex, summary) 
wilcox.test(ab.cor_ps ~ sex, data = cor_ps) #P=0.91

#age: P=0.82
cor.test(cor_ps$ab.cor_ps, cor_ps$age, method = "spearman", exact = FALSE)

#wood smoke
tapply(cor_ps$ab.cor_ps, cor_ps$wood, summary)
wilcox.test(ab.cor_ps ~ wood, data = cor_ps)  #P=0.14

#season
tapply(cor_ps$ab.cor_ps, cor_ps$season, summary) 
wilcox.test(ab.cor_ps ~ season, data = cor_ps)  #P=0.98

#abx
tapply(cor_ps$ab.cor_ps, cor_ps$abx, summary) 
wilcox.test(ab.cor_ps ~ abx, data = cor_ps) #Trend towards higher rel abund in absence of abx exposure (P=0.08)

#current URI
tapply(cor_ps$ab.cor_ps, cor_ps$uri_cur2, summary) 
wilcox.test(ab.cor_ps ~ uri_cur2, data = cor_ps)  #Current URI a/w lower rel abund (P=0.04)

#recent URI
tapply(cor_ps$ab.cor_ps, cor_ps$uri_rec2, summary) 
wilcox.test(ab.cor_ps ~ uri_rec2, data = cor_ps)  #Recent URI a/w lower rel abund (P=0.03)

#Maternal age?? (strong a/w HIV status): P=0.98
cor.test(cor_ps$ab.cor_ps, cor_ps$mat_age2, method = "spearman", exact = FALSE)