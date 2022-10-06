#Impact of HIV infection or exposure on the nasopharyngeal microbiome of children in Botswana
#Sub-analysis with Kaiju output
#Compiled May 26, 2022 by Sweta Patel

#***GENERAL ORDER OF STEPS***
#Within HIV+ kids: compare virally suppressed to not, immune suppressed to not
#1. Alpha diversity analyses
#2. Transform data (CLR) and do Beta-diversity analyses (PCoA plots, PERMANOVA)
#3. Generate relative abundance plots using NON-transformed data
    #are rel abundances of 10 most common species a/w specific characteristics of the kids with HIV? 
#4. Filter on mean relative abundance of 50 and remove NOS species --> ddx abundance analyses using Maaslin and FILTERED data
    #set Q value to 0.1-0.2 and control for age --> REMOVED JULY 2022
#5. UPDATE JULY 2022: compare immune competent HEI to HEU, HUU (composition, KW on 10 most abundant species)
#6. UPDATE AUG 2022: compare immune compromised HEI to HEU, HUU (composition, KW on 10 most abundant species)
#7. UPDATE SEPT 2022: see if possible to run negative binomial models to see if TMP-SMX, abx, vl, cd4 are associated with rel abundance of coryne, d pigrum, s. aureus

#For HIV-HEU sibling pairs:
#1. Compare overall composition using paired permanova (pairwise adonis)
#2. Compare rel abundance of 10 most abundant species using paired Wilcoxon

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

set.seed(1234)

setwd("OneDrive - Duke University/Fogarty coding/Sequencing/")

pruned <- readRDS("phy.pruned05092022.rds")
#confirm the taxtable, etc looks good:
test <- as.data.frame(as(tax_table(pruned),"matrix"),stringsAsFactors=FALSE)  #success!
meta_prune <- data.frame(sample_data(pruned))
remove(test)

#Subsetting out HIV+ kids
chi_prune = subset_samples(pruned, subject == "child")
hiv = subset_samples(chi_prune, hiv == "Infected")
chi_hiv <- data.frame(sample_data(hiv))

#Subsetting out sibling-HIV+ pairs (we have 10): B-09, B-16, B-17, B-20, B-23, B-25, B-32, B-35, B-47, B-48
  #Note: missing child data for B-47, so will exclude that sibling from our analysis
sib_pair = subset_samples(pruned, sample_id == "B.09.CHI" | sample_id == "B.09.SIB" | sample_id == "B.16.CHI" | 
                            sample_id == "B.16.SIB" | sample_id == "B.17.CHI" | sample_id == "B.17.SIB" | sample_id == "B.20.CHI" |
                            sample_id == "B.20.SIB" | sample_id == "B.23.CHI" | sample_id == "B.23.SIB" | sample_id == "B.25.CHI" |
                            sample_id == "B.25.SIB" | sample_id == "B.32.CHI" | sample_id == "B.32.SIB" | sample_id == "B.35.CHI" |
                            sample_id == "B.35.SIB" | sample_id == "B.48.CHI" | sample_id == "B.48.SIB")
sib <- data.frame(sample_data(sib_pair))

remove(pruned, meta_prune, chi_prune)

#Don't need sibling variables in our chi datasets
chi_hiv <- within(chi_hiv, rm(sib_weight, sib_sex, sib_height, sib_muac, sib_race, sib_bfeed_1,
                                  sib_bfeed_2, sib_uri_recent, sib_uri_current, sib_pcr_1, sib_pcr_2,
                                  sib_bcg, sib_hepb, sib_clinic, sib_clinic_dx___5, sib_clinic_dx_oth,
                                  sib_hosp, sib_meds, sib_meds_name, sibmo, sibyr, sib_vl1_weeks, sib_vl_weeks2,
                                  sib_med_days, sib_dpt, sib_pcv, sib_rota, sib_polio, sib_measles, sib_age))

##############
#NEW VARIABLES
##############
#How many kids receiving TMP-SMX? NOTE: we did not include these kids in the abx category. Need to re-do analyses with this included?
table(chi_hiv$child_tmpsmx) #15 kids on TMP-SMX
#manual inspection: 10 kids on TMP-SMX without other abx exposure

#Need to make a binary variable for tmpsmx use
chi_hiv$tmpsmx <- NA
chi_hiv$tmpsmx [chi_hiv$child_tmpsmx == 1] <- "Yes"
chi_hiv$tmpsmx[chi_hiv$child_tmpsmx == 2] <- "No"
table(chi_hiv$tmpsmx)

#Don't need new variable, but how many kids had received abx for acute infxn?
table(chi_hiv$abx)

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

#make new variable for immune competent CD4 or not, with CD4 percentage >=25% = immune competent
  #d/w MK re: potentially adding a cutoff at 15% for immunosuppression, but we only have 2 kids with CD4% < 15 so will stay with 25% binary for now
chi_cd4$cd4nL <- NA
chi_cd4$cd4nL [chi_cd4$newest_cd4perc >=25] <- "1"
chi_cd4$cd4nL [chi_cd4$newest_cd4perc < 25] <- "0"
table(chi_cd4$cd4nL)  #11 kids with CD4% < 25

#REPEAT the above steps to get date associated with most recent CD4
#CD4 collection date
chi_cd4$Var1 <- chi_cd4$chi_cd4_weeks1
chi_cd4$Var2 <- chi_cd4$chi_cd4_weeks2

chi_cd4 <- chi_cd4 %>% rowwise() %>%
  mutate(newest_cd4date = get(paste0('Var', newest_cd4)))   #now have info for when cd4 were checked in relation to enrollment

summary(chi_cd4$newest_cd4date)

chi_cd4 <- within(chi_cd4, rm(Var1, Var2))

#In case you need to go back and redo everything with only those children who have BOTH viral loads and CD4s:
test <- chi_viral[c("sample_id", "newest_vlnum", "vl_suppr", "newest_vldate")]
hiv_labs <- merge(chi_cd4, test, by="sample_id")
remove(test)

##################################
#DEMOGRAPHICS BY VIRAL SUPPRESSION
##################################
#will look at HIV specific variables (particularly viral suppression, CD4 perc) and general vars
#ex if you are not virally suppressed, is that a/w SES? 

#Age
hist(chi_viral$age)
summary(chi_viral$age)
tapply(chi_viral$age, chi_viral$vl_suppr, summary)
wilcox.test(age ~ vl_suppr, data = chi_viral) #P=0.19 with warning: cannot compute exact p-value with ties
t.test(age ~ vl_suppr, data = chi_viral)  #P=0.17

#Sex
table(chi_viral$sex)  #22 girls
table(chi_viral$sex, chi_viral$vl_suppr)
prop.table(table(chi_viral$sex, chi_viral$vl_suppr), 2)
summary(table(chi_viral$sex, chi_viral$vl_suppr)) #P=0.74

#Maternal age
hist(chi_viral$mat_age2)  #looks quite normal, but will still present median and IQR for consistency
summary(chi_viral$mat_age2)
tapply(chi_viral$mat_age2, chi_viral$vl_suppr, summary)
wilcox.test(mat_age2 ~ vl_suppr, data = chi_viral) #P=0.15
t.test(mat_age2 ~ vl_suppr, data = chi_viral) #P=0.11

#Season
table(chi_viral$season)  
table(chi_viral$season, chi_viral$vl_suppr)
prop.table(table(chi_viral$season, chi_viral$vl_suppr), 2)
summary(table(chi_viral$season, chi_viral$vl_suppr))  #P=0.44

#Maternal education
table(chi_viral$mat_educ2)  
table(chi_viral$mat_educ2, chi_viral$vl_suppr)
prop.table(table(chi_viral$mat_educ2, chi_viral$vl_suppr), 2)
summary(table(chi_viral$mat_educ2, chi_viral$vl_suppr))  
fisher.test(chi_viral$mat_educ2, chi_viral$vl_suppr) #P=0.34

#Electricity
table(chi_viral$elec)  #33 kids have electricity
table(chi_viral$elec, chi_viral$vl_suppr)
prop.table(table(chi_viral$elec, chi_viral$vl_suppr), 2)
summary(table(chi_viral$elec, chi_viral$vl_suppr))
fisher.test(chi_viral$elec, chi_viral$vl_suppr) #P=1.0

#Wood
table(chi_viral$wood)  #25 kids live in wood-burning households
table(chi_viral$wood, chi_viral$vl_suppr)
prop.table(table(chi_viral$wood, chi_viral$vl_suppr), 2)
summary(table(chi_viral$wood, chi_viral$vl_suppr))
fisher.test(chi_viral$wood, chi_viral$vl_suppr) #P=0.18

#Num household members
hist(chi_viral$hhsize)
summary(chi_viral$hhsize)
tapply(chi_viral$hhsize, chi_viral$vl_suppr, summary)
wilcox.test(hhsize ~ vl_suppr, data = chi_viral) #P=0.97 with warning: cannot compute exact p-value with ties
t.test(hhsize ~ vl_suppr, data = chi_viral) #P=0.68

#antibiotics
table(chi_viral$abx)  #6 kids received abx
table(chi_viral$abx, chi_viral$vl_suppr)
prop.table(table(chi_viral$abx, chi_viral$vl_suppr), 2)
fisher.test(chi_viral$abx, chi_viral$vl_suppr) #P=0.39

#Clinic in past 3 months
table(chi_viral$clinic)  #11 kids with clinic visit
table(chi_viral$clinic, chi_viral$vl_suppr)
prop.table(table(chi_viral$clinic, chi_viral$vl_suppr), 2)
fisher.test(chi_viral$clinic, chi_viral$vl_suppr) #P=0.46

#URI in past 1 month
table(chi_viral$uri_rec2)  #18 kids with recent URI
table(chi_viral$uri_rec2, chi_viral$vl_suppr)
prop.table(table(chi_viral$uri_rec2, chi_viral$vl_suppr), 2)
summary(table(chi_viral$uri_rec2, chi_viral$vl_suppr)) #P=0.92

#URI currently (don't need both this and uri recent in table 1; need to choose)
table(chi_viral$uri_cur2)  #16 kids with current URI
table(chi_viral$uri_cur2, chi_viral$vl_suppr)
prop.table(table(chi_viral$uri_cur2, chi_viral$vl_suppr), 2)
summary(table(chi_viral$uri_cur2, chi_viral$vl_suppr)) #P=0.75

#Hospital in past 3 months -->> 4 children hospitalized, all 4 HIV+
#*MENTION IN MANUSCRIPT
table(chi_viral$hosp)  #3 kids hospitalized who had viral load data
table(chi_viral$hosp, chi_viral$vl_suppr) #2/3 kids not virally suppressed
prop.table(table(chi_viral$hosp, chi_viral$vl_suppr), 2)
fisher.test(chi_viral$hosp, chi_viral$vl_suppr) #P=0.26 (but numbers are SO SMALL)

#PCV doses
table(chi_viral$pcv)  
table(chi_viral$pcv, chi_viral$vl_suppr)
prop.table(table(chi_viral$pcv, chi_viral$vl_suppr), 2)
fisher.test(chi_viral$pcv, chi_viral$vl_suppr) #P=1.0

#HiB
table(chi_viral$dpt)  
table(chi_viral$dpt, chi_viral$vl_suppr)
prop.table(table(chi_viral$dpt, chi_viral$vl_suppr), 2)
fisher.test(chi_viral$dpt, chi_viral$vl_suppr) #P=1.0

#TMP-SMX
table(chi_viral$tmpsmx)  
table(chi_viral$tmpsmx, chi_viral$vl_suppr)
prop.table(table(chi_viral$tmpsmx, chi_viral$vl_suppr), 2)
summary(table(chi_viral$tmpsmx, chi_viral$vl_suppr))  
fisher.test(chi_viral$tmpsmx, chi_viral$vl_suppr) #P=0.17

##########################################
#DEMOGRAPHICS BY CD4 STATUS (< vs >/= 25%)  
##########################################
#Age
hist(chi_cd4$age)
summary(chi_cd4$age)
tapply(chi_cd4$age, chi_cd4$cd4nL, summary)
wilcox.test(age ~ cd4nL, data = chi_cd4) #P=0.10 with warning: cannot compute exact p-value with ties
t.test(age ~ cd4nL, data = chi_cd4)  #P=0.10

#Sex
table(chi_cd4$sex)  #23 girls
table(chi_cd4$sex, chi_cd4$cd4nL)
prop.table(table(chi_cd4$sex, chi_cd4$cd4nL), 2)
summary(table(chi_cd4$sex, chi_cd4$cd4nL)) #P=0.99

#Maternal age
hist(chi_cd4$mat_age2)
summary(chi_cd4$mat_age2)
tapply(chi_cd4$mat_age2, chi_cd4$cd4nL, summary)
wilcox.test(mat_age2 ~ cd4nL, data = chi_cd4) #P=0.31
t.test(mat_age2 ~ cd4nL, data = chi_cd4)  #P=0.60

#Season
table(chi_cd4$season) 
table(chi_cd4$season, chi_cd4$cd4nL)
prop.table(table(chi_cd4$season, chi_cd4$cd4nL), 2)
fisher.test(chi_cd4$season, chi_cd4$cd4nL)  #P=1.0

#Maternal education
table(chi_cd4$mat_educ2) 
table(chi_cd4$mat_educ2, chi_cd4$cd4nL)
prop.table(table(chi_cd4$mat_educ2, chi_cd4$cd4nL), 2)
fisher.test(chi_cd4$mat_educ2, chi_cd4$cd4nL)  #P=1.0

#Electricity
table(chi_cd4$elec) 
table(chi_cd4$elec, chi_cd4$cd4nL)
prop.table(table(chi_cd4$elec, chi_cd4$cd4nL), 2)
fisher.test(chi_cd4$elec, chi_cd4$cd4nL)  #P=0.68

#Wood
table(chi_cd4$wood) #26 kids from wood-burning households
table(chi_cd4$wood, chi_cd4$cd4nL)
prop.table(table(chi_cd4$wood, chi_cd4$cd4nL), 2)
fisher.test(chi_cd4$wood, chi_cd4$cd4nL)  #P=0.16

#Num household members
hist(chi_viral$hhsize)
summary(chi_viral$hhsize)
tapply(chi_viral$hhsize, chi_viral$vl_suppr, summary)
wilcox.test(hhsize ~ vl_suppr, data = chi_viral) #P=0.97 with warning: cannot compute exact p-value with ties
t.test(hhsize ~ vl_suppr, data = chi_viral) #P=0.68

#antibiotics
table(chi_cd4$abx) #7 kids received abx
table(chi_cd4$abx, chi_cd4$cd4nL)
prop.table(table(chi_cd4$abx, chi_cd4$cd4nL), 2)
fisher.test(chi_cd4$abx, chi_cd4$cd4nL)  #P=0.009

#Clinic in past 3 months
table(chi_cd4$clinic) #11 kids with clinic visit
table(chi_cd4$clinic, chi_cd4$cd4nL)
prop.table(table(chi_cd4$clinic, chi_cd4$cd4nL), 2)
fisher.test(chi_cd4$clinic, chi_cd4$cd4nL)  #P=0.12

#URI in past 1 month
table(chi_cd4$uri_rec2) #19 kids with recent URI
table(chi_cd4$uri_rec2, chi_cd4$cd4nL)
prop.table(table(chi_cd4$uri_rec2, chi_cd4$cd4nL), 2)
fisher.test(chi_cd4$uri_rec2, chi_cd4$cd4nL)  #P=0.0008

#URI currently (don't need both this and uri recent in table 1; need to choose)
table(chi_cd4$uri_cur2) #16 kids with current URI
table(chi_cd4$uri_cur2, chi_cd4$cd4nL)
prop.table(table(chi_cd4$uri_cur2, chi_cd4$cd4nL), 2)
fisher.test(chi_cd4$uri_cur2, chi_cd4$cd4nL)  #P=0.07

#Hospital in past 3 months -->> 4 children hospitalized, all 4 HIV+
#*MENTION IN MANUSCRIPT
table(chi_cd4$hosp) 
table(chi_cd4$hosp, chi_cd4$cd4nL)  #3 of 4 kids had CD4 < 25%
prop.table(table(chi_cd4$hosp, chi_cd4$cd4nL), 2)
fisher.test(chi_cd4$hosp, chi_cd4$cd4nL)  #P=0.05 (numbers are SO SMALL though)

#PCV doses
table(chi_cd4$pcv) #35 kids fully vaxed
table(chi_cd4$pcv, chi_cd4$cd4nL)
prop.table(table(chi_cd4$pcv, chi_cd4$cd4nL), 2)
fisher.test(chi_cd4$pcv, chi_cd4$cd4nL)  #P=0.48

#HiB
table(chi_cd4$dpt) #34 kids fully vaxed 
table(chi_cd4$dpt, chi_cd4$cd4nL)
prop.table(table(chi_cd4$dpt, chi_cd4$cd4nL), 2)
fisher.test(chi_cd4$dpt, chi_cd4$cd4nL)  #P=1.0

#TMPSMX use
table(chi_cd4$tmpsmx) #15 kids on TMP-SMX 
table(chi_cd4$tmpsmx, chi_cd4$cd4nL)
prop.table(table(chi_cd4$tmpsmx, chi_cd4$cd4nL), 2)
fisher.test(chi_cd4$tmpsmx, chi_cd4$cd4nL)  #P=0.03; more kids with low CD4 on TMP-SMX

###COME BACK TO: TMP-SMX GUIDELINES IN KIDS

#Summary: no significant association identified b/w viral suppression and SES, clinical factors
  #BUT kids with low CD4: higher percent received abx and TMP-SMX and reported recent URI
    #incorporate abx and rec_uri into maaslin2 model below

################
#ALPHA DIVERSITY
################
#Create dataframe with diversity indices for each specimen using HIV pseq object
diversity1 <- estimate_richness(hiv, measures = c("Shannon", "Simpson", "Chao1", "Observed"))

setDT(diversity1, keep.rownames = TRUE)[]
colnames(diversity1) <- c("sample_id", "Observed", "Chao1", "se.chao1", "Shannon", "Simpson")
rownames(diversity1) <- diversity1$sample_id
diversity <- merge(diversity1, chi_viral, by="sample_id")
div_cd4 <- merge(diversity1, chi_cd4, by="sample_id")
div_all <- merge(diversity1, chi_hiv, by="sample_id")

#
#OBSERVED SPECIES
#
#How many species on average per sample?
hist(diversity$Observed)
summary(diversity$Observed) 
sd(diversity$Observed)
#Median (IQR) # spec: 113 (67, 201)
#Mean (SD): 163 (173)

#what does it look like by TMP-SMX use? NOTE: N=44
tapply(div_all$Observed, div_all$tmpsmx, summary)   
wilcox.test(Observed ~ tmpsmx, data = div_all)  #P=0.42 with warning (ties)

#what does it look like by viral supp? NOTE: N=41
tapply(diversity$Observed, diversity$vl_suppr, summary)   
wilcox.test(Observed ~ vl_suppr, data = diversity)  #P=0.87 with warning (ties)

#
#SHANNON
#

#Is our SDI data normally distributed? 
shapiro.test(diversity$Shannon)  #yes! P=0.42
histogram(diversity$Shannon)  #looks good
qqnorm(diversity$Shannon, pch = 1, frame = FALSE)
qqline(diversity$Shannon, col = "steelblue", lwd = 2)

#Shannon diversity median and IQR for all samples
summary(diversity$Shannon)
sd(diversity$Shannon)
#Median (IQR) SDI: 1.73 (1.44, 2.26)
#Mean (SD) 1.78 (0.65)  

#Viral suppression
tapply(diversity$Shannon, diversity$vl_suppr, summary)
t.test(Shannon ~ vl_suppr, data = diversity)        #P=0.62

#TMP-SMX use
tapply(div_all$Shannon, div_all$tmpsmx, summary)
t.test(Shannon ~ tmpsmx, data = div_all)  #P=0.28

#CD4 (continuous and binary)
summary(div <- lm(Shannon ~ newest_cd4perc, data = div_cd4))  #P=0.09
tapply(div_cd4$Shannon, div_cd4$cd4nL, summary)
t.test(Shannon ~ cd4nL, data = div_cd4)  #P=0.89

#Abx use
tapply(div_all$Shannon, div_all$abx, summary)
t.test(Shannon ~ abx, data = div_all)  #P=0.99

#MULTIVARIABLE MODEL
summary(div_mv <- lm(Shannon ~ vl_suppr + tmpsmx, data = diversity))
  #with abx
summary(div_mv <- lm(Shannon ~ vl_suppr + tmpsmx + abx, data = diversity))

#
#CHAO1
#

#Is our Chao1 data normally distributed? 
shapiro.test(diversity$Chao1)  #no, because p<0.05; also looks skewed on histogram and QQ is wonky
histogram(diversity$Chao1)
qqnorm(diversity$Chao1, pch = 1, frame = FALSE)
qqline(diversity$Chao1, col = "steelblue", lwd = 2)

#What does it look like if we log transform it? Will this allow for parametric testing?
diversity1$Chao1_log <- log(diversity1$Chao1)
#NOW RE-RUN ALL DIVERSITY SUB-DATAFRAMES ABOVE (CD4, VL, ALL)
shapiro.test(diversity$Chao1_log) 
histogram(diversity$Chao1_log)
qqnorm(diversity$Chao1_log, pch = 1, frame = FALSE)
qqline(diversity$Chao1_log, col = "steelblue", lwd = 2)
#Looks better! Can use t.test and linear regression on transformed data

#Log transform for all diversity df (n=44 for all HIV+ kids, n=42 for CD4 kids)
div_all$Chao1_log <- log(div_all$Chao1)
div_cd4$Chao1_log <- log(div_cd4$Chao1)

#Chao1 median and IQR for all samples
summary(diversity$Chao1)
sd(diversity$Chao1)
#Median (IQR) Chao1: 109 (63, 205)
#Mean (SD) 216 (300.4)

#Viral suppression
tapply(diversity$Chao1_log, diversity$vl_suppr, summary)
t.test(Chao1_log ~ vl_suppr, data = diversity)        #P=0.96

#TMP-SMX use (N=44)
tapply(div_all$Chao1_log, div_all$tmpsmx, summary)
t.test(Chao1_log ~ tmpsmx, data = div_all)  #P=0.52

#CD4 (continuous and binary)
summary(div <- lm(Chao1_log ~ newest_cd4perc, data = div_cd4))  #P=0.65
tapply(div_cd4$Chao1_log, div_cd4$cd4nL, summary)
t.test(Chao1_log ~ cd4nL, data = div_cd4)  #P=0.25

#abx (N=44)
tapply(div_all$Chao1_log, div_all$abx, summary)
t.test(Chao1_log ~ abx, data = div_all)  #P=0.92

#MULTIVARIABLE MODEL
summary(div_mv <- lm(Chao1_log ~ vl_suppr + tmpsmx, data = diversity))
  #including abx
summary(div_mv <- lm(Chao1_log ~ vl_suppr + tmpsmx + abx, data = diversity))

#Skipping boxplots
remove(div, div_all, div_cd4, diversity, diversity1)

################
##BETA DIVERSITY
################
#remove kids missing viral load data: B-04, B-07, B-52 and add chi_viral as sample data
vl <- hiv
chi_viral <- as.data.frame(chi_viral)
row.names(chi_viral) <- chi_viral$sample_id
chi_viral <- sample_data(chi_viral)
head(chi_viral)
vl <- merge_phyloseq(vl, chi_viral)
#now see if it worked:
bloop <- data.frame(sample_data(vl))
#not working for some reason on 6/2/22; merged new sample data with pseq FIRST, then removed samples and it worked
#success!
remove(bloop)
vl = subset_samples(vl, sample_id != "B.04.CHI" & sample_id != "B.07.CHI" & sample_id != "B.52.CHI")
chi_viral <- data.frame(sample_data(vl))

#do the same to create a CD4 pseq object to plot
cd4 <- hiv
chi_cd4 <- as.data.frame(chi_cd4)
row.names(chi_cd4) <- chi_cd4$sample_id
chi_cd4 <- sample_data(chi_cd4)
head(chi_cd4)
cd4 <- merge_phyloseq(cd4, chi_cd4)
#now see if it worked:
bloop <- data.frame(sample_data(cd4))
#success!
remove(bloop)
#missing: B-45 and B-52
cd4 = subset_samples(cd4, sample_id != "B.45.CHI" & sample_id != "B.52.CHI")
chi_cd4 <- data.frame(sample_data(cd4))

#COME BACK TO: use same steps to make hiv_labs pseq that has n=40 and both vl and CD4 variables
labs <- hiv
row.names(hiv_labs) <- hiv_labs$sample_id
hiv_labs <- sample_data(hiv_labs)
labs <- merge_phyloseq(labs, hiv_labs)
#now see if it worked:
bloop <- data.frame(sample_data(labs))
#success!
remove(bloop)
#missing: B-45 and B-52
labs = subset_samples(labs, sample_id != "B.04.CHI" & sample_id != "B.07.CHI" & sample_id != "B.45.CHI" & sample_id != "B.52.CHI")
hiv_labs <- data.frame(sample_data(labs))

#Now transform our pseq objects
vl_clr <- transform(vl, "clr")
cd4_clr <- transform(cd4, "clr")
hiv_clr <- transform(hiv, "clr")

#PCoA plots#
#VIRAL SUPPRESSION
clr_pcoa <- ordinate(vl_clr, method = "PCoA", distance = "euclidean")
lines <- c("dashed", "dashed")
clr_vl <- plot_ordination(vl_clr, clr_pcoa, color = "vl_suppr", title = "Euclidean PCoA Plot by Viral Suppression") +
  geom_point(size = 2) + theme_classic() + theme(
    plot.title = element_text(size = 14, hjust = 0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 10),
    legend.position = "right",
    legend.background = element_blank(),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12)) +
  scale_color_manual(values=c('#800000FF','#155F83FF'))+
  stat_ellipse(aes(linetype=vl_suppr), geom="polygon", alpha=0, type="t", level=0.8, size=0.7) +
  xlab("PC1 (13.4%)") + ylab("PC2 (12%)") +
  scale_linetype_manual(values=lines)
png(file="clr_vlsuppr.png",
    width = 5, height = 4, units = 'in', res=600)
clr_vl
dev.off()


#CD4 (binary) using cd4_clr
###AUG 3 UPDATE: change labels: so will need to create new variable and import back into cd4_clr
chi_cd4$cd4cat <- NA
chi_cd4$cd4cat [chi_cd4$cd4nL=="1"] <- "Normal CD4"
chi_cd4$cd4cat [chi_cd4$cd4nL=="0"] <- "Low CD4"
table(chi_cd4$cd4cat)
class(chi_cd4$cd4cat)
chi_cd4$cd4cat <- as.factor(chi_cd4$cd4cat)
chi_cd4$cd4cat <- reorder(chi_cd4$cd4cat, new.order=c("Low CD4", "Normal CD4"))

#Import back into pseq (source: https://jacobrprice.github.io/2017/08/26/phyloseq-to-vegan-and-back.html)
sample_data(cd4_clr) <- as.data.frame(chi_cd4)
#now see if it worked:
bloop <- data.frame(sample_data(cd4_clr))
#success!
remove(bloop)

clr_pcoa <- ordinate(cd4_clr, method = "PCoA", distance = "euclidean")
lines <- c("dashed", "dashed")
clr_cd4 <- plot_ordination(cd4_clr, clr_pcoa, color = "cd4cat") +
  geom_point(size = 2) + theme_classic() + theme(
    plot.title = element_text(size = 14, hjust = 0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 10),
    legend.position = "right",
    legend.background = element_blank(),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12)) +
  scale_color_manual(values=c('#800000FF','#155F83FF'))+
  stat_ellipse(aes(linetype=cd4cat), geom="polygon", alpha=0, type="t", level=0.8, size=0.7) +
  xlab("PC1 (13.4%)") + ylab("PC2 (12%)") +
  scale_linetype_manual(values=lines)
# png(file="clr_cd4nl.png",
#     width = 6, height = 4, units = 'in', res=600)
clr_cd4
p3 <- clr_cd4
p3
# dev.off()

#TMP-SMX use
clr_tmpsmx <- plot_ordination(cd4_clr, clr_pcoa, color = "tmpsmx", title = "Euclidean PCoA Plot by TMP-SMX use") +
  geom_point(size = 2) + theme_classic() + theme(
    plot.title = element_text(size = 14, hjust = 0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 10),
    legend.position = "right",
    legend.background = element_blank(),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12)) +
  scale_color_manual(values=c('#800000FF','#155F83FF'))+
  stat_ellipse(aes(linetype=tmpsmx), geom="polygon", alpha=0, type="t", level=0.8, size=0.7) +
  xlab("PC1 (13.4%)") + ylab("PC2 (12%)") +
  scale_linetype_manual(values=lines)
png(file="clr_tmpsmx.png",
    width = 5, height = 4, units = 'in', res=600)
clr_tmpsmx
dev.off()

#PERMANOVA: viral load & TMP-SMX
otu <- as.data.frame(otu_table(vl_clr, taxa_are_rows=TRUE))
#Need species in columns and observations in rows: transpose using t() function and save as a dataframe
#Source: https://github.com/edamame-course/2015-tutorials/blob/master/final/Ashley_Intro_to_R.md
otu_trans <- data.frame(t(otu))
chi_viral <- data.frame(sample_data(vl_clr))
BC.dist=vegdist(otu_trans, method="euclidean")

#viral suppression
adonis(BC.dist ~ vl_suppr, data = chi_viral, permutations = 1000)
#No significant difference by HIV status among children (P = 0.83)
  #viral suppression + age (per 6/9 MK discussion to account for age in all PERMANOVAs)
adonis2(BC.dist ~ vl_suppr + age, data = chi_viral, permutations = 1000)
#No significant difference by HIV status among children (P = 0.82) or age (P=0.15)

#TMPSMX use (n=41)
adonis(BC.dist ~ tmpsmx, data = chi_viral, permutations = 1000)
#Trend towards difference in overall composition by tmp-smx use but did not reach signif (P = 0.07, R2 = 0.034)
  #TMPSMX + age (per 6/9 MK discussion to account for age in all PERMANOVAs)
adonis(BC.dist ~ tmpsmx + age, data = chi_viral, permutations = 1000)
#No significant difference by TMPSMX among children (P = 0.07) or age (P=0.22)

#PERMANOVA: CD4 (both continuous and binary), TMPSMX, abx
otu <- as.data.frame(otu_table(cd4_clr, taxa_are_rows=TRUE))
#Need species in columns and observations in rows: transpose using t() function and save as a dataframe
#Source: https://github.com/edamame-course/2015-tutorials/blob/master/final/Ashley_Intro_to_R.md
otu_trans <- data.frame(t(otu))
chi_cd4 <- data.frame(sample_data(cd4_clr))
BC.dist=vegdist(otu_trans, method="euclidean")

#CD4 (continuous)
adonis(BC.dist ~ newest_cd4perc, data = chi_cd4, permutations = 1000)
#Trend towards difference in overall composition by CD4 percentage but did not reach signif (P = 0.052, R2 = 0.035)
    #CD4 + age (per 6/9 MK discussion to account for age in all PERMANOVAs)
adonis(BC.dist ~ newest_cd4perc + age, data = chi_cd4, permutations = 1000)
#Trend towards difference in overall composition by CD4 percentage but did not reach signif (P = 0.06, R2 = 0.035); age P=0.4

#CD4 (binary)
adonis(BC.dist ~ cd4nL, data = chi_cd4, permutations = 1000)  #P=0.01, R = 0.042
  #CD4 + age (per 6/9 MK discussion to account for age in all PERMANOVAs)
adonis2(BC.dist ~ cd4nL + age, data = chi_cd4, permutations = 1000) #P=0.01, R = 0.042 for CD4; age P = 0.20

#TMPSMX use (n=42)  #NOTE: at some point adonis got deprecated so had to switch to adonis2 for updates
adonis(BC.dist ~ tmpsmx, data = chi_cd4, permutations = 1000) #P=0.12
    #TMPSMX + age (per 6/9 MK discussion to account for age in all PERMANOVAs)
adonis2(BC.dist ~ tmpsmx + age, data = chi_cd4, permutations = 1000) #P=0.13, R2=0.031

#abx (n=42) with age
adonis2(BC.dist ~ abx + age, data = chi_cd4, permutations = 1000) #P=0.19, R2=0.029

#PERMANOVA: TMP-SMX and abx for all children with HIV (n=44)
otu <- as.data.frame(otu_table(hiv_clr, taxa_are_rows=TRUE))
otu_trans <- data.frame(t(otu))
BC.dist=vegdist(otu_trans, method="euclidean")
  #TMP-SMX
adonis2(BC.dist ~ tmpsmx + age, data = chi_hiv, permutations = 1000) #P=0.07, R2=0.032
  #ABX
adonis2(BC.dist ~ abx + age, data = chi_hiv, permutations = 1000) #P=0.20, R2=0.028

#clean up
remove(hiv_clr, cd4_clr, vl_clr, otu, otu_trans, clr_cd4, clr_pcoa, clr_tmpsmx, clr_vl)
remove(BC.dist, lines)

##############################
#RELATIVE ABUNDANCE PLOTS: CD4
##############################
#Because we see a compositional difference by binary CD4, we will generate relative abundance plots for this pseq and for viral load
chi.rel <- transform_sample_counts(cd4, function(Abundance) Abundance/sum(Abundance))
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
  #merging on Phylum AND Genus gets rid of the phylum.x and phylum.y problem we were having earlier

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
relative_df <- merge(melted_df, abundances, by=c("Genus", "Species"))
#keep getting vector memory exhausted error; need to reboot R so will save melted_df and abundances as csv files
#source for fix: https://stackoverflow.com/questions/51295402/r-on-macos-error-vector-memory-exhausted-limit-reached

relative_df$Species[!(relative_df$Species %in% TOPSpecies)] <- "Other"
sum(relative_df$Abundance) #sums to 42.00003 (and we have n=42 kids...?)
table(relative_df$Species)

#To generate color palette of 11 colors used previously: 
library(RColorBrewer)
library(scales)
chisp2 <- c("#293352", "#800000FF", "#767676FF", "#FFA319FF", "#8A9045FF", "#155F83FF", "#C16622FF", "#4D004B", "#8F3931FF", "#58593FFF", "#350E20FF")
show_col(chisp2)

#CREATING RELATIVE ABUNDANCE PLOT
relative_df$subject_abundances <- (relative_df$Abundance)/(length(unique(relative_df$Sample))) 
#Clarify CD4 cutoff 
relative_df$cd4cat <- NA
relative_df$cd4cat [relative_df$cd4nL==1] <- "Normal CD4"
relative_df$cd4cat [relative_df$cd4nL==0] <- "Low CD4"
table(relative_df$cd4cat)
class(relative_df$cd4cat)
relative_df$cd4cat <- as.factor(chi_cd4$cd4cat)
relative_df$cd4cat <- reorder(chi_cd4$cd4cat, new.order=c("Low CD4", "Normal CD4"))

# relative_df$cd4nL2 <- NA
# relative_df$cd4nL2 [relative_df$cd4nL == 0] <- "CD4 < 25%"
# relative_df$cd4nL2 [relative_df$cd4nL == 1] <- "CD4 >/= 25%"
# table(relative_df$cd4nL2)
# #matches up with table(relative_df$hiv)
# relative_df$hiv2 <- as.factor(relative_df$hiv2)
# relative_df$hiv2 <- reorder(relative_df$hiv2, new.order=c("HUU", "HEU", "HEI"))

relative_df$Species <- factor(relative_df$Species)
relative_df$Species <- reorder(relative_df$Species, new.order=c("Corynebacterium propinquum", "Corynebacterium pseudodiphtheriticum", "Dolosigranulum pigrum",  
                                                                "Haemophilus influenzae", "Micrococcus luteus", "Moraxella catarrhalis", "Moraxella lincolnii", 
                                                                "Moraxella nonliquefaciens", "Staphylococcus aureus", "Streptococcus pneumoniae", "Other"))
subject_plot2 <- ggplot(relative_df[order(relative_df$Species, decreasing = TRUE),], aes(x=cd4cat, y=subject_abundances, fill=Species)) + 
  geom_bar(stat="identity", position="fill") + theme_classic() + theme(legend.position="right", 
                                                                       legend.text=element_text(size=10, face="italic"),
                                                                       axis.title.y = element_text(size=12, margin=margin(0,20,0,0)), axis.text.y = element_text(size=10), 
                                                                       axis.title.x = element_blank(), axis.text.x = element_text(size=10, angle = 45, vjust = 1, hjust = 1, colour = "black")) + 
  scale_fill_manual(values = chisp2) +
  # scale_x_discrete(labels=c("No", "Yes")) +
  # xlab("CD4 Category") + element_text(size=12, margin=margin(10,0,0,0))
  ylab("Mean relative abundance") 

p4 <- plot(subject_plot2)

#NOW: combine all using patchwork vs cowplot (great ref: https://wilkelab.org/cowplot/articles/aligning_plots.html)
#For future reference: p1 and p2 code is in Fog Kaiju 5-26-22.R file, p3 and p4 are here
#After MK feedback, removed x axis titles for B and D, only capitalized "Mean" in y axis titles for B&D, and increased width of A&C (so ratio from 1, 1.65 to 1, 1.55)
fig2 <- plot_grid(p1, p2, p3, p4, labels = "AUTO", nrow = 2, ncol = 2, align = "h", axis = "b", rel_widths = c(1,1.5))

png(file="fig2.png",
    width = 12, height = 9, units = 'in', res = 600)
fig2
dev.off()

# png(file="cd4plot_kaiju_06022022.png",
#     width = 9, height = 5, units = 'in', res = 600)
# plot(subject_plot2)
# dev.off()

#Due to Maaslin2 results below with differences by URI status: what does the rel abundance plot for uri_rec2 look like?
relative_df$Species <- factor(relative_df$Species)
relative_df$Species <- reorder(relative_df$Species, new.order=c("Corynebacterium propinquum", "Corynebacterium pseudodiphtheriticum", "Dolosigranulum pigrum",  
                                                                "Haemophilus influenzae", "Micrococcus luteus", "Moraxella catarrhalis", "Moraxella lincolnii", 
                                                                "Moraxella nonliquefaciens", "Staphylococcus aureus", "Streptococcus pneumoniae", "Other"))
subject_plot2 <- ggplot(relative_df[order(relative_df$Species, decreasing = TRUE),], aes(x=uri_rec2, y=subject_abundances, fill=Species)) + 
  geom_bar(stat="identity", position="fill") + theme_classic() + theme(legend.position="right", 
                                                                       legend.text=element_text(size=12, face="italic"),
                                                                       axis.title.y = element_text(size=14, margin=margin(0,20,0,0)), axis.text.y = element_text(size=12), 
                                                                       axis.title.x = element_text(size=14, margin=margin(10,0,0,0)), axis.text.x = element_text(size=12)) + 
  scale_fill_manual(values = chisp2) +
  # scale_x_discrete(labels=c("No", "Yes")) +
  xlab("Recent URI") + ylab("Relative Abundance") 
png(file="uriplot_kaiju_06062022.png",
    width = 9, height = 5, units = 'in', res = 600)
plot(subject_plot2)
dev.off()

#############################################################################
#RELATIVE ABUNDANCE M-W TESTS: CD4, TMP-SMX, ABX (for children with CD4 meas)
#############################################################################
#Corynebacterium propinquum
cor_p <- filter(melted_df, Species == "Corynebacterium propinquum")
summary(cor_p$Abundance)
hist(cor_p$Abundance)
tapply(cor_p$Abundance, cor_p$cd4nL, summary) #CD4<25% kids have lower median rel abund of C.propinquum
wilcox.test(Abundance ~ cd4nL, data = cor_p, exact=FALSE)  #P=0.004 (warning: ties --> gone wiith exact=FALSE)

  #TMPSMX
tapply(cor_p$Abundance, cor_p$tmpsmx, summary) #tmpsmx exposed kids have lower median rel abund of C.propinquum
wilcox.test(Abundance ~ tmpsmx, data = cor_p, exact=FALSE)  #P=0.11

  #abx
tapply(cor_p$Abundance, cor_p$abx, summary) #abx exposed kids have lower median rel abund of C.propinquum
wilcox.test(Abundance ~ abx, data = cor_p, exact=FALSE)  #P=0.13

#Corynebacterium pseudodiphtheriticum
cor_ps <- filter(melted_df, Species == "Corynebacterium pseudodiphtheriticum")
summary(cor_ps$Abundance)
hist(cor_ps$Abundance)
tapply(cor_ps$Abundance, cor_ps$cd4nL, summary) #CD4<25% kids have lower median rel abund of C.pseudo
wilcox.test(Abundance ~ cd4nL, data = cor_ps, exact=FALSE)  #P=0.007 (warning: ties)

  #TMPSMX
tapply(cor_ps$Abundance, cor_ps$tmpsmx, summary) #tmpsmx exposed kids have lower median rel abund of C. pseudo
wilcox.test(Abundance ~ tmpsmx, data = cor_ps, exact=FALSE)  #P=0.64

  #abx
tapply(cor_ps$Abundance, cor_ps$abx, summary) #abx exposed kids have lower median rel abund of C.pseudo
wilcox.test(Abundance ~ abx, data = cor_ps, exact=FALSE)  #P=0.41

#Dolosigranulum pigrum
dpig <- filter(melted_df, Species == "Dolosigranulum pigrum")
hist(dpig$Abundance)
summary(dpig$Abundance)
tapply(dpig$Abundance, dpig$cd4nL, summary) #CD4<25% kids have lower median rel abund of D.pig
wilcox.test(Abundance ~ cd4nL, data = dpig, exact=FALSE)  #P=0.005 (warning: ties)

  #TMPSMX
tapply(dpig$Abundance, dpig$tmpsmx, summary) #tmpsmx exposed kids have lower median rel abund of D.pig
wilcox.test(Abundance ~ tmpsmx, data = dpig, exact=FALSE)  #P=0.043

  #abx
tapply(dpig$Abundance, dpig$abx, summary) #abx exposed kids have lower median rel abund of C.propinquum
wilcox.test(Abundance ~ abx, data = dpig, exact=FALSE)  #P=0.20

#Haemophilus influenzae
hflu <- filter(melted_df, Species == "Haemophilus influenzae")
hist(hflu$Abundance)  #lots of zeros
summary(hflu$Abundance)
tapply(hflu$Abundance, hflu$cd4nL, summary) #CD4<25% kids have higher median rel abund of H.flu
wilcox.test(Abundance ~ cd4nL, data = hflu, exact=FALSE)  #P=0.22 (warning: ties)

  #TMPSMX
tapply(hflu$Abundance, hflu$tmpsmx, summary) #median roughly same
wilcox.test(Abundance ~ tmpsmx, data = hflu, exact=FALSE)  #P=0.85

  #abx
tapply(hflu$Abundance, hflu$abx, summary) #median roughly same
wilcox.test(Abundance ~ abx, data = hflu, exact=FALSE)  #P=1.0

#Micrococcus luteus
mic <- filter(melted_df, Species == "Micrococcus luteus")
hist(mic$Abundance)  #lots of zeros
summary(mic$Abundance)
tapply(mic$Abundance, mic$cd4nL, summary) #CD4<25% kids have lower median rel abund of M.luteus
wilcox.test(Abundance ~ cd4nL, data = mic, exact=FALSE)  #P=0.16 (warning: ties)

  #TMPSMX
tapply(mic$Abundance, mic$tmpsmx, summary) #tmpsmx exposed kids have lower median rel abund of M. luteus
wilcox.test(Abundance ~ tmpsmx, data = mic, exact=FALSE)  #P=0.14

  #abx
tapply(mic$Abundance, mic$abx, summary) #abx exposed kids have lower median rel abund of M. luteus
wilcox.test(Abundance ~ abx, data = mic, exact=FALSE)  #P=0.61

#Moraxella cattarhalis
mcat <- filter(melted_df, Species == "Moraxella catarrhalis")
hist(mcat$Abundance)  #lots of zeros
summary(mcat$Abundance)
tapply(mcat$Abundance, mcat$cd4nL, summary) #CD4<25% kids have higher median rel abund of M.cat
wilcox.test(Abundance ~ cd4nL, data = mcat, exact=FALSE)  #P=0.49 (warning: ties)

  #TMPSMX
tapply(mcat$Abundance, mcat$tmpsmx, summary) #medians similar
wilcox.test(Abundance ~ tmpsmx, data = mcat, exact=FALSE)  #P=0.71

  #abx
tapply(mcat$Abundance, mcat$abx, summary) #abx exposed kids have higher median rel abund of M.cat
wilcox.test(Abundance ~ abx, data = mcat, exact=FALSE)  #P=0.99

#Moraxella lincolnii
mlin <- filter(melted_df, Species == "Moraxella lincolnii")
hist(mlin$Abundance)  #lots of zeros
summary(mlin$Abundance)
tapply(mlin$Abundance, mlin$cd4nL, summary) #CD4<25% kids have lower median rel abund of M.linc
wilcox.test(Abundance ~ cd4nL, data = mlin, exact=FALSE)  #P=0.22 (warning: ties)

  #TMPSMX
tapply(mlin$Abundance, mlin$tmpsmx, summary) #medians similar
wilcox.test(Abundance ~ tmpsmx, data = mlin, exact=FALSE)  #P=0.43

  #abx
tapply(mlin$Abundance, mlin$abx, summary) #medians similar
wilcox.test(Abundance ~ abx, data = mlin, exact=FALSE)  #P=0.21

#Moraxella nonliquefaciens
mnon <- filter(melted_df, Species == "Moraxella nonliquefaciens")
hist(mnon$Abundance)  #lots of zeros but more evenly distrib than other Moraxellas
summary(mnon$Abundance)
tapply(mnon$Abundance, mnon$cd4nL, summary) #CD4<25% kids have lower median rel abund of M.non
wilcox.test(Abundance ~ cd4nL, data = mnon, exact=FALSE)  #P=0.13 (warning: ties)

  #TMPSMX
tapply(mnon$Abundance, mnon$tmpsmx, summary) #medians similar
wilcox.test(Abundance ~ tmpsmx, data = mnon, exact=FALSE)  #P=0.92

  #abx
tapply(mnon$Abundance, mnon$abx, summary) #abx exposed kids have lower median rel abund of M.non
wilcox.test(Abundance ~ abx, data = mnon, exact=FALSE)  #P=0.19

#Staphylococcus aureus
sa <- filter(melted_df, Species == "Staphylococcus aureus")
hist(sa$Abundance)  #lots of zeros
summary(sa$Abundance)
tapply(sa$Abundance, sa$cd4nL, summary) #CD4<25% kids have higher median rel abund of S.aureus
wilcox.test(Abundance ~ cd4nL, data = sa, exact=FALSE)  #P=0.04 (warning: ties)

  #TMPSMX
tapply(sa$Abundance, sa$tmpsmx, summary) #tmpsmx exposed kids have higher median rel abund of S.aureus (why???)
wilcox.test(Abundance ~ tmpsmx, data = sa, exact=FALSE)  #P=0.001

  #abx
tapply(sa$Abundance, sa$abx, summary) #abx exposed kids have higher median rel abund of S. aureus
wilcox.test(Abundance ~ abx, data = sa, exact=FALSE)  #P=0.04

#Streptococcus pneumoniae
spneum <- filter(melted_df, Species == "Streptococcus pneumoniae")
hist(spneum$Abundance)  #lots of zeros
summary(spneum$Abundance)
tapply(spneum$Abundance, spneum$cd4nL, summary) #CD4<25% kids have higher median rel abund of S.pneumo
wilcox.test(Abundance ~ cd4nL, data = spneum, exact=FALSE)  #P=0.16 (warning: ties)

  #TMPSMX
tapply(spneum$Abundance, spneum$tmpsmx, summary) #tmpsmx exposed kids have higher median rel abund of S.pneumo
wilcox.test(Abundance ~ tmpsmx, data = spneum, exact=FALSE)  #P=0.03

  #abx
tapply(spneum$Abundance, spneum$abx, summary) #medians similar
wilcox.test(Abundance ~ abx, data = spneum, exact=FALSE)  #P=0.54

##################
#CORRELATIONS: CD4 
##################
cor_p <- rename(cor_p, ab.cor_p = Abundance)
cor_ps <- rename(cor_ps, ab.cor_ps = Abundance)
dpig <- rename(dpig, ab.dpig = Abundance)
hflu <- rename(hflu, ab.hflu = Abundance)
mcat <- rename(mcat, ab.mcat = Abundance)
spneum <- rename(spneum, ab.spneum = Abundance)
sa <- rename(sa, ab.sa = Abundance)

#Create simple dataframe with abundance variables, sample_id, cd4, age, and merge on sample_id
corr <- cor_ps[c("sample_id", "age", "pcv", "dpt", "ab.cor_ps", "newest_cd4perc", "cd4nL", "tmpsmx", "abx")]
test <- cor_p[c("sample_id", "ab.cor_p")]
corr <- merge(corr, test, by ="sample_id", sort = TRUE)
test <- dpig[c("sample_id", "ab.dpig")]
corr <- merge(corr, test, by ="sample_id", sort = TRUE)
test <- hflu[c("sample_id", "ab.hflu")]
corr <- merge(corr, test, by ="sample_id", sort = TRUE)
test <- mcat[c("sample_id", "ab.mcat")]
corr <- merge(corr, test, by ="sample_id", sort = TRUE)
test <- spneum[c("sample_id", "ab.spneum")]
corr <- merge(corr, test, by ="sample_id", sort = TRUE)
test <- sa[c("sample_id", "ab.sa")]
corr <- merge(corr, test, by ="sample_id", sort = TRUE)

#save as csv in case we want to come back to it
write.csv(corr, "corrtable_cd4.csv", row.names=F)

#Correlations
#source: https://stats.stackexchange.com/questions/8071/how-to-choose-between-pearson-and-spearman-correlation
#based on the above, will run both tests to look at how they compare (can help us see if the relationship is monotonic and/or linear)

#Does CD4 percentage correlate with Coryne abundance?
cor.test(corr$ab.cor_ps, corr$newest_cd4perc, method = "pearson") #P=0.26, cor = 0.180
cor.test(corr$ab.cor_ps, corr$newest_cd4perc, method = "spearman", exact = FALSE) #P=0.013, cor = 0.38
  #Yes it does!
#Plot:
library(ggtext)
png(file="/Users/swetapatel/OneDrive - Duke University/Fogarty coding/Coryne_cd4.png",
    width = 5, height = 5, units = 'in', res = 600)
ggplot(data = corr, aes(x = newest_cd4perc, y = ab.cor_ps)) + 
  geom_point() + 
  geom_smooth(method = "lm", color="red", size=0.5, se = FALSE) +
  labs(title="Species correlation with CD4%",x= "CD4 Percentage", y = "*C. pseudodiphtheriticum*")+
  geom_text(x=30, y=0.37, label="r = 0.38, P = 0.013") +
  theme_bw() +
  theme(
    plot.title = element_text(size = 14, hjust = 0.5),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12),
    legend.position = "right",
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12),
    axis.title.x = ggtext::element_markdown(),
    axis.title.y = ggtext::element_markdown())
dev.off()

#CD4 and dpig abundance
cor.test(corr$ab.dpig, corr$newest_cd4perc, method = "pearson") #P=0.04, cor = 0.324
cor.test(corr$ab.dpig, corr$newest_cd4perc, method = "spearman", exact = FALSE) #P=0.008, cor = 0.404

png(file="/Users/swetapatel/OneDrive - Duke University/Fogarty coding/Dpig_cd4.png",
    width = 5, height = 5, units = 'in', res = 600)
ggplot(data = corr, aes(x = newest_cd4perc, y = ab.dpig)) + 
  geom_point() + 
  geom_smooth(method = "lm", color="red", size=0.5, se = FALSE) +
  labs(title="Species correlation with CD4%",x= "CD4 Percentage", y = "*D. pigrum*")+
  geom_text(x=20, y=0.39, label="r = 0.404, P = 0.008") +
  theme_bw() +
  theme(
    plot.title = element_text(size = 14, hjust = 0.5),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12),
    legend.position = "right",
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12),
    axis.title.x = ggtext::element_markdown(),
    axis.title.y = ggtext::element_markdown())
dev.off()

#CD4 and S.pneumo abundance
cor.test(corr$ab.spneum, corr$newest_cd4perc, method = "pearson") #P=0.36, cor = -0.144
cor.test(corr$ab.spneum, corr$newest_cd4perc, method = "spearman", exact = FALSE) #P=0.23, cor = -0.188
  #won't plot as not signif

#CD4 and S.aureus abundance
cor.test(corr$ab.sa, corr$newest_cd4perc, method = "pearson") #P=0.17, cor = -0.218
cor.test(corr$ab.sa, corr$newest_cd4perc, method = "spearman", exact = FALSE) #P=0.009, cor = -0.397

png(file="/Users/swetapatel/OneDrive - Duke University/Fogarty coding/Saureus_cd4.png",
    width = 5, height = 5, units = 'in', res = 600)
ggplot(data = corr, aes(x = newest_cd4perc, y = ab.sa)) + 
  geom_point() + 
  geom_smooth(method = "lm", color="red", size=0.5, se = FALSE) +
  labs(title="Species correlation with CD4%",x= "CD4 Percentage", y = "*S. aureus*")+
  geom_text(x=30, y=0.45, label="r = -0.397, P = 0.009") +
  theme_bw() +
  theme(
    plot.title = element_text(size = 14, hjust = 0.5),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12),
    legend.position = "right",
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12),
    axis.title.x = ggtext::element_markdown(),
    axis.title.y = ggtext::element_markdown())
dev.off()

###CORYNE PSEUDO VS STREP PNEUMO
cor.test(corr$ab.cor_ps, corr$ab.spneum, method = "pearson")  #P=0.46, cor = -0.118
cor.test(corr$ab.cor_ps, corr$ab.spneum, method = "spearman", exact = FALSE)  #P=0.39, estimate = -0.136

#need to plot this out and see what it looks like (overall and then by HIV status)
#source: https://bookdown.org/dli/rguide/scatterplots-and-best-fit-lines-two-sets.html#two-scatterplots-using-ggplot2

#Correlations by CD4 status:
corr_lowcd4 <- subset(corr, cd4nL == 0)
corr_nlcd4 <- subset(corr, cd4nL == 1)

cor.test(corr_lowcd4$ab.cor_ps, corr_lowcd4$ab.spneum, method = "pearson")  #P=0.66, cor = -0.148
cor.test(corr_lowcd4$ab.cor_ps, corr_lowcd4$ab.spneum, method = "spearman", exact = FALSE) #P=0.42, cor = 0.27

cor.test(corr_nlcd4$ab.cor_ps, corr_nlcd4$ab.spneum, method = "pearson")  #P=0.49, cor = -0.128
cor.test(corr_nlcd4$ab.cor_ps, corr_nlcd4$ab.spneum, method = "spearman", exact = FALSE) #P=0.47, cor = -0.136

#Plot out
#C. pseudodiphtheriticum vs S. pneumo plot: all participants together
#Source for adding text: http://www.sthda.com/english/wiki/ggplot2-texts-add-text-annotations-to-a-graph-in-r-software
library(ggtext)
# png(file="/Users/swetapatel/OneDrive - Duke University/Fogarty coding/Coryne_SP_all.png",
#     width = 5, height = 5, units = 'in', res = 600)
ggplot(data = corr, aes(x = ab.cor_ps, y = ab.spneum)) + 
  geom_point() + 
  geom_smooth(method = "lm", color="red", size=0.5, se = FALSE) +
  labs(title="Interspecies correlations",x="*C. pseudodiphtheriticum*", y = "*S. pneumoniae*")+
  geom_text(x=0.3, y=0.55, label="r = - 0.136, P = 0.39") +
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
# dev.off()

#What about 2 plots, one for each CD4 category? (Facet wrap)
# png(file="/Users/swetapatel/OneDrive - Duke University/Fogarty coding/Coryne_SP_facet.png",
#     width = 5, height = 5, units = 'in', res = 600)
corr$hiv <- as.factor(corr$hiv)
corr$hiv <- reorder(corr$hiv, new.order=c("Unexposed", "Exposed Uninfected", "Infected"))
ggplot(data = corr, aes(x = ab.cor_ps, y = ab.spneum)) + 
  geom_point() + 
  geom_smooth(method = "lm", color="red", size=0.5, se = FALSE) +
  facet_grid(~cd4nL) +
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
# dev.off()
#manually annotate for now after export 

###CORYNE PSEUDO VS STAPH AUREUS
  #going forward, will look at correlations and plot + corr by CD4 status if significant
cor.test(corr$ab.cor_ps, corr$ab.sa, method = "pearson")  #P=0.39, cor = -0.137
cor.test(corr$ab.cor_ps, corr$ab.sa, method = "spearman", exact = FALSE)  #P=0.01, estimate = -0.394

#signif; what does it look like by CD4 status?
cor.test(corr_lowcd4$ab.cor_ps, corr_lowcd4$ab.sa, method = "pearson")  #P=0.59, cor = -0.182
cor.test(corr_lowcd4$ab.cor_ps, corr_lowcd4$ab.sa, method = "spearman", exact = FALSE) #P=0.052, cor = -0.596

cor.test(corr_nlcd4$ab.cor_ps, corr_nlcd4$ab.sa, method = "pearson")  #P=0.49, cor = -0.128
cor.test(corr_nlcd4$ab.cor_ps, corr_nlcd4$ab.sa, method = "spearman", exact = FALSE) #P=0.313, cor = -0.187

#Plot out
#C. pseudodiphtheriticum vs S. aureus plot
#Source for adding text: http://www.sthda.com/english/wiki/ggplot2-texts-add-text-annotations-to-a-graph-in-r-software
png(file="/Users/swetapatel/OneDrive - Duke University/Fogarty coding/Coryne_SA_cd4.png",
    width = 5, height = 5, units = 'in', res = 600)
ggplot(data = corr, aes(x = ab.cor_ps, y = ab.sa)) + 
  geom_point() + 
  geom_smooth(method = "lm", color="red", size=0.5, se = FALSE) +
  labs(title="Interspecies correlations",x="*C. pseudodiphtheriticum*", y = "*S. aureus*")+
  geom_text(x=0.3, y=0.55, label="r = - 0.394, P = 0.01") +
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

#What about 2 plots, one for each CD4 category? (Facet wrap)
png(file="/Users/swetapatel/OneDrive - Duke University/Fogarty coding/Coryne_SAcd4_facet.png",
    width = 5, height = 5, units = 'in', res = 600)
ggplot(data = corr, aes(x = ab.cor_ps, y = ab.sa)) + 
  geom_point() + 
  geom_smooth(method = "lm", color="red", size=0.5, se = FALSE) +
  facet_grid(~cd4nL) +
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

#CORYNE PSEUDO VS M. CAT
cor.test(corr$ab.cor_ps, corr$ab.mcat, method = "pearson")  #P=0.33, cor = -0.153
cor.test(corr$ab.cor_ps, corr$ab.mcat, method = "spearman", exact = FALSE)  #P=0.60, estimate = -0.083

#CORYNE PSEUDO VS H. FLU
cor.test(corr$ab.cor_ps, corr$ab.hflu, method = "pearson")  #P=0.23, cor = -0.189
cor.test(corr$ab.cor_ps, corr$ab.hflu, method = "spearman", exact = FALSE)  #P=0.15, estimate = -0.224

###############
#MAASLIN 2: CD4
###############
cd4.filter <- filter_taxa(cd4, function(Abundance) mean(Abundance)>=50, TRUE)
ntaxa(cd4.filter) #45

otu <- as.data.frame(otu_table(cd4.filter, taxa_are_rows=TRUE))
#need to remove the unidentified species from here: row 10, 25, 29, 41
#source: https://stackoverflow.com/questions/12328056/how-do-i-delete-rows-in-a-data-frame
#corresponding respectively to: Corynebacterium sp. KPL1859, Streptococcus sp. SK643, Paenibacillus sp. P22, Prevotellaceae bacterium Marseille-P2826
otu <- otu[-c(10, 25, 29, 41),]
#for some reason code to turn chi_cd4 from sample_data to data.frame is not working
chi_cd4 <- as.data.frame(chi_cd4)
  #recreate dataframe from cd4 pseq object
chi_cd4 = data.frame(sample_data(cd4))

#Univariable with CD4 (binary and continuous)
fit_data2 = Maaslin2(
  input_data = otu, 
  input_metadata = chi_cd4, 
  min_prevalence = 0.2,
  max_significance = 0.2,
  output = "maaslin_uv_kaiju_cd4binary", 
  fixed_effects = c("cd4nL"))

fit_data2 = Maaslin2(
  input_data = otu, 
  input_metadata = chi_cd4, 
  min_prevalence = 0.2,
  max_significance = 0.2,
  output = "maaslin_uv_kaiju_cd4cont", 
  fixed_effects = c("newest_cd4perc"))

#Multivariable incorporating age, TMP-SMX use
fit_data2 = Maaslin2(
  input_data = otu, 
  input_metadata = chi_cd4, 
  min_prevalence = 0.2,
  max_significance = 0.2,
  output = "maaslin_mv_kaiju_cd4", 
  fixed_effects = c("cd4nL", "tmpsmx", "age"))

#we saw differences by CD4 status in terms of abx and recent_uri; what does model look like with these?
fit_data2 = Maaslin2(
  input_data = otu, 
  input_metadata = chi_cd4, 
  min_prevalence = 0.2,
  max_significance = 0.2,
  output = "maaslin_mv2_kaiju_cd4", 
  fixed_effects = c("cd4nL", "abx", "uri_rec2"))
#No significant results

#What about with CD4 as a continuous variable?
fit_data2 = Maaslin2(
  input_data = otu, 
  input_metadata = chi_cd4, 
  min_prevalence = 0.2,
  max_significance = 0.2,
  output = "maaslin_mv3_kaiju_cd4", 
  fixed_effects = c("newest_cd4perc", "abx", "uri_rec2"))
#no significant results

fit_data2 = Maaslin2(
  input_data = otu, 
  input_metadata = chi_cd4, 
  min_prevalence = 0.2,
  max_significance = 0.2,
  output = "maaslin_uv_kaiju_cd4_abx", 
  fixed_effects = c("abx"))
#No significant results

fit_data2 = Maaslin2(
  input_data = otu, 
  input_metadata = chi_cd4, 
  min_prevalence = 0.2,
  max_significance = 0.2,
  output = "maaslin_uv_kaiju_cd4_urirec", 
  fixed_effects = c("uri_rec2"))
#URI recently also a/w lower abundance of Coryne pseudo

#what about uri_rec2 + cd4 together?
fit_data2 = Maaslin2(
  input_data = otu, 
  input_metadata = chi_cd4, 
  min_prevalence = 0.2,
  max_significance = 0.2,
  output = "maaslin_mv4_kaiju_cd4", 
  fixed_effects = c("cd4nL", "uri_rec2"))
#Moraxella a/w uri; seems like cd4 and uri cancel each other out...

fit_data2 = Maaslin2(
  input_data = otu, 
  input_metadata = chi_cd4, 
  min_prevalence = 0.2,
  max_significance = 0.2,
  output = "maaslin_mv5_kaiju_cd4", 
  fixed_effects = c("newest_cd4perc", "uri_rec2"))
#CD4 continuous a/w Coryne propinquum, uri_rec a/w Moraxella

#clean up environment to redo relative abundance analyses by viral load
remove(chi.rel, abundances, cor_ps, corr, corr_lowcd4, corr_nlcd4, dpig, fit_data2, genus_abundances, genus_df,
       hflu, mcat, melted_df, otu, phyla_abundances, phylum_df, relative_df, sa, species_abundances, 
       species_df, spneum, subject_plot2, test)

#####################################
#RELATIVE ABUNDANCE PLOTS: VIRAL LOAD
#####################################
#For completeness, we will generate relative abundance plots for viral load too
chi.rel <- transform_sample_counts(vl, function(Abundance) Abundance/sum(Abundance))
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
#merging on Phylum AND Genus gets rid of the phylum.x and phylum.y problem we were having earlier

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
relative_df <- merge(melted_df, abundances, by=c("Genus", "Species"))

relative_df$Species[!(relative_df$Species %in% TOPSpecies)] <- "Other"
sum(relative_df$Abundance) #sums to 41.0003 (and we have n=41 kids...?)
table(relative_df$Species)

#To generate color palette of 11 colors used previously: 
library(RColorBrewer)
library(scales)
chisp2 <- c("#293352", "#800000FF", "#767676FF", "#FFA319FF", "#8A9045FF", "#155F83FF", "#C16622FF", "#4D004B", "#8F3931FF", "#58593FFF", "#350E20FF")
show_col(chisp2)

#CREATING RELATIVE ABUNDANCE PLOT
relative_df$subject_abundances <- (relative_df$Abundance)/(length(unique(relative_df$Sample))) 
#Clarify VL cutoff 
relative_df$vl_suppr2 <- NA
relative_df$vl_suppr2 [relative_df$vl_suppr == 0] <- "No"
relative_df$vl_suppr2 [relative_df$vl_suppr == 1] <- "Yes"
table(relative_df$vl_suppr2)
table(relative_df$vl_suppr)
#matches up with table(relative_df$vl_suppr)
# relative_df$hiv2 <- as.factor(relative_df$hiv2)
# relative_df$hiv2 <- reorder(relative_df$hiv2, new.order=c("HUU", "HEU", "HEI"))
relative_df$Species <- factor(relative_df$Species)
relative_df$Species <- reorder(relative_df$Species, new.order=c("Corynebacterium propinquum", "Corynebacterium pseudodiphtheriticum", "Dolosigranulum pigrum",  
                                                                "Haemophilus influenzae", "Micrococcus luteus", "Moraxella catarrhalis", "Moraxella lincolnii", 
                                                                "Moraxella nonliquefaciens", "Staphylococcus aureus", "Streptococcus pneumoniae", "Other"))
subject_plot2 <- ggplot(relative_df[order(relative_df$Species, decreasing = TRUE),], aes(x=vl_suppr2, y=subject_abundances, fill=Species)) + 
  geom_bar(stat="identity", position="fill") + theme_classic() + theme(legend.position="right", 
                                                                       legend.text=element_text(size=12, face="italic"),
                                                                       axis.title.y = element_text(size=14, margin=margin(0,20,0,0)), axis.text.y = element_text(size=12), 
                                                                       axis.title.x = element_text(size=14, margin=margin(10,0,0,0)), axis.text.x = element_text(size=12)) + 
  scale_fill_manual(values = chisp2) +
  # scale_x_discrete(labels=c("No", "Yes")) +
  xlab("Viral Suppression") + ylab("Relative Abundance") 
png(file="vLplot_kaiju_06082022.png",
    width = 9, height = 5, units = 'in', res = 600)
plot(subject_plot2)
dev.off()

#########################################
#RELATIVE ABUNDANCE M-W TESTS: VIRAL LOAD
#########################################
#Corynebacterium propinquum
cor_p <- filter(melted_df, Species == "Corynebacterium propinquum")
summary(cor_p$Abundance)
hist(cor_p$Abundance)
tapply(cor_p$Abundance, cor_p$vl_suppr, summary)
wilcox.test(Abundance ~ vl_suppr, data = cor_p)  #P=0.61 (warning: ties)

#Corynebacterium pseudodiphtheriticum
cor_ps <- filter(melted_df, Species == "Corynebacterium pseudodiphtheriticum")
summary(cor_ps$Abundance)
hist(cor_ps$Abundance)
tapply(cor_ps$Abundance, cor_ps$vl_suppr, summary) 
wilcox.test(Abundance ~ vl_suppr, data = cor_ps)  #P=0.29 (warning: ties)

#Dolosigranulum pigrum
dpig <- filter(melted_df, Species == "Dolosigranulum pigrum")
hist(dpig$Abundance)
summary(dpig$Abundance)
tapply(dpig$Abundance, dpig$vl_suppr, summary) 
wilcox.test(Abundance ~ vl_suppr, data = dpig)  #P=0.40 (warning: ties)

#Haemophilus influenzae
hflu <- filter(melted_df, Species == "Haemophilus influenzae")
hist(hflu$Abundance)  #lots of zeros
summary(hflu$Abundance)
tapply(hflu$Abundance, hflu$vl_suppr, summary)
wilcox.test(Abundance ~ vl_suppr, data = hflu)  #P=0.72 (warning: ties)

#Micrococcus luteus
mic <- filter(melted_df, Species == "Micrococcus luteus")
hist(mic$Abundance)  #lots of zeros
summary(mic$Abundance)
tapply(mic$Abundance, mic$vl_suppr, summary)
wilcox.test(Abundance ~ vl_suppr, data = mic)  #P=0.77 (warning: ties)

#Moraxella cattarhalis
mcat <- filter(melted_df, Species == "Moraxella catarrhalis")
hist(mcat$Abundance)  #lots of zeros
summary(mcat$Abundance)
tapply(mcat$Abundance, mcat$vl_suppr, summary)
wilcox.test(Abundance ~ vl_suppr, data = mcat)  #P=0.46 (warning: ties)

#Moraxella lincolnii
mlin <- filter(melted_df, Species == "Moraxella lincolnii")
hist(mlin$Abundance)  #lots of zeros
summary(mlin$Abundance)
tapply(mlin$Abundance, mlin$vl_suppr, summary)
wilcox.test(Abundance ~ vl_suppr, data = mlin)  #P=0.31 (warning: ties)

#Moraxella nonliquefaciens
mnon <- filter(melted_df, Species == "Moraxella nonliquefaciens")
hist(mnon$Abundance)  #lots of zeros but more evenly distrib than other Moraxellas
summary(mnon$Abundance)
tapply(mnon$Abundance, mnon$vl_suppr, summary)
wilcox.test(Abundance ~ vl_suppr, data = mnon)  #P=0.53 (warning: ties)

#Staphylococcus aureus
sa <- filter(melted_df, Species == "Staphylococcus aureus")
hist(sa$Abundance)  #lots of zeros
summary(sa$Abundance)
tapply(sa$Abundance, sa$vl_suppr, summary) #Medians of both groups are 0
wilcox.test(Abundance ~ vl_suppr, data = sa)  #P=0.08 (warning: ties)

#Streptococcus pneumoniae
spneum <- filter(melted_df, Species == "Streptococcus pneumoniae")
hist(spneum$Abundance)  #lots of zeros
summary(spneum$Abundance)
tapply(spneum$Abundance, spneum$vl_suppr, summary)
wilcox.test(Abundance ~ vl_suppr, data = spneum)  #P=0.96 (warning: ties)

######################
#MAASLIN 2: VIRAL LOAD
######################
vl.filter <- filter_taxa(vl, function(Abundance) mean(Abundance)>=50, TRUE)
ntaxa(vl.filter) #46

otu <- as.data.frame(otu_table(vl.filter, taxa_are_rows=TRUE))
#need to remove the unidentified species from here: row 11, 26, 30, 42
#source: https://stackoverflow.com/questions/12328056/how-do-i-delete-rows-in-a-data-frame
#corresponding respectively to: Corynebacterium sp. KPL1859, Streptococcus sp. SK643, Paenibacillus sp. P22, Prevotellaceae bacterium Marseille-P2826
otu <- otu[-c(11, 26, 30, 42),] #42 species left
#turn chi_viral from sample data into dataframe (by pulling it out of pseq object again)
chi_viral = data.frame(sample_data(vl))

#univariable binary
fit_data2 = Maaslin2(
  input_data = otu, 
  input_metadata = chi_viral, 
  min_prevalence = 0.2,
  max_significance = 0.2,
  output = "maaslin_uv_kaiju_vl_binary", 
  fixed_effects = c("vl_suppr"))
#no signficant results

#univariable with VL as continuous variable
fit_data2 = Maaslin2(
  input_data = otu, 
  input_metadata = chi_viral, 
  min_prevalence = 0.2,
  max_significance = 0.2,
  output = "maaslin_uv_kaiju_vl_cont", 
  fixed_effects = c("newest_vlnum"))
#vL a/w Haemophilus C1 and parahaemolyticus and Neisseria lactamica

#multivariable with binary vl
fit_data2 = Maaslin2(
  input_data = otu, 
  input_metadata = chi_viral, 
  min_prevalence = 0.2,
  max_significance = 0.2,
  output = "maaslin_mv_kaiju_vl", 
  fixed_effects = c("vl_suppr", "tmpsmx", "age"))
#no significant results

#multivariable with continuous vl
fit_data2 = Maaslin2(
  input_data = otu, 
  input_metadata = chi_viral, 
  min_prevalence = 0.2,
  max_significance = 0.2,
  output = "maaslin_mv_kaiju_vl_cont", 
  fixed_effects = c("newest_vlnum", "tmpsmx", "age"))
#viral load a/w increasing H. parahaemolyticus, TMP-SMX a/w increased Staph aureus, Age a/w M. lincolnii

#################################
#HIV-HEU SIBLING PAIR COMPARISONS
#################################
#1. Compare age, sex, and clinical characteristics using paired Wilcoxon
#2. Compare overall composition using age-adjusted permanova (?pairwise adonis)
#3. Compare rel abundance of C. pseudo using paired Wilcoxon

#we want to compare pairwise age, sex, season, abx, recent URI, PCV, and HIB
#pull these variables out separately for children and the siblings, standardize variable names, merge, and then do a test
sib2 <- sib[c("sample_id", "study_id", "age", "sex", "season", "subject", "abx", "uri_rec2", "pcv", "dpt", "sib_age", "sib_sex",
              "sib_uri_recent", "sib_meds", "sib_meds_name", "sib_dpt", "sib_pcv", "sibmo")]

#FIRST: create new variables for sibling abx, season, and uri_rec2
#abx: cleaned up sib med names in separate metadata file --> "sib_abx 8-10-22"csv file
#GO BACK AND CLEAN IN REDCAP ONCE ACCESS RE-ESTABLISHED
abx <- read.csv("~/Library/CloudStorage/OneDrive-DukeUniversity/Fogarty coding/sib_abx 08-10-22.csv")
#need to add NAs? 
#drop current med_names variable from metadata and then merge with abx df by study_id
#standardize study_id
toString(abx$study_id)
abx$study_id <- gsub("-", ".", abx$study_id)
sib2$sib_meds_name <- NULL
sib2$sib_meds <- NULL
sib2 <- merge(sib2, abx, by ="study_id", sort = TRUE)
table(sib2$sib_meds_name)
sib2$sib_abx <- "No"
sib2$sib_abx [sib2$sib_meds_name == "amoxicillin"] <- "Yes"
table(sib2$sib_abx)

#season
#For consistency with MK: Rainy season = Nov to March. Dry season = April to Oct
sib2$sib_seas <- "99"
sib2$sib_seas [sib2$sibmo == 1 | sib2$sibmo == 2 | sib2$sibmo == 3 | sib2$sibmo == 11 | sib2$sibmo == 12] <- "Rainy"
sib2$sib_seas [sib2$sibmo == 4 | sib2$sibmo == 5 | sib2$sibmo == 6 | sib2$sibmo == 7 | sib2$sibmo == 8 | sib2$sibmo == 9 | sib2$sibmo == 10] <- "Dry"
table(sib2$sib_seas)

#recent uri
table(sib2$sib_uri_recent)
sib2$sib_uri_rec2 <- NA
sib2$sib_uri_rec2 [sib2$sib_uri_recent == 1] <- "Yes"
sib2$sib_uri_rec2 [sib2$sib_uri_recent == 2] <- "No"
table(sib2$sib_uri_rec2)

#NEXT: pull the variables out separately, rename the sibling variables so standardized with chi vars, then merge back together
sib2 <- within(sib2, rm(sib_uri_recent, sibmo, sib_meds, sib_meds_name))
sib_chi <- subset(sib2, sib2$subject == "child")
sib_chi <- sib_chi[c("sample_id", "study_id", "age", "sex", "season", "subject", "abx", "uri_rec2", "pcv", "dpt")]
sib_sib <- subset(sib2, sib2$subject == "sib")
sib_sib <- sib_sib[c("sample_id", "study_id", "subject", "sib_age", "sib_sex",
              "sib_uri_rec2", "sib_abx", "sib_dpt", "sib_pcv", "sib_seas")]
#rename sib variables
sib_sib <- rename(sib_sib, age = "sib_age")
sib_sib <- rename(sib_sib, sex = "sib_sex")
sib_sib <- rename(sib_sib, uri_rec2 = "sib_uri_rec2")
sib_sib <- rename(sib_sib, abx = "sib_abx")
sib_sib <- rename(sib_sib, dpt = "sib_dpt")
sib_sib <- rename(sib_sib, pcv = "sib_pcv")
sib_sib <- rename(sib_sib, season = "sib_seas")

#combine: use rbind? 
sib3 <- rbind(sib_chi, sib_sib)
#save as csv for easier reference
write.csv(sib3, 'sibpairs.csv')
sib3 <- read.csv("sibpairs.csv")

#now compare demographic characteristics by group AND pairwise (?)
  #paired comparisons for numeric data: wilcoxon paired
  #paired comparisons for 2x2 tables: McNemar test: https://rpubs.com/kaz_yos/mcnemar
library(exact2x2)
#Age
hist(sib3$age)
summary(sib3$age)
tapply(sib3$age, sib3$subject, summary)
wilcox.test(age ~ subject, data = sib3, exact=FALSE) #P=0.57
  #paired
wilcox.test(age ~ subject, data=sib3, paired=TRUE, exact=FALSE) #P=0.59

#Sex
table(sib3$sex)  #11 girls
table(sib3$sex, sib3$subject)
prop.table(table(sib3$sex, sib3$subject), 2)
fisher.test(sib3$sex, sib3$subject) #P=1.0
  #paired
mcnemar.test(sib3$sex, sib3$subject)  #P=0.72

#Season
table(sib3$season)  
table(sib3$season, sib3$subject)
prop.table(table(sib3$season, sib3$subject), 2)
fisher.test(sib3$season, sib3$subject)  #P=1.0
  #paired
mcnemar.test(sib3$season, sib3$subject) #P=1.0

#antibiotics
table(sib3$abx)  #3 kids received abx
table(sib3$abx, sib3$subject)
prop.table(table(sib3$abx, sib3$subject), 2)
fisher.test(sib3$abx, sib3$subject) #P=1.0
  #paired
mcnemar.test(sib3$abx, sib3$subject) #P=0.08

#URI in past 1 month
table(sib3$uri_rec2)  #11 kids with recent URI
table(sib3$uri_rec2, sib3$subject)
prop.table(table(sib3$uri_rec2, sib3$subject), 2)
fisher.test(sib3$uri_rec2, sib3$subject)
  #paired
mcnemar.test(sib3$uri_rec2, sib3$subject) #P=0.75

#PCV doses
table(sib3$pcv)  
table(sib3$pcv, sib3$subject)
prop.table(table(sib3$pcv, sib3$subject), 2)
fisher.test(sib3$pcv, sib3$subject) #P=1.0
  #paired
wilcox.test(pcv ~ subject, data=sib3, paired=TRUE, exact=FALSE) #P=0.85

#HiB
table(sib3$dpt)  
table(sib3$dpt, sib3$subject)
prop.table(table(sib3$dpt, sib3$subject), 2)
fisher.test(sib3$dpt, sib3$subject) #P=1.0
  #paired
wilcox.test(dpt ~ subject, data=sib3, paired=TRUE, exact=FALSE) #P=0.85

sib_clr <- transform(sib_pair, "clr")

#PCoA plots#
clr_pcoa <- ordinate(sib_clr, method = "PCoA", distance = "euclidean")
lines <- c("dashed", "dashed")
clr_sib <- plot_ordination(sib_clr, clr_pcoa, color = "subject", title = "Euclidean PCoA Plot by Subject Classification") +
  geom_point(size = 2) + theme_classic() + theme(
    plot.title = element_text(size = 14, hjust = 0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 10),
    legend.position = "right",
    legend.background = element_blank(),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12)) +
  scale_color_manual(values=c('#800000FF','#155F83FF'))+
  stat_ellipse(aes(linetype=subject), geom="polygon", alpha=0, type="t", level=0.8, size=0.7) +
  xlab("PC1 (18.1%)") + ylab("PC2 (17.7%)") +
  scale_linetype_manual(values=lines)
png(file="clr_sib.png",
    width = 5, height = 4, units = 'in', res=600)
clr_sib
dev.off()

#What if we try to do it by pairs?
#Set same color palate first (but with 9 colors)
library(RColorBrewer)
library(scales)
chisp3 <- c("#800000FF", "#767676FF", "#FFA319FF", "#8A9045FF", "#155F83FF", "#C16622FF", "#4D004B", "#8F3931FF", "#58593FFF")
show_col(chisp3)
# lines <- c("dashed", "dashed", "dashed", "dashed", "dashed", "dashed", "dashed", "dashed", "dashed")
  #too few points to make an ellipse around each pair

clr_sib2 <- plot_ordination(sib_clr, clr_pcoa, color = "study_id") +
  geom_point(size = 2) + scale_color_manual(values = chisp3) + theme_classic() + theme(
    plot.title = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 10),
    legend.position = "right",
    legend.background = element_blank(),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12)) + 
  xlab("PC1 (18.1%)") + ylab("PC2 (17.7%)") 
  # scale_color_manual(values=c('#800000FF','#155F83FF')) +
  # stat_ellipse(aes(linetype=subject), geom="polygon", alpha=0, type="t", level=0.8, size=0.7) +
  # xlab("PC1 (18.1%)") + ylab("PC2 (17.7%)") +
  # scale_linetype_manual(values=lines)

png(file="clr_sibs.png",
    width = 6, height = 4, units = 'in', res=600)
clr_sib2
dev.off()

#Follow steps like standard permanova
otu <- as.data.frame(otu_table(sib_clr, taxa_are_rows=TRUE))
#Need species in columns and observations in rows: transpose using t() function and save as a dataframe
#Source: https://github.com/edamame-course/2015-tutorials/blob/master/final/Ashley_Intro_to_R.md
otu_trans <- data.frame(t(otu))
BC.dist=vegdist(otu_trans, method="euclidean")

#below options did not work, so trying different code
  #source: https://rdrr.io/github/Jtrachsel/funfuns/man/pairwise.adonis.html
#turn our dataset into a vector?
#we need some way to differentiate the pairs so we end up with 9 tests
  #did not work: using subject (1 test), study_id (tests all IDs against all other IDs), sample_id(error)
#what if we make a new variable called pair? and assign values 1-9 and see if that works
sib3$pair <- NA
sib3$pair [sib3$study_id == "B.09"] <- 9
sib3$pair [sib3$study_id == "B.16"] <- 16
sib3$pair [sib3$study_id == "B.17"] <- 17
sib3$pair [sib3$study_id == "B.20"] <- 20
sib3$pair [sib3$study_id == "B.23"] <- 23
sib3$pair [sib3$study_id == "B.25"] <- 25
sib3$pair [sib3$study_id == "B.32"] <- 32
sib3$pair [sib3$study_id == "B.35"] <- 35
sib3$pair [sib3$study_id == "B.48"] <- 48


test <- sib3[['pair']]
pairwise.adonis(
  otu_trans,
  test,
  sim.method = "euclidean",
  p.adjust.m = "BH")

#Source for pairwise adonis: https://github.com/pmartinezarbizu/pairwiseAdonis
#let's give it a shot
# library(devtools)
# install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)
#pairwise adonis
# pairwise.adonis(BC.dist ~ study_id, data = sib) #error message: unused arguments (data=sib)
# pairwise.adonis(BC.dist, sib$subject) #1 result of child vs sib; not what we want
# pairwise.adonis2(BC.dist, sib3$pair, strata = "study_id") 
# pairwise.adonis2(BC.dist ~ subject, data = sib3, strata = 'pair')
pairwise.adonis2(BC.dist ~ subject, data = sib3, strata = 'study_id')
#adjust for FDR with BH method
pairwise.adonis2(BC.dist ~ subject, data = sib3, strata = 'study_id', p.adjust.m = "BH")
  #same result, no error message.
  ###THIS CODE (ABOVE) WORKS
adonis2(BC.dist ~ subject, data = sib3, permutations = 999) #no difference between children and siblings overall (NOT pairwise)

# #try with strata
# pairwise.adonis(BC.dist ~ subject, data = sib, strata = "study_id")
# #not working
# pairwise.adonis(BC.dist, sib$study_id)

# #try different package
# install.packages("remotes")
# remotes::install_github("ECGen/ComGenR")  #can't get the library to load
# pair.permanova(BC.dist, sib$study_id, nits = 999)
# 
# pairwiseAdonis::pairwise.adonis2(BC.dist ~ subject, data = sib, strata = "study_id", nperm = 999)
# pairwiseAdonis::pairwise.adonis2(otu_trans, sib, sim.function = "vegdist", sim.method = "euclidean",
#                                 p.adjust.m = "bonferroni", reduce = NULL, perm=999)

#I'm hitting a dead end; skip for now and move on to comparing rel abundance of C.pseudo
sib.rel <- transform_sample_counts(sib_pair, function(Abundance) Abundance/sum(Abundance))
sample_sums(sib.rel)  # This is a sanity check to make sure that relative abundance was calculated for each sample prior to pooling data
#All samples add to 1
melted_df <- psmelt(sib.rel)

#UPDATE 8/30/22: for completeness and consistency, are there differences in any other of the top 10 species?
# Create dataframes with overall relative abundances of phyla and genera
melted_df$Phylum <- as.character(melted_df$Phylum)
phyla_abundances <- aggregate(melted_df$Abundance, by=list(Phylum=melted_df$Phylum), FUN=sum)
phyla_abundances$x <- (phyla_abundances$x)/(nsamples(sib.rel))
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

relative_df$Species[!(relative_df$Species %in% TOPSpecies)] <- "Other"
sum(relative_df$Abundance)#sums to 18; n=18 so tracks?
table(relative_df$Species) #We have 11 levels for the Species variable, which are consistent with our top 10 + other

#Create figure of most common species by subject (child vs sib)
#To generate color palette of 11 colors used previously: 
library(RColorBrewer)
library(scales)
chisp2 <- c("#293352", "#800000FF", "#767676FF", "#FFA319FF", "#8A9045FF", "#155F83FF", "#C16622FF", "#4D004B", "#8F3931FF", "#58593FFF", "#350E20FF")
show_col(chisp2)

#Creating relative abundance figure: siblings
#CREATING RELATIVE ABUNDANCE PLOT
#source for rotating x labels: https://stackoverflow.com/questions/1330989/rotating-and-spacing-axis-labels-in-ggplot2
relative_df$subject_abundances <- (relative_df$Abundance)/(length(unique(relative_df$Sample))) 
relative_df$Species <- factor(relative_df$Species)
relative_df$Species <- reorder(relative_df$Species, new.order=c("Corynebacterium propinquum", "Corynebacterium pseudodiphtheriticum", "Dolosigranulum pigrum",  
                                                                "Haemophilus influenzae", "Moraxella catarrhalis", "Moraxella lincolnii", "Moraxella nonliquefaciens", 
                                                                "Staphylococcus aureus", "Streptococcus mitis", "Streptococcus pneumoniae", "Other"))
subject_plot2 <- ggplot(relative_df[order(relative_df$Species, decreasing = TRUE),], aes(x=subject, y=subject_abundances, fill=Species)) + 
  geom_bar(stat="identity", position="fill") + theme_classic() + theme(legend.position="right", 
                                                                       legend.text=element_text(size=10, face="italic"),
                                                                       axis.title.y = element_text(size=12, margin=margin(0,20,0,0)), axis.text.y = element_text(size=10), 
                                                                       axis.title.x = element_blank(), axis.text.x = element_text(size=10, angle = 45, vjust = 1, hjust = 1, colour = "black")) + 
  scale_fill_manual(values = chisp2) +
  # scale_x_discrete(labels=c("No", "Yes")) +
  xlab("Subject classification") + ylab("Mean relative abundance") 

#element_text(size=12, margin=margin(10,0,0,0))

sib2 <- plot(subject_plot2)

png(file="sib_genplot_08302022.png",
    width = 9, height = 5, units = 'in', res = 600)
sib2
dev.off()

#Coryne pseudodiphtheriticum
cor_ps <- filter(melted_df, Species == "Corynebacterium pseudodiphtheriticum")
summary(cor_ps$Abundance)
hist(cor_ps$Abundance)
tapply(cor_ps$Abundance, cor_ps$subject, summary) #HIV+ kids have lower med rel abundance than sibs
wilcox.test(Abundance ~ subject, data=cor_ps, paired=TRUE) #P=0.008 BUT not sure this is the correct code
#Try again with only study_id, abundance, and subject in a dataframe and see what results are
test <- cor_ps[c("Sample", "study_id", "Abundance", "subject")]
wilcox.test(Abundance ~ subject, data=test, paired=TRUE)  #same result

#D. pigrum
dpig <- filter(melted_df, Species == "Dolosigranulum pigrum")
summary(dpig$Abundance)
hist(dpig$Abundance)
tapply(dpig$Abundance, dpig$subject, summary) #HIV+ kids have lower med rel abundance than sibs
wilcox.test(Abundance ~ subject, data=dpig, paired=TRUE) #P=0.03

#REMAINING TOP 10 GENERA
#Corynebacterium propinquum
cor_p <- filter(melted_df, Species == "Corynebacterium propinquum")
summary(cor_p$Abundance)
hist(cor_p$Abundance)
tapply(cor_p$Abundance, cor_p$subject, summary)
wilcox.test(Abundance ~ subject, data=cor_p, paired=TRUE) #P=0.008

#Haemophilus influenzae
hflu <- filter(melted_df, Species == "Haemophilus influenzae")
hist(hflu$Abundance)  #lots of zeros
summary(hflu$Abundance)
tapply(hflu$Abundance, hflu$subject, summary) #limited by zeros (all medians = 0)
wilcox.test(Abundance ~ subject, data=hflu, paired=TRUE)  #P=0.06 (warning: cannot compute exact p-value with zeroes)

#Moraxella cattarhalis
mcat <- filter(melted_df, Species == "Moraxella catarrhalis")
hist(mcat$Abundance)  #lots of zeros
summary(mcat$Abundance)
tapply(mcat$Abundance, mcat$subject, summary) 
wilcox.test(Abundance ~ subject, data=mcat, paired=TRUE) #P=0.29 (warning: cannot compute exact p-value with zeroes)

#Moraxella lincolnii
mlin <- filter(melted_df, Species == "Moraxella lincolnii")
hist(mlin$Abundance)  #lots of zeros
summary(mlin$Abundance)
tapply(mlin$Abundance, mlin$subject, summary) #More abundant in children with HIV
wilcox.test(Abundance ~ subject, data=mlin, paired=TRUE)  #P=0.004

#Moraxella nonliquefaciens
mnon <- filter(melted_df, Species == "Moraxella nonliquefaciens")
hist(mnon$Abundance)  #lots of zeros but more evenly distrib than other Moraxellas
summary(mnon$Abundance)
tapply(mnon$Abundance, mnon$subject, summary)
wilcox.test(Abundance ~ subject, data=mnon, paired=TRUE) #P=1.0 (warning: cannot compute exact p-value with zeroes)

#Staphylococcus aureus
sa <- filter(melted_df, Species == "Staphylococcus aureus")
hist(sa$Abundance)  #lots of zeros
summary(sa$Abundance)
tapply(sa$Abundance, sa$subject, summary) #limited by zeros (all medians = 0)
wilcox.test(Abundance ~ subject, data=sa, paired=TRUE) #P=0.06 (warning: cannot compute exact p-value with zeroes)

#Streptococcus mitis
sm <- filter(melted_df, Species == "Streptococcus mitis")
hist(sm$Abundance)  #lots of zeros
summary(sm$Abundance)
tapply(sm$Abundance, sm$subject, summary)
wilcox.test(Abundance ~ subject, data=sm, paired=TRUE) #P=0.55 (warning: cannot compute exact p-value with zeroes)

#Streptococcus pneumoniae
spneum <- filter(melted_df, Species == "Streptococcus pneumoniae")
hist(spneum$Abundance)  #lots of zeros
summary(spneum$Abundance)
tapply(spneum$Abundance, spneum$subject, summary)
wilcox.test(Abundance ~ subject, data=spneum, paired=TRUE) #P=0.02 (warning: cannot compute exact p-value with zeroes)

##############################################
#ANALYSIS REMOVING HEI CHILDREN WITH CD4 < 25%
##############################################
#1. Compare overall composition using age-adjusted permanova
#2. Rel abundance plots
#3. Maaslin2

#now need to remove the specific kids with HIV that have low CD4 (N=11)
test <- chi_cd4[c("sample_id", "cd4nL")]
#remove B.01, B.03, B.04, B.05, B.09, B.11, B.18, B.25, B.27, B.32, B.46
  #also need to remove the kids missing CD4: B.45 and B.52
#so 50 HUU kids + 44 HEI kids - 2 kids without CD4 - 11 kids with low CD4 = 81
hucd4 = subset_samples(chi_prune, sample_id != "B.01.CHI" & sample_id != "B.03.CHI" & sample_id != "B.04.CHI" &
                         sample_id != "B.05.CHI" & sample_id != "B.09.CHI" & sample_id != "B.11.CHI" &
                         sample_id != "B.18.CHI" & sample_id != "B.25.CHI" & sample_id != "B.27.CHI" &
                         sample_id != "B.32.CHI" & sample_id != "B.45.CHI" & sample_id != "B.46.CHI" &
                         sample_id != "B.52.CHI")
hu <- data.frame(sample_data(hucd4))
table(hu$hiv) #adds up: 50 HUU kids + 49 HEU kids + 31 kids with HIV (44 orig - 13 removed above)

#SECOND: PCOA plots and PERMANOVA
#Transform our pseq object
hucd4_clr <- transform(hucd4, "clr")

#PCoA plot#
clr_pcoa <- ordinate(hucd4_clr, method = "PCoA", distance = "euclidean")
lines <- c("dashed", "dashed", "dashed")
clr_eiuu <- plot_ordination(hucd4_clr, clr_pcoa, color = "hiv", title = "Euclidean PCoA Plot by HIV Status_IC only") +
  geom_point(size = 2) + theme_classic() + theme(
    plot.title = element_text(size = 14, hjust = 0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 10),
    legend.position = "right",
    legend.background = element_blank(),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12)) +
  scale_color_manual(values=c('#800000FF','#155F83FF', '#FFA319FF'))+
  stat_ellipse(aes(linetype=hiv), geom="polygon", alpha=0, type="t", level=0.8, size=0.7) +
  xlab("PC1 (14.9%)") + ylab("PC2 (9.1%)") +
  scale_linetype_manual(values=lines)
png(file="clr_hiv_IC only.png",
    width = 6, height = 4, units = 'in', res=600)
clr_eiuu
dev.off()

#PERMANOVA
otu <- as.data.frame(otu_table(hucd4_clr, taxa_are_rows=TRUE))
#Need species in columns and observations in rows: transpose using t() function and save as a dataframe
#Source: https://github.com/edamame-course/2015-tutorials/blob/master/final/Ashley_Intro_to_R.md
otu_trans <- data.frame(t(otu))
BC.dist=vegdist(otu_trans, method="euclidean")

#hiv + age (per 6/9 MK discussion to account for age in all PERMANOVAs)
adonis(BC.dist ~ hiv + age, data = hu, permutations = 1000)
#No significant difference by HIV status among children (P = 0.18), but age remains significant (P=0.02, R2 = 0.014)
#Now getting message that adonis will be deprecated and I need to use adonis2
adonis2(BC.dist ~ hiv + age, data = hu, permutations = 1000)
#No significant difference by HIV status among children (P = 0.15, R2=0.018), but age remains significant (P=0.012, R2 = 0.014)
#check with age first given issues with order of variables in adonis; not signif different
adonis2(BC.dist ~ age + hiv, data = hu, permutations = 1000)

#clean up
remove(hucd4_clr, test, otu, otu_trans, clr_eiuu, clr_pcoa, bloop, BC.dist, lines)

#RELATIVE ABUNDANCE PLOTS
chi.rel <- transform_sample_counts(hucd4, function(Abundance) Abundance/sum(Abundance))
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
#merging on Phylum AND Genus gets rid of the phylum.x and phylum.y problem we were having earlier

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
relative_df <- merge(melted_df, abundances, by=c("Genus", "Species"))
#keep getting vector memory exhausted error; need to reboot R so will save melted_df and abundances as csv files
#source for fix: https://stackoverflow.com/questions/51295402/r-on-macos-error-vector-memory-exhausted-limit-reached

relative_df$Species[!(relative_df$Species %in% TOPSpecies)] <- "Other"
sum(relative_df$Abundance) #sums to 130.0008 (and we have n=130 kids...?)
table(relative_df$Species)

#To generate color palette of 11 colors used previously: 
library(RColorBrewer)
library(scales)
chisp2 <- c("#293352", "#800000FF", "#767676FF", "#FFA319FF", "#8A9045FF", "#155F83FF", "#C16622FF", "#4D004B", "#8F3931FF", "#58593FFF", "#350E20FF")
show_col(chisp2)

#CREATING RELATIVE ABUNDANCE PLOT
relative_df$subject_abundances <- (relative_df$Abundance)/(length(unique(relative_df$Sample))) 
#abbreviate HIV status
relative_df$hiv2 <- NA
relative_df$hiv2 [relative_df$hiv == "Infected"] <- "HEI"
relative_df$hiv2 [relative_df$hiv == "Exposed Uninfected"] <- "HEU"
relative_df$hiv2 [relative_df$hiv == "Unexposed"] <- "HUU"
table(relative_df$hiv2)
#matches up with table(relative_df$hiv)
relative_df$hiv2 <- as.factor(relative_df$hiv2)
relative_df$hiv2 <- reorder(relative_df$hiv2, new.order=c("HUU", "HEU", "HEI"))
relative_df$Species <- factor(relative_df$Species)
relative_df$Species <- reorder(relative_df$Species, new.order=c("Corynebacterium propinquum", "Corynebacterium pseudodiphtheriticum", "Dolosigranulum pigrum",  
                                                                "Haemophilus influenzae", "Micrococcus luteus", "Moraxella catarrhalis", "Moraxella lincolnii", 
                                                                "Moraxella nonliquefaciens", "Staphylococcus aureus", "Streptococcus pneumoniae", "Other"))
subject_plot2 <- ggplot(relative_df[order(relative_df$Species, decreasing = TRUE),], aes(x=hiv2, y=subject_abundances, fill=Species)) + 
  geom_bar(stat="identity", position="fill") + theme_classic() + theme(legend.position="right", 
                                                                       legend.text=element_text(size=12, face="italic"),
                                                                       axis.title.y = element_text(size=14, margin=margin(0,20,0,0)), axis.text.y = element_text(size=12), 
                                                                       axis.title.x = element_text(size=14, margin=margin(10,0,0,0)), axis.text.x = element_text(size=12)) + 
  scale_fill_manual(values = chisp2) +
  # scale_x_discrete(labels=c("No", "Yes")) +
  xlab("HIV Status") + ylab("Relative Abundance") 
png(file="species_IConly_07052022.png",
    width = 9, height = 5, units = 'in', res = 600)
plot(subject_plot2)
dev.off()

#MAASLIN 2: HEI (immunocompetent only), HEU, HUU
hu.filter <- filter_taxa(hucd4, function(Abundance) mean(Abundance)>=50, TRUE)
ntaxa(hu.filter) #51

otu <- as.data.frame(otu_table(hu.filter, taxa_are_rows=TRUE))
#need to remove the unidentified species from here: row 10, 27, 44, 45, 47
#source: https://stackoverflow.com/questions/12328056/how-do-i-delete-rows-in-a-data-frame
#corresponding respectively to: Corynebacterium sp. KPL1859, Streptococcus sp. SK643, Alkalibacterium sp. 20, uncultured Clostridium sp., Nocardioides sp. Root122
otu <- otu[-c(10, 27, 44, 45, 47),]

#Univariable
fit_data2 = Maaslin2(
  input_data = otu, 
  input_metadata = hu, 
  min_prevalence = 0.2,
  max_significance = 0.2,
  output = "maaslin_immunocompetent_UV", 
  fixed_effects = c("hiv"),
  reference = c("hiv,Unexposed"))
#No significant assoc

fit_data2 = Maaslin2(
  input_data = otu, 
  input_metadata = hu, 
  min_prevalence = 0.2,
  max_significance = 0.2,
  output = "maaslin_immunocompetent_mv", 
  fixed_effects = c("hiv", "season", "age", "abx", "wood", "uri_rec2"),
  reference = c("hiv,Unexposed"))
#Signif assoc with age and season, but not HIV

###K-W tests: HEI (IC only), HEU, HUU
#Corynebacterium propinquum
cor_p <- filter(melted_df, Species == "Corynebacterium propinquum")
summary(cor_p$Abundance)
hist(cor_p$Abundance)
tapply(cor_p$Abundance, cor_p$hiv, summary) #medians similar
kruskal.test(cor_p$Abundance, cor_p$hiv)  #P=0.59

#Corynebacterium pseudodiphtheriticum
cor_ps <- filter(melted_df, Species == "Corynebacterium pseudodiphtheriticum")
summary(cor_ps$Abundance)
hist(cor_ps$Abundance)
tapply(cor_ps$Abundance, cor_ps$hiv, summary) #medians similar
kruskal.test(cor_ps$Abundance, cor_ps$hiv)  #P=0.12

#Dolosigranulum pigrum
dpig <- filter(melted_df, Species == "Dolosigranulum pigrum")
hist(dpig$Abundance)
summary(dpig$Abundance)
tapply(dpig$Abundance, dpig$hiv, summary) #medians similar
kruskal.test(dpig$Abundance, dpig$hiv)  #P=0.40

#Haemophilus influenzae
hflu <- filter(melted_df, Species == "Haemophilus influenzae")
hist(hflu$Abundance)  #lots of zeros
summary(hflu$Abundance)
tapply(hflu$Abundance, hflu$hiv, summary) #medians similar
kruskal.test(hflu$Abundance, hflu$hiv)  #P=0.39

#Micrococcus luteus
mic <- filter(melted_df, Species == "Micrococcus luteus")
hist(mic$Abundance)  #lots of zeros
summary(mic$Abundance)
tapply(mic$Abundance, mic$hiv, summary) #lowest median among children with HIV
kruskal.test(mic$Abundance, mic$hiv)  #P=0.17

#Moraxella cattarhalis
mcat <- filter(melted_df, Species == "Moraxella catarrhalis")
hist(mcat$Abundance)  #lots of zeros
summary(mcat$Abundance)
tapply(mcat$Abundance, mcat$hiv, summary) #highest median abundance among HUU
kruskal.test(mcat$Abundance, mcat$hiv)  #P=0.44

#Moraxella lincolnii
mlin <- filter(melted_df, Species == "Moraxella lincolnii")
hist(mlin$Abundance)  #lots of zeros
summary(mlin$Abundance)
tapply(mlin$Abundance, mlin$hiv, summary) #medians similar
kruskal.test(mlin$Abundance, mlin$hiv)  #P=0.07

#Moraxella nonliquefaciens
mnon <- filter(melted_df, Species == "Moraxella nonliquefaciens")
hist(mnon$Abundance)
summary(mnon$Abundance)
tapply(mnon$Abundance, mnon$hiv, summary) #highest med abundance among children with HIV
kruskal.test(mnon$Abundance, mnon$hiv)  #P=0.24

#Staphylococcus aureus
sa <- filter(melted_df, Species == "Staphylococcus aureus")
hist(sa$Abundance)  #lots of zeros
summary(sa$Abundance)
tapply(sa$Abundance, sa$hiv, summary) #limited by zeros (all medians = 0)
kruskal.test(sa$Abundance, sa$hiv)  #P=0.03 but hard to comment on significance

#Streptococcus pneumoniae
spneum <- filter(melted_df, Species == "Streptococcus pneumoniae")
hist(spneum$Abundance)  #lots of zeros
summary(spneum$Abundance)
tapply(spneum$Abundance, spneum$hiv, summary)   #lowest abundance among children with HIV
kruskal.test(spneum$Abundance, spneum$hiv)  #P=0.37

###############################################
#ANALYSIS REMOVING HEI CHILDREN WITH CD4 >= 25%
###############################################
#1. Compare overall composition using age-adjusted permanova
#2. K-W tests

#now need to remove the specific kids with HIV that have high CD4 (N=31)
test <- chi_cd4[c("sample_id", "cd4nL")]
#remove B.06, B.07, B.12, B.13, B.14, B.15, B.16, B.17, B.19, B.20, B.22, B.23, 
  #B.26, B.28, B.29, B.30, B.31, B.33, B.34, B.35, B.36, B.37, B.38, B.39, B.40, B.41,
  #B.44, B.48, B.49, B.50, B.51
#also need to remove the kids missing CD4: B.45 and B.52
#so 50 HUU kids + 49 HEU kids + 44 HEI kids - 2 kids without CD4 - 31 kids with NORMAL CD4 = 110 kids
hlowcd4 = subset_samples(chi_prune, sample_id != "B.06.CHI" & sample_id != "B.07.CHI" & sample_id != "B.12.CHI" &
              sample_id != "B.13.CHI" & sample_id != "B.14.CHI" & sample_id != "B.15.CHI" &
              sample_id != "B.16.CHI" & sample_id != "B.17.CHI" & sample_id != "B.19.CHI" &
              sample_id != "B.20.CHI" & sample_id != "B.22.CHI" & sample_id != "B.23.CHI" &
              sample_id != "B.26.CHI" & sample_id != "B.28.CHI" & sample_id != "B.29.CHI" &
              sample_id != "B.30.CHI" & sample_id != "B.31.CHI" & sample_id != "B.33.CHI" &
              sample_id != "B.34.CHI" & sample_id != "B.35.CHI" & sample_id != "B.36.CHI" &
              sample_id != "B.37.CHI" & sample_id != "B.38.CHI" & sample_id != "B.39.CHI" &
              sample_id != "B.40.CHI" & sample_id != "B.41.CHI" & sample_id != "B.44.CHI" &
              sample_id != "B.45.CHI" & sample_id != "B.48.CHI" & sample_id != "B.49.CHI" &
              sample_id != "B.50.CHI" & sample_id != "B.51.CHI" & sample_id != "B.52.CHI")
hlow <- data.frame(sample_data(hlowcd4))
table(hlow$hiv) #adds up: 50 HUU kids + 49 HEU kids + 11 kids with HIV (44 orig - 33 removed above)

#SECOND: PERMANOVA +/- PCoA plot
#Transform our pseq object
hlowcd4_clr <- transform(hlowcd4, "clr")

#PERMANOVA
otu <- as.data.frame(otu_table(hlowcd4_clr, taxa_are_rows=TRUE))
#Need species in columns and observations in rows: transpose using t() function and save as a dataframe
#Source: https://github.com/edamame-course/2015-tutorials/blob/master/final/Ashley_Intro_to_R.md
otu_trans <- data.frame(t(otu))
BC.dist=vegdist(otu_trans, method="euclidean")

#hiv + age (per 6/9 MK discussion to account for age in all PERMANOVAs)
adonis2(BC.dist ~ hiv + age, data = hlow, permutations = 1000)
#Significant difference by HIV status among children (P = 0.003, R2 = 0.032) and age (P=0.03, R2 = 0.015)
#check with age first given issues with order of variables in adonis; not signif different
adonis2(BC.dist ~ age + hiv, data = hlow, permutations = 1000)

#PCoA plot#
table(hlow$hiv)
hlow$hiv2 <- NA
hlow$hiv2 [hlow$hiv == "Infected"] <- "Immunocompromised children with HIV"
hlow$hiv2 [hlow$hiv == "Exposed Uninfected"] <- "HEU children"
hlow$hiv2 [hlow$hiv == "Unexposed"] <- "HUU children"
table(hlow$hiv2)
class(hlow$hiv2)
hlow$hiv2 <- as.factor(hlow$hiv2)
hlow$hiv2 <- reorder(hlow$hiv2, new.order=c("Immunocompromised children with HIV", "HEU children", "HUU children"))

#Import back into pseq (source: https://jacobrprice.github.io/2017/08/26/phyloseq-to-vegan-and-back.html)
sample_data(hlowcd4_clr) <- as.data.frame(hlow)
#now see if it worked:
bloop <- data.frame(sample_data(hlowcd4_clr))
#success!
remove(bloop)

clr_pcoa <- ordinate(hlowcd4_clr, method = "PCoA", distance = "euclidean")
lines <- c("dashed", "dashed", "dashed")

p5 <- plot_ordination(hlowcd4_clr, clr_pcoa, color = "hiv2") +
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
  xlab("PC1 (15.4%)") + ylab("PC2 (9.2%)") +
  scale_linetype_manual(values= c("dashed", "dashed", "dashed"))

png(file="clr_lowcd4_HEU_HUU.png",
    width = 5, height = 4, units = 'in', res=600)
p5
dev.off()

#clean up
remove(clr_pcoa, hlowcd4_clr, otu, otu_trans, p5, BC.dist, lines)

#CREATE MELTED_DF FOR K-W TESTS
chi.rel <- transform_sample_counts(hlowcd4, function(Abundance) Abundance/sum(Abundance))
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

###K-W tests: HEI (low CD4 only), HEU, HUU
#Corynebacterium propinquum
cor_p <- filter(melted_df, Species == "Corynebacterium propinquum")
summary(cor_p$Abundance)
hist(cor_p$Abundance)
tapply(cor_p$Abundance, cor_p$hiv, summary) #median lowest in low CD4
kruskal.test(cor_p$Abundance, cor_p$hiv)  #P=0.001

#Corynebacterium pseudodiphtheriticum
cor_ps <- filter(melted_df, Species == "Corynebacterium pseudodiphtheriticum")
summary(cor_ps$Abundance)
hist(cor_ps$Abundance)
tapply(cor_ps$Abundance, cor_ps$hiv, summary) #median lowest in low CD4
kruskal.test(cor_ps$Abundance, cor_ps$hiv)  #P=0.002

#Dolosigranulum pigrum
dpig <- filter(melted_df, Species == "Dolosigranulum pigrum")
hist(dpig$Abundance)
summary(dpig$Abundance)
tapply(dpig$Abundance, dpig$hiv, summary) #median lowest in low CD4
kruskal.test(dpig$Abundance, dpig$hiv)  #P=0.03

#Haemophilus influenzae
hflu <- filter(melted_df, Species == "Haemophilus influenzae")
hist(hflu$Abundance)  #lots of zeros
summary(hflu$Abundance)
tapply(hflu$Abundance, hflu$hiv, summary) #median highest in low CD4
kruskal.test(hflu$Abundance, hflu$hiv)  #P=0.32

#Micrococcus luteus
mic <- filter(melted_df, Species == "Micrococcus luteus")
hist(mic$Abundance)  #lots of zeros
summary(mic$Abundance)
tapply(mic$Abundance, mic$hiv, summary) #lowest median among low CD4
kruskal.test(mic$Abundance, mic$hiv)  #P=0.04

#Moraxella cattarhalis
mcat <- filter(melted_df, Species == "Moraxella catarrhalis")
hist(mcat$Abundance)  #lots of zeros
summary(mcat$Abundance)
tapply(mcat$Abundance, mcat$hiv, summary) #highest median abundance among low CD4
kruskal.test(mcat$Abundance, mcat$hiv)  #P=0.39

#Moraxella lincolnii
mlin <- filter(melted_df, Species == "Moraxella lincolnii")
hist(mlin$Abundance)  #lots of zeros
summary(mlin$Abundance)
tapply(mlin$Abundance, mlin$hiv, summary) #medians similar
kruskal.test(mlin$Abundance, mlin$hiv)  #P=0.21

#Moraxella nonliquefaciens
mnon <- filter(melted_df, Species == "Moraxella nonliquefaciens")
hist(mnon$Abundance)
summary(mnon$Abundance)
tapply(mnon$Abundance, mnon$hiv, summary) #highest med abundance among HUU
kruskal.test(mnon$Abundance, mnon$hiv)  #P=0.29

#Staphylococcus aureus
sa <- filter(melted_df, Species == "Staphylococcus aureus")
hist(sa$Abundance)  #lots of zeros
summary(sa$Abundance)
tapply(sa$Abundance, sa$hiv, summary) #limited by zeros (low cd4 has highest median of 0.0001)
kruskal.test(sa$Abundance, sa$hiv)  #P=0.01

#Streptococcus pneumoniae
spneum <- filter(melted_df, Species == "Streptococcus pneumoniae")
hist(spneum$Abundance)  #lots of zeros
summary(spneum$Abundance)
tapply(spneum$Abundance, spneum$hiv, summary)   #highest abundance among children with HIV
kruskal.test(spneum$Abundance, spneum$hiv)  #P=0.23

#RELATIVE ABUNDANCE PLOTS (skipped for now 8-31-22)
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
#merging on Phylum AND Genus gets rid of the phylum.x and phylum.y problem we were having earlier

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
relative_df <- merge(melted_df, abundances, by=c("Genus", "Species"))
#keep getting vector memory exhausted error; need to reboot R so will save melted_df and abundances as csv files
#source for fix: https://stackoverflow.com/questions/51295402/r-on-macos-error-vector-memory-exhausted-limit-reached

relative_df$Species[!(relative_df$Species %in% TOPSpecies)] <- "Other"
sum(relative_df$Abundance) #sums to 130.0008 (and we have n=130 kids...?)
table(relative_df$Species)

#To generate color palette of 11 colors used previously: 
library(RColorBrewer)
library(scales)
chisp2 <- c("#293352", "#800000FF", "#767676FF", "#FFA319FF", "#8A9045FF", "#155F83FF", "#C16622FF", "#4D004B", "#8F3931FF", "#58593FFF", "#350E20FF")
show_col(chisp2)

#CREATING RELATIVE ABUNDANCE PLOT
relative_df$subject_abundances <- (relative_df$Abundance)/(length(unique(relative_df$Sample))) 
#abbreviate HIV status
relative_df$hiv2 <- NA
relative_df$hiv2 [relative_df$hiv == "Infected"] <- "HEI"
relative_df$hiv2 [relative_df$hiv == "Exposed Uninfected"] <- "HEU"
relative_df$hiv2 [relative_df$hiv == "Unexposed"] <- "HUU"
table(relative_df$hiv2)
#matches up with table(relative_df$hiv)
relative_df$hiv2 <- as.factor(relative_df$hiv2)
relative_df$hiv2 <- reorder(relative_df$hiv2, new.order=c("HUU", "HEU", "HEI"))
relative_df$Species <- factor(relative_df$Species)
relative_df$Species <- reorder(relative_df$Species, new.order=c("Corynebacterium propinquum", "Corynebacterium pseudodiphtheriticum", "Dolosigranulum pigrum",  
                                                                "Haemophilus influenzae", "Micrococcus luteus", "Moraxella catarrhalis", "Moraxella lincolnii", 
                                                                "Moraxella nonliquefaciens", "Staphylococcus aureus", "Streptococcus pneumoniae", "Other"))
subject_plot2 <- ggplot(relative_df[order(relative_df$Species, decreasing = TRUE),], aes(x=hiv2, y=subject_abundances, fill=Species)) + 
  geom_bar(stat="identity", position="fill") + theme_classic() + theme(legend.position="right", 
                                                                       legend.text=element_text(size=12, face="italic"),
                                                                       axis.title.y = element_text(size=14, margin=margin(0,20,0,0)), axis.text.y = element_text(size=12), 
                                                                       axis.title.x = element_text(size=14, margin=margin(10,0,0,0)), axis.text.x = element_text(size=12)) + 
  scale_fill_manual(values = chisp2) +
  # scale_x_discrete(labels=c("No", "Yes")) +
  xlab("HIV Status") + ylab("Relative Abundance") 
png(file="species_IConly_07052022.png",
    width = 9, height = 5, units = 'in', res = 600)
plot(subject_plot2)
dev.off()

remove(chi_hiv, chi.rel, cor_p, cor_ps, dpig, hflu, hiv, hlow, hlowcd4, mcat, 
       melted_df, mic, mlin, mnon, sa, spneum)

#######################################################################
#NEGATIVE BINOMIAL MODELS FOR ABUNDANCES OF SPECIFIC SPECIES ABUNDANCES
#######################################################################
#Using labs phyloseq object; n=40 but all children have CD4 and viral load measurements
#1. Create melted_df and then species-specific dataframes
#2. Merge abundance data for all species and save as a csv file
#3. Convert abundances to integers and attempt neg binomial models
  #starting with D. pigrum; if works, will run for C. pseudo, C. propinquum, S.aureus, S.pneumo

#MELTED_DF
chi.rel <- transform_sample_counts(labs, function(Abundance) Abundance/sum(Abundance))
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

#Species-specific dataframes
cor_p <- filter(melted_df, Species == "Corynebacterium propinquum")
cor_ps <- filter(melted_df, Species == "Corynebacterium pseudodiphtheriticum")
dpig <- filter(melted_df, Species == "Dolosigranulum pigrum")
hflu <- filter(melted_df, Species == "Haemophilus influenzae")
mic <- filter(melted_df, Species == "Micrococcus luteus")
mcat <- filter(melted_df, Species == "Moraxella catarrhalis")
mlin <- filter(melted_df, Species == "Moraxella lincolnii")
mnon <- filter(melted_df, Species == "Moraxella nonliquefaciens")
sa <- filter(melted_df, Species == "Staphylococcus aureus")
spneum <- filter(melted_df, Species == "Streptococcus pneumoniae")

#Dataframe containing relevant covariables + abundances for the above species
cor_p <- rename(cor_p, ab.cor_p = Abundance)
cor_ps <- rename(cor_ps, ab.cor_ps = Abundance)
dpig <- rename(dpig, ab.dpig = Abundance)
hflu <- rename(hflu, ab.hflu = Abundance)
mic <- rename(mic, ab.mic = Abundance)
mcat <- rename(mcat, ab.mcat = Abundance)
mlin <- rename(mlin, ab.mlin = Abundance)
mnon <- rename(mnon, ab.mnon = Abundance)
spneum <- rename(spneum, ab.spneum = Abundance)
sa <- rename(sa, ab.sa = Abundance)

#Create simple dataframe with abundance variables, sample_id, cd4, age, and merge on sample_id
corr <- cor_ps[c("sample_id", "age", "pcv", "dpt", "ab.cor_ps", "newest_cd4perc", "cd4nL", "tmpsmx", "abx", "vl_suppr")]
test <- cor_p[c("sample_id", "ab.cor_p")]
corr <- merge(corr, test, by ="sample_id", sort = TRUE)
test <- dpig[c("sample_id", "ab.dpig")]
corr <- merge(corr, test, by ="sample_id", sort = TRUE)
test <- hflu[c("sample_id", "ab.hflu")]
corr <- merge(corr, test, by ="sample_id", sort = TRUE)
test <- mcat[c("sample_id", "ab.mcat")]
corr <- merge(corr, test, by ="sample_id", sort = TRUE)
test <- spneum[c("sample_id", "ab.spneum")]
corr <- merge(corr, test, by ="sample_id", sort = TRUE)
test <- sa[c("sample_id", "ab.sa")]
corr <- merge(corr, test, by ="sample_id", sort = TRUE)
test <- mic[c("sample_id", "ab.mic")]
corr <- merge(corr, test, by ="sample_id", sort = TRUE)
test <- mlin[c("sample_id", "ab.mlin")]
corr <- merge(corr, test, by ="sample_id", sort = TRUE)
test <- mnon[c("sample_id", "ab.mnon")]
corr <- merge(corr, test, by ="sample_id", sort = TRUE)

#save as csv in case we want to come back to it
write.csv(corr, "data_negbin.csv", row.names=F)
negbin <- read.csv("data_negbin.csv")

#RUNNING ZERO-INFLATED NEGATIVE BINOMIAL MODELS
#From RY teaching in RSV mb analysis: Negative binomial is meant for count data --> so multiply your relative abundance by 100 and round to the nearest integer
#Use mutate function to multiply, and there are rounding functions in R (and dplyr has one) --> RY sent link
#Interpretation of neg binomial models --> will get a p-value and an incident rate ratio
#If IRR was 1.2 with CI not crossing 1 --> 20% increase in the relative abundance in genus A in people who had low CD4 compared to nL CD4
#KNOW YOUR REFERENCE GROUP
#Model will show you what your overdispersion constant is --> look at dispersion parameter to help you decide if you should be using a neg binomial or a zero-inflated neg binomial
#UCLA zero-inflated neg bin tutorial: https://stats.oarc.ucla.edu/r/dae/zinb/
#another tutorial: https://fukamilab.github.io/BIO202/04-C-zero-data.html
#Looking at our data, we will multiply all abundances by 1000 and then round

###NOTE: Also created RMD file for negative binomial models: HIV MB Zinf NB RMD (for MK review)

library(pscl)
#Question: which covariable to include in the logit portion of the model that tells you what process the zeros are a/w?
#CD4? because low or normal CD4 might influence whether a species abundance is 0 or not?
class(negbin$cd4nL) #integer --> need to convert to factor and relevel so normal CD4 is the ref
class(negbin$vl_suppr) #integer --> need to convert to factor and relevel so viral suppr is the ref

negbin$cd4nL <- as.factor(negbin$cd4nL)
negbin$vl_suppr <- as.factor(negbin$vl_suppr)
negbin$cd4nL <- relevel(negbin$cd4nL, ref = "1")
negbin$vl_suppr <- relevel(negbin$vl_suppr, ref = "1")

#D. pigrum was associated with both TMP-SMX and CD4, so we will use this species as our test run
#Creating integers for our independent variable
negbin$ab.dpig.new <- negbin$ab.dpig*1000
negbin$ab.dpig.new <- round(negbin$ab.dpig.new, digits = 0)

summary(m1 <- zeroinfl(ab.dpig.new ~ cd4nL + vl_suppr + tmpsmx + abx | cd4nL,
               data = negbin, dist = "negbin"))
AIC(m1) #487.68
#nothing is significantly associated; does this change if we change the var in the logit portion? No
summary(m1 <- zeroinfl(ab.dpig.new ~ cd4nL + vl_suppr + tmpsmx + abx | tmpsmx,
                       data = negbin, dist = "negbin"))

#is the CD4 percentage as a continuous variable associated? No
summary(m1 <- zeroinfl(ab.dpig.new ~ newest_cd4perc + vl_suppr + tmpsmx + abx | cd4nL,
                       data = negbin, dist = "negbin"))
AIC(m1) #488.08

#run standard negative binomial models for each organism and see which AIC is smaller
#CD4 binary
summary(m2 <- glm.nb(ab.dpig.new ~ cd4nL + vl_suppr + tmpsmx + abx, data = negbin))
#AIC 497.3

#CD4 continuous
summary(m2 <- glm.nb(ab.dpig.new ~ newest_cd4perc + vl_suppr + tmpsmx + abx, data = negbin))
#AIC 496.28

#C. propinquum
negbin$ab.cor_p.new <- negbin$ab.cor_p*1000
negbin$ab.cor_p.new <- round(negbin$ab.cor_p.new, digits = 0)
#CD4 as binary variable: low CD4 (p<0.0001) and TMP-SMX exposure (p=0.03) associated with lower abundance
  #detectable viral load a/w more C. propinquum? Positive estimate; weird. p=0.03
summary(m1 <- zeroinfl(ab.cor_p.new ~ cd4nL + vl_suppr + tmpsmx + abx | cd4nL,
                       data = negbin, dist = "negbin"))
AIC(m1) #347.9
#CD4 as continuous variable: no associations. Interesting.
summary(m1 <- zeroinfl(ab.cor_p.new ~ newest_cd4perc + vl_suppr + tmpsmx + abx | cd4nL,
                       data = negbin, dist = "negbin"))
AIC(m1) #354.1

#run standard negative binomial models for each organism and see which AIC is smaller
#CD4 binary
summary(m2 <- glm.nb(ab.cor_p.new ~ cd4nL + vl_suppr + tmpsmx + abx, data = negbin))
#AIC 344.73

#CD4 continuous
summary(m2 <- glm.nb(ab.cor_p.new ~ newest_cd4perc + vl_suppr + tmpsmx + abx, data = negbin))
#AIC 350 BUT ALGORITHM DID NOT CONVERGE

#C. pseudodiphtheriticum
negbin$ab.cor_ps.new <- negbin$ab.cor_ps*1000
negbin$ab.cor_ps.new <- round(negbin$ab.cor_ps.new, digits = 0)
#CD4 as binary variable: detectable viral load (p=0.004) associated with lower abundance of C. pseudo
summary(m1 <- zeroinfl(ab.cor_ps.new ~ cd4nL + vl_suppr + tmpsmx + abx | cd4nL,
                       data = negbin, dist = "negbin"))
AIC(m1) #343.6
#CD4 as continuous variable: negative assoc with cd4% (p=0.05) and negative estimate for detectable vL (p=0.001)
summary(m1 <- zeroinfl(ab.cor_ps.new ~ newest_cd4perc + vl_suppr + tmpsmx + abx | cd4nL,
                       data = negbin, dist = "negbin"))
AIC(m1)  #341.4

#run standard negative binomial models for each organism and see which AIC is smaller
#CD4 binary
summary(m2 <- glm.nb(ab.cor_ps.new ~ cd4nL + vl_suppr + tmpsmx + abx, data = negbin))
#AIC 355 BUT ALGORITHM DID NOT CONVERGE

#CD4 continuous
summary(m2 <- glm.nb(ab.cor_ps.new ~ newest_cd4perc + vl_suppr + tmpsmx + abx, data = negbin))
#AIC 354 BUT ALGORITHM DID NOT CONVERGE

#S. aureus: abundances are so low, so will multiply by 10K here
negbin$ab.sa.new <- negbin$ab.sa*10000
negbin$ab.sa.new <- round(negbin$ab.sa.new, digits = 0)
#CD4 as binary variable: TMP-SMX associated with higher abundance of S. aureus BUT error message:
  #error: In sqrt(diag(object$vcov)) : NaNs produced --> because of the large number of zeros?
summary(m1 <- zeroinfl(ab.sa.new ~ cd4nL + vl_suppr + tmpsmx + abx | cd4nL,
                       data = negbin, dist = "negbin"))
AIC(m1) #162.6
#CD4 as continuous variable: positive assoc with cd4% (p=0.003) and abx (p=0.002)
  #no error message with this model
summary(m1 <- zeroinfl(ab.sa.new ~ newest_cd4perc + vl_suppr + tmpsmx + abx | cd4nL,
                       data = negbin, dist = "negbin"))
AIC(m1) #170.4

#run standard negative binomial models for each organism and see which AIC is smaller
#CD4 binary
summary(m2 <- glm.nb(ab.sa.new ~ cd4nL + vl_suppr + tmpsmx + abx, data = negbin))
#CD4 continuous
summary(m2 <- glm.nb(ab.sa.new ~ newest_cd4perc + vl_suppr + tmpsmx + abx, data = negbin))
#NOTE: neither model works; error message about valid set of coefficients; suspect 2/2 high #0s?

#S. pneumoniae
negbin$ab.spneum.new <- negbin$ab.spneum*1000
negbin$ab.spneum.new <- round(negbin$ab.spneum.new, digits = 0)
#CD4 as binary variable: no associations with S. pneumo abundance
summary(m1 <- zeroinfl(ab.spneum.new ~ cd4nL + vl_suppr + tmpsmx + abx | cd4nL,
                       data = negbin, dist = "negbin"))
AIC(m1) #320.6
#CD4 as continuous variable: no associations
summary(m1 <- zeroinfl(ab.spneum.new ~ newest_cd4perc + vl_suppr + tmpsmx + abx | cd4nL,
                       data = negbin, dist = "negbin"))
AIC(m1) #319.6

#run standard negative binomial models for each organism and see which AIC is smaller
#CD4 binary
summary(m2 <- glm.nb(ab.spneum.new ~ cd4nL + vl_suppr + tmpsmx + abx, data = negbin))
#AIC 322 BUT ALGORITHM DID NOT CONVERGE

#CD4 continuous
summary(m2 <- glm.nb(ab.spneum.new ~ newest_cd4perc + vl_suppr + tmpsmx + abx, data = negbin))
#AIC 318 BUT ALGORITHM DID NOT CONVERGE

###################################
#TABLE S1. TAXONOMY OF CONTAMINANTS
###################################
#For manuscript table S1, will create table just of taxonomy and export
contam <- read.csv("Kaiju_contam_freq_0.1.csv")
contam_sp <- contam[c("Species")]
#Now split the one column with all tax data into several columns
contam_sp <- contam_sp %>% 
  dplyr::rename("tax" = "Species") %>% 
  dplyr::select(tax) %>%
  separate(tax, sep = "\\;", remove = FALSE, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))
contam_sp$tax <- NULL
#export as csv
write.csv(contam_sp, "contam_tableS1.csv", row.names=F)