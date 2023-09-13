#Fogarty phyloseq object creation using Kaiju pipeline results 
#Compiled Feb 16, 2022 by Sweta Patel
#UPDATED with new Kaiju pipeline results from MK on Feb 22, 2022
#Used code from MK (2/16/22 email) + my prior phyloseq code 
#UPDATED to run decontam package on April 12, 2022
#UPDATED to run decontam package and only remove taxa by frequency method on May 6, 2022
#Lancet Microbe Revisions: September 12, 2023

#***LOADING PACKAGES***
library(plyr)
library(dplyr)
library(ggplot2)
library(gplots)
library(phyloseq)
library(tidyr)
library(vegan)
library(metagenomeSeq)
library(httr)
library(gridExtra)
library(data.table)
library(readxl)
library(gsubfn)
library(decontam)

set.seed(1234)

setwd("OneDrive - Duke University/Fogarty coding/Sequencing/")

#***LOADING DATA***
#Cleaned metadata
meta <- read.csv("/Users/swetapatel/OneDrive - Duke University/Fogarty coding/Sequencing/metadata_clean_10-7-21.csv")

#reads
reads <- read.csv("/Users/swetapatel/OneDrive - Duke University/Fogarty coding/Sequencing/HIV_NGS_kneaddata_counts.csv")

#taxonomy
orgs <- read.csv("/Users/swetapatel/OneDrive - Duke University/Fogarty coding/Sequencing/HIV_NGS_kaiju.csv")

#DNA concentration
conc <- read.csv("/Users/swetapatel/OneDrive - Duke University/Fogarty coding/Sequencing/6968_LIB qubit.csv")
# library(readxl)
# conc <- read_excel("6968_LIB qubit.xlsx", 
#                    col_types = c("skip", "text", "numeric"))
View(conc)
#for some reason have a few blank rows at the end; delete
conc <- conc[-c(312, 313, 314, 315, 316), ] 

#Need to standardize sample_id names and reformat data so we can put it into pseq obj
names(orgs)[names(orgs)=="file"] <- "sample_id"
orgs$sample_id <- gsub("/work/msk37/HIV_NGS/6_kaiju/","", as.factor(orgs$sample_id))
orgs$sample_id <- gsub("_kaiju_mem.out","", as.factor(orgs$sample_id))
kaiju <- subset(orgs, !is.na(orgs$taxon_id))
kaiju <- subset(kaiju, reads!=0)
kaiju <- kaiju[,c("sample_id","reads","taxon_name")]
kaiju <- reshape(kaiju, idvar="sample_id", timevar="taxon_name", direction="wide")
kaiju[is.na(kaiju)] <- 0
row.names(kaiju) <- kaiju$sample_id
kaiju <- kaiju[,-1]
kaiju <- data.frame(t(kaiju))
rownames(kaiju) <- gsub("reads.","", row.names(kaiju))
write.csv(kaiju, "HIV_NGS_kaiju_formatted.csv")
kaiju <- read.csv("HIV_NGS_kaiju_formatted.csv", sep=",", stringsAsFactors = T)
row.names(kaiju) <- kaiju[,1]
kaiju <- kaiju[,-1]
kaiju[is.na(kaiju)] <- 0
# kaiju$ext_neg_ctrl <- NULL    #need neg control
kaiju$zymo_pos_ctrl <- NULL
kaiju_ids <- colnames(kaiju)

######################################
#PHYLOSEQ OBJECT CREATION FOR DECONTAM
######################################
otudata <- otu_table(kaiju, taxa_are_rows = TRUE)
kaiju_rows <- nrow(otudata)
taxmat = matrix(sample(letters, kaiju_rows, replace = TRUE), nrow = nrow(otudata), ncol = 7)
rownames(taxmat) <- rownames(otudata)
colnames(taxmat) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
taxmat <- tax_table(taxmat)
phy.kaiju <- phyloseq(otudata, taxmat)

sample_names(phy.kaiju)

#Need to add in sample metadata: sample id, control vs sample, and DNA concentration
  #DNA concentration --> lib qubit values
conc$Sample.Name <- gsubfn("-", ".", conc$Sample.Name)

#rename Sample variable to sample_id to facilitate merging with metadata file
conc <- rename(conc, sample_id = 'Sample.Name')
conc <- rename(conc, qubit = 'Library.Qubit..ng.ul.')
conc$Sample <- NULL

#Add variable identifying true samples vs neg control
conc$class <- NA
conc$class [conc$sample_id == "ext_neg_ctrl"] <- "Control"
conc$class [conc$sample_id != "ext_neg_ctrl"] <- "Sample"
table(conc$class)

#Now merge with phy.kaiju
rownames(conc)<- conc$sample_id
sample_data(phy.kaiju) <- conc
head(sample_data(phy.kaiju))
sample_names(phy.kaiju)

# #REMOVE pos control; don't think will need for decontam
phy.kaiju = subset_samples(phy.kaiju, sample_id != "zymo_pos_ctrl")

#SEPT 2023: save phy.kaiju in case we need to run more contamination frequencies, etc
saveRDS(phy.kaiju, file = "/Users/swetapatel/OneDrive - Duke University/Fogarty coding/Sequencing/phy.kaiju_decon09132023.rds", ascii = FALSE, version = NULL,
        compress = FALSE, refhook = NULL)

#########
#DECONTAM
#########
#Source: https://benjjneb.github.io/decontam/vignettes/decontam_intro.html

#Step 1: inspect library sizes
df <- as.data.frame(sample_data(phy.kaiju)) # Put sample_data into a ggplot-friendly data.frame
df$LibrarySize <- sample_sums(phy.kaiju)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
ggplot(data=df, aes(x=Index, y=LibrarySize, color=class)) + geom_point()

#Step 2: identify contaminants with frequency approach
contamdf.freq <- isContaminant(phy.kaiju, method="frequency", conc="qubit")
  #As we did not specify the threshold, the default value of threshold = 0.1 was used
    #Thus $contaminant=TRUE if $p < 0.1
table(contamdf.freq$contaminant)  #53 sequences classified as contaminants
head(which(contamdf.freq$contaminant))  #499, 1017, 1068, 1162, 1483, 1527
#what species are classified as contaminants?
confreq <- subset(contamdf.freq, contamdf.freq$contaminant=="TRUE")
setDT(confreq, keep.rownames = "Species")
table(confreq$Species)
#how does a plot of the most common taxa compare to the 499th?
plot_frequency(phy.kaiju, taxa_names(phy.kaiju)[c(1,499)], conc="qubit") + 
  xlab("Paired reads")
  #Taxa 1 is def more frequent at higher concentrations, taxa 499 is less frequent at higher concentrations
#check 3 more random ASVs
set.seed(100)
plot_frequency(phy.kaiju, taxa_names(phy.kaiju)[sample(which(contamdf.freq$contaminant),3)], conc="qubit") +
  xlab("Paired reads")
#save confreq df as a csv file so we can compare the species identified compared to our original "either" approach
write.csv(confreq, "Kaiju_contam_freq_0.1.csv")

###LANCET MICROBE REVISIONS###
#Reviewer requesting the following: "I believe there is a way to generate percentage removed based upon filtering as a contaminant vs. not contaminant -
#[cont'd] I would like to see those bar plots. You can include them as supplemental figures."
#Does the reviewer want barplots by different frequency thresholds as in this paper? https://journals.asm.org/doi/10.1128/msystems.00290-19
#We can generate a single barplot for our frequency threshold of 0.1
contamdf.freq$thresh <- 0.1
bar_0.1 <- ggplot(contamdf.freq, aes(x=1, fill=contaminant)) + geom_bar(stat= "count", position = "fill") +
  scale_fill_manual(values = c('#800000FF','#155F83FF')) +
  # scale_y_continuous(labels = waiver()) +  
  # scale_x_discrete(labels=xlab) +
  labs(title=" ", x="Threshold = 0.1", y="Proportion of taxa (n=8265)", fill = "Contaminant")+
  theme_classic() +
  theme(
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title = element_text(size = 12))

#We can see how many taxa get kicked out with 0.3 and 0.5 and try to combine all of these into one single dataframe and barplot
contamdf.freq.3 <- isContaminant(phy.kaiju, method="frequency", conc="qubit", threshold=0.3)
table(contamdf.freq.3$contaminant) #208 contaminants here, compared to 53 for a threshold of 0.1
  #if reviewer is concerned about picking up rare taxa, then using a higher frequency is not the way to go.
contamdf.freq.3$thresh <- 0.3

#given the reviewer's concern, what if we use a threshold of 0.05?
contamdf.freq.05 <- isContaminant(phy.kaiju, method="frequency", conc="qubit", threshold=0.05)
table(contamdf.freq.05$contaminant) #29 contaminants here, compared to 53 for a threshold of 0.1
#if reviewer is concerned about picking up rare taxa, then using a higher frequency is not the way to go.
contamdf.freq.05$thresh <- 0.05

#make row names into a column, and then merge on that column (?); ideally will end up with 8265x3 = 24,795 rows
contamdf.freq$taxa <- rownames(contamdf.freq)
contamdf.freq.3$taxa <- rownames(contamdf.freq.3)
contamdf.freq.05$taxa <- rownames(contamdf.freq.05)
# test <- merge(contamdf.freq, contamdf.freq.05, by="taxa")

#tricky to merge because contaminant may be True for one threshold and false for another; we'll create 3 barplots and put in panel figure with cowplot
bar_0.3 <- ggplot(contamdf.freq.3, aes(x=1, fill=contaminant)) + geom_bar(stat= "count", position = "fill") +
  scale_fill_manual(values = c('#800000FF','#155F83FF')) +
  # scale_y_continuous(labels = waiver()) +  
  # scale_x_discrete(labels=xlab) +
  labs(title=" ", x="Threshold = 0.3", y="Proportion of taxa (n=8265)", fill = "Contaminant")+
  theme_classic() +
  theme(
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title = element_text(size = 12))

bar_0.05 <- ggplot(contamdf.freq.05, aes(x=1, fill=contaminant)) + geom_bar(stat= "count", position = "fill") +
  scale_fill_manual(values = c('#800000FF','#155F83FF')) +
  # scale_y_continuous(labels = waiver()) +  
  # scale_x_discrete(labels=xlab) +
  labs(title=" ", x="Threshold = 0.05", y="Proportion of taxa (n=8265)", fill = "Contaminant")+
  theme_classic() +
  theme(
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title = element_text(size = 12))

#now combine into one figure and save for supplemental file
library(cowplot)
bar_decontam <- plot_grid(bar_0.05, bar_0.1, bar_0.3, labels = "AUTO",  nrow = 1, align = "h", axis = "b")
bar_decontam

png(file="/Users/swetapatel/OneDrive - Duke University/Fogarty coding/decontam_barplots.png", width = 9, height = 5, units = 'in', res = 600)
bar_decontam
dev.off()

#Step 3: identify contaminants using prevalence approach [SKIPPING 6 MAY 2022]
sample_data(phy.kaiju)$is.neg <- sample_data(phy.kaiju)$class == "Control"
contamdf.prev <- isContaminant(phy.kaiju, method="prevalence", neg="is.neg")
table(contamdf.prev$contaminant)  #This says no samples are contaminants; big discrepancy
  #what if we use a threshold of 0.5?
  #(to identify as contaminants all sequences that are more prevalent in negative controls than in positive samples)
contamdf.prev05 <- isContaminant(phy.kaiju, method="prevalence", neg="is.neg", threshold=0.5)
table(contamdf.prev05$contaminant)  #identifies 1 contaminant
head(which(contamdf.prev05$contaminant))  #row 21 --> cutibacterium acnes (the only species identified in our neg control)

#According to tutorial and the ?isContaminant section in R, can use method = 'either' to
  #ID contaminants identified by either freq or prevalence; try it:
contamdf.eith <- isContaminant(phy.kaiju, method="either", conc="qubit", neg="is.neg")
table(contamdf.eith$contaminant)
#what species are classified as contaminants?
coneith <- subset(contamdf.eith, contamdf.eith$contaminant=="TRUE")
setDT(coneith, keep.rownames = "Species")
#save coneith df as a csv file so we can compare the species identified compared to our original "either" approach
write.csv(coneith, "Kaiju_contam_eith.csv")

#Directly compared the species listed in each csv file on 5/6
  #Because the same threshold value was used for freq and prev contam identification (0.1), the species contained in both CSV files are identical
  #C. acnes was only thrown out once we got to a prev threshold of 0.5 using the prev method above and was not thrown out in the either method
  #Thus, we do NOT need to make a new phyloseq object :)

#Now prune taxa; can prune all 53 species identified using frequency method
  #have 8265 taxa; after pruning contaminants, should end up with 8212 taxa
kaiju.noncontam <- prune_taxa(!contamdf.freq$contaminant, phy.kaiju)  #8212 taxa!

remove(conc, contamdf.freq, df, phy.kaiju)
######################################################
#PHYLOSEQ OBJECT CREATION: ADDING METADATA TO DECONTAM 
######################################################
###METADATA TABLE###
meta$X <- NULL
#standardize study_id to include "." instead of "-"
toString(meta$study_id)
meta$study_id <- gsub("-", ".", meta$study_id)

#now need to add sample_id variable to the metadata
#probably faster way to do this, but will make a dataframe that has our sample_ids and then merge that with metadata
test <- kaiju
colnames(test)

mat = matrix(ncol=0, nrow=309)
test2=data.frame(mat)

test2$sample_id <- colnames(test)
test2 <- as.data.frame(test2)
colnames(test2) <- c("sample_id")
remove(test)
test2$study_id <- test2$sample_id
toString(test2$study_id)
test2$study_id <- gsub(".MOM", "", test2$study_id)
test2$study_id <- gsub(".CHI", "", test2$study_id)
test2$study_id <- gsub(".SIB", "", test2$study_id)

test3 <- merge(meta, test2, by = 'study_id', all.x=TRUE)
meta2 <- test3

#make sample names row names 
rownames(meta2)<- meta2$sample_id

# #and delete the first column because it is now redundant
# meta2<- meta2 %>% 
#   select(-sample_id)

metadata_kaiju <- subset(meta2, sample_id %in% kaiju_ids)
rownames(metadata_kaiju) <- metadata_kaiju$sample_id
sample_data(kaiju.noncontam) <- metadata_kaiju
#success!
splitted_names <- strsplit(taxa_names(kaiju.noncontam),"\\;")
splitted_names_length <- lengths(splitted_names)
kaiju.noncontam <- prune_taxa(taxa_names(kaiju.noncontam)[splitted_names_length==7], kaiju.noncontam)
tax_dt <-  data.frame(name = gsub("[a-zA-Z]__","",taxa_names(kaiju.noncontam)))
suppressWarnings(taxmat <- separate(data = tax_dt, col = name, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = "\\;"))
rownames(taxmat) <- taxa_names(kaiju.noncontam)
tax_table(kaiju.noncontam) <- tax_table(as.matrix(taxmat))
taxa_names(kaiju.noncontam) <- as.character(tax_table(kaiju.noncontam)[,"Species"])
remove(splitted_names, splitted_names_length, tax_dt, taxmat, kaiju_ids, kaiju_rows)

#check if everything looks good
sample_sums(kaiju.noncontam) #count data!

sample_names(kaiju.noncontam) #sample names

rank_names(kaiju.noncontam)        #taxa levels  
## [1] "Kingdom"  "Phylum"  "Class"   "Order"   "Family"  "Genus"   "Species"

sample_variables(kaiju.noncontam)  #all of our metadata variables

otu_table(kaiju.noncontam)[1:3, 1:2]
taxa_names(kaiju.noncontam)[1:5]

#to start getting the feel of the data, let´s check what had the highest abundance
max(kaiju.noncontam@otu_table)    #811453

#Remove non-bacterial species from pseq object
nsamples(kaiju.noncontam)  #308 samples
ntaxa(kaiju.noncontam)     #8211 taxa
phy.fog <- subset_taxa(kaiju.noncontam, Kingdom=="Bacteria")
ntaxa(phy.fog)          #7998 bacterial taxa
summary(sample_sums(phy.fog))
# phy.fog <- transform_sample_counts(phy.fog, function(Species) Species/sum(Species))
# summary(sample_sums(phy.fog))   #all sums to 1 (presumably what we want)
saveRDS(phy.fog, file = "phy.kaiju.decontam.04192022.rds", ascii = FALSE, version = NULL,
        compress = FALSE, refhook = NULL)