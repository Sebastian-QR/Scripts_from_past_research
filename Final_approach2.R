#### READ SUMMARIZARION ####
source("http://bioconductor.org/biocLite.R")
biocLite("Rsubread")

library(Rsubread)
aligned_files <- c("mapped_data_ensembl/1251.1G_ensembl_Aligned.out.bam",
                   "mapped_data_ensembl/1251.1S_ensembl_Aligned.out.bam",
                   "mapped_data_ensembl/1251.2G_ensembl_Aligned.out.bam",
                   "mapped_data_ensembl/1251.2S_ensembl_Aligned.out.bam",
                   "mapped_data_ensembl/658.1G_ensembl_Aligned.out.bam",
                   "mapped_data_ensembl/658.1S_ensembl_Aligned.out.bam",
                   "mapped_data_ensembl/1452.1G_ensembl_Aligned.out.bam",
                   "mapped_data_ensembl/1452.1S_ensembl_Aligned.out.bam",
                   "mapped_data_ensembl/1452.2G_ensembl_Aligned.out.bam",
                   "mapped_data_ensembl/1452.2S_ensembl_Aligned.out.bam")
data_new<-featureCounts(aligned_files,
                        annot.ext="ref_genome_ensembl/Sheep_genome_annotations_ensembl.gtf",
                        isGTFAnnotationFile=TRUE,
                        nthreads=16)

#Write table of the raw and filtered counts to a file, for backup.
colnames(data_new$counts) <- c("1251.1G","1251.1S","1251.2G","1251.2S","658.1G","658.1S","1452.1G","1452.1S","1452.2G","1452.2S")
write.csv(data_new$counts, file="DE_Analysis(new_approach)/counts_new_raw.csv")
write.csv(data_new$counts[rowSums(data_new$counts)>10,], file="DE_Analysis(new_approach)/counts_new.csv")

#Import list of candidate genes from file
candidate_genes<-read.table(file="DE_Analysis(new_approach)/Candidate_genes.txt", header=F, fill=T, col.names="Gene")




##### DE ANALYSIS #####

## NORMALIZATION STEPS edgeR ##

#Get rid of outliers
data_final <- as.data.frame(data_new$counts[,c(1:4,7:10)])

#Get rid of low count genes (less than 10 counts in each sample)
data_final<-subset(data_final, data_final$`1251.1G` >= 10)
data_final<-subset(data_final, data_final$`1251.1S` >= 10)
data_final<-subset(data_final, data_final$`1251.2G` >= 10)
data_final<-subset(data_final, data_final$`1251.2S` >= 10)
data_final<-subset(data_final, data_final$`1452.1G` >= 10)
data_final<-subset(data_final, data_final$`1452.1S` >= 10)
data_final<-subset(data_final, data_final$`1452.2G` >= 10)
data_final<-subset(data_final, data_final$`1452.2S` >= 10)
head(data_final)

#Import experiment matrix and create factors
experiment_design <- read.table("experiment_design_sq.txt",header = T, sep = "\t")
rownames(experiment_design) <- experiment_design$SampleID
samples_names <- as.character(experiment_design$SampleID)
sample_type <- factor(experiment_design$Sample_type)
sheep_number <- factor(experiment_design$Source_animal)

#Create DGEList object
library(edgeR)
DGEList.edgeR <- DGEList(counts=data_final, group = sample_type)
DGEList.edgeR <- calcNormFactors(DGEList.edgeR)
DGEList.edgeR <- estimateCommonDisp(DGEList.edgeR,verbose = T)
DGEList.edgeR <- estimateTagwiseDisp(DGEList.edgeR)
DGEList.tgw <- exactTest(DGEList.edgeR)

#Apply the exact test
DGEList.edgeR.tgw <- exactTest(DGEList.edgeR)
FDR.table <- topTags(DGEList.edgeR.tgw,n=11167)
FDR.table <- FDR.table$table
FDR <- decideTestsDGE(DGEList.tgw, adjust.method="BH", p.value=0.05)
summary(FDR)

#Apply unconditional exact test
library(edgeRun)
un.DGEList.tgw <- UCexactTest(DGEList.edgeR,upper=50000)
topTags(un.DGEList.tgw, n=20)

un.FDR.p_all <- un.DGEList.tgw$table
#subset genes obtained with unconditional test by pvalue
un.FDR.p_0.05 <- subset(un.FDR.p_all, un.FDR.p_all$PValue <= 0.05)
un.FDR.p_0.01 <- subset(un.FDR.p_all, un.FDR.p_all$PValue <= 0.01)

###Write results for both edgeR and edgeRun in table

#edgeR, p<0.05
write.table(normalized_significant, file = "DE_Analysis(new_approach)/edgeR_0.05.txt", 
            quote = F, sep = "\t", row.names = T)

#edgeR, p< 0.01
write.table(normalized_significant_0.01, file = "DE_Analysis(new_approach)/edgeR_0.01.txt", 
            quote = F, sep = "\t", row.names = T)

#edgeRun, p<0.05
write.table(un.FDR.p_0.05, file = "DE_Analysis(new_approach)/edgeRun_0.05.txt", 
            quote = F, sep = "\t", row.names = T)

#edgeRun, p<0.01
write.table(un.FDR.p_0.01, file = "DE_Analysis(new_approach)/edgeRun_0.01.txt", 
            quote = F, sep = "\t", row.names = T)


un.FDR.table <- topTags(un.DGEList.tgw,n=11167)
un.FDR.table <- un.FDR.table$table
un.FDR <- decideTestsDGE(un.DGEList.tgw, adjust.method="BH", p.value=0.05)
summary(un.FDR)

#Save data in file
write.table(DGEList.edgeR$pseudo.counts, file = "DE_Analysis(new_approach)/normalised_data_edgeR(pseudoCounts).txt", 
            quote = F, sep = "\t", row.names = T)



FDR.table.significant <- subset(FDR.table, FDR.table$PValue <= 0.05)

#Save data in file
write.table(DGEList.edgeR.tgw$table, file = "DE_Analysis(new_approach)/normalised_data_edgeR(DGEList.edgeR.tgw$table).txt", 
            quote = F, sep = "\t", row.names = T)

#Normalized data into new variable
normalized_data<-un.FDR.p_0.05


#Try to create a reference list from the annotation file including the names and symbols 
#to match with the list of interest afterwards
library(rtracklayer)
gtf <- import("ref_genome_ensembl/Sheep_genome_annotations_ensembl.gtf")
gtf2 <- subset(gtf, type == "gene")
gtf_df <- as.data.frame(gtf2)[,c("gene_id","gene_name")]

#merge gtf data frame and normalized data
normalized_data[,4]<-c(1:1963)
normalized_data[,4]<-rownames(normalized_data)

merged_normalized <- merge(normalized_data, gtf_df, by.x="V4", by.y="gene_id")
write.table(merged_normalized, file = "DE_Analysis(new_approach)/merged_p.0.05.txt", 
            quote = F, sep = "\t", row.names = T)
merged_normalized_final <- merge(merged_normalized, candidate_genes, by.x="V4", by.y="gene_id")
write.table(merged_normalized_final, file = "DE_Analysis(new_approach)/merged_p.0.05_inliterature.txt", 
            quote = F, sep = "\t", row.names = T)

