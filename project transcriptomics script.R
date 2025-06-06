library(Rsubread)
buildindex(
  basename = 'ref_human',
  reference = 'Homo_sapiens.GRCh38.dna.toplevel.fa.gz',
  memory = 4000,
  indexSplit = TRUE)
setwd("C:/Users/xelis/OneDrive/transcriptomics project")
getwd()
align.controle1 <- align(index = "ref_human", readfile1 = "SRR4785819_1_subset_controle1.fastq", readfile2 = "SRR4785819_2_subset_controle1.fastq", output_file = "controle1.BAM")
align.controle2 <- align(index = "ref_human", readfile1 = "SRR4785820_1_subset_controle2.fastq", readfile2 = "SRR4785820_2_subset_controle2.fastq", output_file = "controle2.BAM")
align.controle3 <- align(index = "ref_human", readfile1 = "SRR4785828_1_subset_controle3.fastq", readfile2 = "SRR4785828_2_subset_controle3.fastq", output_file = "controle3.BAM")
align.controle4 <- align(index = "ref_human", readfile1 = "SRR4785831_1_subset_controle4.fastq", readfile2 = "SRR4785831_2_subset_controle4.fastq", output_file = "controle4.BAM")
align.RA1 <- align(index = "ref_human", readfile1 = "SRR4785979_1_subset_RA1.fastq", readfile2 = "SRR4785979_2_subset_RA1.fastq", output_file = "RA1.BAM")
align.RA2 <- align(index = "ref_human", readfile1 = "SRR4785980_1_subset_RA2.fastq", readfile2 = "SRR4785980_2_subset_RA2.fastq", output_file = "RA2.BAM")
align.RA3 <- align(index = "ref_human", readfile1 = "SRR4785986_1_subset_RA3.fastq", readfile2 = "SRR4785986_2_subset_RA3.fastq", output_file = "RA3.BAM")
align.RA4 <- align(index = "ref_human", readfile1 = "SRR4785988_1_subset_RA4.fastq", readfile2 = "SRR4785988_2_subset_RA4.fastq", output_file = "RA4.BAM")
library(Rsamtools)
samples.project <- c('controle1', 'controle2', 'controle3', 'controle4', 'RA1', 'RA2', 'RA3', 'RA4')
lapply(samples.project, function(s) {sortBam(file = paste0(s, '.BAM'), destination = paste0(s, '.sorted'))})
lapply(samples.project, function(s) {indexBam(file = paste0(s, '.sorted.bam'))})

#count-matrix
library(readr)
library(dplyr)
library(Rsamtools)
library(Rsubread)
allsamples.project <- c("controle1.BAM", "controle2.BAM", "controle3.BAM", "controle4.BAM", "RA1.BAM", "RA2.BAM", "RA3.BAM", "RA4.BAM")
count_matrix.project <- featureCounts(
  files = allsamples.project,
  annot.ext = "Homo_sapiens.GRCh38.114.gtf.gz",
  isPairedEnd = TRUE,
  isGTFAnnotationFile = TRUE,
  GTF.attrType = "gene_id",
  useMetaFeatures = TRUE)
head(count_matrix.project$annotation)
head(count_matrix.project$counts)
View(count_matrix.project$counts)
str(count_matrix.project)
counts.project <- count_matrix.project$counts
colnames(counts.project) <- c("controle1", "controle2", "controle3", "controle4", "RA1", "RA2", "RA3", "RA4")
write.csv(counts.project, "bewerkt_countmatrix_project.csv")

#statistiek en analyse
counts.project <- read.table("count_matrix.txt", row.names = 1)
treatment.project <- c("controle", "controle", "controle", "controle", "RA", "RA", "RA", "RA")
treatment_table_project <- data.frame(treatment.project)
rownames(treatment_table_project) <- c('SRR4785819', 'SRR4785820', 'SRR4785828', 'SRR4785831', 'SRR4785979', 'SRR4785980', 'SRR4785986', 'SRR4785988')
dds.project <- DESeqDataSetFromMatrix(countData = round(counts.project),
                              colData = treatment_table_project,
                              design = ~ treatment.project)
dds.project <- DESeq(dds.project)
resultaten.project <- results(dds.project)
write.table(resultaten.project, file = 'ResultatenProject.csv', row.names = TRUE, col.names = TRUE)
sum(resultaten.project$padj < 0.05 & resultaten.project$log2FoldChange > 1, na.rm = TRUE)
sum(resultaten.project$padj < 0.05 & resultaten.project$log2FoldChange < -1, na.rm = TRUE)
hoogste_fold_change <- resultaten.project[order(resultaten.project$log2FoldChange, decreasing = TRUE), ]
laagste_fold_change <- resultaten.project[order(resultaten.project$log2FoldChange, decreasing = FALSE), ]
laagste_p_waarde <- resultaten.project[order(resultaten.project$padj, decreasing = FALSE), ]
head(laagste_p_waarde)
EnhancedVolcano(resultaten.project,
                lab = rownames(resultaten.project),
                x = 'log2FoldChange',
                y = 'padj')
dev.copy(png, 'VolcanoplotProject2.png', 
         width = 8,
         height = 10,
         units = 'in',
         res = 500)
dev.off()

#GO analyse
all <- rownames(resultaten.project)
resultaten.project <- as.data.frame(resultaten.project)
joo <- resultaten.project %>%
  filter(padj < 0.05) %>%
  filter(log2FoldChange >= 6 | log2FoldChange <= -6)
joo <- rownames (joo)
BiocManager::install("goseq")
BiocManager::install("geneLenDataBase")
BiocManager::install("org.Dm.eg.db")
library(goseq)
library(geneLenDataBase)
library(org.Dm.eg.db)
gene.vector=as.integer(all%in%joo)
names(gene.vector)=all
head(gene.vector)
tail(gene.vector)
gene.vector
pwf=nullp(gene.vector,"hg19", "geneSymbol")
GO.wall=goseq(pwf,"hg19","geneSymbol")
class(GO.wall)
head(GO.wall)
nrow(GO.wall)
enriched.GO=GO.wall$category[GO.wall$over_represented_pvalue<.05]
class(enriched.GO)
head(enriched.GO)
length(enriched.GO)
library(GO.db)
capture.output(for(go in enriched.GO[1:258]) { print(GOTERM[[go]])
  cat("--------------------------------------\n")
}
, file="SigGo.txt")
library(ggplot2)
library(tidyverse)
BiocManager::install("tidyverse")

# Filter top 15 significante oververtegenwoordigde GO-termen
GO.sig <- GO.wall %>%
  filter(over_represented_pvalue < 0.05) %>%
  arrange(over_represented_pvalue) %>%
  slice(1:15)

# Zet GO-termen op volgorde in de y-as
GO.sig$term <- factor(GO.sig$term, levels = rev(GO.sig$term))

ggplot(GO.sig, aes(x = -log10(over_represented_pvalue), y = term)) +
  geom_point(aes(size = numDEInCat, color = ontology), alpha = 0.8) +
  scale_size_continuous(name = "Aantal genen") +
  scale_color_brewer(palette = "Set1", name = "Ontologie") +
  labs(
    title = "Top 15 verrijkte GO-termen",
    x = expression(-log[10](p-waarde)),
    y = "GO-term"
  ) +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 10))

resultaten[1] <- NULL
resultaten[2:5] <- NULL
pathview(
  gene.data = resultaten.project,
  pathway.id = "hsa04612",  
  species = "hsa",          
  gene.idtype = "SYMBOL",     
  limit = list(gene = 5)    
)
