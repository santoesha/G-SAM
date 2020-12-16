#setwd("~/project/G-SAM")

# ---- load packages ----
library(dplyr)
library(tidyverse)
library(cowplot)

# ---- load data ----
source("scripts/R/gsam_metadata.R")

sign_genes <- read.delim("data-local/gene_signatures_ESTIMATE.txt")

counts <- read.csv("output/nc_gsam_deseq2.csv") %>%
              separate(X, c("ensembl", "geneID","chrom_loc"),sep="\\|") %>%
              dplyr::select(-c("ensembl","chrom_loc")) 

counts[which(counts$geneID=="CCN4"),"geneID"] <- "WISP1"
counts[which(counts$geneID=="ADGRA2"),"geneID"] <- "GPR124"
counts[which(counts$geneID=="PLPPR4"),"geneID"] <- "LPPR4"
counts[which(counts$geneID=="NME8"),"geneID"] <- "TXNDC3"
counts[which(counts$geneID=="TENM4"),"geneID"] <- "ODZ4"
counts[which(counts$geneID=="FYB1"),"geneID"] <- "FYB"

counts <- counts %>%
          dplyr::filter(geneID %in%sign_genes$V1) %>%
          column_to_rownames(var="geneID") %>%
          t() %>%
          as.data.frame(.)

#write.csv(counts,file="output/estimate_tumor_pur.csv")
estimate_score <- read.csv("data-local/GlioVis_tumor_purity.csv") 
estimate_score$donor_ID <- substr(estimate_score$Sample,1,4)
estimate_score$pid <- substr(estimate_score$Sample,1,3)
estimate_score$resection <-  ifelse(rownames(estimate_score)%in%grep('[1]',estimate_score$Sample),"r1","r2")

subtypes <- read.csv("data-local/GlioVis_subtypes.csv") 

# ---- ESTIMATE output analysis ----
m <- gsam.rna.metadata
lowq_samples <- m[m$blacklist.pca == T & m$batch == "single",2] %>% substr(1,3)
idh1_mutated <- gsam.cnv.metadata %>%
  dplyr::filter(IDH1 == "Stable" | IDH1 == "Lost" | IDH1 == "Gained") %>%
  dplyr::pull(pid)

estimate_score <- read.csv("data-local/GlioVis_tumor_purity.csv") %>%
  dplyr::filter(!Sample%in%c("EBP1","FAH2","KAE1.new","KAC1.new","AAC1","AAC2")) %>%
  dplyr::left_join(.,m[,c("pid","sample","resection")],by=c("Sample"="sample")) %>%
  dplyr::filter(!pid%in%lowq_samples) %>%
  dplyr::filter(!pid%in%idh1_mutated) %>%
  dplyr::left_join(.,subtypes[,c("Sample","majority.call")])

ggboxplot(estimate_score,
          x = "majority.call",
          y = "ESTIMATEScore",
          color = "resection",
          palette = "lancet",
          xlab = "Subtype",
          ylab = "Tumor Purity Score (ESTIMATE)",
          notch = T,
          add = "jitter",
          shape = "resection") +
  stat_compare_means(aes(group = resection))
#ggsave("output/tumor_purity_boxplot_subtypes.pdf")

ggboxplot(estimate_score,
          x = "majority.call",
          y = "ESTIMATEScore",
          color = "majority.call",
          palette = "lancet",
          xlab = "Subtype",
          ylab = "Tumor Purity Score (ESTIMATE)",
          notch = T,
          add = "jitter",
          shape = "majority.call") +
  stat_compare_means(comparisons=list(c("Classical", "Mesenchymal"),
                                      c("Classical", "Proneural"),
                                      c("Proneural", "Mesenchymal"))) +
  stat_compare_means(label.y = 36)
#ggsave("output/tumor_purity_boxplot_between_subtypes.pdf")

ggboxplot(estimate_score,
          x = "resection",
          y = "ESTIMATEScore",
          color = "resection",
          palette = "lancet",
          xlab = "Resection",
          ylab = "Tumor Purity Score (ESTIMATE)",
          notch = T,
          add = "jitter",
          shape = "resection",
) +
  stat_compare_means(method="t.test") 


#ggsave("output/tumor_purity_boxplot.pdf")
#wilcox.test(ESTIMATEScore ~ resection, data = estimate_score)

ggboxplot(estimate_score,
          x = "resection",
          y = "StromalScore",
          color = "resection",
          palette = "lancet",
          xlab = "Resection",
          ylab = "Stromal Score",
          notch = T,
          add = "jitter",
          shape = "resection",
) +
  stat_compare_means(method="wilcox.test") 
#ggsave("output/Stromal_Score_Boxplot.pdf",width=6,height=5)

ggboxplot(estimate_score,
          x = "resection",
          y = "ImmuneScore",
          color = "resection",
          palette = "lancet",
          xlab = "Resection",
          ylab = "Immune Score",
          notch = T,
          add = "jitter",
          shape = "resection",
) +
  stat_compare_means(method="wilcox.test") 
#ggsave("output/Immune_Score_Boxplot.pdf",width=6,height=5)


tmp <- read.delim('/home/yhoogstrate/projects/gsam/output/tables/cnv/tumor.percentage.estimate.txt',sep=" ")

