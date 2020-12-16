#setwd("~/project/G-SAM")

# --- load packages ----
library(tidyverse)
library(dplyr)
library(matrixStats)
library(cluster)
library(dendextend)

# --- load data ----
source("scripts/R/gsam_metadata.R")
e.vst <- read.csv("data-local/vst_counts_gsam_paired.csv",row.names = c(1))
m <- gsam.rna.metadata
source("scripts/R/job_gg_theme.R")

#sum(grepl("1", colnames(e.vst)))
#sum(grepl("2", colnames(e.vst)))

# ---- PCA ----
# PCA + plot
ntop <- 500
variances <- rowVars(as.matrix(e.vst))
select <- order(variances, decreasing = TRUE)[seq_len(min(ntop, length(variances)))]
high_variance_genes <- e.vst[select,]

pc <- as.data.frame(prcomp(t(high_variance_genes))$x) %>% 
  dplyr::mutate(sample=rownames(.)) %>% 
  dplyr::left_join(.,m,by="sample") %>%
  dplyr::select(sample,PC1,PC2,pid,resection) 

ggplot(pc, aes(x=PC1, y=PC2, label=pid)) +
  geom_point(aes(col=resection)) +
  job_gg_theme
#ggsave("output/pca_resection_col.pdf",width=1106/72* 0.8,height=672/72 * 0.8)

