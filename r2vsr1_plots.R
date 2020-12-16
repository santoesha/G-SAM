#setwd("~/project/G-SAM")

# ---- load packages ----
library(broom)
library(MASS)
library(ggpubr)
library(pastecs)

# ---- load: data ----
source("scripts/R/job_gg_theme.R")
source("scripts/R/ligands.R")
source("scripts/R/gsam_metadata.R")
source("scripts/R/expression_matrix.R")
source("scripts/R/snv.R")
source("scripts/R/gsam_rnaseq_egfrviii_expression.R")

lowq_samples <- gsam.rna.metadata[gsam.rna.metadata$blacklist.pca == T & gsam.rna.metadata$batch == "single",2] %>%
                substr(1,3)

vst_counts <- read.csv("data-local/vst_counts_gsam.csv",row.names = c(1))

vst_counts_pairs <- vst_counts %>% 
                    dplyr::select(-c("EBP1","FAH2","KAC1.new","KAE1.new","AAC1","AAC2")) %>% 
                    dplyr::select(-all_of(colnames(.)[substr(colnames(.),1,3)%in%lowq_samples])) 

ligands <- c("BTC","EREG","AREG","HBEGF","EPGN","TGFA","EGF")

dge_deseq <- read.csv("output/DESeq2_Results_r1vsr2.csv")

# ---- EGFRvIII ----
gsam.viii.rnaseq$vIII_mutated <- ifelse(gsam.viii.rnaseq$egfrviii.pct>10,"Yes","No")
gsam.viii.rnaseq$sample <- gsub("-",".",gsam.viii.rnaseq$sample)
gsam.viii.rnaseq$vIII_mutated[is.na(gsam.viii.rnaseq$vIII_mutated)] <- "No"

# ---- EGFR SNV ----
vaf <- vaf[grepl("R108|A289|G598",vaf$egfr_snv),]

# ---- MEOX2 ----
m <- inner_join(gsam.rna.metadata,gsam.patient.metadata,by=c("pid"="studyID"))

MEOX2_exp <- as.data.frame(t(vst_counts[grep("MEOX2",rownames(vst_counts)),])) %>%
             dplyr::mutate(sample=rownames(.)) %>%
             dplyr::inner_join(.,m[,c("sample","resection","pid","EGFR")],by=c("sample"="sample")) 
  
colnames(MEOX2_exp)[1] <- "MEOX2"

tmp1 <- MEOX2_exp %>% 
       dplyr::group_by(pid) %>% 
       dplyr::mutate(foldChange = log(MEOX2/first(MEOX2))) %>%
       dplyr::filter(resection=="r2")
tmp1$x <- order(tmp1$foldChange)
tmp1 <- tmp1[tmp1$x,]
tmp1$x <- order(tmp1$foldChange)

# pvalue_meox2 <- dge_deseq %>% 
#                 dplyr::filter(X=="MEOX2") %>% 
#                 dplyr::pull(padj) %>% 
#                 round(.,3)

ggplot(tmp1, aes(x=x, y=foldChange, label=sample)) + 
  geom_bar(
    data = subset(tmp1, !is.na(tmp1$foldChange)),
    aes(x=x, y=foldChange),stat="identity",width=0.7) + 
  geom_point(
    data = subset(tmp1, !is.na(tmp1$foldChange)), aes(col=EGFR)
  ) +
  #annotate("text", x=0, y=0.3, label="padj =") +
  #annotate("text", x=7.5, y=0.3, label=pvalue_meox2) +
  labs(y = "MEOX2 log fold change (r2 vs  r1)",
       x = "GSAM patient (ordered on difference in RNA-seq)") + 
  job_gg_theme 
# ---- SOCS2 ----
m <- inner_join(gsam.rna.metadata,gsam.patient.metadata,by=c("pid"="studyID"))

SOCS2_exp <- as.data.frame(t(vst_counts[grep("ENSG00000120833.14_7",rownames(vst_counts)),])) %>%
  dplyr::mutate(sample=rownames(.)) %>%
  dplyr::inner_join(.,m[,c("sample","resection","pid","EGFR")],by=c("sample"="sample"))
  
colnames(SOCS2_exp)[1] <- "SOCS2"

tmp1 <- SOCS2_exp %>% 
  dplyr::group_by(pid) %>% 
  dplyr::mutate(foldChange = log(SOCS2/first(SOCS2))) %>%
  dplyr::filter(resection=="r2")
tmp1$x <- order(tmp1$foldChange)
tmp1 <- tmp1[tmp1$x,]
tmp1$x <- order(tmp1$foldChange)

#pvalue_socs2 <- dge_deseq %>% 
  #dplyr::filter(X=="SOCS2") %>% 
  #dplyr::pull(padj) %>% 
  #round(.,6)

ggplot(tmp1, aes(x=x, y=foldChange, label=sample)) + 
  geom_bar(
    data = subset(tmp1, !is.na(tmp1$foldChange)),
    aes(x=x, y=foldChange),stat="identity",width=0.7) + 
  geom_point(
    data = subset(tmp1, !is.na(tmp1$foldChange)), aes(col=EGFR)
  ) +
  #annotate("text", x=0, y=0.25, label="padj =") +
  #annotate("text", x=8, y=0.25, label=pvalue_socs2) +
  labs(y = "SOCS2 log fold change (r2 vs  r1)",
       x = "GSAM patient (ordered on difference in RNA-seq)") + 
  job_gg_theme 



# ---- CDK4 ----
m <- inner_join(gsam.rna.metadata,gsam.patient.metadata,by=c("pid"="studyID"))

CDK4_exp <- as.data.frame(t(vst_counts[grep("ENSG00000135446.17_6",rownames(vst_counts)),])) %>%
  dplyr::mutate(sample=rownames(.)) %>%
  dplyr::inner_join(.,m[,c("sample","resection","pid","EGFR")],by=c("sample"="sample"))

colnames(CDK4_exp)[1] <- "CDK4"

tmp1 <- CDK4_exp %>% 
  dplyr::group_by(pid) %>% 
  dplyr::mutate(foldChange = log(CDK4/first(CDK4))) %>%
  dplyr::filter(resection=="r2")
tmp1$x <- order(tmp1$foldChange)
tmp1 <- tmp1[tmp1$x,]
tmp1$x <- order(tmp1$foldChange)

pvalue_cdk4 <- dge_deseq %>% 
  dplyr::filter(X=="CDK4") %>% 
  dplyr::pull(padj) %>% 
  round(.,5)

ggplot(tmp1, aes(x=x, y=foldChange, label=sample)) + 
  geom_bar(
    data = subset(tmp1, !is.na(tmp1$foldChange)),
    aes(x=x, y=foldChange),stat="identity",width=0.7) + 
  geom_point(
    data = subset(tmp1, !is.na(tmp1$foldChange)), aes(col=EGFR)
  ) +
  annotate("text", x=0, y=0.18, label="padj =") +
  annotate("text", x=9, y=0.18, label=pvalue_cdk4) +
  labs(y = "CDK4 log fold change (r2 vs  r1)",
       x = "GSAM patient (ordered on difference in RNA-seq)") + 
  job_gg_theme 
# ---- EGFR ----
m <- inner_join(gsam.rna.metadata,gsam.patient.metadata,by=c("pid"="studyID"))

EGFR_exp <- as.data.frame(t(vst_counts[grep("ENSG00000146648.18_5",rownames(vst_counts)),])) %>%
  dplyr::mutate(sample=rownames(.)) %>%
  dplyr::inner_join(.,m[,c("sample","resection","pid","EGFR")],by=c("sample"="sample"))

colnames(EGFR_exp)[1] <- "EGFR_ex"

tmp1 <- EGFR_exp %>% 
  dplyr::group_by(pid) %>% 
  dplyr::mutate(foldChange = log(EGFR/first(EGFR))) %>%
  dplyr::filter(resection=="r2")
tmp1$x <- order(tmp1$foldChange)
tmp1 <- tmp1[tmp1$x,]
tmp1$x <- order(tmp1$foldChange)

pvalue_egfr <- dge_deseq %>% 
  dplyr::filter(X=="EGFR") %>% 
  dplyr::pull(padj) %>% 
  round(.,4)

ggplot(tmp1, aes(x=x, y=foldChange, label=sample)) + 
  geom_bar(
    data = subset(tmp1, !is.na(tmp1$foldChange)),
    aes(x=x, y=foldChange),stat="identity",width=0.7) + 
  geom_point(
    data = subset(tmp1, !is.na(tmp1$foldChange)), aes(col=EGFR)
  ) +
  annotate("text", x=0, y=0.6, label="padj =") +
  annotate("text", x=8, y=0.6, label=pvalue_egfr) +
  labs(y = "EGFR log fold change (r2 vs  r1)",
       x = "GSAM patient (ordered on difference in RNA-seq)") +
  job_gg_theme 

# ---- EGF ----
m <- inner_join(gsam.rna.metadata,gsam.patient.metadata,by=c("pid"="studyID"))

EGF_exp <- as.data.frame(t(vst_counts_pairs[grep("ENSG00000138798.12_4",rownames(vst_counts_pairs)),])) %>%
  dplyr::mutate(sample=rownames(.)) %>%
  dplyr::inner_join(.,m[,c("sample","resection","pid","EGFR")],by=c("sample"="sample"))

colnames(EGF_exp)[1] <- "EGF"

tmp1 <- EGF_exp %>% 
  dplyr::group_by(pid) %>% 
  dplyr::mutate(foldChange = log(EGF/first(EGF))) %>%
  dplyr::filter(resection=="r2")
tmp1$x <- order(tmp1$foldChange)
tmp1 <- tmp1[tmp1$x,]
tmp1$x <- order(tmp1$foldChange)

pvalue_egf <- dge_deseq %>% 
  dplyr::filter(X=="EGF") %>% 
  dplyr::pull(padj) %>% 
  round(.,2)

ggplot(tmp1, aes(x=x, y=foldChange, label=sample)) + 
  geom_bar(
    data = subset(tmp1, !is.na(tmp1$foldChange)),
    aes(x=x, y=foldChange),stat="identity",width=0.7) + 
  geom_point(
    data = subset(tmp1, !is.na(tmp1$foldChange)), aes(col=EGFR)
  ) +
  labs(y = "EGF log fold change (r2 vs  r1)",
       x = "GSAM patient (ordered on difference in RNA-seq)") +
  #annotate("text", x=0, y=0.25, label="padj =") +
  #annotate("text", x=6.5, y=0.25, label=pvalue_egf) +
  job_gg_theme 

# ---- BTC ----
m <- inner_join(gsam.rna.metadata,gsam.patient.metadata,by=c("pid"="studyID"))

BTC_exp <- as.data.frame(t(vst_counts_pairs[grep("BTC",rownames(vst_counts_pairs)),])) %>%
  dplyr::mutate(sample=rownames(.)) %>%
  dplyr::inner_join(.,m[,c("sample","resection","pid","EGFR")],by=c("sample"="sample"))

colnames(BTC_exp)[1] <- "BTC"

tmp1 <- BTC_exp %>% 
  dplyr::group_by(pid) %>% 
  dplyr::mutate(foldChange = log(BTC/first(BTC))) %>%
  dplyr::filter(resection=="r2")
tmp1$x <- order(tmp1$foldChange)
tmp1 <- tmp1[tmp1$x,]
tmp1$x <- order(tmp1$foldChange)

pvalue_btc <- dge_deseq %>% 
              dplyr::filter(X=="BTC") %>% 
              dplyr::pull(padj) %>% 
              round(.,2)

ggplot(tmp1, aes(x=x, y=foldChange, label=sample)) + 
  geom_bar(
    data = subset(tmp1, !is.na(tmp1$foldChange)),
    aes(x=x, y=foldChange),stat="identity",width=0.7) + 
  geom_point(
    data = subset(tmp1, !is.na(tmp1$foldChange)), aes(col=EGFR)
  ) +
  annotate("text", x=0, y=0.40, label="padj =") +
  annotate("text", x=7, y=0.40, label=pvalue_btc) +
  labs(y = "BTC log fold change (r2 vs  r1)",
       x = "GSAM patient (ordered on difference in RNA-seq)") +
  job_gg_theme 

#ggsave("output/BTC_r1vsr2.pdf")

# ---- TGFA ----
m <- inner_join(gsam.rna.metadata,gsam.patient.metadata,by=c("pid"="studyID"))

TGFA_exp <- as.data.frame(t(vst_counts[grep("ENSG00000163235.16_5",rownames(vst_counts)),])) %>%
  dplyr::mutate(sample=rownames(.)) %>%
  dplyr::inner_join(.,m[,c("sample","resection","pid","EGFR")],by=c("sample"="sample"))

colnames(TGFA_exp)[1] <- "TGFA"

tmp1 <- TGFA_exp %>% 
  dplyr::group_by(pid) %>% 
  dplyr::mutate(foldChange = log(TGFA/first(TGFA))) %>%
  dplyr::filter(resection=="r2")
tmp1$x <- order(tmp1$foldChange)
tmp1 <- tmp1[tmp1$x,]
tmp1$x <- order(tmp1$foldChange)

# pvalue_tgfa <- dge_deseq %>% 
#   dplyr::filter(X=="TGFA") %>% 
#   dplyr::pull(padj) %>% 
#   round(.,8)

ggplot(tmp1, aes(x=x, y=foldChange, label=sample)) + 
  geom_bar(
    data = subset(tmp1, !is.na(tmp1$foldChange)),
    aes(x=x, y=foldChange),stat="identity",width=0.7) + 
  geom_point(
    data = subset(tmp1, !is.na(tmp1$foldChange)), aes(col=EGFR)
  ) +
  # annotate("text", x=0, y=0.45, label="padj =") +
  # annotate("text", x=9, y=0.45, label=pvalue_tgfa) +
  labs(y = "TGFA log fold change (r2 vs  r1)",
       x = "GSAM patient (ordered on difference in RNA-seq)") +
  job_gg_theme 

# ---- EREG ----
m <- inner_join(gsam.rna.metadata,gsam.patient.metadata,by=c("pid"="studyID"))

EREG_exp <- as.data.frame(t(vst_counts[grep("EREG",rownames(vst_counts)),])) %>%
  dplyr::mutate(sample=rownames(.)) %>%
  dplyr::inner_join(.,m[,c("sample","resection","pid","EGFR")],by=c("sample"="sample"))

colnames(EREG_exp)[1] <- "EREG"

tmp1 <- EREG_exp %>% 
  dplyr::group_by(pid) %>% 
  dplyr::mutate(foldChange = log(EREG/first(EREG))) %>%
  dplyr::filter(resection=="r2")
tmp1$x <- order(tmp1$foldChange)
tmp1 <- tmp1[tmp1$x,]
tmp1$x <- order(tmp1$foldChange)

# pvalue_ereg <- dge_deseq %>% 
#   dplyr::filter(X=="EREG") %>% 
#   dplyr::pull(padj) %>% 
#   round(.,2)

ggplot(tmp1, aes(x=x, y=foldChange, label=sample)) + 
  geom_bar(
    data = subset(tmp1, !is.na(tmp1$foldChange)),
    aes(x=x, y=foldChange),stat="identity",width=0.7) + 
  geom_point(
    data = subset(tmp1, !is.na(tmp1$foldChange)), aes(col=EGFR)
  ) +
  # annotate("text", x=0, y=0.80, label="padj =") +
  # annotate("text", x=6.8, y=0.80, label=pvalue_ereg) +
  labs(y = "EREG log fold change (r2 vs  r1)",
       x = "GSAM patient (ordered on difference in RNA-seq)") +
  job_gg_theme 

# ---- AREG ----
m <- inner_join(gsam.rna.metadata,gsam.patient.metadata,by=c("pid"="studyID"))

AREG_exp <- as.data.frame(t(vst_counts[grep("ENSG00000109321.11_4",rownames(vst_counts)),])) %>%
  dplyr::mutate(sample=rownames(.)) %>%
  dplyr::inner_join(.,m[,c("sample","resection","pid","EGFR")],by=c("sample"="sample"))

colnames(AREG_exp)[1] <- "AREG"

tmp1 <- AREG_exp %>% 
  dplyr::group_by(pid) %>% 
  dplyr::mutate(foldChange = log(AREG/first(AREG))) %>%
  dplyr::filter(resection=="r2")
tmp1$x <- order(tmp1$foldChange)
tmp1 <- tmp1[tmp1$x,]
tmp1$x <- order(tmp1$foldChange)

pvalue_areg <- dge_deseq %>% 
  dplyr::filter(X=="AREG") %>% 
  dplyr::pull(padj) %>% 
  round(.,4)

ggplot(tmp1, aes(x=x, y=foldChange, label=sample)) + 
  geom_bar(
    data = subset(tmp1, !is.na(tmp1$foldChange)),
    aes(x=x, y=foldChange),stat="identity",width=0.7) + 
  geom_point(
    data = subset(tmp1, !is.na(tmp1$foldChange)), aes(col=EGFR)
  ) +
  # annotate("text", x=0, y=0.70, label="padj =") +
  # annotate("text", x=8, y=0.70, label=pvalue_areg) +
  labs(y = "AREG log fold change (r2 vs  r1)",
       x = "GSAM patient (ordered on difference in RNA-seq)") +
  job_gg_theme 

# ---- HBEGF ----
m <- inner_join(gsam.rna.metadata,gsam.patient.metadata,by=c("pid"="studyID"))

HBEGF_exp <- as.data.frame(t(vst_counts[grep("HBEGF",rownames(vst_counts)),])) %>%
  dplyr::mutate(sample=rownames(.)) %>%
  dplyr::inner_join(.,m[,c("sample","resection","pid","EGFR")],by=c("sample"="sample"))

colnames(HBEGF_exp)[1] <- "HBEGF"

tmp1 <- HBEGF_exp %>% 
  dplyr::group_by(pid) %>% 
  dplyr::mutate(foldChange = log(HBEGF/first(HBEGF))) %>%
  dplyr::filter(resection=="r2")
tmp1$x <- order(tmp1$foldChange)
tmp1 <- tmp1[tmp1$x,]
tmp1$x <- order(tmp1$foldChange)

# pvalue_hbegf <- dge_deseq %>% 
#   dplyr::filter(X=="HBEGF") %>% 
#   dplyr::pull(padj) %>% 
#   round(.,2)

ggplot(tmp1, aes(x=x, y=foldChange, label=sample)) + 
  geom_bar(
    data = subset(tmp1, !is.na(tmp1$foldChange)),
    aes(x=x, y=foldChange),stat="identity",width=0.7) + 
  geom_point(
    data = subset(tmp1, !is.na(tmp1$foldChange)), aes(col=EGFR)
  ) +
  # annotate("text", x=0, y=0.33, label="padj =") +
  # annotate("text", x=7, y=0.33, label=pvalue_hbegf) +
  labs(y = "HBEGF log fold change (r2 vs  r1)",
       x = "GSAM patient (ordered on difference in RNA-seq)") +
  job_gg_theme 

# ---- Ligands pairs ----
BTC_exp_paired <- data.frame(t(vst_counts_pairs[grep("BTC",rownames(vst_counts_pairs)),])) %>%
  dplyr::mutate(sample=rownames(.)) %>%
  dplyr::inner_join(.,m[,c("sample","resection","pid","EGFR")],by=c("sample"="sample")) %>%
  dplyr::rename(BTC=colnames(.)[1]) %>%
  dplyr::group_by(pid) %>%
  dplyr::mutate(foldChange_BTC = log(BTC/first(BTC)))

EREG_exp_paired <- data.frame(t(vst_counts_pairs[grep("EREG",rownames(vst_counts_pairs)),])) %>%
  dplyr::mutate(sample=rownames(.)) %>%
  dplyr::inner_join(.,m[,c("sample","resection","pid","EGFR")],by=c("sample"="sample")) %>%
  dplyr::rename(EREG=colnames(.)[1]) %>%
  dplyr::group_by(pid) %>%
  dplyr::mutate(foldChange_EREG = log(EREG/first(EREG)))

AREG_exp_paired <- data.frame(t(vst_counts_pairs[grep("ENSG00000109321.11_4",rownames(vst_counts_pairs)),])) %>%
  dplyr::mutate(sample=rownames(.)) %>%
  dplyr::inner_join(.,m[,c("sample","resection","pid","EGFR")],by=c("sample"="sample")) %>%
  dplyr::rename(AREG=colnames(.)[1]) %>%
  dplyr::group_by(pid) %>%
  dplyr::mutate(foldChange_AREG = log(AREG/first(AREG)))

HBEGF_exp_paired <- data.frame(t(vst_counts_pairs[grep("HBEGF",rownames(vst_counts_pairs)),])) %>%
  dplyr::mutate(sample=rownames(.)) %>%
  dplyr::inner_join(.,m[,c("sample","resection","pid","EGFR")],by=c("sample"="sample")) %>%
  dplyr::rename(HBEGF=colnames(.)[1]) %>%
  dplyr::group_by(pid) %>%
  dplyr::mutate(foldChange_HBEGF = log(HBEGF/first(HBEGF)))

TGFA_exp_paired <- data.frame(t(vst_counts_pairs[grep("ENSG00000163235.16_5",rownames(vst_counts_pairs)),])) %>%
  dplyr::mutate(sample=rownames(.)) %>%
  dplyr::inner_join(.,m[,c("sample","resection","pid","EGFR")],by=c("sample"="sample")) %>%
  dplyr::rename(TGFA=colnames(.)[1]) %>%
  dplyr::group_by(pid) %>%
  dplyr::mutate(foldChange_TGFA = log(TGFA/first(TGFA)))

EGF_exp_paired <- data.frame(t(vst_counts_pairs[grep("ENSG00000138798.12_4",rownames(vst_counts_pairs)),])) %>%
  dplyr::mutate(sample=rownames(.)) %>%
  dplyr::inner_join(.,m[,c("sample","resection","pid","EGFR")],by=c("sample"="sample")) %>%
  dplyr::rename(EGF=colnames(.)[1]) %>%
  dplyr::group_by(pid) %>%
  dplyr::mutate(foldChange_EGF = log(EGF/first(EGF)))

EGFR_exp_paired <- data.frame(t(vst_counts_pairs[grep("ENSG00000146648.18_5",rownames(vst_counts_pairs)),])) %>%
  dplyr::mutate(sample=rownames(.)) %>%
  dplyr::inner_join(.,m[,c("sample","resection","pid","EGFR")],by=c("sample"="sample")) %>%
  dplyr::rename(EGFR_exp=colnames(.)[1]) %>%
  dplyr::group_by(pid) %>%
  dplyr::mutate(foldChange_EGFR = log(EGFR_exp/first(EGFR_exp)))

# ---- Ligands combined ----
#Paired plot
ligands_combined <- EGFR_exp_paired %>%
                    dplyr::inner_join(.,EREG_exp_paired[,-c(3:5)],by=c("sample"="sample")) %>%
                    dplyr::inner_join(.,BTC_exp_paired[,-c(3:5)],by=c("sample"="sample")) %>%
                    dplyr::inner_join(.,HBEGF_exp_paired[,-c(3:5)],by=c("sample"="sample")) %>%
                    dplyr::inner_join(.,TGFA_exp_paired[,-c(3:5)],by=c("sample"="sample")) %>%
                    dplyr::inner_join(.,EGF_exp_paired[,-c(3:5)],by=c("sample"="sample")) %>%
                    dplyr::inner_join(.,AREG_exp_paired[,-c(3:5)],by=c("sample"="sample")) %>%
                    dplyr::filter(resection=="r2")
                    
ggplot(ligands_combined,aes(x=reorder(sample,foldChange_EGFR),y=foldChange_EGFR)) +
  job_gg_theme +
  geom_point(pch=20) +
  geom_point(mapping = aes(x = reorder(sample, foldChange_EGFR), y=foldChange_BTC,col="BTC")) +
  geom_point(mapping = aes(x = reorder(sample, foldChange_EGFR), y=foldChange_EREG,col="EREG")) +
  geom_point(mapping = aes(x = reorder(sample, foldChange_EGFR), y=foldChange_AREG,col="AREG")) +
  geom_point(mapping = aes(x = reorder(sample, foldChange_EGFR), y=foldChange_HBEGF,col="HBEGF")) +
  geom_point(mapping = aes(x = reorder(sample, foldChange_EGFR), y=foldChange_TGFA,col="TGFA")) +
  geom_point(mapping = aes(x = reorder(sample, foldChange_EGFR), y=foldChange_EGF,col="EGF")) 

model <- lm(foldChange_EGFR ~ foldChange_BTC + foldChange_EREG:foldChange_AREG + foldChange_HBEGF + foldChange_EGF + foldChange_TGFA,data=ligands_combined)
varexplained <- model %>% tidy()
gvmodel <- gvlma(model) 
summary(gvmodel)

# ---- r2 vs r1: EGFRvIII and EGFR SNV ---- 
egfrviii <- gsam.viii.rnaseq %>%
            dplyr::mutate(pid = substr(sample,1,3)) %>%
            dplyr::filter(!is.na(egfrviii.pct))

duplicated <- egfrviii %>% filter(duplicated(pid)) %>% pull(pid)

egfrviii <- egfrviii %>%
            dplyr::filter(pid %in% duplicated) %>%
            dplyr::filter(!pid == "AAC") %>%
            dplyr::mutate(resection = ifelse(substr(sample,4,4)==1,"r1","r2")) %>%
            dplyr::mutate(egfrviii.pct.r1 = ifelse(resection == "r1",egfrviii.pct,NA)) %>%
            dplyr::mutate(egfrviii.pct.r2 = ifelse(resection == "r2",egfrviii.pct,NA))

#length(egfrviii$sample %in% available_vcf_data) --> vcf files lijken compleet?

x <- egfrviii %>%
     dplyr::select(egfrviii.pct.r1,pid) %>%
     dplyr::filter(!is.na(egfrviii.pct.r1)) 

y <- egfrviii %>%
     dplyr::select(egfrviii.pct.r2,pid) %>%
     dplyr::filter(!is.na(egfrviii.pct.r2))

plot_viii <- merge(x,y,by="pid")
plot_viii <- plot_viii %>%
             dplyr::left_join(.,vaf[,c("pid","trans")],by="pid")

plot_viii[is.na(plot_viii$trans),"trans"] <- "none"

ggplot(plot_viii, aes(x=egfrviii.pct.r1,y=egfrviii.pct.r2)) +
  geom_smooth(na.rm=TRUE,
              method="lm",
              se=F) +
  geom_point(na.rm = T) +
  geom_point(aes(col=trans)) +
  xlim(0, 100) +
  ylim(0, 100) +
  labs(x = "% EGFRvIII: resection 1",
       y = "% EGFRvIII: resection 2") +
  job_gg_theme +
  labs(color='Change in EGFR SNV status') 

cor.test(y$egfrviii.pct.r2,x$egfrviii.pct.r1,method="spearman")

tmp <- egfrviii %>%
       dplyr::group_by(pid) %>%
       dplyr::mutate(stable = ifelse(sum(vIII_mutated == "Yes") == 2,"Stable","")) 
  
#ggsave("output/egfr_viii_percentage_r1vs42.pdf")

# ---- r2 vs r1: EGFR amplification ----
egfr_ampli <- gsam.cnv.metadata %>%
              dplyr::select(pid,donor_ID,EGFR,cnEGFRp,cnEGFRr,cnStatusEGFRs) %>%
              dplyr::filter(!pid%in%c("AAC","BAP","CDG","FAE")) %>% #NA values / AAC normal cells were sequenced
              dplyr::mutate(resection = ifelse(substr(donor_ID,4,4)==1,"r1","r2"))

d <- egfr_ampli %>% 
     dplyr::filter(resection == "r1") 

wilcox.test(d$cnEGFRp,d$cnEGFRr,paired=T)

# ---- r2 vs r1: EGFR expression ----
m <- gsam.rna.metadata
lowq_samples <- m[m$blacklist.pca == T & m$batch == "single",2] %>% substr(1,3) #RNA blacklist, remove matching pair 

egfr_expr <- vst_counts[grepl("ENSG00000146648",rownames(vst_counts)),] %>%
             t(.) %>%
             as.data.frame(.) %>%
             dplyr::mutate(sample = rownames(.)) %>%
             dplyr::filter(!substr(sample,1,3)%in%lowq_samples) %>%
             dplyr::filter(!sample %in% c("EBP1","FAH2","KAE1.new","KAC1.new","AAC1","AAC2"))

colnames(egfr_expr)[1] <- "egfr_vst"
egfr_expr$sample <- substr(egfr_expr$sample,1,4)

egfr_vst <- dplyr::inner_join(egfr_ampli,egfr_expr,by=c("donor_ID"="sample")) 

#wilcox.test(egfr_vst ~ resection, data=egfr_vst)

egfr_vst_amplified_stable <- egfr_vst %>%
                             dplyr::filter(cnStatusEGFRs %in% c("Stable"))

#wilcox.test(egfr_vst ~ resection, data=egfr_vst_amplified_stable)

x <- egfr_expr %>%
  dplyr::filter(resection == "r1") %>%
  dplyr::select(egfr_vst,pid) %>%
  dplyr::filter(!is.na(egfr_vst)) %>%
  dplyr::filter(!pid %in% c("ECH","GAL","CBE")) %>% #geen DNA data beschikbaar van sommige samples, dus paren missen 
  distinct()

y <- egfr_vst %>%
  dplyr::filter(resection == "r2") %>%
  dplyr::select(egfr_vst,pid) %>%
  dplyr::filter(!is.na(egfr_vst)) %>%
  dplyr::filter(!pid %in% c("ECH","GAL","CBE")) %>% #geen DNA data beschikbaar van sommige samples, dus paren missen
  distinct()

plot_egfr <- inner_join(x,y,by="pid")
plot_egfr <- plot_egfr %>%
             dplyr::left_join(.,vaf[,c("pid","trans")],by="pid") 

plot_egfr[is.na(plot_egfr$trans),"trans"] <- "none"

ggplot(plot_egfr, aes(x=egfr_vst.x,y=egfr_vst.y)) +
  geom_smooth(na.rm=TRUE,
              method="lm",
              se=F) +
  geom_point(na.rm = T) +
  geom_point(aes(col=trans)) +
  xlim(5, 18) +
  ylim(5, 20) +
  labs(x = "EGFR VST expression: resection 1",
       y = "EGFR VST expression: resection 2") +
  job_gg_theme +
  labs(color='Change in EGFR SNV status') +
  stat_regline_equation(label.y = 18, aes(label = ..rr.label..))
#ggsave("output/EGFR_vst_r1vsr2_ECSNV.pdf")

# ---- r2 vs r1: BTC ----
x <- BTC_exp %>%
  dplyr::filter(resection == "r1") %>%
  dplyr::mutate(pid = substr(sample,1,3)) %>%
  dplyr::select(BTC,pid) 

y <- BTC_exp %>%
  dplyr::filter(resection == "r2") %>%
  dplyr::mutate(pid = substr(sample,1,3)) %>%
  dplyr::select(BTC,pid) 

plot_BTC_exp <- inner_join(x,y,by="pid")
plot_BTC_exp <- plot_BTC_exp %>%
                dplyr::left_join(.,vaf[,c("pid","trans")],by="pid") 

plot_BTC_exp[is.na(plot_BTC_exp$trans),"trans"] <- "none"

ggplot(plot_BTC_exp, aes(x=BTC.x,y=BTC.y)) +
  geom_smooth(na.rm=TRUE,
              method="lm",
              se=F) +
  geom_point(na.rm = T) +
  geom_point(aes(col=trans)) +
  xlim(0, 10) +
  ylim(0, 10) +
  labs(x = "EGFR VST expression: resection 1",
       y = "EGFR VST expression: resection 2") +
  job_gg_theme +
  labs(color='Change in EGFR SNV status') +
  stat_regline_equation(label.y = 8, aes(label = ..rr.label..))