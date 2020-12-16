# ---- initialization setup ----

#setwd("~/project/G-SAM")

# ---- load packages ----
library(EnhancedVolcano)
library(DESeq2)
library(ggplot2)
library(tidyverse)
library(gprofiler2)
library(pheatmap)
library(RColorBrewer)
library(gplots)
library(PCAtools)

# ---- load: data ----
source("scripts/R/job_gg_theme.R")
source("scripts/R/youri_gg_theme.R")
source("scripts/R/ligands.R")
source("scripts/R/gsam_metadata.R")
source("scripts/R/expression_matrix.R")
source("scripts/R/snv.R")
source("scripts/R/gsam_rnaseq_egfrviii_expression.R")

rm(gencode.31)

selected_mut <- c("G598V","G598A","A289V","A289T","A289D","R108K","R108G")
vaf$selected_mut <-  ifelse(rownames(vaf)%in%grep(paste(selected_mut,collapse="|"),vaf$egfr_snv),"Yes",NA)

c <- dplyr::left_join(gsam.rna.metadata,gsam.patient.metadata,by=c("pid"="studyID")) # c is a terrible name as there is a function named `c()`
gsam.cnv.metadata <- dplyr::left_join(gsam.cnv.metadata,vaf[,c(2,3,4)],by.x="donor_ID",by.y="donor_ID",all.x=T)
gsam.cnv.metadata$egfr_mut <- ifelse(is.na(gsam.cnv.metadata$egfr_snv),0,1)

# ---- vIII ----

gsam.viii.rnaseq$vIII_mutated <- ifelse(gsam.viii.rnaseq$egfrviii.pct>10,"Yes","No")
gsam.viii.rnaseq$sample <- gsub("-",".",gsam.viii.rnaseq$sample)
gsam.viii.rnaseq$vIII_mutated[is.na(gsam.viii.rnaseq$vIII_mutated)] <- "No"
c <- dplyr::left_join(c,gsam.viii.rnaseq[,c("sample","vIII_mutated","egfrviii.pct")],by=c("sample"="sample"))


# ---- remove old samples from 1st batch + blacklist ----

m <- gsam.rna.metadata
expr <- expression_matrix_full_new

#replicates <- gsam.rna.metadata$sample[grepl("replicate",gsam.rna.metadata$sample)] %>% gsub(".replicate","",.)
#replicate_pairs <- gsam.rna.metadata %>% dplyr::select(sample,wt.reads.v3,vIII.reads.v3,vIII.percentage) %>% dplyr::filter(str_detect(sample, "CAO$*.|FAB$*.|GAS$*."))
#replicate_pairs$reads <- replicate_pairs$wt.reads.v3 + replicate_pairs$vIII.reads.v3 #exclude CAO1.replicate, GAS2.replicate and FAB2

blacklist <- c("KAC2.new","CAO1.replicate","GAS2.replicate","FAB2","AAC2") #to be updated after DNA contamination calculations

#blacklist.pca <- gsam.rna.metadata %>% dplyr::filter(blacklist.pca == T) %>% dplyr::pull(sample)
#blacklist <- c(blacklist, blacklist.pca) %>% unique()

expr <- expr %>% 
     dplyr::select(-all_of(blacklist)&-all_of(gsam.rna.metadata[gsam.rna.metadata$batch == "old",2])) %>%
     dplyr::filter(!rownames(.)%in%rownames(.)[grepl("chrM",rownames(.))]) %>%
     dplyr::select(which(colSums(.)>200000)) 


IDH1_mutated <- gsam.cnv.metadata %>% 
  dplyr::filter(IDH1 == "Stable" | IDH1 == "Lost" | IDH1 == "Gained") %>% 
  dplyr::select(pid) %>%
  dplyr::left_join(.,m[,c("pid","sample")],by=c("pid"="pid")) %>%
  dplyr::filter(sample%in%colnames(expr)) %>%
  dplyr::distinct(sample) %>% 
  dplyr::pull(sample)

expr <- expr %>% dplyr::select(-all_of(IDH1_mutated))


rm(blacklist)



# ---- DESeq2: random conditions----


# @ comment van youri
# Hoezo was dit plaatje een probleem?
# De MA is niet biassed met random condities toch?


lowq_samples <- m[m$blacklist.pca == T & m$batch == "single",2] 
e <- expr %>% dplyr::select(!which(colnames(.)%in%lowq_samples))

cond <- as.factor(runif(ncol(e)) > 0.5)

dds <- DESeqDataSetFromMatrix(e, data.frame(cond), ~cond)  %>% 
       estimateSizeFactors(.)
nc <- counts(dds, normalized=TRUE)
filter <- rowSums(nc >= 10) >= 5 #normalized count of at least 10 in five or more samples
res <- dds[filter,] %>% DESeq(.,parallel = T) %>% results(.)

DESeq2::plotMA(res,ylim=c(-2,2))


#e.vst <- assay(vst(dds, blind=TRUE))
#write.csv(e.vst,"data-local/vst_counts_gsam.csv")

# ---- Unsupervised analysis----
#e.vst <- read.csv("data-local/vst_counts_gsam.csv",row.names = c(1))

#PCAtools package
metadata <- c %>% 
            dplyr::select(sample,resection,gender.corrected,age,survivalDays) %>%
            dplyr::filter(sample%in%colnames(e.vst))
rownames(metadata) <- metadata$sample
all(colnames(e.vst) == rownames(metadata))

p <- pca(e.vst, metadata = metadata, removeVar = 0.1) #remove lower 10% of variables based on variance

biplot(p,colby="resection")

plotloadings(p, labSize = 3,  
             legendPosition = 'bottom',
             legendLabSize = 10,  
             gridlines.major = F,
             gridlines.minor = F)

eigencorplot(p,
             metavars = c('gender.corrected','age','resection','survivalDays'))

# PCA + plot
ntop <- 500
variances <- rowVars(as.matrix(e.vst))
select <- order(variances, decreasing = TRUE)[seq_len(min(ntop, length(variances)))]
high_variance_genes <- e.vst[select,]

pc <- as.data.frame(prcomp(t(high_variance_genes))$x) %>% 
  dplyr::mutate(sample=rownames(.)) %>% 
  dplyr::inner_join(.,c,by="sample") %>%
  dplyr::select(sample,PC1,PC2,PC3,PC4,percentage.G,percentage.C,resection,gender.corrected,survivalDays) %>%
  dplyr::mutate(gc_bias=percentage.G+percentage.C)

ggplot(data=pc,mapping=aes(x=PC1,y=PC2)) +
  geom_point(aes(colour = resection)) +
  job_gg_theme +
  labs(colour="resection")
#ggsave("output/pc_plot_resection.pdf",width=7,height=5)

ggplot(data=pc,mapping=aes(x=PC1,y=PC2)) +
  geom_point(aes(colour = gc_bias)) +
  scale_colour_gradient(low = "blue", high = "red")  +
  geom_text_repel(data=subset(pc,gc_bias>60),mapping=aes(label=sample),colour="grey30",nudge_x = -3) +
  job_gg_theme +
  labs(colour="%GC")

#pvclust
new_samples <- e.vst[,grep("new",colnames(e.vst))]
pvclust <- pvclust::parPvclust(cl=NULL,new_samples,method.dist = "euclidean", method.hclust = "ward.D2",nboot=10,iseed=NULL)
plot(pvclust,print.num = FALSE,print.pv = FALSE,cex=1)
pvclust::pvrect(pvclust, alpha=0.95)

#Heatmap
data <- e.vst
data$mad <- apply(data, 1, mad)  
data <- data %>% dplyr::filter(mad > 1) %>% dplyr::select(.,-mad)

data_scaled <- t(scale(t(data)))
hypermutated <- c("GAR","CAO","JAC","JAF","JAG","FAG","FAN","BAB","BAD","BAK","AQA","AAT","EAZ") #from paper

meta <- c %>%
  dplyr::filter(sample%in%colnames(data_scaled)) %>%
  dplyr::mutate(hp=ifelse(pid%in%hypermutated & resection == "r2","Yes","No")) %>%
  dplyr::mutate(donor_ID=substr(sample,1,4))
meta$cnStatusEGFRs <- as.factor(meta$cnStatusEGFRs)
meta <- merge(meta,vaf[,c(2,3,4)],by.x="donor_ID",by.y="donor_ID",all.x=T)
meta$selected_mut <- ifelse(is.na(meta$selected_mut),"No","Yes")
rownames(meta) <- meta$sample

annotation <- data.frame(resection = meta$resection,
                         hypermutated = meta$hp,
                         EGFR = meta$cnStatusEGFRs,
                         vIII_mutated = meta$vIII_mutated,
                         A289_G598_R108 = meta$selected_mut,
                         row.names = row.names(meta))

annotation_colors <- list(
                     resection = c(r1="cornflowerblue",r2="maroon"),
                     EGFR=c(Stable = "coral",Lost="red",Gained="darkgreen",Wildtype="gold4"),
                     hypermutated = c(Yes="black",No="grey70"),
                     A289_G598_R108 = c(Yes="palegreen",No="mediumpurple"),
                     vIII_mutated = c(Yes="pink",No="blue"))
#pdf("output/unsupervised_heatmap_r1r2.pdf",width=3*7.67 / 1.45, height = 6.83 / 1.45 * 2)
pheatmap(data,
         show_rownames = F,
         scale = 'row',
         clustering_method = "complete",
         annotation = annotation, 
         annotation_colors = annotation_colors,
         color=colorRampPalette(c("navy", "white", "red"))(50),
         show_colnames = F)
#dev.off()

# ---- DESeq2: resection as condition (no correction)----

m <- gsam.rna.metadata
lowq_samples <- m[m$blacklist.pca == T & m$batch == "single",2] %>% substr(1,3)
e <- expr %>% 
     dplyr::select(-c("EBP1","FAH2","KAE1.new","AAC1")) %>% 
     dplyr::select(-all_of(colnames(.)[substr(colnames(.),1,3)%in%lowq_samples])) #no pair

# make sure order metadata aligns with expression data (!)
stopifnot(sum(colnames(e) == m[match(colnames(e),m$sample),]$sample) == ncol(e))

cond <- factor(m[match(colnames(e),m$sample),]$resection, levels=c('r1','r2') ) # ensure r1 is first named level

colnames(e) == m[match(colnames(e),m$sample),]$sample

dds <- DESeqDataSetFromMatrix(e, DataFrame(cond), ~cond) %>%
       estimateSizeFactors(.)

nc <- counts(dds, normalized=TRUE)
filter <- rowSums(nc >= 10) >= 5 #normalized count of at least 10 in five or more samples
res <- dds[filter,] %>% DESeq(.,parallel = T) %>% results(.)

DESeq2::plotMA(res, alpha = 0.01, main='DESeq2: Resection 1 vs 2',ylim=c(-3,3))

resOrdered <- res[order(res$padj),]
rownames(resOrdered) <- gsub("\\|chr.*","",rownames(resOrdered))
rownames(resOrdered) <- gsub(".*\\|","",rownames(resOrdered))

EnhancedVolcano(title = NULL,
                subtitle = NULL,
                resOrdered,
                lab = rownames(resOrdered),
                axisLabSize = 12,
                legendPosition = "right",
                legendLabSize = 9,
                x = 'log2FoldChange',
                y = 'padj',  
                pCutoff = 0.01,
                FCcutoff = 1.5,
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                col = c("black", "forestgreen", "blue", "red"),
                xlim = c(-2,3), 
                captionLabSize = 9)

dtable <- as.data.frame(resOrdered) %>% dplyr::filter(abs(log2FoldChange) > 1.5 & padj < 0.01)
#Adjusted P-value cutoff: 0.01, fold change cutoff = 1.5

#write.csv(as.data.frame(resOrdered),"output/DESeq2_Results_r1vsr2.csv")


# ---- + 2a: volcano [verhaak subtypes ] -----

# "NEFL, GABRA1, SYT1 and SLC12A5"
# FBXO3 GABRB2 SNCG, MBP


# @todo stop in ./scripts/R/gbm.subtypes.R
genes.neural <- c(
  #"NEFL", "GABRA1", "SYT1", "SLC12A5", "FBXO3", "GABRB2", "SNCG", "MBP"
  "TTPA","SIRT5","CASQ1","AKR7A3","MRPL49","GUK1","VSX1","NDRG2","PPP2R5A","RND1","ZNF323","LYRM1","SEPW1","USP33","ANKRD46","SPAST","PRPSAP2","PDE6D","ORC4L","SCHIP1","NSL1","CRBN","CRYZL1","ACYP2","MGST3","PEX19","MDH1","ATP5L","TSNAX","MAT2B","YPEL5","TCEAL1","CALM2","ATP5F1","COX5B","PEX11B","IMPA1","TTC1","GABARAPL2","NDUFS3","FBXO3","CCDC121","CRYL1","SNX11","GABRB2","SERPINI1","KCNK1","SNCG","CPNE6","KCNJ3","GRM1","VIP","HPCAL4","HPCA","CRYM","CCK","GPR22","CHN1","CA4","ADD3","CAMK2G","NTSR2","AGXT2L1","EDG1","MYBPC1","PPP1R1A","FEZF2","LOC201229","SLCO1A2","DHRS9","FLJ22655","THTPA","SLC30A10","ANXA3","FXYD1","PARP8","MAGEH1","SLC16A7","PGBD5","TPM3","CDC42","PSCD1","ENPP4","MMD","FAM49B","ARRB1","ROGDI","GRM3","SGK3","PLCL1","BEST1","KIAA1598","ENPP2","EDIL3","SLC31A2","EVI2A","TMEM144","CLCA4","MBP","UGT8","BCAS1","SH3GL3","REPS2","BASP1","DYNC1I1","SH3GL2","FUT9","ANKS1B","CDR1","ATRNL1","GNAI1","AGTPBP1","EPB41L3","FHIT","NANOS1","TCEAL2","PPFIA2","ANXA7","UROS","PPA1","SAR1A","CUTC","MSRB2","HPRT1","ACSL4","POPDC3","MGC72080","SEPP1","MORF4L2"
) # test

# wang - https://ars.els-cdn.com/content/image/1-s2.0-S1535610817302532-mmc2.xlsx
genes.mesenchymal <- c("ARPC1B","S100A11","CTSC","GLIPR1","NNMT","VDR","RGS2","CTSB","TGFBI","PLAUR","LY96","BCL3","TNFAIP8","IER3","PRSS23","IL7R","RAB27A","RUNX1","P4HA2","CYP1B1","BACE2","ACPP","FTL","SLPI","RAC2","RARRES1","SYNGR2","THBS1","IL6","CAV1","PI3","CDCP1","ITGB1","LOX","CD72","COL1A2","ANPEP","MMP7","SPAG4","BNC2","NDRG1","CNN2","LUM","PTGS2","COL3A1","COL5A1","SDC1","COL1A1","GPRC5A","COL15A1")
genes.proneural <- c("HN1","RAB33A","HDAC2","MYT1","MTSS1","HOXD3","GPR17","PTTG1","KLRC3","HRASLS","TCP1","NPPA","PFDN2","CA10","EPHB1","UGT8","PAK7","SLC1A1","NARF","DCTN3","SMPD3","ZNF804A","RASL11B","MYB","PDGFRA","ERBB3","CLGN","SOX10","BCL11A","NMU","ZNF643","CDKN1C","JPH3","PCDHA9","IL1RAPL1","MAST1","VIPR2","SIM2","BAMBI","PKMYT1","PLCB4","SLC17A6","KLRK1","CENPJ","NHLH1","GABRB3","KLRC4","KCNK3","GRID2","DACH1")
genes.classical <- c("PTPRA","ELOVL2","MLC1","SOX9","ARNTL","DENND2A","BBS1","ABLIM1","PAX6","ZHX3","USP8","PLCG1","CDH4","RASGRP1","ACSBG1","CST3","BCKDHB","LHFP","VAV3","ACSL3","EYA2","SEPT11","SLC4A4","SLC20A2","C14orf159","CTNND1","ZFHX4","SPRY2","ZNF45","NCOA1","PLCE1","DTNA","POLRMT","SALL1","TYK2","TJP1","MEOX2","FGFR3","STXBP3","GRIK1","GATM","UPF1","NPEPL1","KIAA0494","RBCK1","PHKB","SLC3A2","PPARGC1A","PNPLA6","MYO5C")



plt <- resOrdered %>%
  data.frame(stringsAsFactors = F) %>%
  dplyr::filter(!is.na(padj) & !is.na(log2FoldChange)) %>%
  tibble::rownames_to_column('gene.id') %>%
  dplyr::mutate(neural = gene.id %in% genes.neural) %>%
  dplyr::mutate(proneural = gene.id %in% genes.proneural) %>%
  dplyr::mutate(mesenchymal = gene.id %in% genes.mesenchymal) %>%
  dplyr::mutate(classical = gene.id %in% genes.classical) %>%
  dplyr::mutate(log10padj = -log10(padj)) %>%
  dplyr::mutate(class = as.factor(case_when(
    neural == T ~ "neural",
    proneural == T ~ "proneural",
    mesenchymal == T ~ "mesenchymal",
    classical == T ~ "classical",
    neural == F & proneural == F &  mesenchymal == F & classical == F ~ "-"))
  )


ggplot(plt, aes(x= log2FoldChange, y = log10padj , label = gene.id, col=class)) +
  geom_point(size=0.4, col="gray60") + 
  geom_point(size=1.2, pch=1, data=subset(plt, class != "-")) +
  youri_gg_theme + 
  scale_color_manual(values = c('neural' = 'red',
                                "proneural"  = 'green',
                                '-'='white',
                                "mesenchymal"  = 'blue',
                                'classical' = 'purple')
                     , guide = guide_legend(title = NULL, title.position = 'top', title.hjust = 0.5, ncol = 4, keywidth = 0.75, keyheight = 0.75)) 



ggsave("output/figures/de_volcano_subgroups.png", width=8, height = 4)


# ---- + 2b: volcano [gravendeel subtypes ] -----

genes.X9  <- readRDS(file = "tmp/genes.X9.Rds")
genes.X16 <- readRDS(file = "tmp/genes.X16.Rds")
genes.X17 <- readRDS(file = "tmp/genes.X17.Rds")
genes.X18 <- readRDS(file = "tmp/genes.X18.Rds")
genes.X22 <- readRDS(file = "tmp/genes.X22.Rds")
genes.X23 <- readRDS(file = "tmp/genes.X23.Rds")




plt <- resOrdered %>%
  data.frame(stringsAsFactors = F) %>%
  dplyr::filter(!is.na(padj) & !is.na(log2FoldChange)) %>%
  tibble::rownames_to_column('gene.id') %>%
  dplyr::mutate(X9 = gene.id %in% genes.X9) %>%
  dplyr::mutate(X16 = gene.id %in% genes.X16) %>%
  dplyr::mutate(X17 = gene.id %in% genes.X17) %>%
  dplyr::mutate(X18 = gene.id %in% genes.X18) %>%
  dplyr::mutate(X22 = gene.id %in% genes.X22) %>%
  dplyr::mutate(X23 = gene.id %in% genes.X23) %>%
  
  dplyr::mutate(X9.up = gene.id %in% genes.X9.up) %>%
  dplyr::mutate(X16.up = gene.id %in% genes.X16.up) %>%
  dplyr::mutate(X17.up = gene.id %in% genes.X17.up) %>%
  dplyr::mutate(X18.up = gene.id %in% genes.X18.up) %>%
  dplyr::mutate(X22.up = gene.id %in% genes.X22.up) %>%
  dplyr::mutate(X23.up = gene.id %in% genes.X23.up) %>%
  
  dplyr::mutate(X9.down = gene.id %in% genes.X9.down) %>%
  dplyr::mutate(X16.down = gene.id %in% genes.X16.down) %>%
  dplyr::mutate(X17.down = gene.id %in% genes.X17.down) %>%
  dplyr::mutate(X18.down = gene.id %in% genes.X18.down) %>%
  dplyr::mutate(X22.down = gene.id %in% genes.X22.down) %>%
  dplyr::mutate(X23.down = gene.id %in% genes.X23.down) %>%
  
  dplyr::mutate(log10padj = -log10(padj)) %>%
  dplyr::mutate(class = as.factor(case_when(
    X9 == T ~ "G-9",
    X16 == T ~ "G-16",
    X17 == T ~ "G-17",
    X18 == T ~ "G-18",
    X22 == T ~ "G-22" ,
    X23 == T ~ "G-23" ,
    X9 == F & X16 == F & X17 == F & X18 == F & X22 == F & X23 == F ~ "-"))
  ) %>%
  dplyr::mutate(class = as.factor(case_when(
    X22.up == T ~ "G-9 (up)",
    X22.down == T ~ "G-9 (down)",
    X22 == F ~ "-"))
  )

# X9 lijkt neural


ggplot(plt, aes(x= log2FoldChange, y = log10padj , label = gene.id, col=class)) +
  geom_point(size=0.4, col="gray60") + 
  geom_point(size=1.2, pch=1, data=subset(plt, class != "-")) +
  youri_gg_theme + 
  scale_color_manual(values = c( 'red', 'darkgreen', 'blue', 'brown' , 'purple' , 'black')
                   , guide = guide_legend(title = NULL, title.position = 'top', title.hjust = 0.5, ncol = 4, keywidth = 0.75, keyheight = 0.75)) 





# ---- + 2c: supervised clustering ----

# - welke genen? padj < 0.01 & abS(lfc) > 1.5 ?
# - correctie voor geslacht van vst waardes voor clustering

gs.de <- resOrdered %>%
  data.frame(stringsAsFactors = F) %>%
  dplyr::filter(padj < 0.01 & abs(log2FoldChange) > 1.5) %>%
  tibble::rownames_to_column('gid') %>%
  dplyr::pull('gid')


if(!exists('e.vst')) {
  e.vst <- read.csv("/home/sghisai/project/G-SAM/data-local/vst_counts_gsam.csv",stringsAsFactors = F)
}

stopifnot(colnames(e) %in% colnames(e.vst))
e.vst.cluster <- e.vst %>%
  dplyr::rename(gid = X) %>%
  dplyr::mutate(hgnc.symbol = gsub('^[^|]+\\|','', gid) ) %>%
  dplyr::mutate(hgnc.symbol = gsub('\\|.+$','', hgnc.symbol ) ) %>%
  dplyr::select(c('gid', 'hgnc.symbol',colnames(e))) %>%
  dplyr::filter(hgnc.symbol %in% gs.de) %>%
  dplyr::mutate(sumVar = rowSums(dplyr::select(., contains("1") | contains("2")))) %>% # some genes appear twice such as SNORD73, only pick the one w/ highest count
  dplyr::arrange(- sumVar ) %>%
  dplyr::mutate(sumVar = NULL) %>%
  dplyr::filter(!duplicated(hgnc.symbol)) %>%
  dplyr::mutate(gid = NULL) %>%
  tibble::column_to_rownames('hgnc.symbol')


# apply gender correction
meta <- data.frame(sid = colnames(e.vst.cluster)) %>%
  dplyr::mutate(pid = gsub("^(...).+$","\\1", sid) ) %>%
  dplyr::left_join(gsam.patient.metadata, c('pid' = 'studyID'))


e.vst.cluster.n <- limma::removeBatchEffect(e.vst.cluster, meta$gender.corrected)


annotation_col <- meta %>%
  dplyr::select(c('sid', 'gender.corrected')) %>%
  dplyr::mutate(res = as.factor( grepl ("2", sid) ) ) %>%
  dplyr::mutate(gender.corrected = NULL) %>%
  tibble::column_to_rownames('sid')


annotation_row <- data.frame(gid = rownames(e.vst.cluster)) %>%
  dplyr::mutate(neural = as.factor(gid %in% genes.neural)) %>%
  #dplyr::mutate(mesenchymal = as.factor(gid %in% genes.mesenchymal)) %>%
  #dplyr::mutate(proneural = as.factor(gid %in% genes.proneural)) %>%
  #dplyr::mutate(classical = as.factor(gid %in% genes.classical)) %>%
  tibble::column_to_rownames('gid')


#c <- cor(t(scale(e.vst.cluster.n)), method="pearson")
c <- cor(t(e.vst.cluster.n), method="pearson")
d <- as.dist(1 - c) 
clust_row <- hclust(d)
plot(clust_row)




c <- cor(t(e.vst.cluster.n), method="kendall")
d <- as.dist(1 - c) 
clust_row <- hclust(d)
plot(clust_row)


#c <- cor(e.vst.cluster.n, method="spearman")


hm <- pheatmap::pheatmap(e.vst.cluster.n,
                   scale="row",
                   #cluste_rows = clust_row,
                   annotation_col = annotation_col,
                   annotation_row = annotation_row,
                   clustering_distance_rows = "correlation",
                   clustering_distance_cols = "manhattan",
                   fontsize_row = 2,
                   fontsize_col = 2,
                   filename = "output/figures/supervised_heatmap.png",
                   cutree_rows = 3 ,
                   height=8,
                   width=22)




# ---- DESeq2: resection as condition (correct for gender)----

m <- gsam.rna.metadata
lowq_samples <- m[m$blacklist.pca == T & m$batch == "single",2] %>% substr(1,3)
e <- expr %>% 
  dplyr::select(-c("EBP1","FAH2","KAE1.new")) %>% 
  dplyr::select(-all_of(colnames(.)[substr(colnames(.),1,3)%in%lowq_samples])) #no pair

# make sure order metadata aligns with expression data (!)
stopifnot(sum(colnames(e) == m[match(colnames(e),m$sample),]$sample) == ncol(e))

resection <- as.factor(m[match(colnames(e),m$sample),]$resection)
gender <- as.factor(c[match(colnames(e),c$sample),]$gender.corrected) 
colnames(e) == c[match(colnames(e),c$sample),]$sample

dds <- DESeqDataSetFromMatrix(e, data.frame(resection,gender), ~gender+resection) %>% 
       estimateSizeFactors(.) 
nc <- counts(dds, normalized=TRUE)
filter <- rowSums(nc >= 10) >= 5 #normalized count of at least 10 in five or more samples
res <- dds[filter,] %>% DESeq(.,parallel = T) %>% results(.)

#Batch correction visualisation
#e.vst <- assay(vst(dds, blind=T)) --> for plots (PC3vsPC4) run PCA lines in "Unsupervised analysis section"
#vst <- vst(dds, blind=T)
#mat <- assay(vst) %>% limma::removeBatchEffect(.,vst$gender) #--> for plots (PC3vsPC4) run PCA lines in "Unsupervised analysis section"

DESeq2::plotMA(res, main='DESeq2: Resection 1 vs 2')

resOrdered <- res[order(res$padj),]
rownames(resOrdered) <- gsub("\\|chr.*","",rownames(resOrdered))
rownames(resOrdered) <- gsub(".*\\|","",rownames(resOrdered))

pvalues <- as.data.frame(resOrdered)


#Adjusted pvalue: 0.01, log fold change cutoff = 1.5
EnhancedVolcano(title = NULL,
                subtitle = NULL,
                resOrdered,
                lab = rownames(resOrdered),
                selectLab = genes.neural,
                axisLabSize = 12,
                legendPosition = "bottom",
                legendLabSize = 9,
                x = 'log2FoldChange',
                y = 'padj',  
                pCutoff = 0.01,
                FCcutoff = 1.5,
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                legend=c('Not significant','Log2 FC','p-value',
                               'Adj. P-value & Log2 FC'),
                col = c("black", "forestgreen", "blue", "red"),
                xlim = c(-2,3), 
                captionLabSize = 9)