# ---- initialization setup ----

#setwd("~/project/G-SAM")

# ---- load packages ----
library(igraph)
library(tidyverse)
library(survival)
library(survminer)
library(igraph)
library(ggsci)

# ---- load: data ----
source("scripts/R/job_gg_theme.R")
source("scripts/R/ligands.R")
source("scripts/R/gsam_metadata.R")
source("scripts/R/expression_matrix.R")
source("scripts/R/snv.R")
source("scripts/R/gsam_rnaseq_egfrviii_expression.R")

rm(gencode.31)
c <- inner_join(gsam.rna.metadata,gsam.patient.metadata,by=c("pid"="studyID"))
vst_counts <- read.csv("data-local/vst_counts_gsam.csv",row.names = c(1))
m <- gsam.rna.metadata
gene_signature <- read.csv("data-local/50genesignature.csv",sep=";",header=F)

#vIII
gsam.viii.rnaseq$vIII_mutated <- ifelse(gsam.viii.rnaseq$egfrviii.pct>10,"Yes","No")
gsam.viii.rnaseq$sample <- gsub("-",".",gsam.viii.rnaseq$sample)
gsam.viii.rnaseq$vIII_mutated[is.na(gsam.viii.rnaseq$vIII_mutated)] <- "No"
c <- dplyr::left_join(c,gsam.viii.rnaseq[,c("sample","vIII_mutated")],by=c("sample"="sample"))

# ---- classification GBM subtype (50 gene signature from Wang Q. et al. 2017) ----

# vst_counts$gene <- rownames(vst_counts)
# vst_counts$gene <- gsub("\\|chr.*","",rownames(vst_counts))
# vst_counts$gene <- gsub(".*\\|","",vst_counts$gene)
# vst_counts <- vst_counts %>% dplyr::filter(!duplicated(gene))
# rownames(vst_counts) <- vst_counts$gene
# vst_counts$gene <- NULL
#
# 
# #Adjust names based on aliases in gene signature list
# rownames(vst_counts)[which(rownames(vst_counts)=="JPT1")] <- "HN1"
# rownames(vst_counts)[which(rownames(vst_counts)=="PLAAT1")] <- "HRASLS"
# rownames(vst_counts)[which(rownames(vst_counts)=="PAK5")] <- "PAK7"
# rownames(vst_counts)[which(rownames(vst_counts)=="ZFP69B")] <- "ZNF643"
# rownames(vst_counts)[which(rownames(vst_counts)=="LHFPL6")] <- "LHFP"
# rownames(vst_counts)[which(rownames(vst_counts)=="DGLUCY")] <- "C14orf159"
# rownames(vst_counts)[which(rownames(vst_counts)=="EFCAB14")] <- "KIAA0494"
# 
# length(intersect(gene_signature$V1,rownames(vst_counts)))
# vst_counts <- vst_counts[rownames(vst_counts)%in%gene_signature$V1,]
# 
# vst_counts <- t(vst_counts)
# write.csv(vst_counts,"data-local/expression_sign_genes.csv")

# ---- subtypes from GlioViz ----
lowq_samples <- m[m$blacklist.pca == T & m$batch == "single",2] %>% substr(1,3)
idh1_mutated <- gsam.cnv.metadata %>% 
                dplyr::filter(IDH1 == "Stable" | IDH1 == "Lost" | IDH1 == "Gained") %>% 
                dplyr::pull(pid) 

subtypes <- read.csv("data-local/GlioVis_subtypes.csv") 
subtypes <- subtypes %>% 
            dplyr::filter(!Sample%in%c("EBP1","FAH2","KAE1.new","KAC1.new","AAC1","AAC2")) %>%
            dplyr::inner_join(.,c[,c("sample","resection","pid","survivalDays","survivalFromSecondSurgeryDays","progressionFreeDays","survival.event")],by=c("Sample"="sample")) %>%
            dplyr::filter(!pid%in%lowq_samples) %>%
            dplyr::filter(!pid%in%idh1_mutated) %>%
            dplyr::group_by(pid) %>% 
            dplyr::mutate(transition=length(unique(majority.call))) %>%
            dplyr::ungroup(pid) %>%
            dplyr::mutate(majority.call=as.character(majority.call)) %>%
            dplyr::mutate(cl_unchanged=ifelse(transition==1 & majority.call == "Classical" & resection == "r2",T,"")) %>%
            dplyr::mutate(mes_unchanged=ifelse(transition==1 & majority.call == "Mesenchymal" & resection == "r2",T,"")) %>%
            dplyr::mutate(pn_unchanged=ifelse(transition==1 & majority.call == "Proneural" & resection == "r2",T,"")) %>%
            dplyr::mutate(to_subtype=ifelse(transition==2 & resection == "r2",majority.call,"")) 

no_majority_call <- subtypes %>% 
                    dplyr::mutate(no.majority.call=ifelse(svm_call == knn_call | svm_call == gsea_call | knn_call == gsea_call,F,T)) %>%
                    dplyr::filter(no.majority.call == T) %>%
                    dplyr::pull(pid)
subtypes <- subtypes %>%
            dplyr::filter(!pid %in% no_majority_call)

to_mes <- subtypes %>% dplyr::filter(to_subtype=="Mesenchymal") %>% dplyr::pull(pid)
to_cl <- subtypes %>% dplyr::filter(to_subtype=="Classical") %>% dplyr::pull(pid)
to_pn <- subtypes %>% dplyr::filter(to_subtype=="Proneural") %>% dplyr::pull(pid)

CL <- c(sum(subtypes$majority.call=="Classical"&subtypes$resection=="r1"&subtypes$transition==1), #cl-to-cl,
               sum(subtypes$majority.call=="Classical"&subtypes$resection=="r1"&subtypes$pid%in%to_mes), #cl-to-mes,
               sum(subtypes$majority.call=="Classical"&subtypes$resection=="r1"&subtypes$pid%in%to_pn)) #cl-to-pn

MES <- c(sum(subtypes$majority.call=="Mesenchymal"&subtypes$resection=="r1"&subtypes$pid%in%to_cl), #mes-to-cl,
         sum(subtypes$majority.call=="Mesenchymal"&subtypes$resection=="r1"&subtypes$transition==1), #mes-to-mes,
         sum(subtypes$majority.call=="Mesenchymal"&subtypes$resection=="r1"&subtypes$pid%in%to_pn)) #mes-to-pn,

PN <- c(sum(subtypes$majority.call=="Proneural"&subtypes$resection=="r1"&subtypes$pid%in%to_cl), #pn-to-cl
        sum(subtypes$majority.call=="Proneural"&subtypes$resection=="r1"&subtypes$pid%in%to_mes), #pn-to-mes
        sum(subtypes$majority.call=="Proneural"&subtypes$resection=="r1"&subtypes$transition==1)) #pn-to-pn

# sum(CL[1],PN[3],MES[2])/sum(CL,PN,MES)
# CL[1]/sum(CL) * 100
# PN[3]/sum(PN) * 100
# MES[2]/sum(MES) * 100


#Real values
links <- data.frame(from=c(rep(c("CL"),3),rep(c("MES"),3),rep(c("PN"),3)), to=c(rep(c("CL","MES","PN"),3)),
                    weight=c(CL,MES,PN))

links <- links %>% dplyr::group_by(from) %>% dplyr::mutate(thickness=round(weight/sum(weight),2))

nodes <- data.frame(id=c("CL","MES","PN"),
                    subtype=c("Classical","Mesenchymal","Proneural"),
                    size=c(sum(subtypes$majority.call=="Classical"&subtypes$resection=="r1"),
                           sum(subtypes$majority.call=="Mesenchymal"&subtypes$resection=="r1"),
                           sum(subtypes$majority.call=="Proneural"&subtypes$resection=="r1")))
net <- graph_from_data_frame(d=links, vertices=nodes, directed=T) 
net <- simplify(net, remove.multiple = F, remove.loops = F) 
E(net)$width <- E(net)$thickness*6
E(net)$size <- E(net)$size*0.75

plot(net, 
     edge.arrow.size=.4,
     diag=T,
     edge.label = E(net)$weight,
     arrow.mode=3,
     edge.curved=.2,
     edge.color="gray70",
     vertex.color="coral",
     vertex.label.font=2,
     edge.label.font=1,
     edge.label.cex=1.3,
     edge.label.color="blue",
     )

#Probabilities
links <- data.frame(from=c(rep(c("CL"),3),rep(c("MES"),3),rep(c("PN"),3)), to=c(rep(c("CL","MES","PN"),3)),
                    weight=c(round(CL/sum(subtypes$majority.call=="Classical"&subtypes$resection=="r1"),2),
                             round(MES/sum(subtypes$majority.call=="Mesenchymal"&subtypes$resection=="r1"),2),
                             round(PN/sum(subtypes$majority.call=="Proneural"&subtypes$resection=="r1"),2)))

#write.csv(links,"output/transition_matrix_subtypes.csv")

# ---- Kaplan Meier ----
pn_change <- subtypes %>% 
             dplyr::filter(resection=="r1" & majority.call == "Proneural" & transition == 2) %>%
             dplyr::pull(pid)

surv_data <- subtypes %>%
             dplyr::mutate(pn_switch = ifelse(pid%in%pn_change,1,0)) %>%
             dplyr::mutate(donor_ID = substr(Sample,1,4)) %>% 
             dplyr::left_join(.,vaf[,-c(1,7)],by=c("donor_ID"="donor_ID")) %>%
             dplyr::left_join(.,c[,c("sample","EGFR","vIII_mutated")],by=c("Sample"="sample"))

surv_data$selected_mut <- ifelse(surv_data$A289 == T | surv_data$G598 == T | surv_data$R108 == T,1,NA)
surv_data$selected_mut <- ifelse(is.na(surv_data$selected_mut),0,1)

surv_data$egfr_alteration <- ifelse(surv_data$selected_mut == 1 | surv_data$vIII_mutated == "Yes" ,1,0)
surv_data$transition <- ifelse(surv_data$transition == 1,0,1)


discard <- apply(surv_data[,c("survivalFromSecondSurgeryDays","progressionFreeDays","survival.event")], 1, function(x) any( is.na(x) ))
surv_data <- surv_data[!discard,]

surv_data_r1 <- surv_data %>%
             #dplyr::filter(EGFR == "Gained" | EGFR == "Stable") %>%
             #dplyr::filter(majority.call == "Classical") %>%
             dplyr::filter(resection == "r1")

surv <- Surv(time=surv_data_r1$progressionFreeDays, event = surv_data_r1$survival.event)
fit <- survfit(surv ~ majority.call, data = surv_data_r1)
ggsurvplot(fit, 
           data = surv_data_r1,
           pval = TRUE, # show p-value of log-rank test.
           conf.int.alpha = 0.125,
           palette = "uchicago",
           xlab = "Progression-free Survival (days)")
#ggsave("output/ggsurv_subtype_r1.pdf")
#ggsave("output/ggsurv_egfr_alteration_r1_within_ampli.pdf")

surv_data_r2 <- surv_data %>%
                #dplyr::filter(EGFR == "Gained" | EGFR == "Stable") %>%
                #dplyr::filter(majority.call == "Classical") %>%
                dplyr::filter(resection == "r2")

surv2 <- Surv(time=surv_data_r2$survivalFromSecondSurgeryDays, event = surv_data_r2$survival.event)
fit2 <- survfit(surv2 ~ transition, data = surv_data_r2)
ggsurvplot(fit2, 
           data = surv_data_r2,
           pval = TRUE, # show p-value of log-rank test.
           conf.int.alpha = 0.125,
           palette = "uchicago",
           xlab = "Post-progression survival (days)")
#ggsave("output/ggsurv_subtype_r2.pdf") #~majority.call
#ggsave("output/ggsurv_egfr_alteration_r2_within_ampli.pdf") ~egfr_alteration
#ggsave("output/ggsurv_subtype_switch.pdf") #~transition

# ---- Heatmap ----
annotation_genes <- gene_signature
rownames(annotation_genes) <- annotation_genes$V1
annotation_genes$V1 <- NULL
colnames(annotation_genes) <- 'GBM subtype'

annotation_samples <- subtypes[,c(1,6)]
annotation_samples$majority.call <- factor(annotation_samples$majority.call,levels=c("Mesenchymal","Proneural","Classical"))
annotation_samples <- annotation_samples %>% dplyr::arrange(.,majority.call) %>% as.data.frame(.)
rownames(annotation_samples) <- annotation_samples$Sample
annotation_samples$Sample <- NULL

vst_counts_heatmap <- read.csv("data-local/expression_sign_genes.csv",row.names = c(1))
vst_counts_heatmap <- vst_counts_heatmap[which(rownames(vst_counts_heatmap)%in%subtypes$Sample),]
vst_counts_ordered <- vst_counts_heatmap[match(rownames(annotation_samples),rownames(vst_counts_heatmap)),
                                         match(rownames(annotation_genes),colnames(vst_counts_heatmap))]

#pdf("output/heatmap_subtypes_genesignature.pdf",width=3*7.67 / 1.45, height = 6.83 / 1.45 * 2)
pheatmap(t(vst_counts_ordered),
         color=colorRampPalette(c("navy", "white", "red"))(50),
         cluster_cols = F,
         annotation_col = annotation_samples,
         cluster_rows = F,
         annotation_row = annotation_genes,
         scale = 'row',
         show_colnames = F,
         annotation_names_row = F,
         annotation_names_col = F,
         fontsize_row = 4,
         gaps_col = c(84,137),
         gaps_row = c(50,100))

# ---- PCA plot ----
vst_pca <- as.data.frame(vst_counts_heatmap)

pc <- as.data.frame(prcomp(vst_pca)$x) %>%
      dplyr::mutate(Sample=rownames(.)) %>%   
      dplyr::left_join(.,subtypes[,c("Sample","majority.call")],by="Sample") 

ggplot(data=pc,mapping=aes(x=PC1,y=PC2)) +
        geom_point(aes(colour = majority.call)) +
        job_gg_theme +
        labs(colour="GBM subtype") 
# ggsave("output/pca_plot_subtypes.pdf")
