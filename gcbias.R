#!/usr/bin/env R 

setwd("~/project/G-SAM")

# ---- load libs ----

library(tidyverse)

# ---- load data ----

source('scripts/R/gsam_metadata.R')
source('scripts/R/job_gg_theme.R')
source('scripts/R/youri_gg_theme.R')

gsam.gc.content <- gsam.rna.metadata %>% dplyr::select(batch,sid,percentage.A,percentage.G,percentage.C,percentage.T,RMSE)

# mediane waarden in niet outlier samples:
# a     c     t     g
# 25.11 24.61 25.61 24.67

# , alpha=batch
gsam.gc.content$gc_bias <- gsam.gc.content$percentage.G + gsam.gc.content$percentage.C   

batch = factor(batch, levels=c("single","old","new"))

ggplot(data=gsam.gc.content,aes(x=reorder(sid,gc_bias),y=gc_bias,group=batch,colour=batch)) +
  geom_line(lwd=1) +
  scale_color_manual(values = c("red","forestgreen","skyblue1")) +
  geom_text_repel(data=subset(gsam.gc.content, batch == "old" & gc_bias > 55),
                  aes(label=sid),
                  col="black",
                  size=2.5,
                  nudge_x=2) +
  geom_text_repel(data=subset(gsam.gc.content, batch == "single" & gc_bias > 55),
                  aes(label=sid),
                  col="black",
                  size=2.5,
                  nudge_x=2) +
  job_gg_theme +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  facet_grid(~batch , drop=F, scales="free_x", space = "free_x") +
  labs(x="GSAM RNA-Seq sample", y="%GC",col="batch") 

