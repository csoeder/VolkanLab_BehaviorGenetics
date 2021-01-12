#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
# 1: counts file in
# 2: prefix for deseq out
# 3: name of contrast


library("DGCA")
library("dplyr")
library("yaml")
library("readr")
library("tidyr")
library("stringr")
library("magrittr")

xpr_in <- args[1]
nombre <- args[2]
cor_pre <- args[3]




#	build the sample DF
trammel <- read_yaml("config.yaml")
data_sets.df <- plyr::ldply(trammel$data_sets, data.frame)
data_sets.df$name <- as.factor(data_sets.df$name)
data_sets.df$day<- as.factor(data_sets.df$day)
data_sets.df$subgroups<- as.factor(data_sets.df$subgroups)
data_sets.df$rep<- as.factor(data_sets.df$rep)
data_sets.df$housing<- as.factor(data_sets.df$housing)
data_sets.df$genotype<- as.factor(data_sets.df$genotype)
data_sets.df$tissue<- as.factor(data_sets.df$tissue)

#	build the correlation experiment
corr.df <- plyr::ldply(trammel$correlations, data.frame) %>% filter(name == nombre )

# prep the sample metadata
sbgrp <- (corr.df %>% select(filt))[1,] %>% as.character
#coldata <- data_sets.df %>% as.data.frame() %>% ungroup()%>% filter(subgroups == sbgrp)
#coldata <- data_sets.df %>% ungroup() %>% group_by(name) %>%  mutate(genotype = paste0(as.character(genotype), collapse = "," ))  %>% filter(subgroups == sbgrp)
coldata <- data_sets.df%>% group_by(name) %>%  mutate(genotype = paste0(unique(as.character(genotype)), collapse = "." )) %>% ungroup()  %>% filter(subgroups == sbgrp)  %>% as.data.frame() 

rownames(coldata) <- coldata$name
dsgn.df <- coldata %>% select(name, corr.df$vars %>% as.character() )%>% mutate(present =1) %>% spread(key=genotype, value=present, fill=0) 
rownames(dsgn.df) <- dsgn.df$name
dsgn.df %<>% select(-c("name"))
dsgn.mat <- dsgn.df %>% as.matrix()




#xpr_in <- "expression/rotund.vs_dm6main.dm6_genes.mapspliceMulti.MpBC.rpkm"

#xpr_in <- "expression/rotund.vs_dm6main.dm6_genes.mapspliceMulti.MpBC.rpkm"


#load and format the count data
xpr.df <- read_delim(xpr_in, "\t", escape_double = FALSE, trim_ws = TRUE)
xpr.df %<>% as.data.frame()
rownames(xpr.df) <- xpr.df$Geneid
xpr.df %<>% select(rownames(dsgn.df)) 

ddcor_res <- ddcorAll(inputMat = xpr.df, design = dsgn.mat, compare = c("wt", "rn"),adjust = "none", heatmapPlot = FALSE, nPerm = 0)

#Write to TSV
write_delim(ddcor_res, paste(cor_pre, "corr","tsv", sep="."), delim = "\t")


