#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
# 1: counts file in
# 2: prefix for deseq out
# 3: name of contrast

source("utils/DEXSeq/Subread_to_DEXSeq/load_SubreadOutput.R")


library("dplyr")
library("yaml")
library("readr")
library("tidyr")
library("magrittr")
library("stringr")

counts_in <- args[1] #dxCounts/wtHousing.vs_dm6main.dm6_genes.mapspliceMulti.M.dxsq.counts
dsq_pre <- args[2] #diff_exon_use/dex_wtHousing/wtHousing.vs_dm6main.dm6_genes.mapspliceMulti.M
dsqx_gtf <- args[3] #utils/annotations/DEXSeq/dm6_genes.dxsqReady.gtf
nombre <- args[4] #dex_wtHousing
mPhatic <- args[5] #"chr3R_FBgn0004652-"
#mPhatic <- "chr3R_FBgn0004652-"
#use example:
#	Rscript scripts/deSeqer.R {input.fc_in} diff_exon_use/{wildcards.contrast}/{wildcards.group}.vs_{wildcards.ref_genome}.{wildcards.annot}.{wildcards.aligner}.{wildcards.flag} {wildcards.contrast} {gene_id}
#	Rscript scripts/deSeqer.R dxCounts/wtHousing.vs_dm6main.dm6_genes.mapspliceMulti.M.dxsq.counts diff_exon_use/dex_wtHousing/wtHousing.vs_dm6main.dm6_genes.mapspliceMulti.M utils/annotations/DEXSeq/dm6_genes.dxsqReady.gtf dex_wtHousing "_FBgn0004652-"


#counts_in <- "counts/DEXSeq/wtHousing.vs_dm6main.dm6_genes.mapspliceMulti.M.dxsq.counts"
#dsqx_gtf <- "utils/annotations/DEXSeq/dm6_genes.dxsqReady.gtf"

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
print("metadata frame built")

#	build the contrast experiment
contrasts.df <- plyr::ldply(trammel$xcontrasts, data.frame) %>% filter(name == nombre )
print("contrast defined")

# prep the sample metadata
sbgrp <- (contrasts.df %>% select(filt))[1,] %>% as.character
fact <- (contrasts.df %>% select(vars))[1,] %>% as.character
coldata.df <- data_sets.df%>% group_by(name) %>%  mutate(genotype = paste0(unique(as.character(genotype)), collapse = "." )) %>% ungroup()  %>% filter(subgroups == sbgrp) %>% as.data.frame() 
row.names(coldata.df) <- coldata.df$name
coldata.df %<>% dplyr::select(fact)
names(coldata.df) <- "condition"
print("experimental metadata selecteted")



# Run DEXSeq
#desgn <- "~ sample + exon + condition:exon"
desgn <- as.character(contrasts.df$design[1])
dxd.fc = DEXSeqDataSetFromFeatureCounts(counts_in,flattenedfile = dsqx_gtf, sampleData =coldata.df, design= eval(parse(text=desgn)))
print("DEXSeq Dataset Built")

dxd.fc = estimateSizeFactors( dxd.fc )
print("size factors estimated")

dxd.fc = estimateDispersions( dxd.fc )
print("dispersion estimated")


png(filename = paste(dsq_pre,".dispersionPlot.png", sep=""), height=500, width=800)
plotDispEsts( dxd.fc , main = paste("DEXSeq Dispersion Plot: ",nombre," Contrast", sep=""))
dev.off()

dxd.fc = testForDEU( dxd.fc )
print("DEU tested for")


dxd.fc = estimateExonFoldChanges( dxd.fc, fitExpToVar="condition")
print("Exon fold changes estimated")


dxd.fc.results = DEXSeqResults( dxd.fc )
print("Results gathered")


png(filename = paste(dsq_pre,".",gsub("-","",gsub("\\+","",gsub("_","",mPhatic))),".differentialExonUsePlot.png", sep=""), height=500, width=800)
plotDEXSeq( dxd.fc.results, mPhatic, expression=FALSE, splicing=TRUE, legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2, FDR=0.01 )
dev.off()

png(filename = paste(dsq_pre,".",gsub("-","",gsub("\\+","",gsub("_","",mPhatic))),".differentialExonExpressionPlot.png", sep=""), height=500, width=800)
plotDEXSeq( dxd.fc.results, mPhatic, expression=TRUE, splicing=FALSE, legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 , FDR=0.01)
dev.off()


dxd.fc.results.df <- dxd.fc.results %>% as_tibble()


#Write to TSV
#dxd.fc.results.df <- dxd.fc.results.df %>% select(geneid, everything())
write_delim(dxd.fc.results.df %>% mutate(transcripts=as.character(transcripts)), paste(dsq_pre, "de", sep="."), delim = "\t")
print("Results written")

print("................... DONE!!!")


