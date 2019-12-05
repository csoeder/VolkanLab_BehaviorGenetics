#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
# 1: counts file in
# 2: prefix for deseq out
# 3: name of contrast


library("dplyr")
#library("yaml")
library("readr")
library("tidyr")
library("magrittr")

#library("stringr")

counts_in <- args[1] # counts from featureCount
#counts_in <- "counts/all.vs_dm6main.dm6_genes.mapspliceMulti.MpBC.counts"
aln_meta <- args[2] # metadata for alignments (to get million mapped reads)
#aln_meta <- "summaries/alignments.vs_dm6main.mapspliceMulti.summary"
exp_out <- args[3] # output file to write expression
#exp_out <- "expression.potato"

counts.df <- read_delim(counts_in, "\t", escape_double = FALSE, trim_ws = TRUE)

RPK.df <- counts.df %>% select(-c(Chr, Start, End, Strand)) %>% mutate(Length = Length/1000) # convert bases to kb
#turn counts into per-kb rates
#https://stackoverflow.com/a/48693978
RPK.df[,-(1:2)] %<>% sapply(`/`, RPK.df[,2])
RPK.df.gath <- RPK.df  %>% gather(key = "sample", value = "count", -one_of("Geneid", "Length"))


#collect the million reads mapped (proper pairs, given the current featureCounts settings)
aln_meta.df <- read_delim(aln_meta, "\t", escape_double = FALSE, trim_ws = TRUE,col_names = F)
names(aln_meta.df)<- c("ref_genome","aligner","sample","measure","value")

aln_meta.df %<>% filter(measure == "properly_paired_count") %<>% select(-c("measure", "aligner", "ref_genome"))


RPKM.df.gath <- inner_join(RPK.df.gath, aln_meta.df, by =c("sample"="sample"))

RPKM.df.gath %<>% mutate(MMR = value/1000000, expression = count/MMR )  %>% select(c("Geneid", "sample", "expression"))

RPKM.df.sprud <- RPKM.df.gath %>% spread(key = "sample", value = "expression")

write_delim(RPKM.df.sprud, exp_out, delim = "\t", col_names = T)


