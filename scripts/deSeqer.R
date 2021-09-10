#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
# 1: counts file in
# 2: prefix for deseq out
# 3: name of contrast


library("DESeq2")
library("dplyr")
library("yaml")
library("readr")
library("tidyr")
library("stringr")

counts_in <- args[1]
dsq_pre <- args[2]
nombre <- args[3]
config <- args[4]

#	build the sample DF
#trammel <- read_yaml("config.yaml")
trammel <- read_yaml(config)
data_sets.df <- plyr::ldply(trammel$data_sets, data.frame)
data_sets.df$name <- as.factor(data_sets.df$name)
data_sets.df$day<- as.factor(data_sets.df$day)
data_sets.df$subgroups<- as.factor(data_sets.df$subgroups)
data_sets.df$rep<- as.factor(data_sets.df$rep)
data_sets.df$housing<- as.factor(data_sets.df$housing)
data_sets.df$genotype<- as.factor(data_sets.df$genotype)
data_sets.df$tissue<- as.factor(data_sets.df$tissue)

#	build the contrast experiment
contrasts.df <- plyr::ldply(trammel$contrasts, data.frame) %>% filter(name == nombre )

# prep the sample metadata
sbgrp <- (contrasts.df %>% select(filt))[1,] %>% as.character
#coldata <- data_sets.df %>% as.data.frame() %>% ungroup()%>% filter(subgroups == sbgrp)
#coldata <- data_sets.df %>% ungroup() %>% group_by(name) %>%  mutate(genotype = paste0(as.character(genotype), collapse = "," ))  %>% filter(subgroups == sbgrp)
coldata <- data_sets.df%>% group_by(name) %>%  mutate(genotype = paste0(unique(as.character(genotype)), collapse = "." )) %>% ungroup()  %>% filter(subgroups == sbgrp) %>% unique() %>% as.data.frame() 

rownames(coldata) <- coldata$name
coldata <- coldata %>% select(contrasts.df$vars %>% as.character() )

#load and format the count data
counts.df <- read_delim(counts_in, "\t", escape_double = FALSE, trim_ws = TRUE)
counts.df.gath <- counts.df %>% select(-c(Chr, Start, End, Strand, Length)) %>% gather(key = "sample", value = "count", -one_of("Geneid"))
counts.df.sprud <- counts.df.gath %>% spread(key = sample, value = count) %>% as.data.frame() 
rownames(counts.df.sprud) <- counts.df.sprud$Geneid
counts.df.sprud <- counts.df.sprud %>%  select(rownames(coldata))

# convert counts, metadata, design into DESeq data object
#counts.dds <- DESeqDataSetFromMatrix(countData = counts.df.sprud, colData = coldata, design = ~ housing) 
desgn <- as.character(contrasts.df$design[1])
counts.dds <- DESeqDataSetFromMatrix(countData = counts.df.sprud, colData = coldata, design = eval(parse(text=desgn)))

#prefilter to remove low-count rows
keep <- rowSums(counts(counts.dds)) >= 10
counts.dds <- counts.dds[keep,]

#relevel 
#counts.dds$condition <- relevel(counts.dds$condition, ref = "untreated")
for (i in seq(nrow(contrasts.df))){
	var <- contrasts.df[i,]$vars %>% as.character()
	ruf <-contrasts.df[i,]$ref %>% as.character()
	counts.dds[[var]] <- relevel(counts.dds[[var]], ref = ruf)
}

#Run DESeq
counts.dds.dsq <- DESeq(counts.dds)

# run the vanilla DE
counts.dds.res <- results(counts.dds.dsq) %>% as.data.frame()
counts.dds.res$geneid <- rownames(counts.dds.res)

# Effect-size shrinkage
coff <- resultsNames(counts.dds.dsq)[2]
counts.dds.res.shrunk <- lfcShrink(counts.dds.dsq, coef=coff, type="apeglm") %>% as.data.frame()
#explore shrinkage types??
counts.dds.res.shrunk$geneid <- rownames(counts.dds.res.shrunk)

# scale by gene length
counts.dds.res<-inner_join(counts.dds.res, counts.df%>% select(Geneid, Length), by=c("geneid"="Geneid")) %>% mutate(expression = baseMean/Length) %>% select(-c(Length))
counts.dds.res.shrunk<-inner_join(counts.dds.res.shrunk, counts.df%>% select(Geneid, Length), by=c("geneid"="Geneid")) %>% mutate(expression = baseMean/Length) %>% select(-c(Length))

#Write to TSV
counts.dds.res.full <- inner_join(counts.dds.res %>%  as.data.frame(), counts.dds.res.shrunk %>%  as.data.frame(), by =c("geneid"="geneid"), suffix=c(".fresh",".shrunk")) %>% select(geneid, everything())
write_delim(counts.dds.res.full, paste(dsq_pre,"de",sep="."), delim = "\t")


###		Now run the itemized DE

#get the factors in the design
desgn.splt <- (desgn %>% str_replace_all("~","")%>% str_replace_all(" ","") %>% strsplit("+", fixed=TRUE))[[1]]

#first, collect the geneids and the basemeans: just run the first factor and its first nonref level
fact <- desgn.splt[1]
ref_lev <- (contrasts.df %>% filter(vars == fact))$ref
all_levs<-levels(counts.dds.dsq[[fact]])
alt_levs<- all_levs[all_levs != ref_lev]

basic.df <- results(counts.dds.dsq, contrast=c(fact, as.character(ref_lev), alt_levs[1])) %>% as.data.frame() %>% select(c("baseMean"))
basic.df$geneid <- rownames(basic.df)
basic.df <- select(basic.df, c("geneid", "baseMean"))

#for each factor in factors....

dds.dsq.byFactor.df <- select(basic.df, c("geneid") )


for (fact in desgn.splt){

	ref_lev <- as.character((contrasts.df %>% filter(vars == fact))$ref)
	all_levs<-levels(counts.dds.dsq[[fact]])
	alt_levs<- all_levs[all_levs != ref_lev]


#running basemean by level:
#https://support.bioconductor.org/p/63567/

	#	# gather baseMeans for each level in factor.levels:
	factor.baseMeanPerLvl <- sapply( levels(counts.dds.dsq[[fact]]), function(lvl) rowMeans( counts(counts.dds.dsq,normalized=TRUE)[,counts.dds.dsq[[fact]] == lvl, drop=F] ) ) %>% as.data.frame()

	names(factor.baseMeanPerLvl ) <- paste0(paste("baseMean.",fact,".", sep=""), names(factor.baseMeanPerLvl))
	factor.baseMeanPerLvl$geneid <- rownames(factor.baseMeanPerLvl)
	factor.baseMeanPerLvl <-inner_join(factor.baseMeanPerLvl, counts.df%>% select(Geneid, Length), by=c("geneid"="Geneid"))

	for (levi in all_levs){

		factor.baseMeanPerLvl$x <- factor.baseMeanPerLvl[[paste("baseMean",fact,levi,sep=".")]]
		factor.baseMeanPerLvl[[paste("expression",fact,levi,sep=".")]] <- factor.baseMeanPerLvl$x/factor.baseMeanPerLvl$Length

	}

	factor.baseMeanPerLvl <- factor.baseMeanPerLvl %>% select(-c("x", "Length"))

	dds.dsq.byFactor.df <- full_join(dds.dsq.byFactor.df, factor.baseMeanPerLvl, by = c("geneid"="geneid"))

#	#for each alt_level in factor.levels:
	for (aLev in alt_levs){

		#dds.res_factor_level <- results(counts.dds.dsq, contrast=c(fact, ref_lev, aLev))
		dds.res_factor_level.shrunk <- lfcShrink(counts.dds.dsq, coef=paste(fact, aLev, "vs", ref_lev, sep = "_"), type="apeglm") %>% as.data.frame()

		dds.res_factor_level.shrunk$geneid <- rownames(dds.res_factor_level.shrunk)
		dds.res_factor_level.shrunk<-inner_join(dds.res_factor_level.shrunk, counts.df%>% select(Geneid, Length), by=c("geneid"="Geneid")) %>% mutate(expression = baseMean/Length) %>% select(-c(Length))
		rownames(dds.res_factor_level.shrunk) <- dds.res_factor_level.shrunk$geneid
		dds.res_factor_level.shrunk <- dds.res_factor_level.shrunk %>% select(-c(geneid))
		names(dds.res_factor_level.shrunk ) <- paste0(names(dds.res_factor_level.shrunk),paste(".",fact,".vs_",aLev,".apeglm", sep=""))
		dds.res_factor_level.shrunk$geneid <- rownames(dds.res_factor_level.shrunk)

		dds.dsq.byFactor.df <- full_join(dds.dsq.byFactor.df, dds.res_factor_level.shrunk, by = c("geneid"="geneid"))
#		counts.dds.res_factor_level <- results(counts.dds.dsq, contrast=c(factor1, ref_level, alt_level1)) 
#		counts.dds.res_factor_level$geneid <- rownames(counts.dds.res_factor_level)

	}

}

#Write to TSV
dds.dsq.byFactor.df <- dds.dsq.byFactor.df %>% select(geneid, everything())
write_delim(dds.dsq.byFactor.df, paste(dsq_pre, "itemized.de", sep="."), delim = "\t")


