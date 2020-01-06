#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library("readr")
library("biomaRt")
library("org.Dm.eg.db")
library("topGO")
library("dplyr")
library("yaml")
library("stringr")
library("magrittr")
data("geneList")

deSeq_file <- args[1]
nombre <- args[2]
go_out <- args[3]
#nombre <- "hausWtVsMut"


sig_thresh <- 0.01

#hausWtVsMut_vs_dm6main_dm6_genes_mapspliceMulti_MpBC_itemized <- read_delim("diff_expr/hausWtVsMut/hausWtVsMut.vs_dm6main.dm6_genes.mapspliceMulti.MpBC.itemized.de", "\t", escape_double = FALSE, trim_ws = TRUE)
DESeq.itemized<- read_delim(deSeq_file, "\t", escape_double = FALSE, trim_ws = TRUE)

trammel <- read_yaml("config.yaml")
#	build the contrast experiment
contrasts.df <- plyr::ldply(trammel$contrasts, data.frame) %>% filter(name == nombre )
#get the factors in the design
desgn <- as.character(contrasts.df$design[1])
desgn.splt <- (desgn %>% str_replace_all("~","")%>% str_replace_all(" ","") %>% strsplit("+", fixed=TRUE))[[1]]

query_genes.list <- DESeq.itemized$geneid %>% unique() %>% as.list()

#ensembl_dm = useMart("ensembl",dataset="dmelanogaster_gene_ensembl")
marty <- useDataset("dmelanogaster_gene_ensembl",  useMart("ensembl",  host = "useast.ensembl.org") )
G_list <- getBM(attributes= c("flybase_gene_id", "ensembl_gene_id", "external_gene_name", "go_id", "name_1006"), mart= marty) %>% mutate(go_name = name_1006) %>% select(-c("name_1006"))


### ### ### ### ### ### ### ### 
### plz fix this
go_bad <- c("GO:0110109","GO:0120176","GO:0120177","GO:0120170","GO:0062023")
bad_boiz <- G_list %>% filter(go_id %in% go_bad) %>% dplyr::select("flybase_gene_id")
bad_boiz <- rbind(bad_boiz, "FBgn0050046", "FBgn0033710",  "FBgn0050203", "FBgn0026721")
#mask the baddies
G_list.unq <- G_list %>% filter(!flybase_gene_id %in% bad_boiz$flybase_gene_id ) %>% dplyr::select("flybase_gene_id") %>% unique()
DESeq.itemized %<>% filter(!(geneid %in% bad_boiz$flybase_gene_id))
### ### ### ### ### ### ### ### 

#present_gene_list <- G_list.unq$flybase_gene_id %in% DESeq.itemized$geneid %>% as.numeric()
#names(present_gene_list) <- G_list.unq$flybase_gene_id

jeanOnt.mega.df <- as.data.frame(c())

#for each factor in factors....

for (fact in desgn.splt){

	ref_lev <- (contrasts.df %>% filter(vars == fact))$ref
	#all_levs<-levels(counts.dds.dsq[[fact]])

	#subset the DE data to the factor
	DESeq.itemized.factor <- DESeq.itemized[,!is.na(str_locate(names(DESeq.itemized), paste(fact,".vs_" ,sep=""))[,1])]
	#pull out the alt levels from the data
	alt_levs<- (str_split_fixed(names(DESeq.itemized.factor), paste(fact,".vs_", sep=""), n=2)[,2] %>% str_split_fixed("\\.", n=2))[,1] %>% unique()
	DESeq.itemized.factor$geneid <- DESeq.itemized$geneid

	#	#for each alt_level in factor.levels:

	for (alt in alt_levs) {
		print(alt)
		#subset the DE data to that level
		DESeq.itemized.level <- DESeq.itemized.factor[,!is.na(str_locate(names(DESeq.itemized.factor), paste(fact,"\\.vs_",alt,"\\." ,sep=""))[,1])]
		DESeq.itemized.level$geneid <- DESeq.itemized.factor$geneid
		#pull out the genes & significance values
		sig_gene_list <- DESeq.itemized.level[paste("padj.",fact,".vs_",alt,".apeglm",sep="")][[1]]# < sig_thresh)[,1] #%>% as.numeric()
		sig_gene_list[is.na(sig_gene_list)] <- 1
		names(sig_gene_list) <- DESeq.itemized.level$geneid

		jeanOnt.df <- as.data.frame(c())

		for (ont in c("MF","BP","CC") ){

			#topDiffGenes uses p cutoff of 0.01
			GOdata <- new("topGOdata", ontology = ont, allGenes = sig_gene_list, geneSel = topDiffGenes, annot = annFUN.org, mapping = "org.Dm.eg.db", ID="Ensembl")
			fishy <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
			kissy <- runTest(GOdata, algorithm = "classic", statistic = "ks")

			goRes.tmp <- GenTable(GOdata, classicFisher =fishy, classicKS = kissy,  topNodes = 50) %>% mutate(ontology = as.factor(ont), classicFisher = as.numeric(classicFisher), classicKS = as.numeric(classicKS)) %>% select(-c(Term), Term) %>% as_tibble() 
			goRes.tmp %<>% left_join(G_list%>% select(c("go_id", "go_name")) %>% unique(), by =c("GO.ID" = "go_id")) %>% select(-c("Term"))
			jeanOnt.df <- rbind(jeanOnt.df, goRes.tmp)
			#store & bind the topGO results. 

		}

		jeanOnt.df$variable <- as.factor(fact)
		jeanOnt.df$alt_level <- as.factor(alt)

		jeanOnt.mega.df <- rbind(jeanOnt.mega.df,jeanOnt.df)

	}

}

write.table(jeanOnt.mega.df, file=args[3], row.names=FALSE, col.names = TRUE, sep = "\t", quote=F)
