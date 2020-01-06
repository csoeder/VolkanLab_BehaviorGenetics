#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library("readr")
library("biomaRt")
library("org.Dm.eg.db")
library("topGO")
library("dplyr")

list_file <- read_delim(args[1], "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
query_genes.list <- list_file$X1 %>% unique() %>% as.list()


#ensembl_dm = useMart("ensembl",dataset="dmelanogaster_gene_ensembl")
marty <- useDataset("dmelanogaster_gene_ensembl",  useMart("ensembl",  host = "useast.ensembl.org") )
G_list <- getBM(attributes= c("flybase_gene_id", "ensembl_gene_id", "entrezgene_id", "go_id"), mart= marty) 

### plz fix this
go_bad <- c("GO:0110109","GO:0120176","GO:0120177","GO:0120170","GO:0062023")
bad_boiz <- G_list %>% filter(go_id %in% go_bad) %>% dplyr::select("flybase_gene_id")
bad_boiz <- rbind(bad_boiz, "FBgn0050046", "FBgn0033710",  "FBgn0050203", "FBgn0026721")
G_list.unq <- G_list %>% filter(!flybase_gene_id %in% bad_boiz$flybase_gene_id ) %>% dplyr::select("flybase_gene_id") %>% unique()
### ### ### ### ### ### ### ### 


gene_list <- factor(as.integer( G_list.unq$flybase_gene_id %in% query_genes.list ))
names(gene_list) <- G_list.unq$flybase_gene_id




jeanOnt.df <- as.data.frame(c())

for (ont in c("MF","BP","CC") ){
	
	GOdata <- new("topGOdata", ontology = ont, allGenes = gene_list, annot = annFUN.org, mapping = "org.Dm.eg.db", ID="Ensembl")
	fishy <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
	goRes.tmp <- GenTable(GOdata, classicFisher =fishy, ranksOf = "classicFisher", topNodes = 50) %>% mutate(ontology = as.factor(ont), classicFisher = as.numeric(classicFisher)) %>% select(-c(Term), Term) %>% as_tibble() 
	jeanOnt.df <- rbind(jeanOnt.df, goRes.tmp)

}

write.table(jeanOnt.df, file=args[2], row.names=FALSE, sep = "\t", quote=F)
