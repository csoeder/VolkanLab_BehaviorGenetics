# Behavioral Genetics in Drosophila: RNA-Seq analysis pipeline

Software used to analyze data for: "Changes in splicing and neuromodulatory gene expression programs in sensory neurons with pheromone signaling and social experience", Deanhardt et al. 2021, bioRxiv. doi: https://doi.org/10.1101/2021.06.18.449021
https://www.biorxiv.org/content/10.1101/2021.06.18.449021v2


## Dependencies

The following were required to run the full workflow:
```bash
module load python/3.6.6 bedtools bedops samtools r/3.6.0 rstudio/1.1.453 bowtie sratoolkit subread
```

as well, the Snakefile includes rules with hard paths that will need changing:
```
fastp_clean_sample_pe	:	fastp
write_report	:	pandoc_path
```


## Data

Not all the data in the config.yaml have been published; so full results summary won't run :( Intermediate results can be generated though.


## Use

The full results summary ("VolkanLab_BehaviorGenetics.05_Apr_2021.pdf") was generated simply by running 

```bash
snakemake
```

This recursively generates the differential expression analysis from raw data, summarizes/vizualizes it in a PDF, and bundles/timestamps the results summary.


The underlying differential expression and differential exon use data from Deanhardt et al. 2021 can be generated:

```bash
snakemake diff_expr/grpWtVs47b/grpWtVs47b.vs_dm6main.dm6_genes.mapspliceMulti.MpBC.itemized.de diff_expr/grpWtVs67d/grpWtVs67d.vs_dm6main.dm6_genes.mapspliceMulti.MpBC.itemized.de diff_expr/grpWtVsFru_smolFru/grpWtVsFru_smolFru.vs_dm6main.dm6_genes.mapspliceMulti.MpBC.itemized.de diff_expr/grpWtVsMut/grpWtVsMut.vs_dm6main.dm6_genes.mapspliceMulti.MpBC.itemized.de 
snakemake diff_exon_use/dex_grpWtVs47b/grpWtVs47b.vs_dm6main.dm6_genes.mapspliceMulti.M.de diff_exon_use/dex_grpWtVs67d/grpWtVs67d.vs_dm6main.dm6_genes.mapspliceMulti.M.de diff_exon_use/dex_grpWtVsFru_smolFru/grpWtVsFru_smolFru.vs_dm6main.dm6_genes.mapspliceMulti.M.de diff_exon_use/dex_grpWtVsMut/grpWtVsMut.vs_dm6main.dm6_genes.mapspliceMulti.M.de
```

The results summary from Deanhardt et al. 2021 can be generated: (untested!)

```bash
snakemake --snakefile Snakefile.Deanhardt2021  --configfile config.Deanhardt2021.yaml 
```






## Code of Note

```html
├── config.yaml
│ └── (configuration file for the pipeline)
├── dev
│ ├── ...
│ └── (files from earlier stages in project development, including an older FAIRE-seq experiment)
├── misc
│ ├── ...
│ └── (mostly metadata re: sequencing)
├── README.md
├── quichen_duan
│ ├── ...
│ └── (scripts for followup analyses from Quichen Duan)
├── scripts
│ ├── bam_summarizer.mapspliced.py
│ ├── collapsedSignalMerger.py
│ ├── correlate.R
│ ├── deSeqer.R
│ │ └── (differential expression testing script)
│ ├── dexSeqe.R
│ │ └── (differential exon use testing script)
│ ├── edgeHog.py
│ ├── expressionOmete.R
│ ├── faire_results.Rmd
│ ├── fastp_reporter.py
│ ├── geneOntologe.R
│ ├── multimap_spreader.py
│ ├── overlapSignalCollapser.py
│ ├── read_Hollower.py
│ ├── references.bib
│ ├── replicateCorrelations.png
│ └── RNAseq_results.Rmd
│   └── (summarize results)
├── Snakefile
│ └── (rules file for the pipeline)
├── utils
│ └── annotations
│   ├── ...
│   └── (includes lists of genes of interest)
└── VolkanLab_BehaviorGenetics.05_Apr_2021.pdf
  └── (pipeline walkthrough and results summary)

```





