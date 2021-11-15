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

### Deanhardt et al. 2021


The results summary from Deanhardt et al. 2021 can be generated: (untested!)

```bash

nohup time snakemake --restart-times 6 --latency-wait 60 --jobs 24 -p --cluster "sbatch --time={params.runtime} -n {params.cores} --mem={params.runmem_gb}G "  --snakefile Snakefile.Deanhardt2021  --configfile config.Deanhardt2021.yaml &

```

This recursively generates the differential expression analysis from raw data, summarizes/vizualizes it in a PDF, and bundles/timestamps the results summary (".pdf")

 The underlying differential expression and differential exon use data from Deanhardt et al. 2021 can be generated:

```bash
snakemake diff_expr/grpWtVs47b/grpWtVs47b.vs_dm6main.dm6_genes.mapspliceMulti.MpBC.itemized.de diff_expr/grpWtVs67d/grpWtVs67d.vs_dm6main.dm6_genes.mapspliceMulti.MpBC.itemized.de diff_expr/grpWtVsFru_smolFru/grpWtVsFru_smolFru.vs_dm6main.dm6_genes.mapspliceMulti.MpBC.itemized.de diff_expr/grpWtVsMut/grpWtVsMut.vs_dm6main.dm6_genes.mapspliceMulti.MpBC.itemized.de 
snakemake diff_exon_use/dex_grpWtVs47b/grpWtVs47b.vs_dm6main.dm6_genes.mapspliceMulti.M.de diff_exon_use/dex_grpWtVs67d/grpWtVs67d.vs_dm6main.dm6_genes.mapspliceMulti.M.de diff_exon_use/dex_grpWtVsFru_smolFru/grpWtVsFru_smolFru.vs_dm6main.dm6_genes.mapspliceMulti.M.de diff_exon_use/dex_grpWtVsMut/grpWtVsMut.vs_dm6main.dm6_genes.mapspliceMulti.M.de
```


### Full Results Summary

The full results summary ("VolkanLab_BehaviorGenetics.05_Apr_2021.pdf") was generated simply by running 

```bash
snakemake --cluster "sbatch --time={params.runtime} -n {params.cores} --mem={params.runmem_gb}G "

```

### Bypassing the Hard Parts

Users may wish to avoid the time and memory consuming steps of read downloading, cleaning, mapping, and counting. The raw counts underlying the results of Deanhardt et al. 2021 are stored in the NCBI GEO ( https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE179213 ); a shortcut script is included to download this data and prepare it such that it's where the Snakemake workflow anticipates. Currently only the "multi" alignment is available. 

The shortcut script is run thus:

```bash
bash utils/shortcut.bash
```

Fruitless replicate #1 was problematic and not used in the final analysis; it is not included in the download by default. To include it, run:

```bash
bash utils/shortcut.bash -f
```

The main time bottleneck is the "rando" alignment, which downsamples but does not entirely remove multimapping reads. Alignement strategy ultimately had little impact on results, so users may wish to omit them by generating the "noRando" PDF (untested!)

```bash

nohup time snakemake --restart-times 6 --latency-wait 60 --jobs 24 -p --cluster "sbatch --time={params.runtime} -n {params.cores} --mem={params.runmem_gb}G "  --snakefile Snakefile.Deanhardt2021  --configfile config.Deanhardt2021.yaml results/VolkanLabBehaviorGenetics.Deanhardt2021.noRando.pdf &

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


