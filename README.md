# Behavioral Genetics in Drosophila: RNA-Seq analysis pipeline

Software used to analyze data for [to be published]

## Dependencies

The following were required to run the full workflow:
```bash
module load python/3.6.6 bedtools bedops samtools r/3.6.0 rstudio/1.1.453 bowtie sratoolkit subread
```

as well, the Snakefile includes rules with hard paths that will need changing:
fastp_clean_sample_pe	:	fastp
write_report	:	pandoc_path


## Data

Not all the data in the config.yaml have been published; so full results summary won't run :( Intermediate results can be generated though.


## Use

The full results summary ("VolkanLab_BehaviorGenetics.05_Apr_2021.pdf") was generated simply by running 

```bash
snakemake
```

This recursively generates the differential expression analysis from raw data, summarizes/vizualizes it in a PDF, and bundles/timestamps the results summary.

## Code of Note


├── config.yaml
│ └── <i>configuration file for the pipeline</i>
├── dev
│ ├── ...
│ └── <i>files from earlier stages in project development, including an older FAIRE-seq experiment</i>
├── misc
│ ├── ...
│ └── <i>mostly metadata re: sequencing</i>
├── README.md
├── scripts
│ ├── bam_summarizer.mapspliced.py
│ ├── collapsedSignalMerger.py
│ ├── correlate.R
│ ├── deSeqer.R
│ │ └── <i>differential expression testing script</i>
│ ├── dexSeqe.R
│ │ └── <i>differential exon use testing script</i>
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
│   └── <i>summarize results</i>
├── Snakefile
│ └── <i>rules file for the pipeline</i>
├── utils
│ └── annotations
│   ├── ...
│   └── <i>includes lists of genes of interest</i>
└── VolkanLab_BehaviorGenetics.05_Apr_2021.pdf
  └── <i>pipeline walkthrough and results summary</i>







