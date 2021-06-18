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

```html
├── config.yaml<br>
│ └── <i>configuration file for the pipeline</i><br>
├── dev<br>
│ ├── ...<br>
│ └── <i>files from earlier stages in project development, including an older FAIRE-seq experiment</i><br>
├── misc<br>
│ ├── ...<br>
│ └── <i>mostly metadata re: sequencing</i><br>
├── README.md<br>
├── scripts<br>
│ ├── bam_summarizer.mapspliced.py<br>
│ ├── collapsedSignalMerger.py<br>
│ ├── correlate.R<br>
│ ├── deSeqer.R<br>
│ │ └── <i>differential expression testing script</i><br>
│ ├── dexSeqe.R<br>
│ │ └── <i>differential exon use testing script</i><br>
│ ├── edgeHog.py<br>
│ ├── expressionOmete.R<br>
│ ├── faire_results.Rmd<br>
│ ├── fastp_reporter.py<br>
│ ├── geneOntologe.R<br>
│ ├── multimap_spreader.py<br>
│ ├── overlapSignalCollapser.py<br>
│ ├── read_Hollower.py<br>
│ ├── references.bib<br>
│ ├── replicateCorrelations.png<br>
│ └── RNAseq_results.Rmd<br>
│   └── <i>summarize results</i><br>
├── Snakefile<br>
│ └── <i>rules file for the pipeline</i><br>
├── utils<br>
│ └── annotations<br>
│   ├── ...<br>
│   └── <i>includes lists of genes of interest</i><br>
└── VolkanLab_BehaviorGenetics.05_Apr_2021.pdf<br>
  └── <i>pipeline walkthrough and results summary</i><br>

```





