configfile: 'config.yaml'
#	module load python/3.6.6 bedtools bedops samtools r/3.5.0 rstudio/1.1.453 bowtie sratoolkit

#from itertools import combinations

sample_by_name = {c['name'] : c for c in config['data_sets']}
ref_genome_by_name = { g['name'] : g for g in config['reference_genomes']}
annotation_by_name = { a['name'] : a for a in config['annotations']}
genelist_by_name = { l['name'] : l for l in config['gene_queries']}
#distance_categories = ["I","II","III","IV","V","VI","VII","VIII","IX","X"]
#distance_breaks = { distance_categories[n] : [config["distance_divisions"][n]['distance'],config["distance_divisions"][n+1]['distance']] for n in range(0, len(config["distance_divisions"])-1) }
# samplist_by_experimental = {c['experimental'] : [] for c in config['data_sets']}
# for c in config['data_sets']:
# 	samplist_by_experimental[c["experimental"]].append(c)

# samples_by_put = { c['put'] : [] for c in config['data_sets'] }
# for c in config['data_sets']:
# 	samples_by_put[c['put']].append(c)

sampname_by_group = {}
for s in sample_by_name.keys():
	subgroup_lst = sample_by_name[s]['subgroups']
	for g in subgroup_lst:
		if g in sampname_by_group.keys():
			sampname_by_group[g].append(s)
		else:
			sampname_by_group[g] = [s]

sra_by_fqpath = {}
for s in sample_by_name.keys():
	if "SRA" in sample_by_name[s].keys():
		sra_by_fqpath[sample_by_name[s]["path"]] = sample_by_name[s]["SRA"]



def return_filename_by_sampname(sampname):
	filenames = []
	if sample_by_name[sampname]['paired']:
		filenames.append(sample_by_name[sampname]['readsfile1'])
		filenames.append(sample_by_name[sampname]['readsfile2'])
	else:
		filenames.append(sample_by_name[sampname]['readsfile'])
	return filenames


def return_file_relpath_by_sampname(wildcards):
	sampname = wildcards.samplename
	pathprefix = sample_by_name[sampname]["path"]
	filesin = return_filename_by_sampname(sampname)
	pathsout = ["".join([pathprefix, fq]) for fq in filesin]
	return pathsout

def partition_experimental_by_put(exp):
	partition = {"in":[], "out":[]}
	for samp in samplist_by_experimental[exp]:
		partition[samp["put"]].append(samp)
	return partition



rule all:
	input: 
		pdf_out="results/VolkanLab_BehaviorGenetics.pdf",
	params:
		runmem_gb=1,
		runtime="0:01:00",
		cores=1,
	run:
		shell(""" mv results/VolkanLab_BehaviorGenetics.pdf results/VolkanLab_BehaviorGenetics.$(date | cut -f 2,3,6 -d " " | awk '{{print$2"_"$1"_"$3}}').pdf """)
		shell(""" tar cf results.$(date | tr -s " " | cut -f 2,3,6 -d " " | awk '{{print$2"_"$1"_"$3}}').tar results/ """)


rule summon_reads_SRA:
	output: 
		reads1='FASTQ/{path}/{prefix}_1.fastq',
		reads2='FASTQ/{path}/{prefix}_2.fastq',
	params:
		runmem_gb=8,
		runtime="3:00:00",
		cores=1,
	run:
		print("%s%s" % tuple(["FASTQ/", wildcards.path]))
		try:
			sra = sra_by_fqpath["%s%s/" % tuple(["FASTQ/", wildcards.path])]
			shell(""" mkdir -p FASTQ/{wildcards.path}/ """)
			shell("""
				fasterq-dump  --split-3 --outdir FASTQ/{wildcards.path}/ {sra}
			""")
		except KeyError:
			raise KeyError("Sample is listed as empirical but no reads found and no SRA to download!" )      		

rule mapsplice_builder:
	output:
		windowed='utils/MapSplice-v2.1.8/mapsplice.py.bak'
	params:
		runmem_gb=8,
		runtime="30:00",
		cores=1,
	run:
		shell("mkdir -p utils/")
		shell(
			""" 
			cd utils;
			wget http://protocols.netlab.uky.edu/~zeng/MapSplice-v2.1.8.zip ;
			unzip MapSplice-v2.1.8.zip;
			2to3 -w MapSplice-v2.1.8/;
			cd MapSplice-v2.1.8/;
			make;
			cd ..;
			rm MapSplice-v2.1.8.zip;
			"""
		)		



rule window_maker:
	output:
		windowed='utils/{ref_genome}_w{window_size}_s{slide_rate}.windows.bed'
	params:
		runmem_gb=64,
		runtime="30:00",
		cores=2,
	run:
		fai_path = ref_genome_by_name[wildcards.ref_genome]['fai'],
		shell("mkdir -p utils")
		shell(
			'bedtools makewindows -w {wildcards.window_size} -s {wildcards.slide_rate} -g {fai_path} -i winnum | bedtools sort -i - > {output.windowed}'
		)

rule annotation_subsetter:
	output:
		annot_out = "utils/annotations/{subname}.bed"
	params:
		runmem_gb=1,
		runtime="5:00",
		cores=1,
	run:
		shell("mkdir -p utils/annotations")
		base_path = annotation_by_name[genelist_by_name[wildcards.subname]["base"]]["bed_path"]
		list_file = genelist_by_name[wildcards.subname]["list_path"]
		shell(""" grep -wFf <( tail -n +2 {list_file} | cut -f 2 ) {base_path} > {output.annot_out} """)

rule reference_genome_reporter:
	input:
		fai_in = lambda wildcards: ref_genome_by_name[wildcards.ref_gen]['fai'],
	output:
		report_out = "summaries/reference_genomes/{ref_gen}.fai.report"
	params:
		runmem_gb=1,
		runtime="5:00",
		cores=1,
	shell:
		"""
		mkdir -p summaries/reference_genomes/
		cat {input.fai_in} | awk '{{sum+=$2}} END {{ print "number_contigs\t",NR; print "number_bases\t",sum}}' | sed -e 's/^/{wildcards.ref_gen}\t/g' > {output.report_out};
		"""

rule demand_reference_genome_summary:
	input:
		refgen_reports = lambda wildcards: expand("summaries/reference_genomes/{ref_gen}.fai.report", ref_gen=ref_genome_by_name.keys())
	output:
		refgen_summary = "summaries/reference_genomes.summary"
	params:
		runmem_gb=1,
		runtime="5:00",
		cores=1,
	shell:
		"cat {input.refgen_reports} > {output.refgen_summary}"


rule annotation_reporter:
	input:
		annot = lambda wildcards: annotation_by_name[wildcards.annot_name]["bed_path"]
	output:
		report_out = "summaries/annotations/{annot_name}.stats"
	params:
		runmem_gb=1,
		runtime="5:00",
		cores=1,
	run:
		shell(""" mkdir -p summaries/annotations/ """)
		shell(""" rm -f {output.report_out} """)
		shell(""" cat {input.annot} | cut -f 1 | grep -v "chr2110000222\|chrmitochondrion\|Un\|rand\|chrrDNA" | sort | uniq -c | tr -s " " | tr " " "\t" | awk '{{print"count\t"$2"\t"$1}}' >> {output.report_out} """)
		shell(""" cat {input.annot} | wc -l | awk '{{print"count\ttotal\t"$0}}' >> {output.report_out} """)
		shell(""" cat {input.annot} | awk '{{print$3-$2;}}' | awk '{{sum+=$1; sumsq+=$1*$1}} END {{ print "size\ttotal\t",sum; print "size\tavg\t",sum/NR; print "size\tstd\t",sqrt(sumsq/NR - (sum/NR)**2)}}'  >> {output.report_out} """)


rule demand_annotation_summaries:
	input:
		refgen_reports = lambda wildcards: expand("summaries/annotations/{ref_ann}.stats", ref_ann= [a["name"] for a in config['annotations'] if not a["derived"] ] ) # annotation_by_name.keys())
	output:
		refann_summary = "summaries/reference_annotations.summary"
	params:
		runmem_gb=1,
		runtime="5:00",
		cores=1,
	run:
		print([a["name"] for a in config['annotations'] if not a["derived"] ])
		shell(""" rm -f {output.refann_summary} """)
		for anne in [a["name"] for a in config['annotations'] if not a["derived"] ]:
			shell(""" cat summaries/annotations/{anne}.stats | awk '{{print"{anne}\t"$0}}' >> {output.refann_summary}""")



rule geneListReporter:
	input:
		list_in = lambda wildcards: genelist_by_name[wildcards.listicle]["list_path"],
		bed_in = lambda wildcards: genelist_by_name[wildcards.listicle]["bed_path"],
	output:
		statsOut = "summaries/geneLists/{listicle}.stats",
	params:
		runmem_gb=1,
		runtime="1:00",
		cores=1,
	run:
		shell(""" rm -f {output.statsOut} """)
		shell(""" echo "count	total	$( tail -n +2 {input.list_in} | wc -l  )" >> {output.statsOut} """)
		annot_path = annotation_by_name[genelist_by_name[wildcards.listicle]["base"]]["bed_path"]
		annot_name = annotation_by_name[genelist_by_name[wildcards.listicle]["base"]]["name"]
		shell(""" echo "count\tannotated\t$(cat {input.bed_in} | wc -l )">> {output.statsOut} """)
		shell(""" cat {input.bed_in} | awk '{{print$3-$2;}}' | awk '{{sum+=$1; sumsq+=$1*$1}} END {{ print "size\ttotal\t",sum; print "size\tavg\t",sum/NR; print "size\tstd\t",sqrt(sumsq/NR - (sum/NR)**2)}}'  >> {output.statsOut} """)


rule reportAllGeneLists:
	input:
		statsIn = expand( "summaries/geneLists/{listicle}.stats", listicle=genelist_by_name.keys()),
	output:
		allStats = "summaries/geneLists.stats",
	params:
		runmem_gb=1,
		runtime="1:00",
		cores=1,
	run:
		for listicle in genelist_by_name.keys():
			annot = annotation_by_name[genelist_by_name[listicle]["base"]]["name"]
			shell(""" cat summaries/geneLists/{listicle}.stats | awk '{{print"{listicle}\t{annot}\t"$0}}' >> {output.allStats} """)



rule fastp_clean_sample_pe:
	input:
		fileIn = lambda wildcards: return_file_relpath_by_sampname(wildcards)
	output:
		fileOut = ["{pathprefix}/{samplename}.clean.R1.fastq","{pathprefix}/{samplename}.clean.R2.fastq"],
		jason = "{pathprefix}/{samplename}.json"
	params:
		runmem_gb=8,
		runtime="6:00:00",
		cores=1,
		#--trim_front1 and -t, --trim_tail1
		#--trim_front2 and -T, --trim_tail2. 
		common_params = "--json {pathprefix}/{samplename}.json",# --html meta/FASTP/{samplename}.html", 
		pe_params = "--detect_adapter_for_pe --correction",
	message:
		"FASTP QA/QC on paired-ended reads ({wildcards.samplename}) in progress.... "
	shell:
		"/nas/longleaf/home/csoeder/modules/fastp/fastp {params.common_params} {params.pe_params} --in1 {input.fileIn[0]} --out1 {output.fileOut[0]} --in2 {input.fileIn[1]} --out2 {output.fileOut[1]}"



rule FASTP_summarizer:
	input: 
		jason = lambda wildcards: expand("{path}{samp}.json", path=sample_by_name[wildcards.samplename]['path'], samp = wildcards.samplename, )
	output:
		jason_pruned = "summaries/FASTP/{samplename}.json.pruned"
	params:
		runmem_gb=1,
		runtime="5:00",
		cores=1,
	message:
		"Summarizing reads for sample ({wildcards.samplename}) .... "	
	shell:
		"""
		cp {input.jason} summaries/FASTP/{wildcards.samplename}.json
		python3 scripts/fastp_reporter.py {input.jason} {output.jason_pruned} -t {wildcards.samplename}
		"""

rule demand_FASTQ_analytics:	#forces a FASTP clean
	input:
		jasons_in = expand("summaries/FASTP/{samplename}.json.pruned", samplename = sample_by_name.keys())
	output:
		summary = "summaries/sequenced_reads.dat"
	params:
		runmem_gb=1,
		runtime="1:00",
		cores=1,
	message:
		"Collecting read summaries for all samples ...."
	shell:
		"cat {input.jasons_in} > {output.summary}"



rule mapsplice2_align_raw:
	input:
		mapsplice = "utils/MapSplice-v2.1.8/mapsplice.py.bak",
		reads_in = lambda wildcards: expand("{path}{sample}.clean.R{pair}.fastq", path=sample_by_name[wildcards.sample]['path'], sample=wildcards.sample, pair = [1,2]),
	output:
		raw_bam = "mapped_reads/mapspliceRaw/{sample}/{sample}.vs_{ref_genome}.mapspliceRaw.sort.bam",
#		bam_out = "mapped_reads/mapsplice/{sample}/{sample}.vs_{ref_genome}.mapsplice.sort.bam",
#		splutflg = "utils/mapsplut.{sample}.vs_{ref_genome}.flg"
	params:
		runmem_gb=32,
		runtime="16:00:00",
		cores=8,
	message:
		"aligning reads from {wildcards.sample} to reference_genome {wildcards.ref_genome} .... "
	run:
		ref_genome_split = ref_genome_by_name[wildcards.ref_genome]['split'],
		ref_genome_bwt = ref_genome_by_name[wildcards.ref_genome]['bowtie'],
		shell(""" mkdir -p mapped_reads/mapspliceRaw/{wildcards.sample}/ summaries/BAMs/""")
		shell(""" 
			python utils/MapSplice-v2.1.8/mapsplice.py -c {ref_genome_split} -x {ref_genome_bwt} -1 {input.reads_in[0]} -2 {input.reads_in[1]} -o mapped_reads/mapspliceRaw/{wildcards.sample} 
			""")
		shell(""" 
			samtools view -Sbh mapped_reads/mapspliceRaw/{wildcards.sample}/alignments.sam | samtools sort - > mapped_reads/mapspliceRaw/{wildcards.sample}/{wildcards.sample}.vs_{wildcards.ref_genome}.mapspliceRaw.sort.bam;
			samtools index {output.raw_bam};
			rm mapped_reads/mapspliceRaw/{wildcards.sample}/alignments.sam;
			""")
#			mv mapped_reads/mapsplice/{wildcards.sample}/stats.txt summaries/BAMs/{wildcards.sample}.vs_{wildcards.ref_genome}.mapsplice.fresh.stats;
#			mv mapped_reads/mapsplice/{wildcards.sample}/{wildcards.sample}.vs_{wildcards.ref_genome}.mapsplice.dupStats summaries/BAMs/{wildcards.sample}.vs_{wildcards.ref_genome}.mapsplice.dupStats;

#			samtools view -Sbh alignments.sam | samtools sort -n | samtools fixmate -m | samtools sort - | samtools markdup {params.dup_flg} -s - {wildcards.sample}.vs_{wildcards.ref_genome}.mapsplice.sort.bam 2> {wildcards.sample}.vs_{wildcards.ref_genome}.mapsplice.dupStats;

#		



rule mapsplice2_align_multi:
	input:
		raw_bam = "mapped_reads/mapspliceRaw/{sample}/{sample}.vs_{ref_genome}.mapspliceRaw.sort.bam",
	output:
		multi_bam = "mapped_reads/mapspliceMulti/{sample}/{sample}.vs_{ref_genome}.mapspliceMulti.sort.bam",
#		bam_out = "mapped_reads/mapsplice/{sample}/{sample}.vs_{ref_genome}.mapsplice.sort.bam",
#		splutflg = "utils/mapsplut.{sample}.vs_{ref_genome}.flg"
	params:
		runmem_gb=8,
		runtime="1:00:00",
		cores=8,
		dup_flg = "-rS ",#remove marked duplicates; remove secondary alignments of marked duplicates
		quality="-q 20 -F 0x0200 -F 0x04 -f 0x0002", # no QC fails, no duplicates, proper pairs only, mapping qual >20
	message:
		"filtering raw {wildcards.sample} alignment to {wildcards.ref_genome} for quality, duplication.... "
	run:
		shell(""" mkdir -p mapped_reads/mapspliceMulti/{wildcards.sample}/ summaries/BAMs/""")
		shell(""" 
			samtools sort -n {input.raw_bam} | samtools fixmate -m - - | samtools sort - | samtools markdup {params.dup_flg} - - | samtools view -bh {params.quality} > {output.multi_bam} ;
			samtools index {output.multi_bam} 
			""")



rule mapsplice2_align_uniq:
	input:
		multi_bam = "mapped_reads/mapspliceMulti/{sample}/{sample}.vs_{ref_genome}.mapspliceMulti.sort.bam",
	output:
		uniq_bam = "mapped_reads/mapspliceUniq/{sample}/{sample}.vs_{ref_genome}.mapspliceUniq.sort.bam",
	params:
		runmem_gb=8,
		runtime="1:00:00",
		cores=8,
	message:
		"filtering raw {wildcards.sample} alignment to {wildcards.ref_genome} for quality, duplication.... "
	run:
		shell(""" mkdir -p mapped_reads/mapspliceUniq/{wildcards.sample}/ summaries/BAMs/""")
		shell(""" 
			samtools view {input.multi_bam} | grep -w "IH:i:1" | cat <( samtools view -H {input.multi_bam} ) - | samtools view -Sbh > {output.uniq_bam};
			samtools index {output.uniq_bam}
			""")
### filtering on uniqueness: MapSplice2 uses the IH:i:1 flag in the SAM specs. 

rule mapsplice2_align_rando:
	input:
		multi_bam = "mapped_reads/mapspliceMulti/{sample}/{sample}.vs_{ref_genome}.mapspliceMulti.sort.bam",
		uniq_bam = "mapped_reads/mapspliceUniq/{sample}/{sample}.vs_{ref_genome}.mapspliceUniq.sort.bam",
	output:
		rando_bam = "mapped_reads/mapspliceRando/{sample}/{sample}.vs_{ref_genome}.mapspliceRando.sort.bam",
	params:
		runmem_gb=16,
		runtime="12:00:00",
		cores=8,
	message:
		"filtering raw {wildcards.sample} alignment to {wildcards.ref_genome} for quality, duplication.... "
	run:
		shell(""" mkdir -p mapped_reads/mapspliceRando/{wildcards.sample}/ summaries/BAMs/""")
		shell(""" 

			samtools sort -n {input.multi_bam} | samtools view | grep -v "IH:i:[01]" > {output.rando_bam}.multi.unsampled

			rm -f {output.rando_bam}.multi.sampled
			for idx in $(cut -f 1 {output.rando_bam}.multi.unsampled | sort | uniq); do 
				grep -w $idx {output.rando_bam}.multi.unsampled | sort  --random-sort | head -n 1 >> {output.rando_bam}.multi.sampled
			done;

			cat <( samtools view -h {input.uniq_bam} ) {output.rando_bam}.multi.sampled | samtools view -Sbh | samtools sort - > {output.rando_bam};
			samtools index {output.rando_bam}
			rm {output.rando_bam}.multi.unsampled {output.rando_bam}.multi.sampled
			""")


rule spliced_alignment_reporter:
	input:
		bam_in = "mapped_reads/{aligner}/{sample}/{sample}.vs_{ref_genome}.{aligner}.sort.bam",
	output:
		report_out = "summaries/BAMs/{sample}.vs_{ref_genome}.{aligner}.summary"
	params:
		runmem_gb=8,
		runtime="4:00:00",
		cores=1,
	message:
		"Collecting metadata for the {wildcards.aligner} alignment of {wildcards.sample} to {wildcards.ref_genome}.... "
	run:
 		ref_genome_idx=ref_genome_by_name[wildcards.ref_genome]['fai']
 		shell(""" which samtools """)
 		shell("samtools idxstats {input.bam_in} > {input.bam_in}.idxstats")
		shell("samtools flagstat {input.bam_in} > {input.bam_in}.flagstat")
		shell("bedtools genomecov -max 1 -i <( bedtools bamtobed -i {input.bam_in} | cut -f 1-3 ) -g {ref_genome_idx} > {input.bam_in}.span.genomcov")
		shell("bedtools genomecov -max 1 -i <( bedtools bamtobed -split -i {input.bam_in} | cut -f 1-3 ) -g {ref_genome_idx} > {input.bam_in}.split.genomcov")
# 		#change the -max flag as needed to set 
		shell(""" samtools depth -a {input.bam_in} | awk '{{sum+=$3; sumsq+=$3*$3}} END {{ print "average_depth\t",sum/NR; print "std_depth\t",sqrt(sumsq/NR - (sum/NR)**2)}}' > {input.bam_in}.dpthStats """)
# 		#https://www.biostars.org/p/5165/
# 		#save the depth file and offload the statistics to the bam_summarizer script?? 
		### filtering on uniqueness: MapSplice2 uses the IH:i:1 flag in the SAM specs. 
		shell("""samtools view -F 4 {input.bam_in} | cut -f 13  | cut -f 3 -d ":" | sort | uniq -c | tr -s " " | awk -F " " '{{print$2"\t"$1}}' > {input.bam_in}.mapmult""")
		shell(""" which samtools """)
		shell(""" samtools sort -n {input.bam_in} | samtools fixmate -m - -  | samtools sort - | samtools markdup -sS - - 1> /dev/null 2> {input.bam_in}.dupe """)# /dev/null 2> {input.bam_in}.dupe;""")
		shell("python3 scripts/bam_summarizer.mapspliced.py -f {input.bam_in}.flagstat -i {input.bam_in}.idxstats -g {input.bam_in}.split.genomcov -G {input.bam_in}.span.genomcov -d {input.bam_in}.dpthStats -D {input.bam_in}.dupe -m {input.bam_in}.mapmult -o {output.report_out} -t {wildcards.sample}")

# cat stats.pass1.txt  | grep -A4 "By mapping uniqueness" | grep -v "By mapping uniqueness" | cut -f 2 | paste <(echo -e "total_reads\nuniquely_mapped\nmultimapped\nunmapped") - 
# cat stats.pass1.txt | grep paired_reads | grep -v fusion | cut -f 1,2
# cat stats.pass1.txt | grep -A 3 "alignment types" | tail -n +2 | head -n 2 | sed -e 's/spliced/spliced_alignments/g'
# cat stats.pass1.txt | grep -A 3 "Splice junctions" | tail -n +2 | head -n 2 | cut -f 1,2
#https://github.com/broadinstitute/rnaseqc
#awk '($6 ~ /N/)' | wc -l

### filtering on uniqueness: MapSplice2 uses the IH:i:1 flag in the SAM specs. 


rule demand_BAM_analytics:
	input:
		bam_reports = lambda wildcards: expand("summaries/BAMs/{sample}.vs_{ref_genome}.{aligner}.summary", sample=sampname_by_group['all'], ref_genome=wildcards.ref_genome, aligner=wildcards.aligner)
	output:
		full_report = "summaries/alignments.vs_{ref_genome}.{aligner}.summary"
	params:
		runmem_gb=1,
		runtime="1:00",
		cores=1,
	message:
		"collecting all alignment metadata.... "
	shell:
		"""cat {input.bam_reports} | awk '{{print "{wildcards.ref_genome}\t{wildcards.aligner}\t"$0}}'> {output.full_report}"""









rule write_report:
	input:
		reference_genome_summary = ["summaries/reference_genomes.summary"],
		reference_annotation_summary = ["summaries/reference_annotations.summary"],
		gene_lists_summary = ["summaries/geneLists.stats"],
		sequenced_reads_summary=["summaries/sequenced_reads.dat"],
		aligned_reads_summary = expand("summaries/alignments.vs_dm6main.{aligner}.summary", aligner=["mapspliceRaw","mapspliceMulti"])#,"mapspliceRando","mapspliceUniq",]),
	output:
		pdf_out="results/VolkanLab_BehaviorGenetics.pdf",
	params:
		runmem_gb=8,
		runtime="1:00:00",
		cores=2,
	message:
		"writing up the results.... "
	run:
		pandoc_path="/nas/longleaf/apps/rstudio/1.0.136/bin/pandoc"
		pwd = subprocess.check_output("pwd",shell=True).decode().rstrip()+"/"
		shell("""mkdir -p results/figures/supp/ results/tables/supp/""")
		shell(""" R -e "setwd('{pwd}');Sys.setenv(RSTUDIO_PANDOC='{pandoc_path}')" -e  "peaDubDee='{pwd}'; rmarkdown::render('scripts/RNAseq_results.Rmd',output_file='{pwd}{output.pdf_out}')"  """)
#		shell(""" tar cf results.tar results/ """)



# rule geneList_ontologist:
# 	input:
# 		gene_list = "features/genelists/closest/{lust_prefix}.list"
# 	output:
# 		list_ontology = "features/ontologies/{lust_prefix}.go"
# 	params:
# 		runmem_gb=8,
# 		runtime="1:00:00",
# 		cores=8,
# 	run:
# 		shell(""" mkdir -p features/ontologies/ """)
# 		shell(""" Rscript scripts/geneOntologer.R {input.gene_list} {output.list_ontology} """)
