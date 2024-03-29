configfile: 'config.yaml'
#	module load python/3.6.6 bedtools bedops samtools r/4.0.3 rstudio/1.1.453 bowtie sratoolkit subread
####	module load python/3.6.6 bedtools bedops samtools r/3.6.0 rstudio/1.1.453 bowtie sratoolkit subread

sample_by_name = {c['name'] : c for c in config['data_sets']}
ref_genome_by_name = { g['name'] : g for g in config['reference_genomes']}
annotation_by_name = { a['name'] : a for a in config['annotations']}
genelist_by_name = { l['name'] : l for l in config['gene_queries']}
contrasts_by_name = {c['name'] : c for c in config['contrasts']}


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
		pdf_out="results/VolkanLabBehaviorGenetics.full.pdf",
	params:
		runmem_gb=1,
		runtime="0:01:00",
		cores=1,
	run:
		shell(""" mkdir -p results/figures/; touch results/figures/null.png; for fig in results/figures/*png; do mv $fig $(echo $fig| rev | cut -f 2- -d . | rev ).$(date +%d_%b_%Y).png; done;  rm results/figures/null.*.png; """)
		shell(""" mkdir -p results/figures/supp/ ; touch results/figures/supp/null.png; for fig in results/figures/supp/*png; do mv $fig $(echo $fig| rev | cut -f 2- -d . | rev ).$(date +%d_%b_%Y).png; done; rm results/figures/supp/null.*.png; """)

		shell(""" mkdir -p results/tables/ ; touch results/tables/null.tmp ; for phial in $(ls -p results/tables/ | grep -v /); do pre=$(echo $phial | rev | cut -f 2- -d . | rev ); suff=$(echo $phial | rev | cut -f 1 -d . | rev ); mv results/tables/$phial results/tables/$pre.$(date +%d_%b_%Y).$suff; done ; rm results/tables/null.*.tmp; """)
		shell(""" mkdir -p results/tables/supp/ ; touch results/tables/supp/null.tmp ; for phial in $(ls -p results/tables/supp/ | grep -v /); do pre=$(echo $phial | rev | cut -f 2- -d . | rev ); suff=$(echo $phial | rev | cut -f 1 -d . | rev ); mv results/tables/supp/$phial results/tables/supp/$pre.$(date +%d_%b_%Y).$suff; done ; rm results/tables/supp/null.*.tmp; """)

		shell(""" mv results/VolkanLabBehaviorGenetics.full.pdf results/VolkanLab_BehaviorGenetics.full.$(date +%d_%b_%Y).pdf """)
		shell(""" tar cf results.full.$(date +%d_%b_%Y).tar results/ """)


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


rule fruIsoids:
	input:
		fruGene_in = "utils/annotations/fru.gtf",
	output:
		fruIsoid_out = "utils/annotations/fru.isoids.gtf",
	params:
		runmem_gb=8,
		runtime="1:00",
		cores=8,
	message:
		"running the edgeHog on Fru .... "
	run:
		shell(""" python scripts/edgeHog.py --gtf_in $(pwd)/utils/annotations/fru.gtf -o utils/annotations/fru -n FBgn0004652 """)


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
	params:
		runmem_gb=16,
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



rule mapsplice2_align_multi:
	input:
		raw_bam = "mapped_reads/mapspliceRaw/{sample}/{sample}.vs_{ref_genome}.mapspliceRaw.sort.bam",
	output:
		multi_bam = "mapped_reads/mapspliceMulti/{sample}/{sample}.vs_{ref_genome}.mapspliceMulti.sort.bam",
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
		"filtering raw {wildcards.sample} alignment to {wildcards.ref_genome} to remove multimappers entirely.... "
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
		runtime="150:00:00",
		cores=8,
	message:
		"Downsampling multimappers from {wildcards.sample} alignment to {wildcards.ref_genome}..."
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


rule alignmentSpliceExtractor:
	input:
		in_bam = "mapped_reads/{aligner}/{sample}/{sample}.vs_{ref_genome}.{aligner}.sort.bam",
	output:
		out_bam = "mapped_reads/{aligner}_SpliceOnly/{sample}/{sample}.vs_{ref_genome}.{aligner}_SpliceOnly.sort.bam",
	params:
		runmem_gb=8,
		runtime="6:00:00",
		cores=8,
	message:
		"Downsampling multimappers from {wildcards.sample} alignment to {wildcards.ref_genome}..."
	run:
		shell(""" mkdir -p mapped_reads/{wildcards.aligner}_SpliceOnly/{wildcards.sample}/ summaries/BAMs/""")
		shell(""" 
			cat <(samtools view -H {input.in_bam}) <(samtools view {input.in_bam}| awk '($6 ~ /N/)' | less) | samtools view -hbS > {output.out_bam}.tmp;
			samtools index {output.out_bam}.tmp;
			python scripts/read_Hollower.py -i {output.out_bam}.tmp -o {output.out_bam}.tmp.sam;
			cat {output.out_bam}.tmp.sam |  sed -e 's/\tN\t/\t*\t/g' | samtools view -hbS > {output.out_bam};
			samtools index {output.out_bam};
			rm {output.out_bam}.tmp {output.out_bam}.tmp.sam;
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
		shell(""" samtools view -F 0x4 {input.bam_in} | cut -f 1 | sort | uniq | wc -l > {input.bam_in}.maptCount """)
		shell(""" samtools sort -n {input.bam_in} | samtools fixmate -m - -  | samtools sort - | samtools markdup -sS - - 1> /dev/null 2> {input.bam_in}.dupe """)# /dev/null 2> {input.bam_in}.dupe;""")
		shell("python3 scripts/bam_summarizer.mapspliced.py -f {input.bam_in}.flagstat --mapped_count {input.bam_in}.maptCount -i {input.bam_in}.idxstats -g {input.bam_in}.split.genomcov -G {input.bam_in}.span.genomcov -d {input.bam_in}.dpthStats -D {input.bam_in}.dupe -m {input.bam_in}.mapmult -o {output.report_out} -t {wildcards.sample}")

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


rule extract_fru_coverage:
	input:
		bam_in = "mapped_reads/{aligner}/{sample}/{sample}.vs_{ref_genome}.{aligner}.sort.bam",
	output:
		bg_out = "summaries/coverage/fru.{aligner}.{sample}.{ref_genome}.bedgraph"
	params:
		runmem_gb=1,
		runtime="1:00",
		cores=1,
		location = "chr3R:18,413,579-18,546,928"
	message:
		"collecting all alignment metadata.... "
	run:
		fai = ref_genome_by_name[wildcards.ref_genome]['fai']
		shell(""" mkdir -p summaries/coverage/ """)
		shell(""" samtools view -hb {input.bam_in} {params.location} | bedtools genomecov -split -bga -ibam - -g {fai} > {output.bg_out} """)



rule summon_fru_coverages:
	input:
		beegees = lambda wildcards: expand("summaries/coverage/fru.mapspliceMulti.{sample}.dm6main.bedgraph",sample= sampname_by_group['hausWtVsMut'] ),
	output:
		bg_out = "utils/frubeds.graphed"
	params:
		runmem_gb=1,
		runtime="1:00",
		cores=1,
		location = "chr3R:18,413,579-18,546,928"
	message:
		"collecting all alignment metadata.... "
	run:
		shell(""" touch {output} """)



rule count_features:
	input:
		alignments_in = lambda wildcards: expand("mapped_reads/{aligner}/{sample}/{sample}.vs_{ref_genome}.{aligner}.sort.bam", aligner = wildcards.aligner, sample=sampname_by_group[wildcards.group], ref_genome = wildcards.ref_genome  ),
		annot_in = lambda wildcards: annotation_by_name[wildcards.annot]["gtf_path"]
	output:
		counted_features = "counts/{group}.vs_{ref_genome}.{annot}.{aligner}.{flags}.counts",
		count_stats = "counts/{group}.vs_{ref_genome}.{annot}.{aligner}.{flags}.counts.summary",
	params:
		runmem_gb=8,
		runtime="64:00:00",
		cores=8,
		fc_params = " -J --verbose ",#count multimappers, record junctions, count paired-end fragments, well-mapped   <- some of this is moved into the "flag" wildcard
		#-O ? -f? 
	message:
		"counting reads from {wildcards.group} aligned to {wildcards.ref_genome} with {wildcards.aligner} overlapping {wildcards.annot} .... "
	run:
		flug_str = wildcards.flags
		flug_lst = [ s for s in flug_str if s != "_"]
		flug_str = " -%s "*len(flug_lst) % tuple(flug_lst)


		sed_suff = ".vs_%s.%s.sort.bam" % tuple([wildcards.ref_genome, wildcards.aligner])		
		sed_cmd =  " sed -e 's/mapped_reads\/[a-zA-Z0-9\/_]*\///g' | sed -e 's/%s//g' " % tuple([sed_suff])

		annot_gtf = annotation_by_name[wildcards.annot]["gtf_path"]
		shell(""" mkdir -p counts summaries/counts/""")
		shell("""
			featureCounts {flug_str} {params.fc_params} -t exon -g gene_id -F GTF -a <(cat {annot_gtf} | awk '{{print"chr"$0}}') -o {output.counted_features}.tmp {input.alignments_in}
			""")
		shell("""
			cat counts/{wildcards.group}.vs_{wildcards.ref_genome}.{wildcards.annot}.{wildcards.aligner}.{wildcards.flags}.counts.tmp | tail -n +2 | {sed_cmd} > {output.counted_features}
			cat counts/{wildcards.group}.vs_{wildcards.ref_genome}.{wildcards.annot}.{wildcards.aligner}.{wildcards.flags}.counts.tmp.jcounts | {sed_cmd} > counts/{wildcards.group}.vs_{wildcards.ref_genome}.{wildcards.annot}.{wildcards.aligner}.{wildcards.flags}.counts.jcounts
			cat counts/{wildcards.group}.vs_{wildcards.ref_genome}.{wildcards.annot}.{wildcards.aligner}.{wildcards.flags}.counts.tmp.summary | {sed_cmd} > counts/{wildcards.group}.vs_{wildcards.ref_genome}.{wildcards.annot}.{wildcards.aligner}.{wildcards.flags}.counts.summary
			rm counts/{wildcards.group}.vs_{wildcards.ref_genome}.{wildcards.annot}.{wildcards.aligner}.{wildcards.flags}.counts.tmp* 
			""")


rule summarize_feature_counts:
	input:
		sum_in = "counts/{group}.vs_{ref_genome}.{annot}.{aligner}.{flag}.counts.summary"
	output:
		count_sum = "summaries/{group}.vs_{ref_genome}.{annot}.{aligner}.{flag}.counts.stat",
	params:
		runmem_gb=1,
		runtime="15:00",
		cores=1,
	message:
		"consolidating read-counting summary of {wildcards.group} aligned to {wildcards.ref_genome} with {wildcards.aligner} overlapping {wildcards.annot} .... "
	run:
		shell("""
			head -n 1 {input.sum_in} > {output.count_sum};
			cat {input.sum_in} | grep -w "Assigned\|Unassigned_NoFeatures\|Unassigned_Ambiguity" >> {output.count_sum};
			""")


rule measure_expression:
	input:
		counted_features = "counts/{group}.vs_{ref_genome}.{annot}.{aligner}.{flags}.counts",
		aln_rprt = "summaries/alignments.vs_{ref_genome}.{aligner}.summary",
	output:
		rpkm = "expression/{group}.vs_{ref_genome}.{annot}.{aligner}.{flags}.rpkm",
		tpm="expression/{group}.vs_{ref_genome}.{annot}.{aligner}.{flags}.tpm",
	params:
		runmem_gb=1,
		runtime="15:00",
		cores=1,
	message:
		"converting counts of {wildcards.group} aligned to {wildcards.ref_genome} with {wildcards.aligner} overlapping {wildcards.annot} to expression values.... "
	run:
		shell(""" mkdir -p expression/ """)
		shell(""" Rscript  scripts/expressionOmete.R {input.counted_features} {input.aln_rprt} expression/{wildcards.group}.vs_{wildcards.ref_genome}.{wildcards.annot}.{wildcards.aligner}.{wildcards.flags} """)


rule summon_expression_measures:
	input:
		rpkm =  lambda wildcards: expand("expression/{group}.vs_{ref_genome}.{annot}.{aligner}.{flags}.rpkm",group = wildcards.group, ref_genome = "dm6main", annot = ["dm6_genes","fru_exons"], aligner=["mapspliceMulti", "mapspliceUniq","mapspliceRando"], flags = ["MpBC","MpBCO"])
	output:
		expFlg = "utils/{group}.expression.flag",
	params:
		runmem_gb=1,
		runtime="1:00",
		cores=1,
	message:
		"converting counts of {wildcards.group} aligned to (dm6main) with (aligner) overlapping (annot) to expression values.... "
	run:
		shell(""" mkdir -p utils/ """)
		shell(""" touch {output.expFlg} """)
## replace this with a data summary script?


rule da_Seeker:
	input:
		fc_in = "counts/all.vs_{ref_genome}.{annot}.{aligner}.{flag}.counts",#, group =["all"])
	output:
		deSq = "diff_expr/{contrast}/{contrast}.vs_{ref_genome}.{annot}.{aligner}.{flag}.de",
		deSq_itemized = "diff_expr/{contrast}/{contrast}.vs_{ref_genome}.{annot}.{aligner}.{flag}.itemized.de",
	params:
		runmem_gb=8,
		runtime="30:00",
		cores=8,
	message:
		"calculating differential expression in {wildcards.contrast} contrast ({wildcards.aligner} to {wildcards.ref_genome} overlapping {wildcards.annot}) .... "
	run:
		#cont_sbgrp = contrasts_by_name[wildcards.contrast]["filt"]
		shell(""" mkdir -p diff_expr/{wildcards.contrast} """)
		shell(""" Rscript scripts/deSeqer.R {input.fc_in} diff_expr/{wildcards.contrast}/{wildcards.contrast}.vs_{wildcards.ref_genome}.{wildcards.annot}.{wildcards.aligner}.{wildcards.flag} {wildcards.contrast} config.yaml """)
		#shell(""" mkdir -p meta/DESeq2_methods/ """)
		#shell(""" mv diff_expr/{wildcards.contrast}/{wildcards.contrast}.vs_{wildcards.ref_genome}.{wildcards.annot}.{wildcards.aligner}.{wildcards.flag}.method meta/DESeq2_methods/ """)



rule de_GOntologer:
	input:
		deSq_itemized = "diff_expr/{contrast}/{contrast}.vs_{ref_genome}.{annot}.{aligner}.{flag}.itemized.de",
	output:
		deSq_go = "gene_ont/{contrast}/{contrast}.vs_{ref_genome}.{annot}.{aligner}.{flag}.go",
	params:
		runmem_gb=8,
		runtime="30:00",
		cores=8,
	message:
		"calculating differential expression in {wildcards.contrast} contrast ({wildcards.aligner} to {wildcards.ref_genome} overlapping {wildcards.annot}) .... "
	run:
		#cont_sbgrp = contrasts_by_name[wildcards.contrast]["filt"]
		shell(""" mkdir -p gene_ont/{wildcards.contrast}/ """)
		shell(""" Rscript scripts/geneOntologe.R {input.deSq_itemized} {wildcards.contrast} {output.deSq_go} config.yaml """)




rule da_Correlator:
	input:
		xpr_in = "expression/{group}.vs_{ref_genome}.{annot}.{aligner}.{flag}.rpkm",#, group =["all"])
	output:
		corr_out = "correlations/{correlation}/{group}.vs_{ref_genome}.{annot}.{aligner}.{flag}.corr.tsv"
	params:
		runmem_gb=64,
		runtime="30:00",
		cores=8,
	message:
		"calculating differential correlation in {wildcards.correlation} contrast ({wildcards.aligner} to {wildcards.ref_genome} overlapping {wildcards.annot}) .... "
	run:
		#cont_sbgrp = contrasts_by_name[wildcards.contrast]["filt"]
		shell(""" mkdir -p correlations/{wildcards.correlation}/ """)
		shell(""" Rscript scripts/correlate.R {input.xpr_in} {wildcards.correlation} correlations/{wildcards.correlation}/{wildcards.group}.vs_{wildcards.ref_genome}.{wildcards.annot}.{wildcards.aligner}.{wildcards.flag} """)


#dexseq_prepper
#~/R/libs/DEXSeq/python_scripts/dexseq_count.py 

rule DEXSeq_builder:
	output:
		scripts = ["utils/DEXSeq/Subread_to_DEXSeq/dexseq_prepare_annotation2.py","utils/DEXSeq/Subread_to_DEXSeq/load_SubreadOutput.R"]
#		scripts=['utils/DEXSeq/dexseq_count.py','utils/DEXSeq/dexseq_prepare_annotation.py']
	params:
		runmem_gb=8,
		runtime="30:00",
		cores=1,
	run:
		shell(
			"""
			mkdir -p utils/DEXSeq/;
			""")

#			dexwd=$(Rscript <(echo "cat(system.file( 'python_scripts', package='DEXSeq' ))" ));
#			2to3 --verbose -n -W --output-dir=utils/DEXSeq/ $dexwd ;
		shell("""
			cd utils/DEXSeq/;
			git clone https://github.com/vivekbhr/Subread_to_DEXSeq;
			""")
		#patch for wonky transcripts:https://stat.ethz.ch/pipermail/bioconductor/2012-June/046494.html
		shell("""
			cp utils/dexseq_prepare_annotation2.patch.py  utils/DEXSeq/Subread_to_DEXSeq/dexseq_prepare_annotation2.py
			""")

rule DEXSeq_annot_prepper:
	input:
		scripts='utils/DEXSeq/Subread_to_DEXSeq/dexseq_prepare_annotation2.py',
		anne = lambda wildcards: annotation_by_name[wildcards.anne]["gtf_path"]
	output:
		nu_gff = "utils/annotations/DEXSeq/{anne}.dxsqReady.gff",
		fc_ready = "utils/annotations/DEXSeq/{anne}.dxsqReady.gtf"
	params:
		runmem_gb=8,
		runtime="10:00",
		cores=1,
		prep_params = " --aggregate no",
	run:
		shell(
			"""	
				mkdir -p utils/annotations/DEXSeq/;
				python utils/DEXSeq/Subread_to_DEXSeq/dexseq_prepare_annotation2.py {params.prep_params} <( cat {input.anne}| awk '{{print"chr"$0}}' ) -f {output.fc_ready} {output.nu_gff}
			""")

rule DEXSeq_annot_counter:#flag = M
	input:
		dxq_ann = "utils/annotations/DEXSeq/{anne}.dxsqReady.gtf",
		alignments_in = lambda wildcards: expand("mapped_reads/{aligner}/{sample}/{sample}.vs_{ref_genome}.{aligner}.sort.bam", aligner = wildcards.aligner, sample=sampname_by_group[wildcards.group], ref_genome = wildcards.ref_genome  ),
	output:
		counted_features = "dxCounts/{group}.vs_{ref_genome}.{anne}.{aligner}.{flags}.dxsq.counts",
		count_stats = "dxCounts/{group}.vs_{ref_genome}.{anne}.{aligner}.{flags}.dxsq.counts.summary",
	params:
		runmem_gb=8,
		runtime="8:00:00",
		cores=8,
		fc_params = " -s 0 -f -O -p -B -C ",#count multimappers, record junctions, count paired-end fragments, well-mapped   <- some of this is moved into the "flag" wildcard
	message:
		"counting reads from {wildcards.group} aligned to {wildcards.ref_genome} with {wildcards.aligner} overlapping the DEXSeq-squashed {wildcards.anne} .... "
	run:
		flug_str = wildcards.flags
		flug_lst = [ s for s in flug_str if s != "_"]
		flug_str = " -%s "*len(flug_lst) % tuple(flug_lst)

		sed_suff = ".vs_%s.%s.sort.bam" % tuple([wildcards.ref_genome, wildcards.aligner])		
		sed_cmd =  " sed -e 's/mapped_reads\/[a-zA-Z0-9\/_]*\///g' | sed -e 's/%s//g' " % tuple([sed_suff])

		shell(""" mkdir -p dxCounts summaries/dxCounts/""")
		shell("""
			featureCounts {flug_str} {params.fc_params} -t exon -g gene_id -F GTF -a {input.dxq_ann} -o {output.counted_features}.tmp {input.alignments_in}
			""")
		shell("""
			cat {output.counted_features}.tmp | tail -n +2 | {sed_cmd} > {output.counted_features}
			cat {output.counted_features}.tmp.summary | {sed_cmd} > {output.count_stats}
			rm {output.counted_features}.tmp* 
			""")
		##cat {output.counted_features}.tmp.jcounts | {sed_cmd} > counts/DEXSeq/{wildcards.group}.vs_{wildcards.ref_genome}.{wildcards.anne}.{wildcards.aligner}.{wildcards.flags}.counts.jcounts



rule deX_seeker:#needs v3.5.0 lol
	input:
		dxq_ann = "utils/annotations/DEXSeq/{anne}.dxsqReady.gtf",
		counted_features = "dxCounts/{group}.vs_{ref_genome}.{anne}.{aligner}.{flags}.dxsq.counts",
	output:
		dexSq ="diff_exon_use/{contrast}/{group}.vs_{ref_genome}.{anne}.{aligner}.{flags}.de",
	params:
		runmem_gb=8,
		runtime="30:00",
		cores=1,
	message:
		"calculating differential expression in {wildcards.contrast} contrast ({wildcards.aligner} to {wildcards.ref_genome} overlapping {wildcards.anne}) .... "
	run:
		#cont_sbgrp = contrasts_by_name[wildcards.contrast]["filt"]
		shell(""" mkdir -p diff_exon_use/{wildcards.contrast}/ """)
		shell("""  Rscript scripts/dexSeqe.R {input.counted_features} diff_exon_use/{wildcards.contrast}/{wildcards.group}.vs_{wildcards.ref_genome}.{wildcards.anne}.{wildcards.aligner}.{wildcards.flags} {input.dxq_ann} {wildcards.contrast} 'chr3R_FBgn0004652-'  """)
#		shell(""" module unload r; module load r/3.5.0; Rscript scripts/dexSeqe.R {input.counted_features} diff_exon_use/{wildcards.contrast}/{wildcards.group}.vs_{wildcards.ref_genome}.{wildcards.anne}.{wildcards.aligner}.{wildcards.flags} {input.dxq_ann} {wildcards.contrast} 'chr3R_FBgn0004652-' ; module unload r; module load r/3.6.0; """)




rule write_report:
	input:
		reference_genome_summary = ["summaries/reference_genomes.summary"],
		reference_annotation_summary = ["summaries/reference_annotations.summary"],
		gene_lists_summary = ["summaries/geneLists.stats"],
		sequenced_reads_summary=["summaries/sequenced_reads.dat"],
		aligned_reads_summary = expand("summaries/alignments.vs_dm6main.{aligner}.summary", aligner=["mapspliceRaw","mapspliceMulti","mapspliceUniq","mapspliceRando"]),
		counting_summary_genes = expand("summaries/{group}.vs_{ref_genome}.{annot}.{aligner}.{flag}.counts.stat", group ="all", ref_genome = "dm6main", annot = ["dm6_genes","fru_exons"], flag=["MpBC","MpBCO"], aligner=["mapspliceMulti","mapspliceUniq","mapspliceRando"]),#"mapspliceRaw",
#		counting_summary_exons = expand("summaries/{group}.vs_{ref_genome}.{annot}.{aligner}.{flag}.counts.stat", group ="all", ref_genome = "dm6main", annot = ["fru_junct", "fru_intron"], flag=["MpBCO"], aligner=["mapspliceMulti_SpliceOnly","mapspliceUniq_SpliceOnly","mapspliceRando_SpliceOnly"]),#"mapspliceRaw",
#		expFlg = "utils/all.expression.flag",
#		diff_exprs = expand("diff_expr/{contrast}/{contrast}.vs_dm6main.{annot}.{aligner}.{flag}.{suff}", contrast = ["grpWtVs47b","grpWtVs67d","grpWtVsFru","grpWtVsMut","hausWtVsMut","wildTypeHousing"], annot = ["dm6_genes",], flag= ["MpBC", "MpBCO"],aligner=["mapspliceMulti","mapspliceUniq","mapspliceRando"], suff = ["de"]),#,"itemized.de"]),
#		diff_item = expand("diff_expr/{contrast}/{contrast}.vs_dm6main.{annot}.{aligner}.{flag}.{suff}", contrast = ["hausWtVsMut", ], annot = ["dm6_genes","fru_exons"], flag= ["MpBC", "MpBCO"],aligner=["mapspliceMulti","mapspliceUniq","mapspliceRando"], suff = ["itemized.de"]),#,""]),
#		fruQC = expand("diff_expr/{contrast}/{contrast}.vs_dm6main.{annot}.{aligner}.{flag}.{suff}", contrast = ["hausWtVsMut_smolFru", "hausWtVsMut_noFru"], annot = ["dm6_genes",], flag= ["MpBC", "MpBCO"],aligner=["mapspliceMulti","mapspliceUniq","mapspliceRando"], suff = ["itemized.de"]),#,""]),
#		diff_item2 = expand("diff_expr/{contrast}/{contrast}.vs_dm6main.{annot}.{aligner}.{flag}.{suff}", contrast = ["2_days_difference", "47b_on_88a", "cantonAmos", "cantonWt", "cantonAmosWt"], annot = ["dm6_genes",], flag= ["MpBC", ],aligner=["mapspliceMulti","mapspliceUniq","mapspliceRando"], suff = ["itemized.de"]),#,""]),
#		diff_fru = expand("diff_expr/{contrast}/{contrast}.vs_dm6main.{annot}.{aligner}.{flag}.{suff}", contrast = ["hausWtVsMut", ], annot = ["fru_junct", "fru_intron", "fru_exons"], flag= [ "MpBCO"],aligner=["mapspliceMulti_SpliceOnly","mapspliceUniq_SpliceOnly","mapspliceRando_SpliceOnly"], suff = ["itemized.de"]),#,""]),
#		gene_onts = expand("gene_ont/{contrast}/{contrast}.vs_dm6main.dm6_genes.mapsplice{aligner}.MpBC.go", aligner = ["Multi","Rando","Uniq"], contrast = ["hausWtVsMut","47b_on_88a","2_days_difference", "cantonAmos", "cantonWt", "cantonAmosWt"]),
#		dxsq_readr = "utils/DEXSeq/Subread_to_DEXSeq/load_SubreadOutput.R",
	output:
		pdf_out="results/VolkanLabBehaviorGenetics.full.pdf",
	params:
		runmem_gb=64,
		runtime="2:00:00",
		cores=8,
	message:
		"writing up the results.... "
	run:
		pandoc_path="/nas/longleaf/apps/rstudio/1.0.136/bin/pandoc"
		pwd = subprocess.check_output("pwd",shell=True).decode().rstrip()+"/"
		shell("""mkdir -p results/figures/supp/ results/tables/supp/""")
		shell(""" R -e "setwd('{pwd}');Sys.setenv(RSTUDIO_PANDOC='{pandoc_path}')" -e  "peaDubDee='{pwd}'; rmarkdown::render('scripts/RNAseq_results.Rmd',output_file='{pwd}{output.pdf_out}')"  """)
#		shell(""" tar cf results.tar results/ """)


