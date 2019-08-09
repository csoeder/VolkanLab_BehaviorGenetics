configfile: 'config.yaml'
#	module load python/3.6.6 fseq bedtools bedops samtools bwa r/3.5.0 rstudio/1.1.453
#PATH=$PATH:/nas/longleaf/home/csoeder/modules/gem_old/GEM_prerelease3/bin:/nas/longleaf/home/csoeder/modules/bffBuilder/bin

#from itertools import combinations

sample_by_name = {c['name'] : c for c in config['data_sets']}
ref_genome_by_name = { g['name'] : g for g in config['reference_genomes']}
annotation_by_name = { a['name'] : a for a in config['annotations']}
genelist_by_name = { l['name'] : l for l in config['gene_queries']}
distance_categories = ["I","II","III","IV","V","VI","VII","VIII","IX","X"]
distance_breaks = { distance_categories[n] : [config["distance_divisions"][n]['distance'],config["distance_divisions"][n+1]['distance']] for n in range(0, len(config["distance_divisions"])-1) }




sampname_by_group = {}
for s in sample_by_name.keys():
	subgroup_lst = sample_by_name[s]['subgroups']
	for g in subgroup_lst:
		if g in sampname_by_group.keys():
			sampname_by_group[g].append(s)
		else:
			sampname_by_group[g] = [s]


rule all:
	input: 
		peaks_in = lambda wildcards: expand("fSeq/{sample}.vs_{ref_genome}.{aligner}.calledPeaks.bed", ref_genome = "dm6", sample = sample_by_name.keys(), aligner = "bwa" )
#	output:
#		dummy = "test.txt"
#	params:
#		runmem_gb=1,
#		runtime="0:01:00",
#		cores=1,
#	shell:
#		"touch {output.dummy}"





def return_file_relpath_by_sampname(wildcards):
	sampname = wildcards.samplename
	pathprefix = sample_by_name[sampname]["path"]
	filesin = sample_by_name[sampname]['readsfile']
	pathsout = ["".join([pathprefix, filesin])]
	return pathsout


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
		report_out = "meta/reference_genomes/{ref_gen}.fai.report"
	params:
		runmem_gb=1,
		runtime="5:00",
		cores=1,
	shell:
		"""
		mkdir -p meta/reference_genomes/
		cat {input.fai_in} | awk '{{sum+=$2}} END {{ print "number_contigs\t",NR; print "number_bases\t",sum}}' | sed -e 's/^/{wildcards.ref_gen}\t/g' > {output.report_out};
		"""

rule demand_reference_genome_summary:
	input:
		refgen_reports = lambda wildcards: expand("meta/reference_genomes/{ref_gen}.fai.report", ref_gen=ref_genome_by_name.keys())
	output:
		refgen_summary = "meta/reference_genomes.summary"
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
		report_out = "meta/annotations/{annot_name}.stats"
	params:
		runmem_gb=1,
		runtime="5:00",
		cores=1,
	run:
		shell(""" mkdir -p meta/annotations/ """)
		shell(""" rm -f {output.report_out} """)
		shell(""" cat {input.annot} | cut -f 1 | grep -v "chr2110000222\|chrmitochondrion\|Un\|rand\|chrrDNA" | sort | uniq -c | tr -s " " | tr " " "\t" | awk '{{print"count\t"$2"\t"$1}}' >> {output.report_out} """)
		shell(""" cat {input.annot} | wc -l | awk '{{print"count\ttotal\t"$0}}' >> {output.report_out} """)
		shell(""" cat {input.annot} | awk '{{print$3-$2;}}' | awk '{{sum+=$1; sumsq+=$1*$1}} END {{ print "size\ttotal\t",sum; print "size\tavg\t",sum/NR; print "size\tstd\t",sqrt(sumsq/NR - (sum/NR)**2)}}'  >> {output.report_out} """)


rule demand_annotation_summaries:
	input:
		refgen_reports = lambda wildcards: expand("meta/annotations/{ref_ann}.stats", ref_ann= [a["name"] for a in config['annotations'] if not a["derived"] ] ) # annotation_by_name.keys())
	output:
		refann_summary = "meta/reference_annotations.summary"
	params:
		runmem_gb=1,
		runtime="5:00",
		cores=1,
	run:
		print([a["name"] for a in config['annotations'] if not a["derived"] ])
		shell(""" rm -f {output.refann_summary} """)
		for anne in [a["name"] for a in config['annotations'] if not a["derived"] ]:
			shell(""" cat meta/annotations/{anne}.stats | awk '{{print"{anne}\t"$0}}' >> {output.refann_summary}""")



rule geneListReporter:
	input:
		list_in = lambda wildcards: genelist_by_name[wildcards.listicle]["list_path"],
		bed_in = lambda wildcards: genelist_by_name[wildcards.listicle]["bed_path"],
	output:
		statsOut = "meta/geneLists/{listicle}.stats",
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
		statsIn = expand( "meta/geneLists/{listicle}.stats", listicle=genelist_by_name.keys()),
	output:
		allStats = "meta/geneLists.stats",
	params:
		runmem_gb=1,
		runtime="1:00",
		cores=1,
	run:
		for listicle in genelist_by_name.keys():
			annot = annotation_by_name[genelist_by_name[listicle]["base"]]["name"]
			shell(""" cat meta/geneLists/{listicle}.stats | awk '{{print"{listicle}\t{annot}\t"$0}}' >> {output.allStats} """)

rule fastp_clean_sample_se:
	input:
		fileIn = lambda wildcards: return_file_relpath_by_sampname(wildcards)
	output:
		fileOut = ["{pathprefix}/{samplename}.clean.fastq"],
		jason = "{pathprefix}/{samplename}.json",
	params:
		runmem_gb=8,
		runtime="3:00:00",
		cores=1,
		se_params = "",
		common_params = "--json {pathprefix}/{samplename}.json",# --html meta/FASTP/{samplename}.html", 
	message:
		"FASTP QA/QC on single-ended reads ({wildcards.samplename}) in progress.... "
	shell:
		"/nas/longleaf/home/csoeder/modules/fastp/fastp {params.common_params} {params.se_params} --in1 {input.fileIn[0]} --out1 {output.fileOut[0]}"


rule FASTP_summarizer:
	input: 
		jason = lambda wildcards: expand("{path}{samp}.json", path=sample_by_name[wildcards.samplename]['path'], samp = wildcards.samplename, )
	output:
		jason_pruned = "meta/FASTP/{samplename}.json.pruned"
	params:
		runmem_gb=1,
		runtime="5:00",
		cores=1,
	message:
		"Summarizing reads for sample ({wildcards.samplename}) .... "	
	shell:
		"""
		cp {input.jason} meta/FASTP/{wildcards.samplename}.json
		python3 scripts/fastp_reporter.py {input.jason} {output.jason_pruned} -t {wildcards.samplename}
		"""

rule demand_FASTQ_analytics:	#forces a FASTP clean
	input:
		jasons_in = expand("meta/FASTP/{samplename}.json.pruned", samplename = sample_by_name.keys())
	output:
		summary = "meta/sequenced_reads.dat"
	params:
		runmem_gb=1,
		runtime="1:00",
		cores=1,
	message:
		"Collecting read summaries for all samples ...."
	shell:
		"cat {input.jasons_in} > {output.summary}"



rule bwa_align:
	input:
		reads_in = lambda wildcards: expand("{path}{sample}.clean.fastq", path=sample_by_name[wildcards.sample]['path'], sample=wildcards.sample),
		ref_genome_file = lambda wildcards: ref_genome_by_name[wildcards.ref_genome]['path'],
	output:
		bam_out = "mapped_reads/{sample}.vs_{ref_genome}.bwa.sort.bam",
	params:
		runmem_gb=32,
		runtime="8:00:00",
		cores=8,
	message:
		"aligning reads from {wildcards.sample} to reference_genome {wildcards.ref_genome} .... "
	run:
		shell("bwa aln {input.ref_genome_file} {input.reads_in} > {input.reads_in}.sai ")
		shell("bwa samse {input.ref_genome_file} {input.reads_in}.sai {input.reads_in} | samtools view -Shb | samtools addreplacerg -r ID:{wildcards.sample} -r SM:{wildcards.sample} - | samtools sort -o {output.bam_out} - ")
		shell("samtools index {output.bam_out}")



rule bam_reporter:
	input:
		bam_in = "mapped_reads/{sample}.vs_{ref_genome}.{aligner}.sort.bam"
	output:
		report_out = "meta/BAMs/{sample}.vs_{ref_genome}.{aligner}.summary"
	params:
		runmem_gb=8,
		runtime="4:00:00",
		cores=1,
	message:
		"Collecting metadata for the {wildcards.aligner} alignment of {wildcards.sample} to {wildcards.ref_genome}.... "
	run:
		ref_genome_idx=ref_genome_by_name[wildcards.ref_genome]['fai']
		shell("samtools idxstats {input.bam_in} > {input.bam_in}.idxstats")
		shell("samtools flagstat {input.bam_in} > {input.bam_in}.flagstat")
		shell("bedtools genomecov -max 1 -ibam {input.bam_in} -g {ref_genome_idx} > {input.bam_in}.genomcov")
		#change the -max flag as needed to set 
		shell("""samtools depth -a {input.bam_in} | awk '{{sum+=$3; sumsq+=$3*$3}} END {{ print "average_depth\t",sum/NR; print "std_depth\t",sqrt(sumsq/NR - (sum/NR)**2)}}' > {input.bam_in}.dpthStats""")
		#https://www.biostars.org/p/5165/
		#save the depth file and offload the statistics to the bam_summarizer script?? 
		shell("python3 scripts/bam_summarizer.py -f {input.bam_in}.flagstat -i {input.bam_in}.idxstats -g {input.bam_in}.genomcov -d {input.bam_in}.dpthStats -o {output.report_out} -t {wildcards.sample}")


rule demand_BAM_analytics:
	input:
		bam_reports = lambda wildcards: expand("meta/BAMs/{sample}.vs_{ref_genome}.{aligner}.summary", sample=sampname_by_group['all'], ref_genome=wildcards.ref_genome, aligner=wildcards.aligner)
	output:
		full_report = "meta/alignments.vs_{ref_genome}.{aligner}.summary"
	params:
		runmem_gb=1,
		runtime="1:00",
		cores=1,
	message:
		"collecting all alignment metadata.... "
	shell:
		"cat {input.bam_reports} > {output.full_report}"



#build virtual environment rule??


rule peakCall_MACS:
	input:
		bam_in = "mapped_reads/{sample}.vs_{ref_genome}.bwa.sort.bam",

	output:
		peaks_out = "MACS/{sample}.vs_{ref_genome}_summits.bed"
	params:
		runmem_gb=8,
		runtime="6:00:00",
		cores=8,
	run:
		shell(""" mkdir -p MACS/ """)
		shell(""" set +u; source py2/venv/bin/activate; module load macs/2.1.2;  macs2 callpeak -t {input.bam_in} -f BAM --outdir MACS/ --name {wildcards.sample}.vs_{wildcards.ref_genome}; deactivate; """ )


rule gem_mappability:
	input:
		refgen = lambda wildcards: ref_genome_by_name[wildcards.ref_genome]["path"],
	output:
		map_out = "utils/gem/{ref_genome}.mappability.wig"
	params:
		runmem_gb=8,
		runtime="6:00:00",
		cores=8,
		read_len=100,
	run:
		shell("mkdir -p utils/gem/")
		shell("gem-indexer -i {input.refgen} -o utils/gem/{wildcards.ref_genome}.indx")
		shell("gem-mappability -I utils/gem/{wildcards.ref_genome}.indx.gem -o utils/gem/{wildcards.ref_genome} -l {params.read_len}")
		shell("gem-2-wig -I utils/gem/{wildcards.ref_genome}.indx.gem -i utils/gem/{wildcards.ref_genome}.mappability -o utils/gem/{wildcards.ref_genome}.mappability")


#fSeq needs a mappability track from GEMtools
#https://sourceforge.net/projects/gemlibrary/files/gem-library/
#PATH=$PATH:/nas/longleaf/home/csoeder/modules/gem_old/GEM_prerelease3/bin
#gem-indexer -i /proj/cdjones_lab/Genomics_Data_Commons/genomes/drosophila_melanogaster/dm6.main.fa -o gemTest 
#gem-mappability -I gemTest.gem -o mappableTest -l 100
#gem-2-wig -I gemTest.gem -i mappableTest.mappability -o mapTest
#bedtools makewindows -w 1 -s 1 -g /proj/cdjones_lab/Genomics_Data_Commons/genomes/drosophila_melanogaster/dm6.main.fa.fai > utils/dm6_wig_base.bed
#echo $header > wigTest.wig
#bedtools intersect -wa -wb -a utils/dm6_wig_base.bed -b <(cat mapTest.wig | wig2bed ) | cut -f 8 >> wigTest.wig 
#~/modules/bffBuilder/bin/bffBuilder wigTest.wig > test.bff

		# shell(
		# 	"""
		# 	for chrom in $( cat {input.ref_indx} | cut -f 1 ); do 
		# 		echo "fixedStep chrom=$chrom start=1 step=1" > utils/gem/{ref_genome}_background/"$chrom".wig ;
		# 		bedtools intersect -wa -wb -a {input.wig_base} -b <(cat utils/gem/{wildcards.ref_genome}.mappability | wig2bed ) | cut -f 8 >> utils/gem/{ref_genome}_background/"$chrom".wig ;
		# 		bffBuilder utils/gem/{ref_genome}_background/"$chrom".wig  > utils/gem/{ref_genome}_background/"$chrom".bff;
		# 		rm utils/gem/{ref_genome}_background/"$chrom".wig ;
		# 	done
		# 	"""
		# 	)


rule wig2bff:
	input:
		map_wig = "utils/gem/{ref_genome}.mappability.wig",
		wig_base = "utils/{ref_genome}_w1_s1.windows.bed",
	output:
		bff_out = "utils/{ref_genome}.background/{chrom}.bff",
	params:
		runmem_gb=16	,
		runtime="6:00:00",
		cores=8,
		read_len=100,
	run:
		shell("""mkdir -p utils/gem/{wildcards.ref_genome}.background/ """)
		shell("""echo "fixedStep chrom={wildcards.chrom} start=1 step=1" > utils/gem/{wildcards.ref_genome}.background/{wildcards.chrom}.wig ; """)
		shell("""bedtools intersect -wa -wb -a <(cat {input.wig_base} | grep -w {wildcards.chrom}) -b <(cat {input.map_wig} | wig2bed ) | cut -f 9 >> utils/gem/{wildcards.ref_genome}.background/{wildcards.chrom}.wig;""")
		shell("""  (cd utils/gem/{wildcards.ref_genome}.background/ && bffBuilder {wildcards.chrom}.wig); """)
		shell(""" mv utils/gem/{wildcards.ref_genome}.background/{wildcards.chrom}.bff {output.bff_out}""" )


rule gangsAllHere:
	input:
		lambda wildcards: expand("utils/{ref_gen}.background/{chrom}.bff", ref_gen = wildcards.ref_gen, chrom=[line.split("\t")[0] for line in open(ref_genome_by_name[wildcards.ref_gen]['fai'], "r").readlines()] )
	output:
		gang_flag = "utils/{ref_gen}.bffBuilt.flg"
	params:
		runmem_gb=1,
		runtime="00:30",
		cores=8,
		read_len=100,
	run:
		shell(""" touch {output} """)	


#http://fureylab.web.unc.edu/software/fseq/
#NPF specs: http://genome.ucsc.edu/FAQ/FAQformat.html#format12
rule peakCall_Fseq:
	input:
		bam_in = "mapped_reads/{sample}.vs_{ref_genome}.{aligner}.sort.bam",
		bkgrnd_ready = "utils/{ref_genome}.bffBuilt.flg",
	output:
		peaks_out = "fSeq/{sample}.vs_{ref_genome}.{aligner}.calledPeaks.bed"
	params:
		runmem_gb=8,
		runtime="1:00:00",
		cores=8,
	run:
		shell(""" mkdir -p fSeq/{wildcards.sample}/ """)
		shell(""" fseq -b utils/{wildcards.ref_genome}.background/ -v <( bedtools bamtobed -i {input.bam_in} ) -of npf -o fSeq/{wildcards.sample}/""" )
		shell(""" cat fSeq/{wildcards.sample}/*.npf | awk '{{print$1,$2,$3,"{wildcards.sample}",$5,$6,$7,$8,$9,$10}}' | tr " " "\t" > {output.peaks_out} """)
		shell(""" rm -rf fSeq/{wildcards.sample}/ """)


rule callAllThePeaks:
	input:
		expand("fSeq/{sample}.vs_dm6.bwa.calledPeaks.bed", sample = sampname_by_group['all'])
	output:
		flug = "utils/fseqDone.flg"
	params:
		runmem_gb=8,
		runtime="1:00:00",
		cores=8,
	run:
		shell(""" touch {output.flug} """)

rule basicPeakStats:
	input:
		peak2Peek = "fSeq/{sample}.vs_dm6.bwa.calledPeaks.bed",
	output:
		peakStatsOut = "meta/peakStats/{sample}.vs_dm6.bwa.calledPeaks.stats",
	params:
		runmem_gb=1,
		runtime="1:00",
		cores=1,
	run:
		shell(""" mkdir -p meta/peakStats/ """)
		shell(""" rm -f {output.peakStatsOut} """)
		shell(""" cat {input.peak2Peek} | cut -f 1 | uniq -c | tr -s " " | tr " " "\t" | awk '{{print"count\t"$2"\t"$1}}' >> {output.peakStatsOut}""" )
		shell(""" echo "count	total	$(wc -l {input.peak2Peek} | cut -f 1 -d " " )" >> {output.peakStatsOut}""" )
		shell(""" cat {input.peak2Peek}| awk '{{print$3-$2}}' | awk '{{sum+=$1; sumsq+=$1*$1}} END {{ print "width\tavg\t",sum/NR; print "width\tstd\t",sqrt(sumsq/NR - (sum/NR)**2)}}' >> {output.peakStatsOut}""" )
		shell(""" cat {input.peak2Peek}| awk '{{sum+=$7; sumsq+=$7*$7}} END {{ print "signal\tavg\t",sum/NR; print "signal\tstd\t",sqrt(sumsq/NR - (sum/NR)**2)}}' >> {output.peakStatsOut}""" )
		shell(""" touch {output}""")

rule summonBasicPeakStats:
	input:
		peakStatsIn = expand("meta/peakStats/{sample}.vs_dm6.bwa.calledPeaks.stats", sample=sampname_by_group['all']),
	output:
		fullPeakStatsOut = "meta/basicPeakStats.vs_dm6.bwa.summary",
	params:
		runmem_gb=1,
		runtime="1:00",
		cores=1,
	run:
		for samp in sampname_by_group['all']:
			shell(""" cat meta/peakStats/{samp}.vs_dm6.bwa.calledPeaks.stats | awk '{{print"{samp}\t"$0}}' >> {output.fullPeakStatsOut} """ )

rule singleDistance:
	input:
		peek1 = "fSeq/{samp1}.vs_dm6.bwa.calledPeaks.bed", 
		peek2 = "fSeq/{samp2}.vs_dm6.bwa.calledPeaks.bed", 
	output:
		distOut = "fSeq/closest/all/{samp1}.to.{samp2}.vs_dm6.bwa.closestPeaks.bed",
	params:
		runmem_gb=8,
		runtime="1:00:00",
		cores=8,
	run:
		shell(""" mkdir -p fSeq/closest/all/ """)
		shell("""  bedtools closest -io -d -D ref -t all -filenames -a <( bedtools sort -i {input.peek1} ) -b <( bedtools sort -i {input.peek2} ) | awk '{{print$0"\t{wildcards.samp2}"}}' > {output.distOut}""")

rule roundRobinDist:
	input:
		distances_in = expand("fSeq/closest/all/{samp1}.to.{samp2}.vs_dm6.bwa.closestPeaks.bed", samp1 = sampname_by_group['all'], samp2 = sampname_by_group['all']),
	output:
		allClose_flag = 'utils/closestCalcd.flg'
	params:
		runmem_gb=1,
		runtime="1:00",
		cores=1,
	run:
		shell(""" touch {output.allClose_flag} """)

rule consolidate_distance_by_group:
	input:
		allClose_flag = 'utils/closestCalcd.flg'
	output:
		grupFlg = 'utils/closest_{group}.flg'
	params:
		runmem_gb=8,
		runtime="1:00:00",
		cores=8,
	run:
		shell(""" mkdir -p fSeq/closest/{wildcards.group} """)
		for samp in sampname_by_group[wildcards.group]:
			for cramp in sampname_by_group[wildcards.group]:
				shell(""" cat fSeq/closest/all/{samp}.to.{cramp}.vs_dm6.bwa.closestPeaks.bed >> fSeq/closest/{wildcards.group}/{samp}.vs_dm6.bwa.closestPeaks.{wildcards.group}.bed """ )
		shell(""" touch {output.grupFlg} """)





rule collapseIntersects_byGroup:
	input:
		uncollapsed_in = lambda wildcards: expand("fSeq/{sample}.vs_dm6.bwa.calledPeaks.bed", sample = [s for s in sampname_by_group[wildcards.group] if (sample_by_name[s]["experimental"] == wildcards.spearmint and sample_by_name[s]["rep"] == "input" )]),
		uncollapsed_out = lambda wildcards: expand("fSeq/{sample}.vs_dm6.bwa.calledPeaks.bed", sample = [s for s in sampname_by_group[wildcards.group] if (sample_by_name[s]["experimental"] == wildcards.spearmint and sample_by_name[s]["rep"] != "input" )]),
	output:
		collapsed_in = "fSeq/collapse/{spearmint}.vs_dm6.bwa.group_{group}.input.signalsCollapsed.bed",
		collapsed_in_stats = "meta/collapseStats/{spearmint}.vs_dm6.bwa.group_{group}.input.signalsCollapsed.stat",
		collapsed_out = "fSeq/collapse/{spearmint}.vs_dm6.bwa.group_{group}.output.signalsCollapsed.bed",
		collapsed_out_stats = "meta/collapseStats/{spearmint}.vs_dm6.bwa.group_{group}.output.signalsCollapsed.stat",
	params:
		runmem_gb=8,
		runtime="1:00:00",
		cores=8,
	run:
		ref_fai = ref_genome_by_name["dm6"]["fai"]
		rep_count = len(input.uncollapsed_in)

		put = "input"
		#*peaksCollapsed.bed: chrom/start/stop/contributorList/signals/maxPkHt/peakWidth
		shell(""" mkdir -p fSeq/collapse/ meta/collapseStats/""")
		shell(""" cat {input.uncollapsed_in} | sort -k1,1 -k2,2n | bedtools merge -i - | bedtools map -b <(cat {input.uncollapsed_in}  | sort -k1,1 -k2,2n | awk '{{print$0"\t"$7*($3-$2)"\t"$3-$2}}' ) -a - -c 1,4,7,10,11,12 -o count,collapse,collapse,collapse,collapse,collapse > fSeq/collapse/{wildcards.spearmint}.vs_dm6.bwa.group_{wildcards.group}.input.peaksCollapsed.bed """ )
		shell(""" paste <(cat fSeq/collapse/{wildcards.spearmint}.vs_dm6.bwa.group_{wildcards.group}.input.peaksCollapsed.bed | cut -f 1-4 ) <(python3 scripts/overlapSignalCollapser.py -i fSeq/collapse/{wildcards.spearmint}.vs_dm6.bwa.group_{wildcards.group}.input.peaksCollapsed.bed  -r {rep_count}) | rev |cut -f 2- | rev > {output.collapsed_in} """)
		shell(""" echo "{put}\tcount\ttotal\t$(cat {output.collapsed_in} | wc -l)" > {output.collapsed_in_stats}""")
		shell(""" cat {output.collapsed_in} | awk '{{print$3-$2}}' | awk '{{sum+=$1; sumsq+=$1*$1}} END {{ print "width\tavg\t",sum/NR; print "width\tstd\t",sqrt(sumsq/NR - (sum/NR)**2)}}' | awk '{{print"{put}\t"$0}}' >> {output.collapsed_in_stats} """)
		shell(""" cat {output.collapsed_in} |cut -f 4 | awk '{{sum+=$1; sumsq+=$1*$1}} END {{ print "contribs\tavg\t",sum/NR; print "contribs\tstd\t",sqrt(sumsq/NR - (sum/NR)**2)}}' | awk '{{print"{put}\t"$0}}' >> {output.collapsed_in_stats} """)
		shell(""" cat {output.collapsed_in} |cut -f 8 | awk '{{sum+=$1; sumsq+=$1*$1}} END {{ print "signal\tavg\t",sum/NR; print "signal\tstd\t",sqrt(sumsq/NR - (sum/NR)**2)}}' | awk '{{print"{put}\t"$0}}' >> {output.collapsed_in_stats} """)

		shell(""" bedtools genomecov -bg -g {ref_fai} -i <( cat {input.uncollapsed_in} | sort -k1,1 -k2,2n ) > fSeq/collapse/{wildcards.spearmint}.vs_dm6.bwa.group_{wildcards.group}.input.subIntervals.bg """)
		shell(""" rm -rf fSeq/collapse/{wildcards.spearmint}.vs_dm6.bwa.group_{wildcards.group}.input.subIntervals.bed """)
		for num_included in range(1,rep_count+1):
			shell(""" bedtools intersect -wb -a  <( cat {input.uncollapsed_in} | sort -k1,1 -k2,2n ) -b <( cat fSeq/collapse/{wildcards.spearmint}.vs_dm6.bwa.group_{wildcards.group}.input.subIntervals.bg | awk '{{if($4=={num_included})print;}}')  >> fSeq/collapse/{wildcards.spearmint}.vs_dm6.bwa.group_{wildcards.group}.input.subIntervals.bed """)
		shell("""cat fSeq/collapse/{wildcards.spearmint}.vs_dm6.bwa.group_{wildcards.group}.input.subIntervals.bed | sort -k1,1 -k2,2n | bedtools merge -i - | bedtools map -b <(cat fSeq/collapse/{wildcards.spearmint}.vs_dm6.bwa.group_{wildcards.group}.input.subIntervals.bed| cut -f 1-7,9,10,14 | sort -k1,1 -k2,2n | awk '{{print$0"\t"$7*($3-$2)"\t"$3-$2}}' ) -a - -c 1,4,7,10,11,12 -o count,collapse,collapse,collapse,collapse,collapse > fSeq/collapse/{wildcards.spearmint}.vs_dm6.bwa.group_{wildcards.group}.input.subIntervals.collapsed.bed """)
		shell(""" touch {output.collapsed_in} """)

		put = "output"
		rep_count = len(input.uncollapsed_out)
		shell(""" cat {input.uncollapsed_out} | sort -k1,1 -k2,2n | bedtools merge -i - | bedtools map -b <(cat {input.uncollapsed_out}  | sort -k1,1 -k2,2n | awk '{{print$0"\t"$7*($3-$2)"\t"$3-$2}}' ) -a - -c 1,4,7,10,11,12 -o count,collapse,collapse,collapse,collapse,collapse > fSeq/collapse/{wildcards.spearmint}.vs_dm6.bwa.group_{wildcards.group}.output.peaksCollapsed.bed """ )
		shell(""" paste <(cat fSeq/collapse/{wildcards.spearmint}.vs_dm6.bwa.group_{wildcards.group}.output.peaksCollapsed.bed | cut -f 1-4 ) <(python3 scripts/overlapSignalCollapser.py -i fSeq/collapse/{wildcards.spearmint}.vs_dm6.bwa.group_{wildcards.group}.output.peaksCollapsed.bed  -r {rep_count}) | rev |cut -f 2- | rev > {output.collapsed_out} """)
		shell(""" echo "{put}\tcount\ttotal\t$(cat {output.collapsed_out} | wc -l)" > {output.collapsed_out_stats}""")
		shell(""" cat {output.collapsed_out} | awk '{{print$3-$2}}' | awk '{{sum+=$1; sumsq+=$1*$1}} END {{ print "width\tavg\t",sum/NR; print "width\tstd\t",sqrt(sumsq/NR - (sum/NR)**2)}}' | awk '{{print"{put}\t"$0}}' >> {output.collapsed_out_stats} """)
		shell(""" cat {output.collapsed_out} |cut -f 4 | awk '{{sum+=$1; sumsq+=$1*$1}} END {{ print "contribs\tavg\t",sum/NR; print "contribs\tstd\t",sqrt(sumsq/NR - (sum/NR)**2)}}' | awk '{{print"{put}\t"$0}}' >> {output.collapsed_out_stats} """)
		shell(""" cat {output.collapsed_out} |cut -f 8 | awk '{{sum+=$1; sumsq+=$1*$1}} END {{ print "signal\tavg\t",sum/NR; print "signal\tstd\t",sqrt(sumsq/NR - (sum/NR)**2)}}' | awk '{{print"{put}\t"$0}}' >> {output.collapsed_out_stats} """)

		shell(""" bedtools genomecov -bg -g {ref_fai} -i <( cat {input.uncollapsed_out} | sort -k1,1 -k2,2n ) > fSeq/collapse/{wildcards.spearmint}.vs_dm6.bwa.group_{wildcards.group}.output.subIntervals.bg """)
		shell(""" rm -rf fSeq/collapse/{wildcards.spearmint}.vs_dm6.bwa.group_{wildcards.group}.output.subIntervals.bed """)
		for num_included in range(1,rep_count+1):
			shell(""" bedtools intersect -wb -a  <( cat {input.uncollapsed_out} | sort -k1,1 -k2,2n ) -b <( cat fSeq/collapse/{wildcards.spearmint}.vs_dm6.bwa.group_{wildcards.group}.output.subIntervals.bg | awk '{{if($4=={num_included})print;}}')  >> fSeq/collapse/{wildcards.spearmint}.vs_dm6.bwa.group_{wildcards.group}.output.subIntervals.bed """)
		shell("""cat fSeq/collapse/{wildcards.spearmint}.vs_dm6.bwa.group_{wildcards.group}.output.subIntervals.bed | sort -k1,1 -k2,2n | bedtools merge -i - | bedtools map -b <(cat fSeq/collapse/{wildcards.spearmint}.vs_dm6.bwa.group_{wildcards.group}.output.subIntervals.bed| cut -f 1-7,9,10,14 | sort -k1,1 -k2,2n | awk '{{print$0"\t"$7*($3-$2)"\t"$3-$2}}' ) -a - -c 1,4,7,10,11,12 -o count,collapse,collapse,collapse,collapse,collapse > fSeq/collapse/{wildcards.spearmint}.vs_dm6.bwa.group_{wildcards.group}.output.subIntervals.collapsed.bed """)
		shell(""" touch {output.collapsed_out} """)
		shell(""" rm fSeq/collapse/{wildcards.spearmint}.vs_dm6.bwa.group_{wildcards.group}.input.subIntervals.bg fSeq/collapse/{wildcards.spearmint}.vs_dm6.bwa.group_{wildcards.group}.output.subIntervals.bg """)


rule collapse_all_intersects:
	input:
		clpsd_in = expand("meta/collapseStats/{spearmint}.vs_dm6.bwa.group_{group}.input.signalsCollapsed.stat", spearmint = list(set([ s["experimental"] for s in config['data_sets'] ])), group=["A","B","C"]),
		clpsd_out = expand("meta/collapseStats/{spearmint}.vs_dm6.bwa.group_{group}.output.signalsCollapsed.stat", spearmint = list(set([ s["experimental"] for s in config['data_sets'] ])), group=["A","B","C"]),
	output:
		clpsd_stat = "meta/collapsedPeakStats.vs_dm6.bwa.summary",
	params:
		runmem_gb=1,
		runtime="1:00",
		cores=1,
	run:
		shell(""" rm -rf {output.clpsd_stat} """)
		for spearmint in list(set([ s["experimental"] for s in config['data_sets'] ])) : 
			for grup in ["A","B","C"]:
				shell(""" cat meta/collapseStats/{spearmint}.vs_dm6.bwa.group_{grup}.input.signalsCollapsed.stat meta/collapseStats/{spearmint}.vs_dm6.bwa.group_{grup}.output.signalsCollapsed.stat | awk '{{print"{spearmint}\t{grup}\t"$0}}' >> {output.clpsd_stat} """)

rule merge_nearby_collapsed:
	input:
		unmerged_in = "fSeq/collapse/{spearmint}.vs_dm6.bwa.group_{group}.input.signalsCollapsed.bed",
		unmerged_out = "fSeq/collapse/{spearmint}.vs_dm6.bwa.group_{group}.output.signalsCollapsed.bed",
	output:
		merged_in = "fSeq/merge/{spearmint}.vs_dm6.bwa.group_{group}.input.signalsMerged.bed",
		merged_in_stats = "meta/mergeStats/{spearmint}.vs_dm6.bwa.group_{group}.input.signalsMerged.stat",
		merged_out = "fSeq/merge/{spearmint}.vs_dm6.bwa.group_{group}.output.signalsMerged.bed",
		merged_out_stats = "meta/mergeStats/{spearmint}.vs_dm6.bwa.group_{group}.output.signalsMerged.stat",
	params:
		runmem_gb=8,
		runtime="1:00:00",
		cores=8,
		merge_dist = 100,
		col2merge = 9, # column containing the input signal. default is 9, the Pessimist
		col2report = 7, # column to report in the *.stats
	run:
		shell(""" mkdir -p fSeq/merge/ meta/mergeStats """)

		put = "input"
		shell(""" cat {input.unmerged_in}  | awk '{{print$0"\t"$3-$2}}' | bedtools merge -d {params.merge_dist} -c 1,4,5,6,7,8,9 -o count,collapse,collapse,collapse,collapse,collapse,collapse -i - > fSeq/merge/{wildcards.spearmint}.vs_dm6.bwa.group_{wildcards.group}.input.peaksMerged.bed """)
		# output: chromStartStop, contribs to merged, contribs to collapsed, collapsed avg, collapsed rescale, collapsed weighted, collapsed pess, collapsed int lengths
		shell(""" paste <(cat fSeq/merge/{wildcards.spearmint}.vs_dm6.bwa.group_{wildcards.group}.input.peaksMerged.bed | cut -f 1-4 ) <(python3 scripts/collapsedSignalMerger.py -i fSeq/merge/{wildcards.spearmint}.vs_dm6.bwa.group_{wildcards.group}.input.peaksMerged.bed  -c {params.col2merge} ) | rev |cut -f 2- | rev > {output.merged_in} """)

		shell(""" echo "{put}\tcount\ttotal\t$(cat {output.merged_in} | wc -l) " > {output.merged_in_stats}""")
		shell(""" cat {output.merged_in} | awk '{{print$3-$2}}' | awk '{{sum+=$1; sumsq+=$1*$1}} END {{ print "width\tavg\t",sum/NR; print "width\tstd\t",sqrt(sumsq/NR - (sum/NR)**2)}}' | awk '{{print"{put}\t"$0}}' >> {output.merged_in_stats} """)
		shell(""" cat {output.merged_in} |cut -f 4 | awk '{{sum+=$1; sumsq+=$1*$1}} END {{ print "contribs\tavg\t",sum/NR; print "contribs\tstd\t",sqrt(sumsq/NR - (sum/NR)**2)}}' | awk '{{print"{put}\t"$0}}' >> {output.merged_in_stats} """)
		shell(""" cat {output.merged_in} |cut -f {params.col2report} | awk '{{sum+=$1; sumsq+=$1*$1}} END {{ print "signal\tavg\t",sum/NR; print "signal\tstd\t",sqrt(sumsq/NR - (sum/NR)**2)}}' | awk '{{print"{put}\t"$0}}' >> {output.merged_in_stats} """)

		put = "output"
		shell(""" cat {input.unmerged_out}  | awk '{{print$0"\t"$3-$2}}' | bedtools merge -d {params.merge_dist} -c 1,4,5,6,7,8,9 -o count,collapse,collapse,collapse,collapse,collapse,collapse -i - > fSeq/merge/{wildcards.spearmint}.vs_dm6.bwa.group_{wildcards.group}.output.peaksMerged.bed """)
		shell(""" paste <(cat fSeq/merge/{wildcards.spearmint}.vs_dm6.bwa.group_{wildcards.group}.output.peaksMerged.bed | cut -f 1-4 ) <(python3 scripts/collapsedSignalMerger.py -i fSeq/merge/{wildcards.spearmint}.vs_dm6.bwa.group_{wildcards.group}.output.peaksMerged.bed  -c {params.col2merge} ) | rev |cut -f 2- | rev > {output.merged_out} """)

		shell(""" echo "{put}\tcount\ttotal\t$(cat {output.merged_out} | wc -l) " > {output.merged_out_stats}""")
		shell(""" cat {output.merged_out} | awk '{{print$3-$2}}' | awk '{{sum+=$1; sumsq+=$1*$1}} END {{ print "width\tavg\t",sum/NR; print "width\tstd\t",sqrt(sumsq/NR - (sum/NR)**2)}}' | awk '{{print"{put}\t"$0}}' >> {output.merged_out_stats} """)
		shell(""" cat {output.merged_out} |cut -f 4 | awk '{{sum+=$1; sumsq+=$1*$1}} END {{ print "contribs\tavg\t",sum/NR; print "contribs\tstd\t",sqrt(sumsq/NR - (sum/NR)**2)}}' | awk '{{print"{put}\t"$0}}' >> {output.merged_out_stats} """)
		shell(""" cat {output.merged_out} |cut -f {params.col2report}  | awk '{{sum+=$1; sumsq+=$1*$1}} END {{ print "signal\tavg\t",sum/NR; print "signal\tstd\t",sqrt(sumsq/NR - (sum/NR)**2)}}' | awk '{{print"{put}\t"$0}}' >> {output.merged_out_stats} """)


rule merge_all_collapsed:
	input:
		mrgd_in = expand("meta/mergeStats/{spearmint}.vs_dm6.bwa.group_{group}.input.signalsMerged.stat", spearmint = list(set([ s["experimental"] for s in config['data_sets'] ])), group=["A","B","C"]),
		mrgd_out = expand("meta/mergeStats/{spearmint}.vs_dm6.bwa.group_{group}.output.signalsMerged.stat", spearmint = list(set([ s["experimental"] for s in config['data_sets'] ])), group=["A","B","C"]),
	output:
		mrgd_stats = "meta/mergedPeakStats.vs_dm6.bwa.summary",
	params:
		runmem_gb=1,
		runtime="1:00",
		cores=1,
	run:
		shell(""" rm -rf {output.mrgd_stats} """)
		for spearmint in list(set([ s["experimental"] for s in config['data_sets'] ])) : 
			for grup in ["A","B","C"]:
				shell(""" cat meta/mergeStats/{spearmint}.vs_dm6.bwa.group_{grup}.input.signalsMerged.stat meta/mergeStats/{spearmint}.vs_dm6.bwa.group_{grup}.output.signalsMerged.stat | awk '{{print"{spearmint}\t{grup}\t"$0}}' >> {output.mrgd_stats} """)






rule write_report:
	input:
		reference_genome_summary = ["meta/reference_genomes.summary"],
		reference_annotation_summary = ["meta/reference_annotations.summary"],
		gene_lists_summary = ["meta/geneLists.stats"],
		sequenced_reads_summary=["meta/sequenced_reads.dat"],
		aligned_reads_summary = ["meta/alignments.vs_dm6.bwa.summary"],
		called_peak_stats = ["meta/basicPeakStats.vs_dm6.bwa.summary"],
		collapsed_peak_stats = ["meta/collapsedPeakStats.vs_dm6.bwa.summary"],
		merged_peak_stats = ["meta/mergedPeakStats.vs_dm6.bwa.summary"],
	output:
		pdf_out="VolkanLab_BehaviorGenetics.pdf"
	params:
		runmem_gb=8,
		runtime="1:00:00",
		cores=2,
	message:
		"writing up the results.... "
	run:
		pandoc_path="/nas/longleaf/apps/rstudio/1.0.136/bin/pandoc"
		pwd = subprocess.check_output("pwd",shell=True).decode().rstrip()+"/"
		shell(""" R -e "setwd('{pwd}');Sys.setenv(RSTUDIO_PANDOC='{pandoc_path}')" -e  "peaDubDee='{pwd}'; rmarkdown::render('scripts/faire_results.Rmd',output_file='{pwd}{output.pdf_out}')"  """)







# rule intersectCollapsed_byGroup:
# 	input:
# 		inputPeaks= "fSeq/{spearmint}.vs_dm6.bwa.inputs.group_{grup}.collapsed.bed",
# 		outputPeaks= "fSeq/{spearmint}.vs_dm6.bwa.outputs.group_{grup}.collapsed.bed",
# 	output:
# 		sharedOverlap = "fSeq/{spearmint}.vs_dm6.bwa.shared.group_{grup}.collapsed.overlap.bed",
# 		sharedOut = "fSeq/{spearmint}.vs_dm6.bwa.shared.group_{grup}.collapsed.bed",
# 		inputExclusive = "fSeq/{spearmint}.vs_dm6.bwa.input_exclusive.group_{grup}.collapsed.bed",
# 		outputExclusive = "fSeq/{spearmint}.vs_dm6.bwa.output_exclusive.group_{grup}.collapsed.bed", 
# 		inputExclusiveReport = "meta/peakStats/{spearmint}.vs_dm6.bwa.input_exclusive.group_{grup}.collapsed.stats",
# 		outputExclusiveReport = "meta/peakStats/{spearmint}.vs_dm6.bwa.output_exclusive.group_{grup}.collapsed.stats",
# 	params:
# 		runmem_gb=8,
# 		runtime="1:00:00",
# 		cores=8,
# 	run:
# 		shell(""" bedtools intersect -wa -wb -a {input.inputPeaks} -b {input.outputPeaks} | bedtools overlap -i stdin -cols 2,3,6,7 > {output.sharedOverlap} """)
# 		shell(""" cat <(cut -f 1-3 {output.sharedOverlap} ) <(cut -f 5-7 {output.sharedOverlap} ) | sort -k1,1 -k2,2n | bedtools merge -i - > {output.sharedOut} """)
# 		shell(""" bedtools intersect -wa -v -a {input.inputPeaks} -b {input.outputPeaks} > {output.inputExclusive} """)
# 		shell(""" bedtools intersect -wa -v -b {input.inputPeaks} -a {input.outputPeaks} > {output.outputExclusive} """)
# 		shell(""" echo "count	total	$(wc -l {output.inputExclusive} | cut -f 1 -d " " )" >> {output.inputExclusiveReport};
# 				cat {output.inputExclusive}| awk '{{print$3-$2}}' | awk '{{sum+=$1; sumsq+=$1*$1}} END {{ print "width\tavg\t",sum/NR; print "width\tstd\t",sqrt(sumsq/NR - (sum/NR)**2)}}' >> {output.inputExclusiveReport};
# 				cat {output.inputExclusive} | awk '{{sum+=$4; sumsq+=$4*$4}} END {{ print "num_merged\tavg\t",sum/NR; print "num_merged\tstd\t",sqrt(sumsq/NR - (sum/NR)**2)}}' >> {output.inputExclusiveReport};
# 			""" )
# 		shell(""" echo "count	total	$(wc -l {output.outputExclusive} | cut -f 1 -d " " )" >> {output.outputExclusiveReport};
# 				cat {output.outputExclusive}| awk '{{print$3-$2}}' | awk '{{sum+=$1; sumsq+=$1*$1}} END {{ print "width\tavg\t",sum/NR; print "width\tstd\t",sqrt(sumsq/NR - (sum/NR)**2)}}' >> {output.outputExclusiveReport};
# 				cat {output.outputExclusive} | awk '{{sum+=$4; sumsq+=$4*$4}} END {{ print "num_merged\tavg\t",sum/NR; print "num_merged\tstd\t",sqrt(sumsq/NR - (sum/NR)**2)}}' >> {output.outputExclusiveReport};
# 			""" )
# ### add shared stats report?

# rule intersectAllCollapsed:
# 	input:
# 		intersects = lambda wildcards: expand("fSeq/{spearmint}.vs_dm6.bwa.{pile}.group_{grup}.collapsed.bed",spearmint = list(set([ s["experimental"] for s in config['data_sets'] ])), pile=["shared","input_exclusive","output_exclusive"], grup = ["A","B","C"] ),
# 		stats = lambda wildcards: expand("meta/peakStats/{spearmint}.vs_dm6.bwa.{pile}.group_{grup}.collapsed.stats",spearmint = list(set([ s["experimental"] for s in config['data_sets'] ])), pile=["input_exclusive","output_exclusive"], grup = ["A","B","C"] ),
# 	output:
# 		collapsedIntersectStats = "meta/collapsedPeaksExclusive.vs_dm6.bwa.stats",
# 	params:
# 		runmem_gb=1,
# 		runtime="1:00",
# 		cores=1,
# 	run:
# 		shell(""" rm -rf {output.collapsedIntersectStats} """)
# 		for spearmint in list(set([ s["experimental"] for s in config['data_sets'] ])):
# 			for grup in ["A","B","C"]:
# 				shell(""" cat meta/peakStats/{spearmint}.vs_dm6.bwa.input_exclusive.group_{grup}.collapsed.stats | awk '{{print "input\t{spearmint}\t{grup}\t"$0}}' >> {output.collapsedIntersectStats} """)
# 				shell(""" cat meta/peakStats/{spearmint}.vs_dm6.bwa.output_exclusive.group_{grup}.collapsed.stats | awk '{{print "output\t{spearmint}\t{grup}\t"$0}}' >> {output.collapsedIntersectStats} """)


# # rule interSample_collapsedComparator:
# # 	input:
# # 		collapsed1 = "fSeq/{prefix1}.collapsed.bed",
# # 		collapsed2 = "fSeq/{prefix2}.collapsed.bed",
# # 	output:

# # 	params:
# # 		runmem_gb=1,
# # 		runtime="1:00",
# # 		cores=1,
# # 	run:



# # #	this currently assumes that there is a single input .......
# # #	add merging??
# # rule intersectExperimental_byGroup:
# # 	input:
# # 		inputPeaks = lambda wildcards: expand("fSeq/{sample}.vs_dm6.bwa.calledPeaks.bed", sample = [s for s in sampname_by_group[wildcards.group] if (sample_by_name[s]["experimental"] == wildcards.spearmint and sample_by_name[s]["rep"] == "input" )]),
# # 		outputPeaks = lambda wildcards: expand("fSeq/{sample}.vs_dm6.bwa.calledPeaks.bed", sample = [s for s in sampname_by_group[wildcards.group] if (sample_by_name[s]["experimental"] == wildcards.spearmint and sample_by_name[s]["rep"] != "input" )]),
# # 	output:
# # 		sharedOut = "fSeq/{spearmint}.vs_dm6.bwa.group_{group}.sharedPeaks.bed",
# # 		sharedCount = "fSeq/{spearmint}.vs_dm6.bwa.group_{group}.sharedPeaks.counts.bed",
# # 		inputExclusive = "fSeq/{spearmint}.vs_dm6.bwa.group_{group}.inputExclusivePeaks.bed",
# # 		outputExclusive = "fSeq/{spearmint}.vs_dm6.bwa.group_{group}.outputExclusivePeaks.bed",
# # 	params:
# # 		runmem_gb=8,
# # 		runtime="1:00:00",
# # 		cores=8,
# # 	run:
# # 		shell("""bedtools intersect -wa -wb -filenames -a {input.inputPeaks} -b {input.outputPeaks} > {output.sharedOut}""")
# # 		shell("""bedtools intersect -c -wa  -a {input.inputPeaks} -b {input.outputPeaks} > {output.sharedCount}""")
# # 		shell("""  bedtools intersect -v -wa -a {input.inputPeaks} -b {input.outputPeaks} > {output.inputExclusive} """)
# # 		shell("""rm -f {output.outputExclusive}""")
# # 		for nom in input.outputPeaks:
# # 			label = nom.split("/")[1].split(".")[0]
# # 			shell(""" bedtools intersect -v -wa -a {nom} -b {input.inputPeaks} |  awk '{{print$0"\t{label}"}}' >> {output.outputExclusive} """)

# # rule intersectEveryone:
# # 	input:
# # 		outs = expand("fSeq/{spearmint}.vs_dm6.bwa.group_{group}.sharedPeaks.bed", spearmint = list(set([ s["experimental"] for s in config['data_sets'] ])), group=["A","B","C"]),
# # 		counts = expand("fSeq/{spearmint}.vs_dm6.bwa.group_{group}.sharedPeaks.counts.bed", spearmint = list(set([ s["experimental"] for s in config['data_sets'] ])), group=["A","B","C"]),
# # 		inEx = expand("fSeq/{spearmint}.vs_dm6.bwa.group_{group}.inputExclusivePeaks.bed", spearmint = list(set([ s["experimental"] for s in config['data_sets'] ])), group=["A","B","C"]),
# # 		outEx = expand("fSeq/{spearmint}.vs_dm6.bwa.group_{group}.outputExclusivePeaks.bed", spearmint = list(set([ s["experimental"] for s in config['data_sets'] ])), group=["A","B","C"]),
# # 	output:
# # 		ntrsct_flag = "utils/everyoneIntersected.flg",
# # 	params:
# # 		runmem_gb=1,
# # 		runtime="1:00",
# # 		cores=1,
# # 	run:
# #		shell(""" touch {output.ntrsct_flag} """)


# rule singleDistance:
# 	input:
# 		peek1 = "fSeq/{samp1}.vs_dm6.bwa.calledPeaks.bed", 
# 		peek2 = "fSeq/{samp2}.vs_dm6.bwa.calledPeaks.bed", 
# 	output:
# 		distOut = "fSeq/closest/all/{samp1}.to.{samp2}.vs_dm6.bwa.closestPeaks.bed",
# 	params:
# 		runmem_gb=8,
# 		runtime="1:00:00",
# 		cores=8,
# 	run:
# 		shell(""" mkdir -p fSeq/closest/all/ """)
# 		shell("""  bedtools closest -io -d -D ref -t all -filenames -a <( bedtools sort -i {input.peek1} ) -b <( bedtools sort -i {input.peek2} ) | awk '{{print$0"\t{wildcards.samp2}"}}' > {output.distOut}""")

# rule roundRobinDist:
# 	input:
# 		distances_in = expand("fSeq/closest/all/{samp1}.to.{samp2}.vs_dm6.bwa.closestPeaks.bed", samp1 = sampname_by_group['all'], samp2 = sampname_by_group['all']),
# 	output:
# 		allClose_flag = 'utils/closestCalcd.flg'
# 	params:
# 		runmem_gb=1,
# 		runtime="1:00",
# 		cores=1,
# 	run:
# 		shell(""" touch {output.allClose_flag} """)

# rule consolidate_distance_by_group:
# 	input:
# 		allClose_flag = 'utils/closestCalcd.flg'
# 	output:
# 		grupFlg = 'utils/closest_{group}.flg'
# 	params:
# 		runmem_gb=8,
# 		runtime="1:00:00",
# 		cores=8,
# 	run:
# 		shell(""" mkdir -p fSeq/closest/{wildcards.group} """)
# 		for samp in sampname_by_group[wildcards.group]:
# 			for cramp in sampname_by_group[wildcards.group]:
# 				shell(""" cat fSeq/closest/all/{samp}.to.{cramp}.vs_dm6.bwa.closestPeaks.bed >> fSeq/closest/{wildcards.group}/{samp}.vs_dm6.bwa.closestPeaks.{wildcards.group}.bed """ )
# 		shell(""" touch {output.grupFlg} """)



# rule annotation_intersect_bed:
# 	input:
# 		bed3 = "fSeq/{bed_prefix}.bed",
# 	output:
# 		full_intersect = "features/intersects/{bed_prefix}.intersect.{annot_name}.full.bed",
# 		#simpleIntersect = "features/intersects/{bed_prefix}.intersect.{annot_name}.compact.bed",
# 		#strictlyIntronic = "features/intersects/{bed_prefix}.intersect.{annot_name}.intronic.bed",
# 		intersect_report = "meta/intersects/{bed_prefix}.intersect.{annot_name}.full.stats",
# 	params:
# 		runmem_gb=8,
# 		runtime="1:00:00",
# 		cores=8,	
# 	run:
# 		shell(""" mkdir -p features/intersects/ meta/intersects/""")
# 		shell(""" rm -f {output.intersect_report} """)
# 		annot_bed = annotation_by_name[wildcards.annot_name]["bed_path"]
# 		shell(""" bedtools intersect -wa -wb -a <(cat {input.bed3} | cut -f 1-3) -b {annot_bed}  > {output.full_intersect} """)
# 		#shell(""" cat {output.fullIntersect} | grep -w gene | cut -f 1-4,7,8,10,12| tr -d '"' | tr -d ";" | tr " " "\t" | cut -f 1-7,9,11 > {output.simpleIntersect} """)
# 		#shell(""" bedtools intersect -v -a {output.simpleIntersect} -b <(cat {output.fullIntersect} | grep -w exon | cut -f 1-3) > {output.strictlyIntronic} """)
# 		shell(""" echo "count\tuniqA\t"$(cat {output.full_intersect} | cut -f 1-3 | uniq | wc -l ) >> {output.intersect_report} """)		
# 		shell(""" echo "count\tuniqB\t"$(cat {output.full_intersect} | cut -f 4-6 | uniq | wc -l ) >> {output.intersect_report} """)
# 		shell(""" echo "max\tAperB\t"$(cat {output.full_intersect} | cut -f 1-3 | uniq -c | tr -s " " | cut -f 2 -d " " | awk 'BEGIN{{a=   0}}{{if ($1>0+a) a=$1}} END{{print a}}' ) >> {output.intersect_report} """)
# 		shell(""" cat {output.full_intersect} | cut -f 1-3 | uniq -c | tr -s " " | cut -f 2 -d " " | awk '{{sum+=$1}} END {{ if (NR + 0 != 0) print "avg\tAperB\t",sum/NR; else print "avg\tAperB\t0";}}'  >> {output.intersect_report}; """)
# 		shell(""" echo "max\tBperA\t"$(cat {output.full_intersect} | cut -f 4-6 | uniq -c | tr -s " " | cut -f 2 -d " " | awk 'BEGIN{{a=   0}}{{if ($1>0+a) a=$1}} END{{print a}}' ) >> {output.intersect_report} """)
# 		shell(""" cat {output.full_intersect} | cut -f 4-6 | uniq -c | tr -s " " | cut -f 2 -d " " | awk '{{sum+=$1}} END {{ if (NR + 0 != 0) print "avg\tBperA\t",sum/NR; else print "avg\tAperB\t0";}}' >> {output.intersect_report}; """)

# #cat /proj/cdjones_lab/Genomics_Data_Commons/annotations/drosophila_melanogaster/dmel-all-r6.13.gtf | awk '{if($3=="gene")print;}' | awk '{print "chr"$1,$4,$5,$7,$10,$12;}' | tr " " "\t" | tr -d ";" | tr -d '"' | awk '{print$1,$2,$3,$5,"0",$4,$6}' | tr " " "\t" > /proj/cdjones_lab/Genomics_Data_Commons/annotations/drosophila_melanogaster/dmel-all-r6.13.bed 


# rule SummarizeAnnotationIntersects:
# 	input:
# 		intersects = lambda wildcards: expand("features/intersects/{spearmint}.vs_dm6.bwa.{pile}.group_{grup}.collapsed.intersect.{annot_name}.full.bed",spearmint = list(set([ s["experimental"] for s in config['data_sets'] ])), pile=["shared","input_exclusive","output_exclusive"], grup = ["A","B","C"], annot_name=wildcards.annot_name ),
# 		stats = lambda wildcards: expand("meta/intersects/{spearmint}.vs_dm6.bwa.{pile}.group_{grup}.collapsed.intersect.{annot_name}.full.stats",spearmint = list(set([ s["experimental"] for s in config['data_sets'] ])), pile=["input_exclusive","output_exclusive","shared"], grup = ["A","B","C"], annot_name=wildcards.annot_name ),
# 	output:
# 		annotIntersectStats = "meta/vs_dm6.bwa.collapsed.intersect.{annot_name}.full.stats",
# 	params:
# 		runmem_gb=1,
# 		runtime="1:00",
# 		cores=1,
# 	run:
# 		shell(""" rm -rf {output.annotIntersectStats} """)
# 		for spearmint in list(set([ s["experimental"] for s in config['data_sets'] ])):
# 			for grup in ["A","B","C"]:
# 				for pile in ["shared","input_exclusive","output_exclusive"]:
# 					shell(""" cat meta/intersects/{spearmint}.vs_dm6.bwa.{pile}.group_{grup}.collapsed.intersect.{wildcards.annot_name}.full.stats | awk '{{print "{pile}\t{spearmint}\t{grup}\t"$0}}' >> {output.annotIntersectStats} """)

# rule countInteresectingPeaks_byAnnotation:
# 	input:
# 		intersectBed = "features/intersects/{bed_prefix}.intersect.{annot_name}.full.bed",
# 	output:
# 		interesectGeneCounts = "features/intersects/{bed_prefix}.intersect.{annot_name}.byGenePeaks.count",
# 	params:
# 		runmem_gb=8,
# 		runtime="1:00:00",
# 		cores=1,
# 		col = 7,
# 	run:
# 		shell(""" cat {input.intersectBed} | cut -f 7 | sort | uniq -c | tr -s " " | tr " " "\t" | awk '{{print $2"\t"$1}}' > {output.interesectGeneCounts} """)


# rule summonPeakCounts_byAnnotation:
# 	input:
# 		counts_in = lambda wildcards: expand("features/intersects/{spearmint}.vs_dm6.bwa.{pile}.group_{grup}.collapsed.intersect.{annot_name}.byGenePeaks.count",spearmint = list(set([ s["experimental"] for s in config['data_sets'] ])), pile=["shared","input_exclusive","output_exclusive"], grup = ["A","B","C"], annot_name=wildcards.annot_name ),
# 	output:
# 		peakCountsByAnn_tbl = "features/intersects/peakCountsByAnn.{annot_name}.tbl",
# 	params:
# 		runmem_gb=1,
# 		runtime="1:00",
# 		cores=1,
# 	run:
# 		#shell(""" touch {output.peakCountsByAnn_flag} """)
# 		shell(""" rm -rf {output.peakCountsByAnn_tbl} """)
# 		for spearmint in list(set([ s["experimental"] for s in config['data_sets'] ])):
# 			for grup in ["A","B","C"]:
# 				for pile in ["shared","input_exclusive","output_exclusive"]:
# 					shell(""" cat features/intersects/{spearmint}.vs_dm6.bwa.{pile}.group_{grup}.collapsed.intersect.{wildcards.annot_name}.byGenePeaks.count | sed -e 's/^/{spearmint}\t{pile}\t{grup}\t/g' >> {output.peakCountsByAnn_tbl} """)		









# rule bed2ann_dist:
# 	input:
# 		bed3 = "fSeq/{bed_prefix}.bed",
# 	output:
# 		upstream_closest = "features/closest/{bed_prefix}.closestUpstreamFeat.{annot_name}.full.bed",
# 		downstream_closest = "features/closest/{bed_prefix}.closestDownstreamFeat.{annot_name}.full.bed",
# 		overall_closest = "features/closest/{bed_prefix}.closestFeatOverall.{annot_name}.full.bed",
# 	params:
# 		runmem_gb=8,
# 		runtime="1:00:00",
# 		cores=8,
# 	run:
# 		#numfields=awk '{print NF}' closestUpstream.neg | sort | uniq # numfields is currently 11,but use this to generalize?
# 		annot_bed = annotation_by_name[wildcards.annot_name]["bed_path"]
# 		shell(""" mkdir -p features/closest/""" )#" meta/intersects/""")
# 		shell(""" bedtools closest -id -d -D ref -a <( cut -f 1-3 {input.bed3} | bedtools sort -i - ) -b <( bedtools sort -i {annot_bed} | awk '{{if($6=="+")print}}') > features/closest/{wildcards.bed_prefix}.closestUpstreamFeat.{wildcards.annot_name}.full.pos """)
# 		shell(""" bedtools closest -iu -d -D ref -a <( cut -f 1-3 {input.bed3} | bedtools sort -i - ) -b <( bedtools sort -i {annot_bed} | awk '{{if($6=="-")print}}') > features/closest/{wildcards.bed_prefix}.closestUpstreamFeat.{wildcards.annot_name}.full.neg """)		
# 		shell(""" bedtools intersect -wa -wb -a features/closest/{wildcards.bed_prefix}.closestUpstreamFeat.{wildcards.annot_name}.full.neg -b features/closest/{wildcards.bed_prefix}.closestUpstreamFeat.{wildcards.annot_name}.full.pos |  awk '{{if($1==$12 && $2==$13 && $3==$14)print;}}' > features/closest/{wildcards.bed_prefix}.closestUpstreamFeat.{wildcards.annot_name}.full.joint """)
# 		shell(""" cat features/closest/{wildcards.bed_prefix}.closestUpstreamFeat.{wildcards.annot_name}.full.joint | awk '{{if($15==".")print;}}' | awk '{{if($4!=".")print;}}' | cut -f 1-11 > features/closest/{wildcards.bed_prefix}.closestUpstreamFeat.{wildcards.annot_name}.full.unsrt """)
# 		shell(""" cat features/closest/{wildcards.bed_prefix}.closestUpstreamFeat.{wildcards.annot_name}.full.joint | awk '{{if($4==".")print;}}' | awk '{{if($15!=".")print;}}' | cut -f 12- >> features/closest/{wildcards.bed_prefix}.closestUpstreamFeat.{wildcards.annot_name}.full.unsrt """)
# 		shell(""" cat features/closest/{wildcards.bed_prefix}.closestUpstreamFeat.{wildcards.annot_name}.full.joint | awk '{{if( $15!="." && $4!="." )print;}}' | awk 'function abs(x){{return ((x < 0.0) ? -x : x)}} {{if(abs($11)<abs($22))print;}}' | cut -f 1-11 >> features/closest/{wildcards.bed_prefix}.closestUpstreamFeat.{wildcards.annot_name}.full.unsrt """)
# 		shell(""" cat features/closest/{wildcards.bed_prefix}.closestUpstreamFeat.{wildcards.annot_name}.full.joint | awk '{{if( $15!="." && $4!="." )print;}}' | awk 'function abs(x){{return ((x < 0.0) ? -x : x)}} {{if(abs($11)>abs($22))print;}}' | cut -f 12- >> features/closest/{wildcards.bed_prefix}.closestUpstreamFeat.{wildcards.annot_name}.full.unsrt """)
# 		shell(""" bedtools sort -i features/closest/{wildcards.bed_prefix}.closestUpstreamFeat.{wildcards.annot_name}.full.unsrt > {output.upstream_closest} """)

# 		shell(""" bedtools closest -id -d -D ref -a <( cut -f 1-3 {input.bed3} | bedtools sort -i -  ) -b <( bedtools sort -i {annot_bed} | awk '{{if($6=="-")print}}') > features/closest/{wildcards.bed_prefix}.closestDownstreamFeat.{wildcards.annot_name}.full.pos """)
# 		shell(""" bedtools closest -iu -d -D ref -a <( cut -f 1-3 {input.bed3} | bedtools sort -i - ) -b <( bedtools sort -i {annot_bed} | awk '{{if($6=="+")print}}') > features/closest/{wildcards.bed_prefix}.closestDownstreamFeat.{wildcards.annot_name}.full.neg """)
# 		shell(""" bedtools intersect -wa -wb -a features/closest/{wildcards.bed_prefix}.closestDownstreamFeat.{wildcards.annot_name}.full.neg -b features/closest/{wildcards.bed_prefix}.closestDownstreamFeat.{wildcards.annot_name}.full.pos |  awk '{{if($1==$12 && $2==$13 && $3==$14)print;}}' > features/closest/{wildcards.bed_prefix}.closestDownstreamFeat.{wildcards.annot_name}.full.joint """)
# 		shell(""" cat features/closest/{wildcards.bed_prefix}.closestDownstreamFeat.{wildcards.annot_name}.full.joint | awk '{{if($15==".")print;}}' | awk '{{if($4!=".")print;}}' | cut -f 1-11 > features/closest/{wildcards.bed_prefix}.closestDownstreamFeat.{wildcards.annot_name}.full.unsrt """)
# 		shell(""" cat features/closest/{wildcards.bed_prefix}.closestDownstreamFeat.{wildcards.annot_name}.full.joint | awk '{{if($4==".")print;}}' | awk '{{if($15!=".")print;}}' | cut -f 12- >> features/closest/{wildcards.bed_prefix}.closestDownstreamFeat.{wildcards.annot_name}.full.unsrt """)
# 		shell(""" cat features/closest/{wildcards.bed_prefix}.closestDownstreamFeat.{wildcards.annot_name}.full.joint | awk '{{if( $15!="." && $4!="." )print;}}' | awk 'function abs(x){{return ((x < 0.0) ? -x : x)}} {{if(abs($11)<abs($22))print;}}' | cut -f 1-11 >> features/closest/{wildcards.bed_prefix}.closestDownstreamFeat.{wildcards.annot_name}.full.unsrt """)
# 		shell(""" cat features/closest/{wildcards.bed_prefix}.closestDownstreamFeat.{wildcards.annot_name}.full.joint | awk '{{if( $15!="." && $4!="." )print;}}' | awk 'function abs(x){{return ((x < 0.0) ? -x : x)}} {{if(abs($11)>abs($22))print;}}' | cut -f 12- >> features/closest/{wildcards.bed_prefix}.closestDownstreamFeat.{wildcards.annot_name}.full.unsrt """)
# 		shell(""" bedtools sort -i features/closest/{wildcards.bed_prefix}.closestDownstreamFeat.{wildcards.annot_name}.full.unsrt > {output.downstream_closest} """)

# 		shell(""" bedtools closest -d -D ref -a <( cut -f 1-3 {input.bed3} ) -b <( bedtools sort -i {annot_bed} ) > {output.overall_closest} """)

# 		shell(""" rm -f features/closest/{wildcards.bed_prefix}.closestDownstreamFeat.{wildcards.annot_name}.full.unsrt features/closest/{wildcards.bed_prefix}.closestDownstreamFeat.{wildcards.annot_name}.full.joint features/closest/{wildcards.bed_prefix}.closestDownstreamFeat.{wildcards.annot_name}.full.pos features/closest/{wildcards.bed_prefix}.closestDownstreamFeat.{wildcards.annot_name}.full.neg """)
# 		shell(""" rm -f features/closest/{wildcards.bed_prefix}.closestUpstreamFeat.{wildcards.annot_name}.full.unsrt features/closest/{wildcards.bed_prefix}.closestUpstreamFeat.{wildcards.annot_name}.full.joint features/closest/{wildcards.bed_prefix}.closestUpstreamFeat.{wildcards.annot_name}.full.pos features/closest/{wildcards.bed_prefix}.closestUpstreamFeat.{wildcards.annot_name}.full.neg """)



# rule summon_bed2annDist_byContrast:
# 	input:
# 		upstream_closest = expand("features/closest/{spearmint}.vs_dm6.bwa.{pile}.group_{grup}.collapsed.closestUpstreamFeat.{annot_name}.full.bed", annot_name = ["dm6_genes"], spearmint = list(set([ s["experimental"] for s in config['data_sets'] ])), pile=["shared","input_exclusive","output_exclusive"], grup = ["A","B","C"] ),
# 		downstream_closest = expand("features/closest/{spearmint}.vs_dm6.bwa.{pile}.group_{grup}.collapsed.closestDownstreamFeat.{annot_name}.full.bed", annot_name = ["dm6_genes"], spearmint = list(set([ s["experimental"] for s in config['data_sets'] ])), pile=["shared","input_exclusive","output_exclusive"], grup = ["A","B","C"] ),
# 		overall_closest = expand("features/closest/{spearmint}.vs_dm6.bwa.{pile}.group_{grup}.collapsed.closestFeatOverall.{annot_name}.full.bed", annot_name = ["dm6_genes"], spearmint = list(set([ s["experimental"] for s in config['data_sets'] ])), pile=["shared","input_exclusive","output_exclusive"], grup = ["A","B","C"] ),
# 	output:
# 		annot_dist_flag = "utils/dist2ann.flg"
# 	params:
# 		runmem_gb=1,
# 		runtime="1:00:00",
# 		cores=1,
# 	run:
# 		shell( """ touch {output.annot_dist_flag} """)


# rule extractGeneList_byAnnotationDist:
# 	input:
# 		annotDist = "features/closest/{spearmint}.vs_dm6.bwa.{pile}.group_{grup}.collapsed.{featType}.{annot_name}.full.bed",
# 		#closestFeatOverall, annot_name = ["dm6_genes"], spearmint = list(set([ s["experimental"] for s in config['data_sets'] ])), pile=["shared","input_exclusive","output_exclusive"], grup = ["A","B","C"] ),
# 	output:
# 		distList = "features/genelists/closest/{spearmint}.vs_dm6.bwa.{pile}.group_{grup}.{featType}.{annot_name}.{distcat}.list",
# 	params:
# 		runmem_gb=8,
# 		runtime="1:00:00",
# 		cores=8,
# 		name_column=7,
# 		dist_column=11,
# 	run:
# 		break1 = distance_breaks[wildcards.distcat][0]
# 		break2 = distance_breaks[wildcards.distcat][1]
# 		shell(""" cat {input.annotDist} | awk '{{if(${params.name_column}!=".")print;}}'| awk 'function abs(x){{return ((x < 0.0) ? -x : x)}} {{if(abs(${params.dist_column})>{break1})print;}}' |  awk 'function abs(x){{return ((x < 0.0) ? -x : x)}} {{if(abs(${params.dist_column})<={break2})print;}}' | cut -f {params.name_column} | sort | uniq > {output.distList} """)


# rule annotationDist_geneList_similarity:
# 	input: #features/genelists/closest/{spearmint}.vs_dm6.bwa.{pile}.group_{grup}.{featType}.{annot_name}.{distcat}.list
# 		closest_lists = lambda wildcards : expand("features/genelists/closest/{spearmint}.vs_dm6.bwa.{pile}.group_{grup}.{closest}.{annot_name}.{distcat}.list",annot_name = ["dm6_genes"], spearmint = list(set([ s["experimental"] for s in config['data_sets'] ])), pile=["shared","input_exclusive","output_exclusive"], grup = ["A","B","C"], distcat = distance_breaks.keys(), closest=wildcards.closest ),
# 	output:
# 		similarity = "features/genelists/similarities/{closest}.geneLists.similarities"
# 	params:
# 		runmem_gb=8,
# 		runtime="1:00:00",
# 		cores=8,
# 		brakes = [0,100,1000,20000, 2000000000000],
# 	run:
# 		shell(""" rm -f {output.similarity} """)
# 		for list1 in input.closest_lists:
# 			for list2 in input.closest_lists:

# 					shell(""" echo {list1} > {output.similarity}.tmp """)
# 					shell(""" echo {list2} >> {output.similarity}.tmp """)

# 					shell(""" cat {list1} | wc -l >> {output.similarity}.tmp """)
# 					shell(""" cat {list2} | wc -l >> {output.similarity}.tmp """)
# 					shell(""" grep -c -wFf {list1} {list2} >> {output.similarity}.tmp || true """)
# 					shell(""" cat {output.similarity}.tmp  | tr "\n" "\t" >> {output.similarity} """)
# 					shell(""" echo "\n" >> {output.similarity} """)
# 		shell(""" cat {output.similarity} | sed -e 's/features\/genelists\/closest\///g' | sed -e 's/.vs_dm6.bwa./\t/g' | sed -e 's/.group_/\t/g' | sed -e 's/.list//g' | tr -s "\n" > {output.similarity}.tmp """)
# 		shell(""" cat {output.similarity}.tmp |tr -s "\n" | tr "." "\t" > {output.similarity} """)
# 		shell(""" rm -f {output.similarity}.tmp """)


# rule annotationDist_geneList_similarities_allStreams:
# 	input:
# 		upstream_sim = "features/genelists/similarities/closestUpstreamFeat.geneLists.similarities",
# 		downstream_sim = "features/genelists/similarities/closestDownstreamFeat.geneLists.similarities",
# 		overall_sim = "features/genelists/similarities/closestFeatOverall.geneLists.similarities",
# 	output:
# 		upstream_byDist = "features/genelists/closestUpstreamFeat.geneLists.similarities.byDist",
# 		downstream_byDist = "features/genelists/closestDownstreamFeat.geneLists.similarities.byDist",
# 		overall_byDist = "features/genelists/closestFeatOverall.geneLists.similarities.byDist",
# 		stream_flag = "utils/allStreamClosest.GeneListSim.flg"
# 	params:
# 		runmem_gb=1,
# 		runtime="1:00",
# 		cores=1,
# 	run:
# 		shell(""" cat {input.upstream_sim} | awk '{{if($6==$12)print;}}' > {output.upstream_byDist} """)
# 		shell(""" cat {input.downstream_sim} | awk '{{if($6==$12)print;}}' > {output.downstream_byDist} """)
# 		shell(""" cat {input.overall_sim} | awk '{{if($6==$12)print;}}' > {output.overall_byDist} """)
# 		shell(""" touch {output.stream_flag} """)


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


# rule summon_all_proximity_GOs:
# 	input:
# 		closest_lists = lambda wildcards : expand("features/ontologies/{spearmint}.vs_dm6.bwa.{pile}.group_{grup}.{closest}.{annot_name}.{distcat}.go", annot_name = ["dm6_genes"], spearmint = list(set([ s["experimental"] for s in config['data_sets'] ])), pile=["shared","input_exclusive","output_exclusive"], grup = ["A","B","C"], distcat = distance_breaks.keys(), closest=wildcards.closest ),
# 	output:
# 		megaGo = "features/ontologies/{closest}.dm6_genes.go"
# 	params:
# 		runmem_gb=1,
# 		runtime="1:00",
# 		cores=1,
# 	run:
# 		annot_name = "dm6_genes"
# 		shell(""" rm -f {output.megaGo} """)
# 		for spearmint in list(set([ s["experimental"] for s in config['data_sets'] ])):
# 			for grup in ["A","B","C"]:
# 				for pile in ["shared","input_exclusive","output_exclusive"]:
# 					for distcat in distance_breaks.keys():
# 						shell(""" cat features/ontologies/{spearmint}.vs_dm6.bwa.{pile}.group_{grup}.{wildcards.closest}.{annot_name}.{distcat}.go | awk '{{print"{spearmint}\t{pile}\t{grup}\t{wildcards.closest}\t{distcat}\t"$0}}' >> {output.megaGo} """)




# rule annotationDist_geneList_ontologies_allStreams:
# 	input:
# 		upstream_go = "features/ontologies/closestUpstreamFeat.dm6_genes.go",
# 		downstream_go = "features/ontologies/closestDownstreamFeat.dm6_genes.go",
# 		overall_go = "features/ontologies/closestFeatOverall.dm6_genes.go",
# 	output:
# 		ont_flag = "utils/allStreamClosest.GeneListOntolog.flg"
# 	params:
# 		runmem_gb=1,
# 		runtime="1:00",
# 		cores=1,
# 	run:
# 		shell(""" touch {output.ont_flag} """)


# rule ann2bed_dist:
# 	input:
# 		bed3 = "fSeq/{bed_prefix}.bed",		
# 	output:
# 		upstream_closest = "features/closest/{annot_name}.closestUpstreamPeak.{bed_prefix}.full.bed",
# 		downstream_closest = "features/closest/{annot_name}.closestDownstreamPeak.{bed_prefix}.full.bed",
# 		overall_closest = "features/closest/{annot_name}.closestPeakOverall.{bed_prefix}.full.bed",
# 	params:
# 		runmem_gb=8,
# 		runtime="1:00:00",
# 		cores=8,
# 	run:
# 		annot_bed = annotation_by_name[wildcards.annot_name]["bed_path"]
# 		shell(""" bedtools closest -d -io -id -D a -a <( bedtools sort -i {annot_bed} ) -b <( cut -f 1-3 {input.bed3} | bedtools sort -i - ) > {output.upstream_closest} """)
# 		shell(""" bedtools closest -d -io -iu -D a -a <( bedtools sort -i {annot_bed} ) -b <( cut -f 1-3 {input.bed3} | bedtools sort -i - ) > {output.downstream_closest} """)
# 		shell(""" bedtools closest -d -io -D a -a <( bedtools sort -i {annot_bed} ) -b <( cut -f 1-3 {input.bed3} | bedtools sort -i -  ) > {output.overall_closest} """)


# rule summon_ann2bedDist_byContrast:
# 	input:
# 		upstream_closest = expand("features/closest/{annot_name}.closestUpstreamPeak.{spearmint}.vs_dm6.bwa.{pile}.group_{grup}.collapsed.full.bed", annot_name = ["dm6_genes","ionotropic","nervSysDev","histoneMod","mating","synapseSig"], spearmint = list(set([ s["experimental"] for s in config['data_sets'] ])), pile=["shared","input_exclusive","output_exclusive"], grup = ["A","B","C"] ),
# 		downstream_closest = expand("features/closest/{annot_name}.closestDownstreamPeak.{spearmint}.vs_dm6.bwa.{pile}.group_{grup}.collapsed.full.bed", annot_name = ["dm6_genes","ionotropic","nervSysDev","histoneMod","mating","synapseSig"], spearmint = list(set([ s["experimental"] for s in config['data_sets'] ])), pile=["shared","input_exclusive","output_exclusive"], grup = ["A","B","C"] ),
# 		overall_closest = expand("features/closest/{annot_name}.closestPeakOverall.{spearmint}.vs_dm6.bwa.{pile}.group_{grup}.collapsed.full.bed", annot_name = ["dm6_genes","ionotropic","nervSysDev","histoneMod","mating","synapseSig"], spearmint = list(set([ s["experimental"] for s in config['data_sets'] ])), pile=["shared","input_exclusive","output_exclusive"], grup = ["A","B","C"] ),
# 	output:
# 		annot_dist_flag = "utils/ann2dist.flg"
# 	params:
# 		runmem_gb=1,
# 		runtime="1:00",
# 		cores=1,
# 	run:
# 		shell( """ touch {output.annot_dist_flag} """)

# rule fseq_scromble:
# 	input:
# 		bed_in = "fSeq/{bed_prfx1}.vs_{ref_gen}.{bed_prfx2}.bed"
# 	output:
# 		bed_out = "fSeq/{bed_prfx1}.vs_{ref_gen,[a-zA-Z0-9]*}.{bed_prfx2}.scramble_{n}.bed"
# 	params:
# 		runmem_gb=8,
# 		runtime="1:00:00",
# 		cores=8,
# 		seeds = [13,69,420,777,12345],
# 		scram = "-chrom -noOverlapping -maxTries 1000000",
# 	run:
# 		reffai = ref_genome_by_name[wildcards.ref_gen]["fai"]
# 		try:

# 			seed = "-seed %s" % tuple([ params.seeds[int(wildcards.n)] ])
# 		except NameError:
# 			seed = ""
# 		except IndexError:
# 			print("warning: unspecified seed index; using pRNG. ")
# 			seed = ""
# 		shell(""" bedtools shuffle {seed} {params.scram} -i {input.bed_in} -g {reffai} > {output.bed_out} """)






# # rule fseq_scromble_summoner:
# # 	input:
# # 		scrombles_in = expand("fSeq/{spearmint}.vs_dm6.bwa.{pile}.group_{grup}.collapsed.scramble_{scram}.bed", spearmint = list(set([ s["experimental"] for s in config['data_sets'] ])), pile=["shared","input_exclusive","output_exclusive"], grup = ["A","B","C"] , scram = [0,1,2,3,4])
# # 	output:
# # 		scrombled_flg = "utils/scromble.flg"
# # 	params:
# # 		runmem_gb=1,
# # 		runtime="1:00",
# # 		cores=1,
# # 	run:
# # 		shell( """ touch {output.scrombled_flg} """)

# rule ann2bedDist_byContrast_scrombleCompresser:
# 	input:
# 		upstream_scromble = lambda wildcards: expand("features/closest/{annot_name}.closestUpstreamPeak.{spearmint}.vs_dm6.bwa.{pile}.group_{grup}.collapsed.scramble_{scram}.full.bed",scram = [0,1,2,3,4], annot_name = wildcards.annot_name, spearmint = wildcards.spearmint, pile = wildcards.pile, grup = wildcards.grup),
# 		downstream_scromble = lambda wildcards: expand("features/closest/{annot_name}.closestDownstreamPeak.{spearmint}.vs_dm6.bwa.{pile}.group_{grup}.collapsed.scramble_{scram}.full.bed",scram = [0,1,2,3,4], annot_name = wildcards.annot_name, spearmint = wildcards.spearmint, pile = wildcards.pile, grup = wildcards.grup),
# 		overall_scromble = lambda wildcards: expand("features/closest/{annot_name}.closestPeakOverall.{spearmint}.vs_dm6.bwa.{pile}.group_{grup}.collapsed.scramble_{scram}.full.bed",scram = [0,1,2,3,4], annot_name = wildcards.annot_name, spearmint = wildcards.spearmint, pile = wildcards.pile, grup = wildcards.grup),
# 	output:
# 		upstream_compress = "features/scrambled/closest/{annot_name}.closestUpstreamPeak.{spearmint}.vs_dm6.bwa.{pile}.group_{grup}.collapsed.full.scrmbl.bed",
# 		downstream_compress = "features/scrambled/closest/{annot_name}.closestDownstreamPeak.{spearmint}.vs_dm6.bwa.{pile}.group_{grup}.collapsed.full.scrmbl.bed",
# 		overall_compress = "features/scrambled/closest/{annot_name}.closestPeakOverall.{spearmint}.vs_dm6.bwa.{pile}.group_{grup}.collapsed.full.scrmbl.bed",
# 	params:
# 		runmem_gb=8,
# 		runtime="1:00:00",
# 		cores=8,
# 	run:
# 		shell("""mkdir -p features/scrambled/closest/""")

# 		up_base = input.upstream_scromble[0]
# 		shell(""" cat {up_base} | rev | cut -f 5- | rev > {output.upstream_compress} """)
# 		for upstr in input.upstream_scromble:
# 			shell(""" paste {output.upstream_compress} <(cat {upstr} | rev | cut -f 1 | rev ) > {output.upstream_compress}.tmp """)
# 			shell(""" mv {output.upstream_compress}.tmp {output.upstream_compress} """)

# 		down_base = input.downstream_scromble[0]
# 		shell(""" cat {down_base} | rev | cut -f 5- | rev > {output.downstream_compress} """)
# 		for dwnstr in input.downstream_scromble:
# 			shell(""" paste {output.downstream_compress} <(cat {dwnstr} | rev | cut -f 1 | rev ) > {output.downstream_compress}.tmp """)
# 			shell(""" mv {output.downstream_compress}.tmp {output.downstream_compress} """)

# 		all_base = input.overall_scromble[0]
# 		shell(""" cat {all_base} | rev | cut -f 5- | rev > {output.overall_compress} """)
# 		for vrl in input.overall_scromble:
# 			shell(""" paste {output.overall_compress} <(cat {vrl} | rev | cut -f 1 | rev ) > {output.overall_compress}.tmp """)
# 			shell(""" mv {output.overall_compress}.tmp {output.overall_compress} """)

#  		##shell(""" rm {input} """)

# rule annot_scromble_intersector:
# 	input:
# 		moveme = "features/intersects/{bed_prefix}.intersect.{annot_name}.full.bed",
# 	output:
# 		moved = "features/scrambled/intersects/{bed_prefix}.intersect.{annot_name}.full.bed",
# 	params:
# 		runmem_gb=1,
# 		runtime="1:00",
# 		cores=1,
# 	run:
# 		shell(""" mkdir -p features/scrambled/intersects/ """)
# 		shell(""" mv {input.moveme} {output.moved} """)

# rule summon_annot_scromble_intersects:
# 	input:
# 		expand("features/scrambled/intersects/{spearmint}.vs_dm6.bwa.{pile}.group_{grup}.collapsed.scramble_{scrm}.intersect.{annot_name}.full.bed", annot_name=["ionotropic", "nervSysDev","histoneMod","mating","synapseSig"], spearmint = list(set([ s["experimental"] for s in config['data_sets'] ])), grup = [ g for g in sampname_by_group.keys() if g != "all" ], pile=["shared","input_exclusive","output_exclusive"], scrm = range(0,len(rules.fseq_scromble.params.seeds)) )
# 	output:
# 		annotScrmblIntersectStats = "meta/scrambledPeaks_intersectWithAnnots.counts"
# 	params:
# 		runmem_gb=1,
# 		runtime="1:00",
# 		cores=1,
# 	run:
# 		shell(""" rm -rf {output.annotScrmblIntersectStats} """)		
# 		for spearmint in list(set([ s["experimental"] for s in config['data_sets'] ])):
# 			for grup in [ g for g in sampname_by_group.keys() if g != "all" ]:
# 				for pile in ["shared","input_exclusive","output_exclusive"]:
# 					for annot_name in ["ionotropic"]:
# 						for scrm in range(0,len(rules.fseq_scromble.params.seeds)):
# 							shell(""" cat features/scrambled/intersects/{spearmint}.vs_dm6.bwa.{pile}.group_{grup}.collapsed.scramble_{scrm}.intersect.{annot_name}.full.bed |cut -f 4- |sort | uniq | wc -l | awk '{{print "{annot_name}\t{pile}\t{spearmint}\t{grup}\t{scrm}\t"$0}}' >> {output.annotScrmblIntersectStats} """)







# # rule ann2bedDist_byContrast_scrombleCompresser:
# # 	input:
# # 		upstream_scromble = lambda wildcards: expand("features/closest/{annot_name}.closestUpstreamPeak.{spearmint}.vs_dm6.bwa.{pile}.group_{grup}.collapsed.scramble_{scram}.full.bed",scram = [0,1,2,3,4], annot_name = wildcards.annot_name, spearmint = wildcards.spearmint, pile = wildcards.pile, grup = wildcards.grup),
# # 		downstream_scromble = lambda wildcards: expand("features/closest/{annot_name}.closestDownstreamPeak.{spearmint}.vs_dm6.bwa.{pile}.group_{grup}.scramble_{scram}.collapsed.full.bed",  scram = [0,1,2,3,4], annot_name = wildcards.annot_name, spearmint = wildcards.spearmint, pile = wildcards.pile, grup = wildcards.grup),
# # 		overall_scromble = lambda wildcards: expand("features/closest/{annot_name}.closestPeakOverall.{spearmint}.vs_dm6.bwa.{pile}.group_{grup}.scramble_{scram}.collapsed.full.bed" , scram = [0,1,2,3,4], annot_name = wildcards.annot_name, spearmint = wildcards.spearmint, pile = wildcards.pile, grup = wildcards.grup),
# # 	output:
# # 		upstream_compress = "features/closest/{annot_name}.closestUpstreamPeak.{spearmint}.vs_dm6.bwa.{pile}.group_{grup}.collapsed.scrambled.bed",
# # 		downstream_compress = "features/closest/{annot_name}.closestDownstreamPeak.{spearmint}.vs_dm6.bwa.{pile}.group_{grup}.collapsed.scrambled.bed",
# # 		overall_compress = "features/closest/{annot_name}.closestPeakOverall.{spearmint}.vs_dm6.bwa.{pile}.group_{grup}.collapsed.scrambled.bed",
# # 	params:
# # 		runmem_gb=8,
# # 		runtime="1:00:00",
# # 		cores=8,
# # 	run:
# # 		shell("""mkdir -p features/closest/scrambled/""")

# # 		up_base = input.upstream_scromble[0]
# # 		shell(""" cat {up_base} | rev | cut -f 5- | rev > {output.upstream_compress} """)
# # 		for upstr in input.upstream_scromble:
# # 			shell(""" paste {output.upstream_compress} <(cat {upstr} | rev | cut -f 1 | rev ) > {output.upstream_compress}.tmp """)
# # 			shell(""" mv {output.upstream_compress}.tmp {output.upstream_compress} """)

# # 		down_base = input.downstream_scromble[0]
# # 		shell(""" cat {down_base} | rev | cut -f 5- | rev > {output.downstream_compress} """)
# # 		for dwnstr in input.downstream_scromble:
# # 			shell(""" paste {output.downstream_compress} <(cat {dwnstr} | rev | cut -f 1 | rev ) > {output.downstream_compress}.tmp """)
# # 			shell(""" mv {output.downstream_compress}.tmp {output.downstream_compress} """)

# # 		all_base = input.overall_scromble[0]
# # 		shell(""" cat {all_base} | rev | cut -f 5- | rev > {output.overall_compress} """)
# # 		for vrl in input.overall_scromble:
# # 			shell(""" paste {output.overall_compress} <(cat {vrl} | rev | cut -f 1 | rev ) > {output.overall_compress}.tmp """)
# # 			shell(""" mv {output.overall_compress}.tmp {output.overall_compress} """)

# # 		shell(""" rm {input} """)





# rule compressAll_ann2bedDist_byContrast_scrombles:
# 	input:
# 		upstream_compress = expand("features/scrambled/closest/{annot_name}.closestUpstreamPeak.{spearmint}.vs_dm6.bwa.{pile}.group_{grup}.collapsed.full.scrmbl.bed", annot_name = ["dm6_genes","ionotropic","nervSysDev","histoneMod","mating","synapseSig"], spearmint = list(set([ s["experimental"] for s in config['data_sets'] ])), pile=["shared","input_exclusive","output_exclusive"], grup = ["A","B","C"] ),
# 		downstream_compress = expand("features/scrambled/closest/{annot_name}.closestDownstreamPeak.{spearmint}.vs_dm6.bwa.{pile}.group_{grup}.collapsed.full.scrmbl.bed", annot_name = ["dm6_genes","ionotropic","nervSysDev","histoneMod","mating","synapseSig"], spearmint = list(set([ s["experimental"] for s in config['data_sets'] ])), pile=["shared","input_exclusive","output_exclusive"], grup = ["A","B","C"] ),
# 		overall_compress = expand("features/scrambled/closest/{annot_name}.closestPeakOverall.{spearmint}.vs_dm6.bwa.{pile}.group_{grup}.collapsed.full.scrmbl.bed", annot_name = ["dm6_genes","ionotropic","nervSysDev","histoneMod","mating","synapseSig"], spearmint = list(set([ s["experimental"] for s in config['data_sets'] ])), pile=["shared","input_exclusive","output_exclusive"], grup = ["A","B","C"] ),
# 	output:
# 		scrombled_flag = "utils/ann2dist.scromble.flg"
# 	params:
# 		runmem_gb=1,
# 		runtime="1:00",
# 		cores=1,
# 	run:
# 		shell( """ touch {output.scrombled_flag} """)




# # 		upstream_compress = expand("features/closest/scrambled/{annot_name}.closestUpstreamPeak.{spearmint}.vs_dm6.bwa.{pile}.group_{grup}.collapsed.full.bed",annot_name = ["dm6_genes","ionotropic"], spearmint = list(set([ s["experimental"] for s in config['data_sets'] ])), pile=["shared","input_exclusive","output_exclusive"], grup = ["A","B","C"] ,),
# # 		downstream_compress =  expand("features/closest/scrambled/{annot_name}.closestDownstreamPeak.{spearmint}.vs_dm6.bwa.{pile}.group_{grup}.collapsed.full.bed",annot_name = ["dm6_genes","ionotropic"], spearmint = list(set([ s["experimental"] for s in config['data_sets'] ])), pile=["shared","input_exclusive","output_exclusive"], grup = ["A","B","C"] ,),
# # 		overall_compress =  expand("features/closest/scrambled/{annot_name}.closestPeakOverall.{spearmint}.vs_dm6.bwa.{pile}.group_{grup}.collapsed.full.bed",annot_name = ["dm6_genes","ionotropic"], spearmint = list(set([ s["experimental"] for s in config['data_sets'] ])), pile=["shared","input_exclusive","output_exclusive"], grup = ["A","B","C"] ,),
# # 	output:
# # 		scrombled_flag = "utils/ann2dist.scromble.flg"
# # 	params:
# # 		runmem_gb=1,
# # 		runtime="1:00",
# # 		cores=1,
# # 	run:
# # 		shell( """ touch {output.scrombled_flag} """)




# #  

# # rule summon_ann2bedDist_byContrast_scromble:
# # 	input:
# # 		upstream_closest = expand("features/closest/{annot_name}.closestUpstreamPeak.{spearmint}.vs_dm6.bwa.{pile}.group_{grup}.collapsed.scramble_{scram}.full.bed", annot_name = ["dm6_genes","ionotropic"], spearmint = list(set([ s["experimental"] for s in config['data_sets'] ])), pile=["shared","input_exclusive","output_exclusive"], grup = ["A","B","C"] , scram = [0,1,2,3,4]),
# # 		downstream_closest = expand("features/closest/{annot_name}.closestDownstreamPeak.{spearmint}.vs_dm6.bwa.{pile}.group_{grup}.scramble_{scram}.collapsed.full.bed", annot_name = ["dm6_genes","ionotropic"], spearmint = list(set([ s["experimental"] for s in config['data_sets'] ])), pile=["shared","input_exclusive","output_exclusive"], grup = ["A","B","C"] , scram = [0,1,2,3,4]),
# # 		overall_closest = expand("features/closest/{annot_name}.closestPeakOverall.{spearmint}.vs_dm6.bwa.{pile}.group_{grup}.scramble_{scram}.collapsed.full.bed", annot_name = ["dm6_genes","ionotropic"], spearmint = list(set([ s["experimental"] for s in config['data_sets'] ])), pile=["shared","input_exclusive","output_exclusive"], grup = ["A","B","C"] , scram = [0,1,2,3,4]),
# # 	output:
# # 		annot_dist_flag = "utils/ann2dist.scram.flg"
# # 	params:
# # 		runmem_gb=1,
# # 		runtime="1:00",
# # 		cores=1,
# # 	run:
# # 		shell( """ touch {output.annot_dist_flag} """)



