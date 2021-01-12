import argparse
import gffutils


parser = argparse.ArgumentParser()
parser.add_argument("-i", "--gtf_in", help="gene model to cronch")
parser.add_argument("-o", "--output_prefix", help="file to write ORFs to")
parser.add_argument("-n", "--gene_name", help="name of gene to cronch ")
#parser.add_argument("-c", "--length_cutoff", help="ignore ORFs with fewer internal (ie, not start or stops) codons", type=int, default=75)
#parser.add_argument("-x", "--extend_across_gaps", help="extend the ORF even if it encounters a gap/N", action="store_true", default=False)
#parser.add_argument("-a", "--include_all_starts", help="report all ORFs when multiple start codons present (else choose longest)", action="store_true", default=False)

args = parser.parse_args()


file_out = args.output_prefix

#path2file="/proj/cdjones_lab/csoeder/VolkanLab_BehaviorGenetics/utils/annotations/fru.gtf"
path2file="/Users/csoeder/Research/VolkanLab_BehaviorGenetics/utils/annotations/fru.gtf"
path2file=args.gtf_in


fn = gffutils.example_filename(path2file)
print(open(fn).read()) 

#db = gffutils.create_db(fn, dbfn='/proj/cdjones_lab/csoeder/VolkanLab_BehaviorGenetics/test.db', force=True, keep_order=True, merge_strategy='merge', sort_attribute_values=True)
#db = gffutils.create_db(fn, dbfn='/proj/cdjones_lab/csoeder/VolkanLab_BehaviorGenetics/test.db', force=True, keep_order=True, merge_strategy='merge', sort_attribute_values=True, id_spec=['ID', 'Name'] ,disable_infer_genes=True)

path2gtf = "%s.db" % tuple([path2file])
#db = gffutils.create_db(fn, dbfn='/Users/csoeder/Research/VolkanLab_BehaviorGenetics/test.db', force=True, keep_order=True, merge_strategy='merge', sort_attribute_values=True, disable_infer_genes=True)
db = gffutils.create_db(fn, dbfn=path2gtf, force=True, keep_order=True, merge_strategy='merge', sort_attribute_values=True, disable_infer_genes=True)

gene_name="FBgn0004652"
gene_name=args.gene_name

gene = db[gene_name]


transcription_edges = []

exon_to_transcript = {}


for transcript in db.children(gene, featuretype='transcript'): 
	transcription_edges.append(transcript.start)
	transcription_edges.append(transcript.stop)
	for exon in db.children(transcript, featuretype='exon'): 
		if exon.start in exon_to_transcript.keys():
			exon_to_transcript[exon.start].append(transcript.id)
		else:
			exon_to_transcript[exon.start] = [transcript.id]
		if exon.stop in exon_to_transcript.keys():
			exon_to_transcript[exon.stop].append(transcript.id)
		else:
			exon_to_transcript[exon.stop] = [transcript.id]

# currently ignoring the *transcript* start/stop sites .... kick all of them out


startStop = []
for transcript in db.children(gene, featuretype='transcript'): 
	startStop.extend([transcript.start, transcript.stop])
startStop = list(set(startStop))

for term in startStop:
	exon_to_transcript.pop(term)



isoid_definitions = {}
isoid_num=0
for iso in exon_to_transcript.values():
	if iso in isoid_definitions.values():
		pass
	else:
		isoid_definitions["isoid_%s"%tuple([isoid_num])] = iso
		isoid_num +=1


isoids = {}
for key in isoid_definitions.keys():
	isoids[key]={"transcripts":isoid_definitions[key],"junctions":[]}
	for junk in exon_to_transcript.keys():
		if exon_to_transcript[junk] == isoid_definitions[key]:
			#print(junk, key, exon_to_transcript[junk])
			isoids[key]["junctions"].append(junk)




phial_name = "%s.isoids.gtf" % tuple([file_out])
phial_out = open(phial_name, "w")


for isoid_name in isoids.keys():
#isoid_name='isoid_0'
	isoid = isoids[isoid_name]
	j_list =isoid['junctions']

	j_list.sort()


	gene_gtf_list = ["edgeHog","gene",min(j_list), max(j_list), 0, gene.strand, ".", "gene_id '%s'; transcript_id '%s';\n" %tuple([isoid_name, isoid_name]) ]
	gene_gtf_str="\t%s"*len(gene_gtf_list) % tuple(gene_gtf_list)
	gene_gtf_str = "%s%s" % tuple([gene.chrom , gene_gtf_str])


	junk_strungs = []
	j_count=0
	for jay in j_list:
		junct_gtf_list = ["edgeHog","exon",jay, jay, 0, gene.strand, ".", "gene_id '%s'; transcript_id '%s'; junction_id j_%s;\n" %tuple([isoid_name, isoid_name, j_count]) ]
		junct_gtf_str="\t%s"*len(junct_gtf_list) % tuple(junct_gtf_list)
		junct_gtf_str = "%s%s" % tuple([gene.chrom , junct_gtf_str])
		junk_strungs.append(junct_gtf_str)
		j_count += 1


	phial_out.write(gene_gtf_str)
	for junk_line in junk_strungs:
		phial_out.write(junk_line)

phial_out.close()


phial_name = "isoids.list"
phial_name = "%s.isoids.list" % tuple([file_out])

phial_out = open(phial_name, "w")

for isoid_name in isoids.keys():
	for trans in isoids[isoid_name]["transcripts"]:
		phial_out.write("'%s'\t%s\n" %tuple([isoid_name, trans]))


phial_out.close()


###	IV look at read counts vs junctions, isoids









