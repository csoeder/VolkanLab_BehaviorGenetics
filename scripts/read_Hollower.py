
import pysam
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--samfile_in", help="file in SAM format to ")
# parser.add_argument("-i", "--idxstat_in", help="samtools idxstat report")
# parser.add_argument("-g", "--split_coverage", help="bedtools genomecov report, with split reads")
# parser.add_argument("-G", "--spanning_coverage", help="bedtools genomecov report, with spanning reads")
# parser.add_argument("-d", "--depthstats_in", help="samtools depth report")
# parser.add_argument("-D", "--dupestats_in", help="samtools markdup stats report")
# parser.add_argument( "-m", "--mapping_multiplicity", help="histogram of mapping multiplicity scraped from the IH:: tags")
# parser.add_argument( "-c", "--mapped_count", help="count of mapped reads")
parser.add_argument("-o", "--samfile_out", help="flatfile summary")
# parser.add_argument("-t", "--tag", help="line-name for the flatfile", default=None)
args = parser.parse_args()

sam_in = args.samfile_in
sammy = pysam.AlignmentFile(sam_in,"rb")  


sam_out = args.samfile_out
samuel = pysam.AlignmentFile(sam_out, "w", template=sammy)


# cigartuples encoding:
# M = 0
# I = 1
# D = 2
# N = 3
# S = 4


for read in sammy.fetch():
	cig = read.cigar
	nu_cig = []
	suffix = []
	if cig[0][0] == 4:
		cig = cig[1:]
	if cig[-1][0] == 4:	#discard the softclips
		cig = cig[:-1]
	front_pad_tally = 0
	while cig[0][0] in [0,1,2]:
		front_pad_tally += cig[0][1]
		cig.pop(0)
	back_pad_tally = 0
	while cig[-1][0] in [0,1,2]:
		back_pad_tally += cig[-1][1]
		cig.pop(-1)

	nu_cig.extend([(3,front_pad_tally-1),(0,1)]) # new cigar starts with 5' end clipped down to the actual splice

	# if cig[0][0] == 0: #take the terminal match blocks
	# 	num_match = cig[0][-1]
	# 	prefix = [(3,num_match-1),(0,1)]
	# 	nu_cig.extend(prefix) #add a bunch of Ns to pad ; ignore 5' and 3' read ends. 
	# 	cig = cig[1:]
	# if cig[-1][0] == 0:
	# 	num_match = cig[-1][1]
	# 	suffix = 
	# 	cig = cig[:-1]

	rill_tally = 0
	mids = []
	for rillo in cig:
		if rillo[0] in [0,1,2]:
			rill_tally += rillo[1] #consolidate exons
		else:
			if rill_tally > 0:
				mids.extend([(0,1),(3,rill_tally-2),(0,1)]) #write hollowed exon
			mids.extend([rillo]) # write intron
			rill_tally=0 # reset tally
	nu_cig.extend(mids)
	nu_cig.extend([(0,1),(3,back_pad_tally-1),]) # end cigar with splice and Ns to 3' end

	read.cigar = nu_cig
	read.seq = "*"
	samuel.write(read)

samuel.close()

