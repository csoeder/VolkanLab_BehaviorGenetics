from argparse import ArgumentParser
from numpy import mean, array


parser = ArgumentParser()
parser.add_argument("-i", "--bed_in", help="overlap bed in")
parser.add_argument("-c", "--signalColumn", help="which column in the input contains the signal to be merged (signal collapser reports many)", default=9)

#parser.add_argument("-o", "--sam_out", help="significance")

# parser.add_argument("-i", "--idxstat_in", help="samtools idxstat report")
# parser.add_argument("-g", "--genomecov_in", help="bedtools genomecov report")
# parser.add_argument("-d", "--depthstats_in", help="samtools depth report")
# #parser.add_argument("stat_in", help="samtools stats report")
# parser.add_argument("-o", "--flat_out", help="flatfile summary")
# parser.add_argument("-t", "--tag", help="line-name for the flatfile", default=None)
args = parser.parse_args()

col = int(args.signalColumn) - 1

bed_file = open(args.bed_in,'r')
bed_lines = bed_file.readlines()
bed_file.close()

#the signal is the mximum Kernel density in the interval, divided by the 
for bed_line in bed_lines:
	bed_line = bed_line.split("\n")[0].split("\t")
	#contributors = list(set([c for c in bed_line[4].split(',')]))
	s_list = [float(s) for s in bed_line[col].split(',')]
	contrib_num = int(bed_line[3])
	metacontrib_nums = [int(m) for m in bed_line[4].split(",")]
#	spare_list = [x for x in bed_line[6].split(',')]
#	k_list = [float(k) for k in bed_line[7].split(',')]
	component_len = [int(l) for l in bed_line[9].split(',')]
	int_len = int(bed_line[2])-int(bed_line[1])

	Savg = mean(s_list)#flat average of the signal strengths

	k_list = array(s_list)*array(component_len)

	Kmean = mean(k_list)
	Srescale = Kmean/int_len # average the max kernel densities, 

	Kweight = sum(array(k_list)*array(component_len))/sum(component_len)
	Sweight = Kweight/int_len # weighted average of K by contributing peak width, then rescale by final peak width


	# missing =int(args.rep_count) - len(contributors)
	# s_list.extend([0]*missing)
	# k_list.extend([0]*missing)
	# component_len.extend([0]*missing)

	# Kpess = sum(array(k_list)*array(component_len))/sum(component_len)
	#Spess = Sweight*(len(contributors)/float(args.rep_count)) # same as Sweight, but now considers replicates with no peak

	data_out = [Savg, Srescale, Sweight]
	print("%s\t"*len(data_out) % tuple(data_out))
