from argparse import ArgumentParser
from numpy import mean


parser = ArgumentParser()
parser.add_argument("-i", "--bed_in", help="overlap bed in")
parser.add_argument("-r", "--rep_count", help="number of replicates total", default=0)

#parser.add_argument("-o", "--sam_out", help="significance")

# parser.add_argument("-i", "--idxstat_in", help="samtools idxstat report")
# parser.add_argument("-g", "--genomecov_in", help="bedtools genomecov report")
# parser.add_argument("-d", "--depthstats_in", help="samtools depth report")
# #parser.add_argument("stat_in", help="samtools stats report")
# parser.add_argument("-o", "--flat_out", help="flatfile summary")
# parser.add_argument("-t", "--tag", help="line-name for the flatfile", default=None)
args = parser.parse_args()


bed_file = open(args.bed_in,'r')
bed_lines = bed_file.readlines()
bed_file.close()


for bed_line in bed_lines:
	bed_line = bed_line.split("\n")[0].split("\t")
	contributors = list(set([c for c in bed_line[3].split(',')]))
	s_list = [float(s) for s in bed_line[4].split(',')]
	spare_list = [x for x in bed_line[5].split(',')]
	z_list = [float(z) for z in bed_line[6].split(',')]
	int_len = int(bed_line[2])-int(bed_line[1])

	Savg = mean(s_list)
	Zstouff = mean(z_list)
	Sstouff = Zstouff/int_len

	missing =int(args.rep_count) - len(contributors)
	z_list.extend([0]*missing)
	Zpess = mean(z_list)
	Spess = Zpess/int_len

	data_out = [Savg, Zstouff, Sstouff, Zpess, Spess]
	print("%s\t"*len(data_out) % tuple(data_out))
