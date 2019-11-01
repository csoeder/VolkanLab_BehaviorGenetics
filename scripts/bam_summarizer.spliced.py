import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--flagstat_in", help="samtools flagstat report")
parser.add_argument("-i", "--idxstat_in", help="samtools idxstat report")
parser.add_argument("-g", "--split_coverage", help="bedtools genomecov report, with split reads")
parser.add_argument("-G", "--spanning_coverage", help="bedtools genomecov report, with spanning reads")
parser.add_argument("-d", "--depthstats_in", help="samtools depth report")
#parser.add_argument("stat_in", help="samtools stats report")
parser.add_argument("-o", "--flat_out", help="flatfile summary")
parser.add_argument("-t", "--tag", help="line-name for the flatfile", default=None)
args = parser.parse_args()


summary_dict={}

flagstat = open(args.flagstat_in, 'r')
flagstat_lines = flagstat.readlines()
flagstat.close()

idxstat = open(args.idxstat_in, 'r')
idxstat_lines = idxstat.readlines()[:-1]
idxstat.close()

span = open(args.spanning_coverage, 'r')
span_lines = span.readlines()
span.close()

split = open(args.split_coverage, 'r')
split_lines = split.readlines()
split.close()

dpth = open(args.depthstats_in, 'r')
dpth_lines = dpth.readlines()
dpth.close()



summary_dict['total_read_count'] = int(flagstat_lines[0].split(" ")[0])
summary_dict['total_mapped_count'] = int(flagstat_lines[4].split(" ")[0])
summary_dict['properly_paired_count'] = int(flagstat_lines[0].split(" ")[0])
#summary_dict['avg_depth'] = sum([float(p.split('\t')[2]) for p in idxstat_lines ])/sum([int(q.split('\t')[1]) for q in idxstat_lines ])
summary_dict['spanned_breadth'] = float(span_lines[-1].split()[-1])
summary_dict['split_breadth'] = float(split_lines[-1].split()[-1])
summary_dict['avg_depth'] = float(dpth_lines[0].split("\t")[1])
summary_dict['std_depth'] = float(dpth_lines[1].split("\t")[1])

phial_out = open(args.flat_out,'w')

keys = ['total_read_count','total_mapped_count', 'properly_paired_count','avg_depth', 'std_depth', 'spanned_breadth', 'split_breadth']

lines2write = [ [k, summary_dict[k]] for k in keys]
if args.tag:
	[ ell.insert(0, args.tag) for ell in lines2write ]

for preline in lines2write:
	field_count = len(preline)
	line = ("%s" + "\t%s"*(field_count-1) + "\n") % tuple(preline)
	phial_out.write(line)

phial_out.close()



