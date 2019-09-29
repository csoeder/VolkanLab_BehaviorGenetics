from argparse import ArgumentParser
import pysam as ps


parser = ArgumentParser()
parser.add_argument("-i", "--sam_in", help="multimaps in")
parser.add_argument("-o", "--sam_out", help="spread multimaps out")

# parser.add_argument("-i", "--idxstat_in", help="samtools idxstat report")
# parser.add_argument("-g", "--genomecov_in", help="bedtools genomecov report")
# parser.add_argument("-d", "--depthstats_in", help="samtools depth report")
# #parser.add_argument("stat_in", help="samtools stats report")
# parser.add_argument("-o", "--flat_out", help="flatfile summary")
# parser.add_argument("-t", "--tag", help="line-name for the flatfile", default=None)
args = parser.parse_args()

sammy = ps.AlignmentFile(args.sam_in, "rb")
sammySprud = ps.AlignmentFile(args.sam_out, "w", template=sammy)

for reed in sammy:
	sammySprud.write(reed)
	for alt in [csv.split(",") for csv in reed.get_tag("XA").split(";")[:-1]]:
		lou_reed = reed
		lou_reed.pos = abs(int(float(alt[1])))
		lou_reed.reference_name= alt[0]
		sammySprud.write(lou_reed)

sammy.close()
sammySprud.close()


