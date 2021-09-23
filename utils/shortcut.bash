

usage()
{
cat << EOF
usage: bash utils/shortcut.bash [-f]
-f    | --full_data		Includes the sample 'FRULEXAFRU440_1', excluded from Deanhardt2021 but required for some legacy uses " 
-h    | --help		Brings up this menu
EOF
}

full=0

while [ "$1" != "" ]; do
    case $1 in
        -f | --full_data )
            shift
            full=1
        ;;
        -h | --help )    usage
            exit
        ;;
        * )              usage
            exit 1
    esac
    shift
done

echo "Will download the counts of reads mapped ('mapspliceMulti' strategy) to dm6 and overlapping the v6.13 gene annotation from flyBase  "
echo "and save to the location expected by the Snakemake workflow. Samples are those underlying Deanhardt et al 2021:"
echo "'Changes in splicing and neuromodulatory gene expression programs in sensory neurons with pheromone signaling and social experience', Deanhardt et al. 2021, bioRxiv. doi: https://doi.org/10.1101/2021.06.18.449021 url: https://www.biorxiv.org/content/10.1101/2021.06.18.449021v2"
echo ;
sleep 2;
if [[ $full -gt 0 ]]
	then echo "the problematic FruLexaFru440 replicate #1 will also be downloaded and merged.";
	echo ;
fi
sleep 5;


mkdir -p counts/tmp/

echo ;echo ;echo ;echo ;echo ;echo ;echo ;echo ;
echo "downloading counts for Deanhardt et al 2021 ......"
echo ;echo ;echo ;echo ;echo ;
wget -O counts/tmp/GSE179213_new_sequence.vs_dm6main.dm6_genes.mapspliceMulti.MpBC.counts.txt.gz https://ftp.ncbi.nlm.nih.gov/geo/series/GSE179nnn/GSE179213/suppl/GSE179213_new_sequence.vs_dm6main.dm6_genes.mapspliceMulti.MpBC.counts.txt.gz
gzip -d counts/tmp/GSE179213_new_sequence.vs_dm6main.dm6_genes.mapspliceMulti.MpBC.counts.txt.gz
sleep 3;
echo ;echo ;echo ;echo ;echo ;echo ;echo ;echo ;echo ;echo ;echo ;echo ;



if [[ $full -gt 0 ]] 
	then echo "downloading problematic Fru...";
	echo ;echo ;echo ;echo ;echo ;
	wget -O counts/tmp/GSE179213_FRULEXAFRU440_1_new_sequence_update.vs_dm6main.dm6_genes.mapspliceMulti.MpBC.counts.txt.gz https://ftp.ncbi.nlm.nih.gov/geo/series/GSE179nnn/GSE179213/suppl/GSE179213_FRULEXAFRU440_1_new_sequence_update.vs_dm6main.dm6_genes.mapspliceMulti.MpBC.counts.txt.gz ;
	gzip -d counts/tmp/GSE179213_FRULEXAFRU440_1_new_sequence_update.vs_dm6main.dm6_genes.mapspliceMulti.MpBC.counts.txt.gz ;
	sleep 3;
	echo ;echo ;echo ;echo ;echo ;echo ;echo ;echo ;echo ;echo ;echo ;echo ;



	echo "joining tables...";
# head -n 1 counts/tmp/GSE179213_new_sequence.vs_dm6main.dm6_genes.mapspliceMulti.MpBC.counts.txt >  counts/tmp/header.1
# head -n 1 counts/tmp/GSE179213_FRULEXAFRU440_1_new_sequence_update.vs_dm6main.dm6_genes.mapspliceMulti.MpBC.counts.txt | rev | cut -f 1 | rev >  counts/tmp/header.2
# paste  counts/tmp/header.1  counts/tmp/header.2 > counts/tmp/header
	paste  <( head -n 1 counts/tmp/GSE179213_new_sequence.vs_dm6main.dm6_genes.mapspliceMulti.MpBC.counts.txt )  <( head -n 1 counts/tmp/GSE179213_FRULEXAFRU440_1_new_sequence_update.vs_dm6main.dm6_genes.mapspliceMulti.MpBC.counts.txt | rev | cut -f 1 | rev ) | tr -d '\r' > counts/tmp/header ; 
	sleep 1;echo ;
# tail -n +2 counts/tmp/GSE179213_new_sequence.vs_dm6main.dm6_genes.mapspliceMulti.MpBC.counts.txt > counts/tmp/body.1
# tail -n +2 counts/tmp/GSE179213_FRULEXAFRU440_1_new_sequence_update.vs_dm6main.dm6_genes.mapspliceMulti.MpBC.counts.txt > counts/tmp/body.2
# paste counts/tmp/body.1 counts/tmp/body.2 > counts/tmp/body
	paste <( tail -n +2 counts/tmp/GSE179213_new_sequence.vs_dm6main.dm6_genes.mapspliceMulti.MpBC.counts.txt )  <( tail -n +2 counts/tmp/GSE179213_FRULEXAFRU440_1_new_sequence_update.vs_dm6main.dm6_genes.mapspliceMulti.MpBC.counts.txt| rev | cut -f 1 | rev ) | tr -d '\r' > counts/tmp/body ;
	echo ;
	echo "writing new data table to 'counts/all.vs_dm6main.dm6_genes.mapspliceMulti.MpBC.counts' ";
	cat counts/tmp/header counts/tmp/body > counts/all.vs_dm6main.dm6_genes.mapspliceMulti.MpBC.counts ;
	sleep 1;echo ;
else echo "moving data table to 'counts/all.vs_dm6main.dm6_genes.mapspliceMulti.MpBC.counts' ";
	cat counts/tmp/GSE179213_new_sequence.vs_dm6main.dm6_genes.mapspliceMulti.MpBC.counts.txt > counts/all.vs_dm6main.dm6_genes.mapspliceMulti.MpBC.counts
	sleep 1;echo ;
fi
echo ;sleep 1;
echo " .... cleaning up......."
rm -rf counts/tmp/
echo ;sleep 1;
echo "!!!!!!!	DONE	!!!!!!!!!"

