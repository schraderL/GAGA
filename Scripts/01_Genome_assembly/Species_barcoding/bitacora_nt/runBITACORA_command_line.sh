#!/bin/bash

### Modified to search nucleotides in a genome sequence

##########################################################
##                                                      ##
##                       BITACORA                       ##
##                                                      ##
##          Bioinformatics tool to assist the           ##
##      comprehensive annotation of gene families       ##
##                                                      ##
##                                                      ##
## Developed by Joel Vizueta                            ##
## Contact: via github or jvizueta@ub.edu               ##
##                                                      ##
##########################################################

VERSION=1.2.1

# Default values for editable parameters
CLEAN=F
EVALUE=1e-5
THREADS=1
ALGORITHM=gemoma
MAXINTRON=15000
GENOMICBLASTP=F
ADDFILTER=T
FILTERLENGTH=30

BITMODE=full
BPATH=''
HPATH=''
SCRIPTDIR='/home/projects/ku_00039/people/joeviz/bitacora_barcoding/bitacora_nt/Scripts'
GEMOMAP='/home/projects/ku_00039/people/joeviz/programs/GeMoMa-1.6.4/GeMoMa-1.6.4.jar'
NAME=Out

function usage {
	echo
	echo "Usage:"
	echo "  $(basename $0)"
	echo "      -m BITACORA runnning mode. Specify 'full', 'genome' or 'protein' (Default=$BITMODE)"
	echo "      -q Folder containing the query database or multiple databases, named as Example1_db.fasta and Example1_db.hmm (Mandatory)"
	echo "      -g Genome fasta file (Mandatory in full and genome mode)"
	echo "      -f GFF file (Mandatory in full mode)"
	echo "      -p Protein fasta file (Mandatory in full and protein mode)"
	echo "      -n Organism Name, i.e. 'Dmel' (Default=$NAME)"
	echo "      -sp PATH to BITACORA Scripts folder (Mandatory)"
	echo "      -gp PATH to GeMoMa jar file (Mandatory in case of using this algorithm to predict novel genes)"
	echo "      -bp PATH to BLAST binaries (Optional if it is already included in PATH)"
	echo "      -hp PATH to HMMER executable (Optional if it is already included in PATH)"
	echo "      -c Clean output files if '-c T'. Specify 'T' or 'F' (Default=$CLEAN)"
	echo "      -e E-value used to filter BLAST and HMMER output (Default=$EVALUE)"
	echo "      -t Number of threads (Default=$THREADS)"
	echo "      -a Algorithm used to predict novel genes. Specify 'gemoma' or 'proximity' (Default=$ALGORITHM)"
	echo "      -i Maximum intron length to join putative exons in the close-proximity algorithm (Default=$MAXINTRON)"
	echo "      -b Specify '-b T' to conduct an additional BLASTP search in addition to HMMER to validate novel genes (Default=$GENOMICBLASTP)"	
	echo "      -r Conduct an additional filtering of the annotations if -r T. Specify 'T' or 'F' (Default=$ADDFILTER)"		
	echo "      -l Minimum length to retain identified genes (Default=$FILTERLENGTH)"											
	echo "      -h Show this help. See BITACORA documentation for further details about each parameter"
	echo ""
	echo "  Example for running BITACORA full mode using GeMoMa algorithm and 4 threads:"
	echo "      ./$(basename $0) -m full -sp /path/to/Scripts -gp /path/to/GeMoMa.jar -n Dmel -g /path/to/genome.fasta -f /path/to/GFF.gff3 -p /path/to/protein.fasta -q /path/to/query_folder -t 4"
	echo ""	
	echo "  Example for running BITACORA genome mode using GeMoMa algorithm:"
	echo "      ./$(basename $0) -m genome -sp /path/to/Scripts -gp /path/to/GeMoMa.jar -n Dmel -g /path/to/genome.fasta -q /path/to/query_folder"
	echo ""	
	echo "  Example for running BITACORA protein mode and 2 threads:"
	echo "      ./$(basename $0) -m protein -sp /path/to/Scripts -n Dmel -p /path/to/protein.fasta -q /path/to/query_folder -t 2"
	echo ""	
	echo ""
	exit 0
}


# Read options

if [ "$#" -lt "1" ]; # at least 1 argument
then
	usage
fi

while [ $# -gt 0 ]; do
	case "$1" in
		-h|-help) usage
				;;
		-m) shift
			BITMODE=$1
			;;				
		-bp) shift
			BPATH=$1
			;;
		-hp) shift
			HPATH=$1
			;;
		-sp) shift
			SCRIPTDIR=$1
			;;
		-gp) shift
			GEMOMAP=$1
			;;
		-n)	shift
			NAME=$1
			;;	
		-g)	shift
			GENOME=$1
			;;	
		-f)	shift
			GFFFILE=$1
			;;	
		-p)	shift
			PROTFILE=$1
			;;	
		-q)	shift
			QUERYDIR=$1
			;;		
		-c)	shift
			CLEAN=$1
			;;							
		-e)	shift
			EVALUE=$1
			;;	
		-t)	shift
			THREADS=$1
			;;	
		-a)	shift
			ALGORITHM=$1
			;;	
		-i)	shift
			MAXINTRON=$1
			;;	
		-b)	shift
			GENOMICBLASTP=$1
			;;	
		-r)	shift
			ADDFILTER=$1
			;;	
		-l)	shift
			FILTERLENGTH=$1
			;;							
		*)	echo 
			echo "ERROR - Invalid option: $1"
			echo
			usage
			;;
	esac
	shift
done


GEMOMA=''
if [ $ALGORITHM == "gemoma" ] ; then
	GEMOMA=T
elif [ $ALGORITHM == "proximity" ] ; then
	GEMOMA=F
else
	echo -e "\nERROR, no recognized algorithm was detected. Please specify 'gemoma' or 'proximity' in -a option\n"
	exit 1;
fi

export PATH=$BPATH:$PATH
export PATH=$HPATH:$PATH


if [ $BITMODE == "nt" ] ; then
	##########################################################
	##                   PIPELINE - CODE                    ##
	##########################################################

	echo -e "\n#######################  Running BITACORA nucleotide mode  #######################";
	echo "BITACORA protein-mode version $VERSION";
	date



	# Run step 1

	perl $SCRIPTDIR/runanalysis_nt_mode.pl $NAME $PROTFILE $QUERYDIR $EVALUE $THREADS 2>>BITACORAstd.err

	ERRORCHECK="$(grep -c 'ERROR' BITACORAstd.err)"

	if [ $ERRORCHECK != 0 ]; then
		cat BITACORAstd.err;
		echo -e "BITACORA died with error\n";
		exit 1;
	fi

	ERRORCHECK="$(grep -c 'Segmentation' BITACORAstd.err)"

	if [ $ERRORCHECK != 0 ]; then
		cat BITACORAstd.err;
		echo -e "BITACORA died with error\n";
		exit 1;
	fi


	# Run additional filtering and clustering
		
	if [ $ADDFILTER == "T" ] ; then
		perl $SCRIPTDIR/runfiltering_protein_mode.pl $NAME $QUERYDIR $FILTERLENGTH 2>>BITACORAstd.err 2>BITACORAstd.err
	fi

	ERRORCHECK="$(grep -c 'ERROR' BITACORAstd.err)"

	if [ $ERRORCHECK != 0 ]; then
		cat BITACORAstd.err;
		echo -e "BITACORA died with error\n";
		exit 1;
	fi


	# Cleaning 

	if [ $CLEAN = "T" ]; then
		perl $SCRIPTDIR/runcleaning_protein_mode.pl $NAME $QUERYDIR
		echo -e "Cleaning output folders\n";
	fi

	rm */*proteins.fasta # delete full scaffold

	rm BITACORAstd.err

	echo -e "BITACORA completed without errors :)";
	date

	exit

else
	echo -e "\nERROR, please specify a valid BITACORA running mode: 'nt' \n";
	exit
fi

exit

