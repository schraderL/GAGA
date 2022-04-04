#!/usr/bin/bash

##################################
# USAGE: bash annotateArrays.sh GAGA-0001_OR.gff3 <max-distance-between-genes> <min-number-of-genes-per-array>
# e.g. bash ../code/ORs/annotateArrays.sh /Users/lukas/sciebo/Projects/GAGA/0.data/Final_OR_annots/GAGA-0001_ABCENTH_clean_OR_renamed_all.gff3 30000 2
# This Script can be used to annotate OR gene tandem arrays from gff3 files. 
##################################

##################################
#Define Functions
##################################


# extract just genes from gff
##################################
getGenes() {
    awk 'BEGIN{FS="\t";OFS="\t"} {if ($3 == "gene") print $0}'
}

# cut away everything after the gene id
##################################
format() {
    perl -pe 's/ID=(.*?);.*/$1/g'
}

#cutoff based on number of genes per array
##################################
cutoff() {
    awk -v cut=$1 'BEGIN{FS="\t";OFS="\t"} {if ($5 >= cut) print $0}'
}

# cut away everything after the gene id
##################################
formatOut() {
    awk 'BEGIN{FS="\t";OFS="\t"} {print $1,$2,$3,"ID=ORArray"NR";genes="$4";genesContained="$5";first="$6";last="$7}'
}

##################################
# Run commands
##################################

cat $1 | getGenes |format  | bedtools merge -d $2 -c 9,9,9,9 -o collapse,count,first,last | cutoff $3 |formatOut

##################################
#USE
##################################
# parallel

#ls /Users/lukas/sciebo/Projects/GAGA/0.data/Final_OR_annots/*_ABCENTH_clean_OR_renamed_all.gff3 |perl -pe 's/.*\///g'|cut -f 1 -d "_" > allIDs.lst
#cat allIDs.lst|parallel --dryrun "bash ../code/ORs/annotateArrays.sh /Users/lukas/sciebo/Projects/GAGA/0.data/Final_OR_annots/{}_ABCENTH_clean_OR_renamed_all.gff3 30000 2 > {}.ORarrays.bed"
#cat allIDs.lst|parallel "bash ../code/ORs/annotateArrays.sh /Users/lukas/sciebo/Projects/GAGA/0.data/Final_OR_annots/{}_ABCENTH_clean_OR_renamed_all.gff3 30000 2 > {}.ORarrays.bed"