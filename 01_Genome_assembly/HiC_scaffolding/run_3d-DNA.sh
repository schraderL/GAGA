#!/bin/bash

if [ $# -lt 1 ];then
        echo "Usage : sh $0 <GAGA_ID>"
        exit
fi

genome=$1


/usr/bin/python /run/media/dell/data/public/software/3d-DNA/juicer-master/misc/generate_site_positions.py MboI $genome.genome.fasta $genome.genome.fasta
perl fastaDeal.pl --attr id:len $genome.genome.fasta > $genome.genome.fasta_size.txt
# Run juicer, Note that the Hi-C read files (in fastq) should be in the same directory in which this script is run. See juicer.sh help: "[topDir]/fastq  - Should contain the fastq files."
sh /run/media/dell/storage1/User/xiongzj/software/juicer_v2.sh -z $genome.genome.fasta -p $genome.genome.fasta_size.txt -y $genome.$genome.fasta_MboI.txt -s MboI -d Juicer_output -D /run/media/dell/storage1/User/xiongzj/software/juicer/bin -T 2
sh /run/media/dell/data/public/software/3d-DNA/3d-dna-master/run-asm-pipeline.sh -m haploid $genome.genome.fasta Juicer_output/aligned/merged_nodups.txt 
