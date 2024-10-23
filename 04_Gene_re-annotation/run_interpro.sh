#!/bin/sh

mkdir -p Interpro_out # Create output dir for interproscan
cd Interpro_out
interproscan.sh -i /home/projects/ku_00039/people/joeviz/GAGA_genomes/GAGA-0001_final_annotation_representative_repfilt.pep.fasta -t p -goterms -iprlookup -d /home/projects/ku_00039/people/joeviz/GAGA_genomes/Interpro_out -cpu 40
    
