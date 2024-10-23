
if [ $# -lt 3 ];then
        echo "Usage : sh $0 <GAGA_ID> <trainingSet> <species>"
        exit
fi

genome=$1
perl scripts/autoAugTrain.pl --genome=$genome.fasta --trainingset=$2 --species=$3  --useexisting 2> auto_training.log
perl scripts/Run_Augustus.pl --species $genome --run T --queue st.q --pro_code P18Z10200N0102 --outdir ./output1 --maxjob 50 $genome.fasta 
perl scripts/Check_GFF.pl -check_cds -mini_cds 150 -cds_ns 10 -outdir ./ $genome.fasta.augustus.gff ./$genome.fasta
perl scripts/cds2aa.pl --check $genome.fasta.augustus.cds > $genome.fasta.augustus.cds.check
perl scripts/fishInWinter.pl -bf table -ff gff -except $genome.fasta.augustus.cds.check $genome.fasta.augustus.gff.check.gff > $genome.fasta.augustus.gff.check.filterIncomplete.gff
perl scripts/filter_gff_gene_lenght.pl --threshold 350 --exons 2 $genome.fasta.augustus.gff.check.filterIncomplete.gff > $genome.fasta.augustus.gff.check.filterIncomplete.gff.filter
perl scripts/filter_gff_gene_lenght.pl --threshold 500 --exons 1 $genome.fasta.augustus.gff.check.filterIncomplete.gff.filter > $genome.fasta.augustus.gff.check.filterIncomplete.gff.filter.gff
perl scripts/trem2GeneName.pl ../Swissprot1/$genome.fasta.augustus.pep.SwissProt.blast /ldfssz1/ST_DIVERSITY/PUB/USER/fangqi/local/database/uniprot/release-2020_05/swissport/uniprot_sprot.Eukaryota.fasta.simple | grep -iP "transpose|transposon|retro-transposon|retrovirus|retrotransposon|reverse transcriptase|transposase|retroviral" > TE.lst
perl scripts/fishInWinter.pl -bf table -ff gff -except TE.lst $genome.fasta.augustus.gff.check.filterIncomplete.gff.filter.gff > $genome.fasta.augustus.gff.check.filterIncomplete.gff.filterShort.filterTE.gff
