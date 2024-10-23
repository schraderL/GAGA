if [ $# -lt 2 ];then
        echo "Usage : sh $0 <GAGA_ID> <ORF.txt>"
        exit
fi

genome=$1
orf=$2


ln -s ./inputGenome/$genome.fasta $genome.fasta

perl scripts/cds_from_bestorf.pl $orf.txt $genome.stringtie.gtf.fasta 
awk '(!$5)' $genome.stringtie.gtf.fasta.orf > $genome.stringtie.gtf.fasta.orf.filter
perl scripts/getGff4orf.pl $genome.stringtie.gtf.fasta.orf.filter $genome.stringtie.gtf $genome.fasta > $genome.start.stop.gff
perl scripts/snap_frame.pl $genome.start.stop.gff > $genome.start.stop.gff.frame.gff
perl scripts/getGene.pl $genome.start.stop.gff.frame.gff $genome.fasta > $genome.start.stop.gff.frame.gff.cds
perl scripts/cds2aa.pl --check $genome.start.stop.gff.frame.gff.cds > $genome.start.stop.gff.frame.gff.cds.check
perl scripts/filter_gff_gene_lenght.pl --threshold 350 --exons 2 $genome.start.stop.gff.frame.gff > $genome.start.stop.gff.frame.gff.filter
perl scripts/filter_gff_gene_lenght.pl --threshold 500 --exons 1 $genome.start.stop.gff.frame.gff.filter > $genome.start.stop.gff.frame.gff.filter.gff
perl scripts/renameID.pl $genome.start.stop.gff.frame.gff.filter.gff $1\_RNAseq > $genome.start.stop.gff.frame.gff.filter.gff.renamed.gff 2> renamed.lst
perl scripts/perfect_gene/clustergff.pl.new.pl $genome.start.stop.gff.frame.gff.filter.gff.renamed.gff multi 1 3 CDS
