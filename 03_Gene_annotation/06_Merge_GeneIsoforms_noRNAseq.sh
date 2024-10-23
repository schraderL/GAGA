if [ $# -lt 1 ];then
        echo "Usage : sh $0 <GAGA_ID> <Tag>"
        exit
fi

genome=$1
tag=$2


grep mR $genome.combined.gff | grep GeMoMa | perl -ne 'if(/ID=([^;]+);/){print "$1\n"}' > Gemoma.lst
perl -ne 'if(/(\w+_\d+)\./){print "$1\n"}' Gemoma.lst > Gemoma.lst.gene

perl /jdfssz1/ST_EARTH/P18Z10200N0102/xiongzj/GAGA/04.Combined_gene/bin/fishInWinter.pl -bf table -ff gff Gemoma.lst.gene ./input/$genome.final_annotation.gff > Gemoma.lst.gene.gff

grep AUGUSTUS $genome.combined.gff > AUGUSTUS.gff
perl /jdfssz1/ST_EARTH/P18Z10200N0102/xiongzj/GAGA/04.Combined_gene/bin/addParent1.pl AUGUSTUS.gff > AUGUSTUS.add.gff
perl /jdfssz1/ST_EARTH/P18Z10200N0102/xiongzj/GAGA/04.Combined_gene/bin/addParent.pl Gemoma.lst.gene.gff > Gemoma.lst.gene.gff.add.gff
cat Gemoma.lst.gene.gff.add.gff AUGUSTUS.add.gff > merge.gff

perl /jdfssz1/ST_EARTH/P18Z10200N0102/xiongzj/GAGA/04.Combined_gene/bin/gff3Merge.pl merge.gff > merge.gff3
grep gene merge.gff3 > gene.lst

perl /jdfssz1/ST_EARTH/P18Z10200N0102/xiongzj/GAGA/04.Combined_gene/bin/renameID.pl gene.lst $tag > $tag.renamed.lst
perl /jdfssz1/ST_EARTH/P18Z10200N0102/xiongzj/GAGA/04.Combined_gene/bin/renameID2.pl $tag.renamed.lst merge.gff3 > merge.renamed.gff3 2>final_renamed.lst
perl ~/bin/fishInWinter.pl -bf gff -ff table $genome.combined.gff final_renamed.lst > $genome.combined.gff.gene.lst
perl ~/bin/fishInWinter.pl -bf table -ff gff -bc 2 $genome.combined.gff.gene.lst merge.renamed.gff3 > merge.renamed.representative.gff3

/jdfssz1/ST_EARTH/P18Z10200N0102/xiongzj/software/genometools-1.6.1/bin/gt gff3 -force -sort -tidy -retainids -o merge.renamed.fixed.gff3 merge.renamed.gff3
/jdfssz1/ST_EARTH/P18Z10200N0102/xiongzj/software/genometools-1.6.1/bin/gt gff3 -force -sort -tidy -retainids -o merge.renamed.fixed.representative.gff3	merge.renamed.representative.gff3

