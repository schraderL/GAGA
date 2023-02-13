if [ $# -lt 1 ];then
        echo "Usage : sh $0 <GAGA_ID> <Tag>"
        exit
fi

genome=$1
tag=$2


grep mR $genome.combined.gff | grep GeMoMa | perl -ne 'if(/ID=([^;]+);/){print "$1\n"}' > Gemoma.lst
perl -ne 'if(/(\w+_\d+)\./){print "$1\n"}' Gemoma.lst > Gemoma.lst.gene

perl /jdfssz1/ST_EARTH/P18Z10200N0102/xiongzj/GAGA/04.Combined_gene/bin/fishInWinter.pl -bf table -ff gff Gemoma.lst.gene ./input/$genome.final_annotation.gff > Gemoma.lst.gene.gff
ln -s /jdfssz1/ST_EARTH/P18Z10200N0102/xiongzj/GAGA/01.RNAseq_mapping/03.ORF/$genome/*.filter.gff.renamed.gff $genome.RNAseq.all.gff
perl ~/bin/FindOverlapAtCDSlevel.pl $genome.combined.gff $genome.RNAseq.all.gff > Overlap.lst4
awk '($10)'  Overlap.lst4 > Overlap.lst4.CDS
awk '{print $1"\t"$2"\t"$11"\t"$12}'  Overlap.lst4.CDS > id.lst
perl -ne '@t=split /\s+/;if(/\d\.\d/){$t[0]=~s/\.\d+$//;print "$t[0]\t$t[1]\t$t[2]\t$t[3]\n"}else {print "$t[0]\t$t[1]\t$t[2]\t$t[3]\n"}' id.lst > id.lst.filter

perl /jdfssz1/ST_EARTH/P18Z10200N0102/xiongzj/GAGA/04.Combined_gene/bin/merge.pl id.lst.filter $genome.RNAseq.all.gff > $genome.RNAseq.all.gff.merge.gff
grep AUGUSTUS $genome.combined.gff > AUGUSTUS.gff
perl /jdfssz1/ST_EARTH/P18Z10200N0102/xiongzj/GAGA/04.Combined_gene/bin/addParent1.pl AUGUSTUS.gff > AUGUSTUS.add.gff
perl /jdfssz1/ST_EARTH/P18Z10200N0102/xiongzj/GAGA/04.Combined_gene/bin/addParent.pl Gemoma.lst.gene.gff > Gemoma.lst.gene.gff.add.gff
cat Gemoma.lst.gene.gff.add.gff $genome.RNAseq.all.gff.merge.gff AUGUSTUS.add.gff  > merge.gff

perl /home/xiongzj/bin/FindOverlapAtCDSlevel.pl merge.gff merge.gff > Overlap.lst5
awk '($1!=$2 && $11==1 && $12==1)' Overlap.lst5 | awk '{print $1"\t"$2"\t"0.3}' > Overlap.lst5.id
perl /jdfssz1/ST_EARTH/P18Z10200N0102/xiongzj/GAGA/04.Combined_gene/bin/hcluster.pl Overlap.lst5.id -verbose > Overlap.lst5.id.cluster 2>log
grep mR $genome.combined.gff | perl -ne 'if(/ID=([^;]+);/){print $1,"\n"}' > $genome.combined.gff.geneID
perl -ne '@t=split /\s+/;shift @t; shift @t;print join("\n",@t),"\n"' Overlap.lst5.id.cluster > Overlap.lst5.id.cluster.all
perl /jdfssz1/ST_EARTH/P18Z10200N0102/xiongzj/GAGA/04.Combined_gene/bin/retain.pl $genome.combined.gff.geneID Overlap.lst5.id.cluster > retain.lst
perl ~/bin/fishInWinter.pl -except retain.lst Overlap.lst5.id.cluster.all > Overlap.lst5.id.cluster.all.remove
perl /jdfssz1/ST_EARTH/P18Z10200N0102/xiongzj/GAGA/04.Combined_gene/bin/iddelete.pl Overlap.lst5.id.cluster.all.remove merge.gff > merge.non-redundant.gff


perl /jdfssz1/ST_EARTH/P18Z10200N0102/xiongzj/GAGA/04.Combined_gene/bin/gff3Merge.pl merge.non-redundant.gff > merge.gff3
grep gene merge.gff3 > gene.lst

perl /jdfssz1/ST_EARTH/P18Z10200N0102/xiongzj/GAGA/04.Combined_gene/bin/renameID.pl gene.lst $tag > $tag.renamed.lst
perl /jdfssz1/ST_EARTH/P18Z10200N0102/xiongzj/GAGA/04.Combined_gene/bin/renameID2.pl $tag.renamed.lst merge.gff3 > merge.renamed.gff3 2>final_renamed.lst
perl ~/bin/fishInWinter.pl -bf gff -ff table $genome.combined.gff final_renamed.lst > $genome.combined.gff.gene.lst
perl ~/bin/fishInWinter.pl -bf table -ff gff -bc 2 $genome.combined.gff.gene.lst merge.renamed.gff3 > merge.renamed.representative.gff3

/jdfssz1/ST_EARTH/P18Z10200N0102/xiongzj/software/genometools-1.6.1/bin/gt gff3 -force -sort -tidy -retainids -o merge.renamed.fixed.gff3 merge.renamed.gff3
/jdfssz1/ST_EARTH/P18Z10200N0102/xiongzj/software/genometools-1.6.1/bin/gt gff3 -force -sort -tidy -retainids -o merge.renamed.fixed.representative.gff3	merge.renamed.representative.gff3

