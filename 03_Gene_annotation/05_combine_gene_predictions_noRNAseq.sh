if [ $# -lt 1 ];then
        echo "Usage : sh $0 <GAGA_ID>"
        exit
fi

genome=$1


ln -s /jdfssz1/ST_EARTH/P18Z10200N0102/xiongzj/GAGA/02.augustus/$genome/output1/*.filterTE.gff $genome.denovo.gff
ln -s ./input/$genome.final_annotation.gff.nr.gff $genome.homolog.gff
perl ~/bin/findOveralpCds.pl $genome.homolog.gff $genome.denovo.gff > Overlap.lst2 
grep CDS Overlap.lst2 > Overlap.lst2.CDS
awk '($7==$13)' Overlap.lst2.CDS > Overlap.lst2.CDS.same
perl ~/bin/fishInWinter.pl -bf table -ff gff -bc 9 -except Overlap.lst2.CDS.same $genome.denovo.gff > $genome.denovo.gff.except.gff

ln -s /jdfssz1/ST_EARTH/P18Z10200N0102/xiongzj/GAGA/02.augustus/$genome/Swissprot1/*.blast Swissprot.blast
perl ~/bin/fishInWinter.pl -bf table -ff gff Swissprot.blast $genome.denovo.gff.except.gff > Swissprot.gff
perl ~/bin/fishInWinter.pl -bf table -ff gff -except  Swissprot.blast $genome.denovo.gff.except.gff > $genome.denovo.gff.except.gff.except.gff

ln -s /jdfssz1/ST_EARTH/P18Z10200N0102/xiongzj/GAGA/02.augustus/$genome/nr/*.nr.blast  Nr.blast
perl ~/bin/fishInWinter.pl -bf table -ff gff Nr.blast $genome.denovo.gff.except.gff.except.gff > Nr.gff
cat $genome.homolog.gff Swissprot.gff Nr.gff >  $genome.combined.gff
