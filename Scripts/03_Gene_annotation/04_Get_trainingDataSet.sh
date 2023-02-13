if [ $# -lt 1 ];then
        echo "Usage : sh $0 <GAGA_ID>"
        exit
fi

genome=$1

ln -s /run/media/dell/storage1/User/xiongzj/GAGA_project/10.Final_genomes/genomes_repeatmask/$genome*
zcat $genome* > genome.fasta
perl /run/media/dell/data/User/xiongzj/bin/cluster/bgigff.pl final_annotation.gff > final_annotation.bgi.gff
mv final_annotation.gff final_annotation.gemoma.gff
ln -s final_annotation.bgi.gff final_annotation.gff
perl /run/media/dell/data/User/xiongzj/bin/cluster/clustergff.pl final_annotation.gff multi 1 3 CDS
perl /run/media/dell/data/User/xiongzj/bin/GACP/getGene.pl final_annotation.gff.nr.gff genome.fasta > final_annotation.gff.nr.gff.cds
perl /run/media/dell/data/User/xiongzj/bin/GACP/fastaDeal.pl --attr id:len:lenwogap final_annotation.gff.nr.gff.cds |awk '($2==$3)' | awk '($3>200)'  > final_annotation.gff.cds.noGap.lst
perl /run/media/dell/data/User/xiongzj/bin/fishInWinter.pl -bf table -ff gff final_annotation.gff.cds.noGap.lst  final_annotation.gff.nr.gff > final_annotation.gff.cds.noGap.lst.gff
grep mR final_annotation.gff.cds.noGap.lst.gff | perl -ne 'if(/ID=(\w+-\d+_\d+\.\d+);.+score=(\d+)/){print $1,"\t",$2,"\n"}' | sort -k2nr | head -1000  > highScore.lst
perl /run/media/dell/data/User/xiongzj/bin/fishInWinter.pl -bf table -ff gff highScore.lst final_annotation.gff.nr.gff > highScore.lst.gff
perl /run/media/dell/data/User/xiongzj/bin/GACP/getGene.pl highScore.lst.gff genome.fasta > highScore.lst.gff.cds
perl /run/media/dell/data/User/xiongzj/bin/GACP/cds2aa.pl --check highScore.lst.gff.cds > highScore.lst.gff.cds.check

