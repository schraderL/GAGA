
if [ $# -lt 2 ];then
        echo "Usage : sh $0 <GAGA_ID> <prefix>"
        exit
fi

genome=$1
prefix=$2

rm -rf nohup.out *.gz
mv $genome.gff3 $genome.gff3.link
awk '{print $0";"}' $genome.gff3.link> $genome.gff3
sh ../../bin/pipeline_test.sh $genome $prefix
