
if [ $# -lt 2 ];then
        echo "Usage : sh $0 <GAGA_ID> <prefix>"
        exit
fi

genome=$1
prefix=$2

cat ../*/*.gff* > annotpipeline.gff3
cat ../../Final_CSP_annotations/$genome*.gff3 ../../Final_GR_annotations_nofragment/$genome*.gff3 ../../Final_IRiGluR_annotations_nofragment/$genome*.gff3 ../../Final_NPC2_annotations/$genome*.gff3 ../../Final_OBP_annotations/$genome*.gff3 ../../Final_OR_annotations_nofragment/$genome*.gff3 ../../Final_PPK_annotations/$genome*.gff3 ../../Final_SNMP_annotations/$genome*.gff3 > OR.gff3

ln -s /run/media/dell/storage1/User/xiongzj/GAGA_project/15.new_anntations/Gene_Re-annotations/Original_annotations_to_replace/$genome\_final_annotation_repfilt.gff3 $genome.gff3

perl bin/extractID.pl annotpipeline.gff3 $prefix
perl bin/Step1_geneReplace.pl $genome.gff3 annotpipeline.gff3 $prefix.mRNA.lst $prefix.gene.lst $prefix > $genome.gff3.replace.gff3

cat OR.gff3 ./$genome.gff3.replace.gff3 > ./$genome.gff3.replace2.gff3
perl bin/FindOverlapAtCDSlevel.pl $genome.gff3.replace2.gff3 $genome.gff3.replace2.gff3 > Overlap.lst.replace
awk '($10>0)' Overlap.lst.replace > Overlap.lst.replace.filter1
awk '($1!=$2)' Overlap.lst.replace.filter1 > Overlap.lst.replace.filter2
perl bin/filter1.pl Overlap.lst.replace.filter2 $prefix > Overlap.lst.replace.filter3
perl bin/Step3_deleteOverlappingGenes.pl $genome.gff3.replace2.gff3 Overlap.lst.replace.filter3 $prefix > $genome.gff3.replace3.gff3
perl bin/addParent.pl $genome.gff3.replace3.gff3 > $genome.gff3.replace4.gff3
perl bin/addCDS.pl $genome.gff3.replace4.gff3 > $genome.gff3.replace5.gff3
perl bin/addTagCDSlength.pl $genome.gff3.replace5.gff3 $prefix > $genome.gff3.replace5.gff3.len 
perl bin/selectIsoform.pl $genome.gff3.replace5.gff3.len > $genome.gff3.replace5.gff3.len.representative
perl bin/filter2.pl $genome.gff3.replace5.gff3.len.representative $genome.gff3.replace5.gff3 > $genome.gff3.replace5.gff3.filter.gff3
perl bin/Step4_selectRepresentative.pl $genome.gff3.replace5.gff3.len.representative $genome.gff3.replace5.gff3.filter.gff3 > $genome.gff3.replace6.gff3
cp $genome.gff3.replace5.gff3.filter.gff3 $genome\_final_annotation_repfilt_addfunc.gff3
cp $genome.gff3.replace6.gff3 $genome\_final_annotation_repfilt_addfunc.representative.gff3
gzip $genome\_final_annotation_repfilt_addfunc.gff3
gzip $genome\_final_annotation_repfilt_addfunc.representative.gff3
