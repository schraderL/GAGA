#!/usr/bin/perl -w
use strict;

my $file = shift;
open IN, "<$file" or die "failed to open $file:$!\n";

while (<IN>) {
	my $id = $_;
	chomp $id;
	if(-e "/run/media/dell/storage1/User/xiongzj/GAGA_project/15.new_anntations/Gene_Re-annotations/Original_annotations_to_replace/$id.gff3") {
		open F, "/run/media/dell/storage1/User/xiongzj/GAGA_project/15.new_anntations/Gene_Re-annotations/Original_annotations_to_replace/$id.gff3" or die "failed to open $id.gff3\n";
	}elsif(-e "/run/media/dell/storage1/User/xiongzj/GAGA_project/15.new_anntations/Gene_Re-annotations/Original_annotations_to_replace/$id\_final_annotation_repfilt.gff3") {
		open F, "/run/media/dell/storage1/User/xiongzj/GAGA_project/15.new_anntations/Gene_Re-annotations/Original_annotations_to_replace/$id\_final_annotation_repfilt.gff3" or die "failed to open $!\n";
	}
	my $prefix = '';
	while (<F>) {
		my @t=split /\t/;
		$prefix = $1 if(/ID=(\w+)_g\d+/);	
	}
	close (F);
	print "cd /run/media/dell/storage1/User/xiongzj/GAGA_project/15.new_anntations/Gene_Re-annotations1/Final_annotations_output_toZijun/$id/; rm -rf CombinedEvidences1; mkdir CombinedEvidences1; cd CombinedEvidences1; nohup sh ../../bin/pipeline.sh $id $prefix & \n";
}
close (IN);
