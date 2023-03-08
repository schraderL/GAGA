#!usr/bin/perl -w
use strict;
use Getopt::Long;
#use Support_Program;

=head1 Name

   get_annot_info.pl

=head1 Author

   ChuanJun  chuanjun@genomics.org.cn
   2010-08-19

=head1 Description

   
=head1 Options
   
   perl get_annot_info.pl [options]
   -nohead   do not show the first instruction line
   -input    input the BLAST result with the m8 option
   -id       input the id file of the database (id	annot_information)
   -out      give the output file's name and the path
   -topmatch integer, to set the top boundary about how many results that one query matched one subjects to dispaly
   -tophit   integer, to set how many subjects for a query to be display
   -eval     float or exponent,to filter the alignments which worse than the E-value cutoff
   -help     output the help information

=head1 Usage

   perl get_annot_info.pl -nohead -input /ifs1/DGE_SR/liwansh/special_item/orange/orange_result/annot/database/All-Unigene/gene-annotation/All-Unigene.fa.blast.kegg.m8.xls -out ../kegg.out.xls -topmatch 1 -tophit 20 -id id.xls

=cut

my ($input,$id,$outfile,$topmatch,$tophit,$help,$eval,$nohead);
GetOptions(
    "input:s"=>\$input,
	"id:s"=>\$id,
	"out:s"=>\$outfile,
	"topmatch:i"=>\$topmatch,
	"tophit:i"=>\$tophit,
	"e:f"=>\$eval,
	"nohead:s"=>\$nohead,
	"help"=>\$help
);
die `pod2text $0` if($help);
die `pod2text $0` unless(defined $input && defined $id && defined $outfile && defined $topmatch && defined $tophit);
my %hash_id;
open INI,"<$id" || die "Can't open the id file $id,maybee you nedd help from zhangfx!  $!";
while(<INI>){
    chomp;
	next if(/^(\s+)$/ || /^#/);
	my @ids = split /\t/;
	$ids[1] = "---" unless(defined $ids[1]);
	$hash_id{$ids[0]} = $ids[1];
}
open INR,"<$input" || die "Can't open the input m8's blast file $input! $!";
open OUT,">$outfile" || die $!;
print OUT "Query_id\tSubject_id\tIdentity\tAlign_length\tMiss_match\tGap\tQuery_start\tQuery_end\tSubject_start\tSubject_end\tE_value\tScore\tSubject_annotation\n" unless(defined $nohead);
my @arrs;
my $pre = " ";
while(<INR>){
    chomp;
	next if(/^(\s+)$/ || /^#/);
	my @data = split (/\t/,$_,3);
	my @evals = split /\t/,$_;
	next if(defined $eval && $evals[10] > $eval);
    if(exists $hash_id{$data[1]}){
		if($pre eq " "){
		    push (@arrs,$data[0],$data[1],$data[2],$hash_id{$data[1]});
#print "1: ". join("\t",@arrs) . "\n";
	    }elsif($pre eq $data[0]){
	        push (@arrs,$data[1],$data[2],$hash_id{$data[1]});
		}else{
		    push @arrs,$tophit,$topmatch;
#print "first info: ". join("\t",@arrs) . "\n";
		    &getInfo(@arrs) if(@arrs > 2);
#print join("\t",@arrs)."\n";
			undef @arrs;
			push @arrs,$data[0],$data[1],$data[2],$hash_id{$data[1]};
			#query_name   subject1_name   alignment_length   subject1_annotation_info   subject2_name    alignment_length    subject2_annotation_info   
		}
		$pre = $data[0];
	}
}
push @arrs,$tophit,$topmatch;
if(@arrs == 2){
	exit 0;
}
&getInfo(@arrs);

close INR;
close OUT;

#have defined tophit and topmatch
sub getInfo{
    my (@arr,$query,$hit,$match);
	my $len = scalar @_;
	for my $ii (0 .. $len-1){
	    push @arr,$_[$ii];
	}
#print "sub info: ".join("\t",@arr)."\n";
	$query = shift @arr;
	$match = pop @arr;
	$hit = pop @arr;
	my $now_hit = 0;
    my %count;
    my $end = ($len - 3)/3;
#print join("\t",@arr)."\n";
    for my $i(0 .. $end - 1){
	    unless($arr[$i*3+2] =~ /unknown|unnamed|hypothetical|predicted|putative/i){
		  my $key = "$query,$arr[$i*3]";
		  next if(exists $count{$key} && $count{$key} >= $match);
	    $count{$key}++;
			$count{$query}++;
			last if($count{$query} > $hit);
		   	print OUT "$query\t$arr[$i*3]\t$arr[$i*3+1]\t$arr[$i*3+2]\n";
	    }
    }
	
	$count{$query} = 0 unless(defined $count{$query});
    if($count{$query} == 0){
        print OUT "$query\t$arr[0]\t$arr[1]\t$arr[2]\n";
    }
	undef @arr;
	undef %count;
}
