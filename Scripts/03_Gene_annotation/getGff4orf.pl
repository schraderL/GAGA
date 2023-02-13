#!/usr/bin/perl
use strict;
use Data::Dumper;

die "perl $0 <orf> <gtf> <seq>  > outfile \n" if @ARGV < 3;

my $Orf=shift;
my $Gtf=shift;
my $seq=shift;

my %seq;
readSeq($seq,\%seq);

#print Dumper \%seq;
#exit;

my %orf;
open IN,$Orf or die $!;
#RNASEQ3.2488.1  2       3217    +       START_ERR;
while(<IN>){
	my ($id,$start,$end,$strand,$error)=(split /\s+/)[0,1,2,3,4];
	next if($start == 0 || $end ==0);
	@{$orf{$id}}=($start,$end,$strand,$error);
}
close IN;

my %gtf;
open IN,$Gtf or die $!;
my $id_tmp = '';
while(<IN>){
	next if /^#/;
	chomp;
	my @a=split /\t+/;
	my $id;
	@a[3,4]=@a[4,3] if($a[3]>$a[4]);
	my $gene_name = $1 if ($a[8] =~ /gene_name \"(\S+)\"/);
#	if($a[2] eq "transcript" && $a[8]=~/transcript_id\s+\"([^\s\"]+)\"/){
#		$id=$1;
#		$a[2]="mRNA";
#		$a[8]="ID=$id;";
#		@{$gtf{$id}{mRNA}}=@a;
	if($a[2] eq "exon" && $a[8]=~/transcript_id\s+\"([^\s\"]+)\"/){
		$id=$1;
		$a[2]='CDS';
		if (defined $gene_name) {
			$a[8]="Parent=$id;gene_name=$gene_name;";
		}else {
			$a[8]="Parent=$id;";
		}
		push @{$gtf{$id}{CDS}},[@a];
		if ($id_tmp ne $id) {
			if (defined $gene_name) {
				$a[8]="ID=$id;gene_name=$gene_name;";
			}else {
				$a[8]="ID=$id;";		
			}
			@{$gtf{$id}{mRNA}} = @a;
		}
		$id_tmp = $id;
	}

}
close IN;

for my $id(keys %orf){
	next unless(exists $gtf{$id});
	my $p_id=$orf{$id};
	my $q_id=$gtf{$id};
	next if ($p_id->[2] eq '.' && $q_id->{mRNA}->[2] eq '.');
	if ($q_id->{mRNA}->[2] eq '.') {
		$q_id->{mRNA}->[2] = $p_id->[2];
	}elsif($p_id->[2] ne $q_id->{mRNA}->[6]) {
		next;
	}
	my @a=sort {$a->[3] <=> $b->[3]} @{$q_id->{CDS}};
	my @b=@{$q_id->{mRNA}};
	my $cds_len=$p_id->[1]-$p_id->[0]+1;

	my @out;
	my $current_len=0;
	my $len=0;
	my $flog=0;
	for(my $i=0;$i<@a;$i++){
		$len+=$a[$i]->[4]-$a[$i]->[3]+1;
		next if($len < $p_id->[0]);
		
		$a[$i]->[3]=$a[$i]->[4]-($len-$p_id->[0]) if($flog==0);
		$flog++;
		$current_len+=$a[$i]->[4]-$a[$i]->[3]+1;
		if($current_len < $cds_len){
			push @out,[@{$a[$i]}];
		}elsif ($current_len == $cds_len){
			push @out,[@{$a[$i]}];
			last;
		}else {
			$a[$i]->[4]=$a[$i]->[4]-($len-$p_id->[1]);
			push @out,[@{$a[$i]}];
			last;
		}
	}
	$b[2] = "mRNA";
	$b[3] = $out[0]->[3];
	$b[4] = $out[-1]->[4];
	print join("\t",@b)."\n" if ($b[3] ne '');

	if (!defined $p_id->[3]) {
		my @utr_5 = @{$out[0]};
		my @utr_3 = @{$out[-1]};
		my $strand = $p_id->[2];
		start_stop(\@utr_5,\@utr_3,$strand,1,1);
	}elsif (defined $p_id->[3]) {
		my @utr_5 = @{$out[0]};
		my @utr_3 = @{$out[-1]};
		my $strand = $p_id->[2];
		my $error = $p_id->[3];
		start_stop_revise(\@utr_5,\@utr_3,$strand,$error);
	}

	for(@out){
		print join("\t",@$_)."\n" if ($b[3] ne '');
	}
}



################
sub readSeq {
	my ($file,$p) = @_;

	open IN, "<$file" or die "failed to open $file: $!\n";
	$/=">"; <IN>; $/="\n";
	while (<IN>) {
		chomp;
		my $id = $1 if (/^(\S+)/);
		$/=">";
		my $seq = <IN>;
		$/="\n";
		$seq=~s/>$//;
		$seq=~s/\n//g;
		$p->{$id} = $seq;
	}
	close (IN);
}

sub start_stop {
	my ($utr5_p,$utr3_p,$strand,$start,$stop) = @_;

	if ($start) {
		my @out = @{$utr5_p};
		$out[2] = 'start_codon';
		$out[3] = $utr5_p->[3];
		$out[4] = $utr5_p->[3]+2;
		my $start_codon = substr($seq{$utr5_p->[0]},$out[3]-1,$out[4]-$out[3]+1);

		if ($strand eq '-') {
			$out[4] = $utr3_p->[4];
			$out[3] = $utr3_p->[4]-2;
			$start_codon = substr($seq{$utr5_p->[0]},$out[3]-1,$out[4]-$out[3]+1);
			$start_codon = reverse $start_codon;
			$start_codon =~ tr/ACGTacgt/TGCATGCA/;
		}

		if ($start_codon eq 'ATG') {
			print join("\t",@out),"\n";
		}
	}
	if ($stop) {
		my @out = @{$utr3_p};

		$out[2] = 'stop_codon';
		$out[4] = $utr3_p->[4];
		$out[3] = $utr3_p->[4]-2;
		my $stop_codon = substr($seq{$utr3_p->[0]},$out[3]-1,$out[4]-$out[3]+1);
		if ($strand eq '-') {
			$out[4] = $utr5_p->[3]+2;
			$out[3] = $utr5_p->[3];
			$stop_codon = substr($seq{$utr3_p->[0]},$out[3]-1,$out[4]-$out[3]+1);
			$stop_codon = reverse $stop_codon;
			$stop_codon =~ tr/ACGTacgt/TGCATGCA/;
		}

		if ($stop_codon =~ /[TAA|TAG|TGA]/) {
			print join("\t",@out),"\n";
		}
	}
}


sub start_stop_revise {
	my ($utr5_p,$utr3_p,$strand,$error) = @_;

	if (!($error =~ /STOP_ERR/)) {
#		print STDERR "matched !\n";
		start_stop($utr5_p,$utr3_p,$strand,0,1);
	}elsif (!($error =~ /START_ERR/)) {
#		print STDERR "matched !\n";
		start_stop($utr5_p,$utr3_p,$strand,1,0);
	}
}
