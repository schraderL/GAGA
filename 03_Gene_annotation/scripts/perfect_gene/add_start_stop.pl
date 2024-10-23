#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;
use File::Basename qw(dirname basename);

if ( @ARGV < 2){
	print "Usage:\n\tperl $0 <genome.fa> <gff> [add_start(amino number)] [add_stop (amino number)]>output_good.gff\n";
	exit;
}

my ($genome_file, $gene_file,$add_start,$add_stop ) = @ARGV;
my (%genome, %gene);

#stop = 'TAA' => 'U', 'TAG' => 'U', 'TGA' => 'U'
#star = 'ATG'

$add_start ||= 0;
$add_stop ||= 0;

open (IN, $gene_file) || die $!;
my $id;
while (<IN>){
	chomp;
	next if (/^#/);
	my @c = split /\s+/, $_;
	if ($c[2] eq 'mRNA'){
		$id = $1 if ($c[8] =~ /ID=(\S+?);/);
		next if ($c[8] =~ /Shift=(\d+)/ && $1 > 0);
	}else{
		push @{$gene{$c[0]}{$id}}, [@c];
	}
}
close IN;

open (IN, $genome_file) || die $!;
$/ = ">"; <IN>; $/ = "\n";
while (<IN>){
	chomp;
	my $chr = $1 if (/^(\S+)/);
	$/=">";
	my $seq = <IN>;
	chomp $seq;
	$/ ="\n";

	next unless ($gene{$chr});
	$seq = uc($seq);
	$seq =~ s/\s//g;

	foreach my $id (keys %{$gene{$chr}}){
 		my @cds = sort {$a->[3] <=> $b->[3]} @{$gene{$chr}{$id}};

		if(extend_start($seq,\@cds) && extend_stop($seq,\@cds)){
			my $out="$chr\t$cds[0][1]\tmRNA\t$cds[0][3]\t$cds[-1][4]\t100\t$cds[0][6]\t.\tID=$id;Shift=0\n";
			foreach my $p (@cds){
				$out .= join "\t", @$p;
				$out .= "\n";
	 		}
			print $out;
		}
	}
}
close IN;

sub RC{
	my $seq = shift;
	$seq = reverse $seq;
	$seq =~ tr/ATCG/TAGC/;
	return $seq;
}

sub extend_start{
	my ($seq,$cds)=@_;
	for(my $i=0;$i<=$add_start;$i++){
		if($$cds[0][6] eq '+'){
			return 0 if ($$cds[0][3]-1-$i*3 < 0 );
			my $star_codon = substr($seq,$$cds[0][3]-1-$i*3,3);
			return 0 if($star_codon=~/N+/);
			return 0 if ($star_codon eq 'TAA' || $star_codon eq 'TAG' || $star_codon eq 'TGA');
			if($star_codon eq 'ATG'){
				$$cds[0][3]-=$i*3;
				return 1;
			}
		}else{
			my $star_codon = &RC(substr($seq,$$cds[-1][4]-3+$i*3,3));
			return 0 if(length($star_codon) != 3 || $star_codon=~/N+/);
			return 0 if ($star_codon eq 'TAA' || $star_codon eq 'TAG' || $star_codon eq 'TGA');
			if($star_codon eq 'ATG'){
				$$cds[-1][4]+=$i*3;
				return 1;
			}
		}
	}
	return 0;
}


sub extend_stop{
	my ($seq,$cds)=@_;
	for(my $i=0;$i<=$add_stop;$i++){
		if($$cds[0][6] eq '+'){
			my $stop_codon=substr($seq, $$cds[-1][4]-3+$i*3, 3);
			return 0 if(length($stop_codon) != 3 || $stop_codon =~ /N+/);
			if($stop_codon eq 'TAA' || $stop_codon eq 'TAG' || $stop_codon eq 'TGA'){
				$$cds[-1][4] += $i*3;
				return 1;
			}
		}else{
			return 0 if ($$cds[0][3]-1-$i*3 < 0 );
			my $stop_codon=&RC(substr($seq,$$cds[0][3]-1-$i*3,3));
			return 0 if ($stop_codon =~ /N+/);
			if($stop_codon eq 'TAA' || $stop_codon eq 'TAG' || $stop_codon eq 'TGA'){
				$$cds[0][3] -= $i*3;
				return 1;
			}
		}
	}
	return 0;
}
