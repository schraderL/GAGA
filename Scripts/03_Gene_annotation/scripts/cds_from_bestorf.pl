#!perl -w
use strict;

# Author: Bob Zhang
# E-mail: zhangbo@genomics.cn
# create: 2010-03-15
# update: 2010-07-28


die "$0 <BestORF.result> <cufflinks.exons.fa> [minimal_pep_len]\n" if @ARGV < 2 ;

my $min_len = $ARGV[2] || 50;

my ($orfFile, $fasta) = @ARGV;

open (CDS, ">$fasta.cds") 		|| die $!;
open (PEP, ">$fasta.pep") 		|| die $!;
open (ORF, ">$fasta.orf")		|| die $!;

my %seqs = read_fasta_to_hash($fasta);
my %orf = best_orf($orfFile);

for my $id (keys %orf){
	my $annot = '';	
	if (exists $orf{$id}{'noorf'} ){
#		print STDERR "$id\t$orf{$id}{'mrna_len'}bp\t0aa\tNO_ORF\n";
		print ORF "$id\t0\t0\t.\tNO_ORF\n";
	}else{
		$annot .= 'START_ERR;' unless $orf{$id}{'cds'} =~ /^ATG/;
		$annot .= 'STOP_ERR;' unless $orf{$id}{'cds'} =~ /(TAG|TGA|TAA)$/;	
		if ($orf{$id}{'aa_len'} >= $min_len){
			print CDS ">$id $annot\n$orf{$id}{'cds'}\n";
			print PEP ">$id $annot\n$orf{$id}{'pep'}\n";			
		}else{
			$annot .= 'SHORT_ORF;';
		}
#		print STDERR ("$id\t$orf{$id}{'mrna_len'}\t$orf{$id}{'aa_len'}\t$annot\n");	
		print ORF "$id\t$orf{$id}{'start'}\t$orf{$id}{'end'}\t$orf{$id}{'strand'}\t$annot\n";		
	}
}

#----------------------------
#============================
sub best_orf{
	open (FH, shift) || die $!;
	local $/ = "//";
	my %bestorf = ();
	while (<FH>){
		chomp;
		my @tmp = split('>');
		my ($id)  = $tmp[0] =~ /\sSeq name:\s+(\S+)\s/;
		next unless exists $seqs{$id};
		my ($seqlen) = $tmp[0] =~ /Length of sequence:\s+(\d+)/;
		
		next unless $id;
		$bestorf{$id}{'mrna_len'} = $seqlen;		
		
		# NO ORF FOUND	
		if (not exists $tmp[1]){
			$bestorf{$id}{'noorf'} = 1;
			next;
		}

		my ($anno, $pep) = split(/\n/, $tmp[1] , 2);
		$pep =~ s/\W+//g;
		$bestorf{$id}{'pep'}= $pep;
	
		my ($start, $end, $aaLen) = $anno =~ /(\d+)\s+\S\s+(\d+)\s+(\d+)\s+aa/;
		my ($strand) = $anno =~ /chain\s+(\S)/;	
	
		$bestorf{$id}{'start'} = $start;
		$bestorf{$id}{'end'} = $end;
		$bestorf{$id}{'aa_len'} = $aaLen;
		$bestorf{$id}{'strand'} = $strand;
		
		my $cds = substr($seqs{$id}, $start - 1, $end - $start + 1);
		
		if ($strand eq '+'){
			my $stop_codon = 'NNN';
			if($bestorf{$id}{'mrna_len'} - $end >= 3){
				$stop_codon = substr($seqs{$id}, $end, 3);
				$bestorf{$id}{'end'} += 3;
			}
			$cds .= $stop_codon;
			
		}elsif($strand eq '-'){
			my $stop_codon = 'NNN';
			if ($start > 3){
				$stop_codon = substr($seqs{$id}, $start-4, 3) ;
				$bestorf{$id}{'start'} -= 3;
			}
			$cds = reverse( $stop_codon . $cds );
			$cds =~ tr/ACGTacgt/TGCATGCA/;			
		}		

		$bestorf{$id}{'cds'} = $cds;
	}
	close FH;
	return %bestorf;
}
#=============================
sub read_fasta_to_hash{
	my %seqHash = ();
	open (FH , shift) || die $!;
	local $/ = "\n>";
	while (<FH>){
		chomp;
		s/^>//;
		my ($tag , $seq) = split (/\n/, $_, 2);
		my ($id)  = $tag =~ /^(\S+)/; 
		$seq =~ s/\W//g;
		$seqHash{$id} = $seq;
	}
	$/ = "\n"; 
	close FH;
	return %seqHash;
}