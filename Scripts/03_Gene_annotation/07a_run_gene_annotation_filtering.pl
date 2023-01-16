#!/usr/bin/perl
use strict;
use warnings;

#########################################################
#####
#
# Filter repeats from annotation file. Script to run for each single species
#
#####

# Usage  # Also edit lines 35, 36 and 38
# perl run_gene_annotation_filtering.pl input.gff input.proteins.fasta input.genome input.genome.hardmasked input.repeat.gff GAGA-ID

# Usage with real data
# perl run_gene_annotation_filtering.pl /home/projects/ku_00039/people/joeviz/GAGA_annotations/Final_files/GAGA-0014/GAGA-0014.annotation.merge.renamed.fixed.representative.gff3 /home/projects/ku_00039/people/joeviz/GAGA_annotations/Final_files/GAGA-0014/GAGA-0014.annotation.merge.renamed.fixed.representative.pep /home/projects/ku_00039/people/joeviz/GAGA_annotations/Repeat_masked_assemblies/GAGA_annotations/GAGA-0014/GAGA-0014_nextpolish_correct_final_dupsrm_filt.softMasked.fasta /home/projects/ku_00039/people/joeviz/GAGA_annotations/Repeat_masked_assemblies/GAGA_annotations/GAGA-0014/GAGA-0014_nextpolish_correct_final_dupsrm_filt.repeatMasked.fasta /home/projects/ku_00039/people/joeviz/GAGA_annotations/Repeat_masked_assemblies/GAGA_annotations/GAGA-0014/GAGA-0014_nextpolish_correct_final_dupsrm_filt.repeats.gff GAGA-0014

#Load modules in qsub script
#module load tools
#module load ngs
#module load genemark-es/4.62
#module load signalp/4.1c
#module load funannotate/1.8.3

## Input variables

my $gff = $ARGV[0];
my $protfile = $ARGV[1];
my $genomefile = $ARGV[2];
my $genomemaskfile = $ARGV[3];
my $repeatgff = $ARGV[4];
my $gagaid = $ARGV[5];
my $fullgff = $ARGV[6];

my $repeatlib = "/home/projects/ku_00039/people/joeviz/GAGA_annotations/Repeat_filtering/RepeatPeps.dmnd";
my $repeatlibfun = "/home/projects/ku_00039/people/joeviz/GAGA_annotations/Monomorium_pharaonis/funannotate/database/repeats.dmnd";

my $gagaidlist = "/home/projects/ku_00039/people/joeviz/GAGA_species_list.txt";

my $overlap = "0.9"; # Percentage to filter proteins located in repeat annotated regions gff

my $name;

## Start

# Read species list and prefix
my %shortname;
open(Gagafile, "<", $gagaidlist);
while(<Gagafile>){
	chomp;
	my $line = $_;
	next if ($line !~ /\S+/);
	$line =~ s/\n//g; $line =~ s/\r//g;
	my @subl = split (/\t/, $line);
	$shortname{$subl[0]} = $subl[3];
}
close Gagafile;
#

# Running diamond blastp with repeat libraries
print "Running $gagaid\n";
#print "Running diamond blastp with repeat libraries: Already run, uncomment lines to run again\n";
system("diamond blastp --sensitive --query $protfile --threads 4 --out repeats.tsv --db $repeatlib --evalue 1e-10 --max-target-seqs 1 --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen");
system("diamond blastp --sensitive --query $protfile --threads 4 --out repeatsfunannot.tsv --db $repeatlibfun --evalue 1e-10 --max-target-seqs 1 --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen");


my %repeatresults;
open (File, "<", "repeats.tsv");
while (<File>){
	chomp;
	my $line = $_;
	next if ($line !~ /\S+/);
	my @subl = split (/\t/, $line);
	$repeatresults{$subl[0]} = "$line";
}
close File;

my %repeatresultsfun;
open (File, "<", "repeatsfunannot.tsv");
while (<File>){
	chomp;
	my $line = $_;
	next if ($line !~ /\S+/);
	my @subl = split (/\t/, $line);
	$repeatresultsfun{$subl[0]} = "$line";
}
close File;


# Identifying repeat proteins based on overlap with repeat annotated GFF

# Creating first a gff with only protein repeats
system ("grep 'protein' $repeatgff > $gagaid\_repeats_proteins.gff");

# Masking fasta based on repeat proteins only
system ("maskFastaFromBed -fi $genomefile -bed $gagaid\_repeats_proteins.gff -fo $gagaid\_genome.protMasked.fasta");

# Getting both proteins with hardmasked genomes: Name will be GAGA-0014_proteinhardmasked.pep.fasta
system ("perl /home/projects/ku_00039/people/joeviz/programs/bitacora/Scripts/Tools/gff2fasta_v3.pl $genomemaskfile $gff $gagaid\_hardmasked");
system ("perl /home/projects/ku_00039/people/joeviz/programs/bitacora/Scripts/Tools/gff2fasta_v3.pl $gagaid\_genome.protMasked.fasta $gff $gagaid\_proteinhardmasked");


# Read protein fasta and identify proteins in masked regions

my %maskedprots;
my %maskedprotsfasta;
open (File, "<", "$gagaid\_proteinhardmasked.pep.fasta");
while (<File>){
	chomp;
	my $line = $_;
	next if ($line !~ /\S+/);
	if ($line =~ />(\S+)/){
		$name = $1;
	} else {
		$maskedprotsfasta{$name} .= $line;
	}
}
close File;

foreach my $id (sort keys %maskedprotsfasta){
	my $seq = $maskedprotsfasta{$id};
	my $nnumber += $seq =~ s/X/X/gi;
	my $seqlen = length($seq);
	my $nperc = $nnumber/$seqlen;
	$maskedprots{$id} = $nperc;
}

# File full masked
my %fullmaskedprots;
my %fullmaskedprotsfasta;
open (File, "<", "$gagaid\_hardmasked.pep.fasta");
while (<File>){
	chomp;
	my $line = $_;
	next if ($line !~ /\S+/);
	if ($line =~ />(\S+)/){
		$name = $1;
	} else {
		$fullmaskedprotsfasta{$name} .= $line;
	}
}
close File;

foreach my $id (sort keys %fullmaskedprotsfasta){
	my $seq = $fullmaskedprotsfasta{$id};
	my $nnumber += $seq =~ s/X/X/gi;
	my $seqlen = length($seq);
	my $nperc = $nnumber/$seqlen;
	$fullmaskedprots{$id} = $nperc;
}


# Print results

my $filteredprots = "";
my $lenfilt = 0;
my $repfilt = 0;
my $rarefilt = 0;

open (Resultsfull, ">", "$gagaid\_proteins_repeats_fulltable.tsv");
print Resultsfull "Gene id\tProtein length\tPercent masked with repeat GFF (0-1)\tPercent masked with only protein repeats GFF(0-1)\tRepeat annotation\tRepeat annotation funannotate db\n";
open (Resultsfilt, ">", "$gagaid\_proteins_repeats_filtered_table.tsv");
print Resultsfilt "Gene id\tProtein length\tPercent masked with repeat GFF(0-1)\tPercent masked with only protein repeats GFF(0-1)\tRepeat annotation\tRepeat annotation funannotate db\n";
open (Resultsnofilt, ">", "$gagaid\_proteins_repeats_retained_table.tsv");
print Resultsnofilt "Gene id\tProtein length\tPercent masked with repeat GFF(0-1)\tPercent masked with only protein repeats GFF(0-1)\tRepeat annotation\tRepeat annotation funannotate db\n";


foreach my $id (sort keys %maskedprotsfasta){

	my $len = length($maskedprotsfasta{$id});
	my $seq = $maskedprotsfasta{$id};

	my $repres = "-NA-";
	if (exists $repeatresults{$id}){
		$repres = $repeatresults{$id};
	}

	my $represfun = "-NA-";
	if (exists $repeatresultsfun{$id}){
		$represfun = $repeatresultsfun{$id};
	}

	print Resultsfull "$id\t$len\t$fullmaskedprots{$id}\t$maskedprots{$id}\t$repres\t$represfun\n";

	# Filter used

	if ($seq =~ /TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT/ || $seq =~ /SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS/ || $seq =~ /DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD/){
		print Resultsfilt "$id\t$len\t$fullmaskedprots{$id}\t$maskedprots{$id}\t$repres\t$represfun\n";
		$filteredprots .= "$id\_";
		$rarefilt++;
#		print ">$id\n$seq\n";

#	} elsif ($len < 50) { # Short length filter: Not used
#		print Resultsfilt "$id\t$len\t$fullmaskedprots{$id}\t$maskedprots{$id}\t$repres\t$represfun\n";
#		$filteredprots .= "$id\_";
#		$lenfilt++;

	} elsif ($repres !~ /-NA-/){
		my @subline = split (/\t/, $repres);

		my $filtro1 = ($subline[12]*2)/3; # Alignment is 2/3 of query protein 
		my $filtro2 = ($subline[13]*0.8); # Alignment is 80% of subject protein

		my $filtro3 = ($subline[12]*1)/2; # Alignment is half of query protein 

		if ($subline[3] >= $filtro1 || $subline[3] >= $filtro2) { # BLAST filtering: alignment covering 2/3 of subject, or 80% of query
			print Resultsfilt "$id\t$len\t$fullmaskedprots{$id}\t$maskedprots{$id}\t$repres\t$represfun\n";		
			$filteredprots .= "$id\_";
			$repfilt++;
		} elsif ($subline[3] >= $filtro3 && $maskedprots{$id} >= $overlap) { # BLAST filtering: alignment covering 50% of query and more than 90% of protein masked from protein repeat GFF 
			print Resultsfilt "$id\t$len\t$fullmaskedprots{$id}\t$maskedprots{$id}\t$repres\t$represfun\n";		
			$filteredprots .= "$id\_";
			$repfilt++;
#		} elsif ($subline[3] >= 10) { # BLAST filtering: All hits (alignment longer than 100 which is all)
#			print Resultsfilt "$id\t$len\t$fullmaskedprots{$id}\t$maskedprots{$id}\t$repres\t$represfun\n";		
#			$filteredprots .= "$id\_";
#			$repfilt++;	
#		}  elsif ($maskedprots{$id} >= $overlap) { # Filtering also genes with 90% of protein masked from protein repeat GFF 
#			print Resultsfilt "$id\t$len\t$fullmaskedprots{$id}\t$maskedprots{$id}\t$repres\t$represfun\n";		
#			$filteredprots .= "$id\_";
#			$repfilt++;
		} else {
			print Resultsnofilt "$id\t$len\t$fullmaskedprots{$id}\t$maskedprots{$id}\t$repres\t$represfun\n";				
		}

	} else {
		print Resultsnofilt "$id\t$len\t$fullmaskedprots{$id}\t$maskedprots{$id}\t$repres\t$represfun\n";		
	}


}
print "\n$lenfilt proteins were filtered as shorter than 50aa, $repfilt were filtered as repeat proteins, and $rarefilt were filtered as low complexity repeat proteins\n\n";

close Resultsnofilt;
close Resultsfull;
close Resultsfilt;


# Filter proteins from GFF and get final fasta

open (FileGFF, "<", "$gff");
open (ResultsGFF, ">", "$gagaid\_final_annotation_representative_repfilt.gff3");
print ResultsGFF "##gff-version 3\n";

while (<FileGFF>){
	chomp;
	my $line = $_;
	next if ($line !~ /\S+/);	
	if ($line =~ /^#/){
		#print ResultsGFF "$line\n";
		next;
	} elsif ($line !~ /\t/) {
		print ResultsGFF "$line\n";
	} else {
		my @subl = split (/\t/, $line);	
		if ($subl[2] =~ /mRNA/){
			my $genename = "";
			my $mrnaname = "";
			if ($subl[8] =~ /ID=([^;]+)/){
				$mrnaname = $1;
			} else {die "ERROR: It does not recognize ID in the GFF3 in: $line\n";}

			if ($mrnaname =~ /(\S+)_i\d+/){
				$genename = $1;
			} else {die "ERROR: It does not recognize ID isoforms in the GFF3 in: $line\n";}

			next if ($filteredprots =~ /$mrnaname\_/); # Exclude filtered proteins
			print ResultsGFF "$subl[0]\t$subl[1]\tgene\t$subl[3]\t$subl[4]\t\.\t$subl[6]\t$subl[7]\tID=$genename\;\n";
			print ResultsGFF "$line\n";
		} else {
			my $genename = "";
			if ($subl[8] =~ /Parent=([^;]+)/){
				$genename = $1;
			} else {die "ERROR: It does not recognize Parent ID in the GFF3 in: $line\n";}
			next if ($filteredprots =~ /$genename\_/); # Exclude filtered proteins
			print ResultsGFF "$line\n";			
		}
	}
}
close FileGFF;
close ResultsGFF;

system ("perl /home/projects/ku_00039/people/joeviz/programs/bitacora/Scripts/Tools/gff2fasta_v3.pl $genomefile $gagaid\_final_annotation_representative_repfilt.gff3 $gagaid\_final_annotation_representative_repfilt");
system ("sed \'s/X\$//\' $gagaid\_final_annotation_representative_repfilt.pep.fasta > $gagaid\_final_annotation_representative_repfilt.pep.fasta.tmp");
system ("mv $gagaid\_final_annotation_representative_repfilt.pep.fasta.tmp $gagaid\_final_annotation_representative_repfilt.pep.fasta");


# Filter proteins from GFF with isoforms and get final fasta

open (FileGFF, "<", "$fullgff");
open (ResultsGFF, ">", "$gagaid\_final_annotation_repfilt.gff3");
print ResultsGFF "##gff-version 3\n";

while (<FileGFF>){
	chomp;
	my $line = $_;
	next if ($line !~ /\S+/);	
	if ($line =~ /^#/){
		#print ResultsGFF "$line\n";
		next;
	} elsif ($line !~ /\t/) {
		print ResultsGFF "$line\n";
	} else {
		my @subl = split (/\t/, $line);	
		if ($subl[2] =~ /mRNA/){
			my $genename = "";
			my $mrnaname = "";
			if ($subl[8] =~ /ID=([^;]+)/){
				$mrnaname = $1;
			} else {die "ERROR: It does not recognize ID in the GFF3 in: $line\n";}

			if ($mrnaname =~ /(\S+)_i\d+/){
				$genename = $1;
			} else {die "ERROR: It does not recognize ID isoforms in the GFF3 in: $line\n";}

			next if ($filteredprots =~ /$genename\_/); # Exclude filtered proteins
			print ResultsGFF "$line\n";
		} elsif ($subl[2] =~ /gene/) {
			my $genename = "";
			if ($subl[8] =~ /ID=([^;]+)/){
				$genename = $1;
			} else {die "ERROR: It does not recognize Parent ID in the GFF3 in: $line\n";}
			next if ($filteredprots =~ /$genename\_/); # Exclude filtered proteins
			print ResultsGFF "$line\n";			
		} else {
			my $genename = "";
			my $mrnaname = "";
			if ($subl[8] =~ /Parent=([^;]+)/){
				$mrnaname = $1;
			} else {die "ERROR: It does not recognize ID in the GFF3 in: $line\n";}

			if ($mrnaname =~ /(\S+)_i\d+/){
				$genename = $1;
			} else {die "ERROR: It does not recognize ID isoforms in the GFF3 in: $line\n";}

			next if ($filteredprots =~ /$genename\_/); # Exclude filtered proteins
			print ResultsGFF "$line\n";			
		}
	}
}
close FileGFF;
close ResultsGFF;

system ("perl /home/projects/ku_00039/people/joeviz/programs/bitacora/Scripts/Tools/gff2fasta_v3.pl $genomefile $gagaid\_final_annotation_repfilt.gff3 $gagaid\_final_annotation_repfilt");
system ("sed \'s/X\$//\' $gagaid\_final_annotation_repfilt.pep.fasta > $gagaid\_final_annotation_repfilt.pep.fasta.tmp");
system ("mv $gagaid\_final_annotation_repfilt.pep.fasta.tmp $gagaid\_final_annotation_repfilt.pep.fasta");

