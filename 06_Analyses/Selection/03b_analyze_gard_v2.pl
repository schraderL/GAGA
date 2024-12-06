#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
my $dirname = dirname(__FILE__);

# Hyphy scripts creates partitions with disrupted codons (no 3 multiple), so my script (v2) creates correct the partitions

# usage: perl analyze_gard_v2.pl

## Change line 70 with the gard output (._gard.output) to remove the dot
## Commented line to run hyphy and get partitions, included now in the script that runs gard

## Load computerome2 modules
=h
module load ngs tools
module load anaconda3/4.4.0
module load fasttree/2.1.11
module load perl/5.24.0
module load openmpi/gcc
module load hyphy/2.5.29
=cut

#my $cdsdir = "/home/projects/ku_00039/people/joeviz/Suz/ortholog_alignments/all/codon_alignments/";
#my $garddir = "/home/projects/ku_00039/people/joeviz/Suz/ortholog_alignments/all/Hyphy/GARD_nozorromasked_normal/gard/";
#my $cdsdir = "/home/projects/ku_00039/people/joeviz/Suz/ortholog_alignments/all_2batch/codon_alignments/";
#my $garddir = "/home/projects/ku_00039/people/joeviz/Suz/ortholog_alignments/all_2batch/Hyphy/GARD/gard/";

my $cdsdir = "/home/projects/ku_00039/people/joeviz/orthology_alignments/all_orthogroups/codon_alignments/codon_alignments_qualfiltered/"; # *cds.nonstop.aln 
my $treedir = "/home/projects/ku_00039/people/joeviz/orthology_trees/all_orthologs/codon_dna_trees/codon_dna_aln_trees/";
my $garddir = "/home/projects/ku_00039/people/joeviz/orthology_trees/all_orthologs/Hyphy_gard/gard/";

my $outdir = "partition_files_v2";
system ("mkdir -p $outdir");

my $aicdif = "100"; # Minimum AIC difference to consider multiple partitions

my ($line, $name);

open (Scores, ">", "GARD_summary_v2.txt");
print Scores "Orthogroup\tRun status\tBest model\tc-AIC best model\tSingle partition model\tc-AIC single partition\tAll GARD models\n";

# Reading Protein fasta

#my %protfasta;
system ("ls $cdsdir\/\*cds.nonstop.aln > tmp_protlist.txt");
#system ("ls $cdsdir\/\*\.treefile > tmp_protlist.txt");

open(File, "<", "tmp_protlist.txt");
while(<File>){
	chomp;
	$line = $_;
	next if ($line !~ /\S+/);	

	my $id = "";
	if ($line =~ /.*\/(\S+)\.cds/){
#	if ($line =~ /.*\/(\S+)\.treefile/){		
		$id = $1;
	} else {
		die "Can't find OG id in $line\n";
	}  

	my $cdsfile = "$line";
	my $treefile = "$treedir\/$id\.treefile";

=h
	my %fasta;
	open(Filef, "<", $line);
	while(<Filef>){
		chomp;
		my $line2 = $_;
		if ($line2 =~ />(\S+)/){
			$name = $1;
		} else {
			$fasta{$name} .= "$line2";
		}
	}
	close Filef;
=cut 

	# Now mask the alignment based on zorro scores

	my $gardfile = "$garddir\/$id\_gard.output";
	my $gardnex = "$garddir\/$id\.best-gard";

	my $status = "None"; # If it does not enter any if (No best model, killed or error, it will be none in the output table)
	my $numpart = 0; my $bestaic = "None"; my $bestaicmodel = "None";
	my $basemodel = "None"; my $basemodelaic = "None";
	my $linecontrol = "0"; my $prevmodel = ""; my $allmodels = "";
	open (Gardout, "<", $gardfile);
	while (<Gardout>){
		chomp;
		my $line2 = $_;
		next if ($line2 !~ /\S+/);
		if ($line2 =~ /estimated parameters/){ # Base model
			if ($line2 =~ /AIC\-c\s+\=\s+(\S+)/){
				$basemodelaic = $1;
				$basemodel = $line2;
				$bestaic = $basemodelaic;
				$bestaicmodel = "single-partition; no breakpoints";
			} else {
				die "It can't read the c-AIC in $line2\nFile $line\n";
			}
			$status = "OK";
		}
		elsif ($line2 =~ /Best sinlge break point location:/){ # 1 break model
			$linecontrol = 1;
			$allmodels .= "$line2";
			$prevmodel = "$line2";
		}		
		elsif ($line2 =~ /Best break point locations:/){ # multiple break model
			$linecontrol = 2;
			$allmodels .= "$line2";
			$prevmodel = "$line2";
		}	
		elsif ($line2 =~ /The alignment is too short to permit c-AIC based model compariso/){ # Error running: alignment too short
			$status = "Error short aln";
		}
		elsif ($line2 =~ /killed/){ # Error running: killed
			$status = "Killed";
		}

		if ($linecontrol > 0){ # 1 or multiple break model
			if ($line2 =~ /c\-AIC\s+\=\s+(\S+)/){
				my $aic = $1;
				$allmodels .= "$line2";
				my $aicdifference = $bestaic - $aic;
				if ($aicdifference > $aicdif){
					$bestaic = $aic;
					$bestaicmodel = $prevmodel;		
					$numpart++;			
				}
			}
		}



	}	
	close Gardout;

#print Scores "Orthogroup\tRun satus\tBest model\tc-AIC best model\tSingle partition model\tc-AIC single partition\n";
	print Scores "$id\t$status\t$bestaicmodel\t$bestaic\t$basemodel\t$basemodelaic\t$allmodels\n";

# Comment not to move the partition files
#=h
	if ($numpart == 0 || $bestaic =~ /None/){
		# Single partition is best
		system ("cp $cdsfile $outdir\/");
		system ("cp $treefile $outdir\/");

	} else {

		# Get partitions
		my @partitions = ();
		if ($bestaicmodel =~ /Best sinlge break point location: (\d+)/){
			push (@partitions, $1);
		} elsif ($bestaicmodel =~ /Best break point locations: (.*)/){
			my $partnum = $1;
			$partnum =~ s/\,//g;
			@partitions = split (/ /, $partnum);

		} else {
			print "Error in $id cannot find the number of partitions to split the alignment\n";
			system ("cp $cdsfile $outdir\/");
			system ("cp $treefile $outdir\/");
			next;			
		}

		# Read fasta of the alignments
		my %fastaln;
		open(Filef, "<", $line);
		while(<Filef>){
			chomp;
			my $line2 = $_;
			if ($line2 =~ />(\S+)/){
				$name = $1;
			} else {
				$fastaln{$name} .= "$line2";
			}
		}
		close Filef;

		my $partnum = scalar (@partitions);
		my $num = "1";
		my $prevpart = "0";
		foreach my $part (sort { $a <=> $b } @partitions){

			open (Resultspart, ">", "$outdir/$id\_parts_SPAN_$num\.cds.nonstop.aln");
			my $totalseqs = "0";

			# Check that the partition is multiple of three
			my $partcheck = ($part-$prevpart)/3;
			if ($partcheck =~ /\.33/){
			        $part = $part + 2;
			} elsif ($partcheck =~ /\.66/) {
			        $part = $part +1;
			} else {
			        $part = $part;
			}
			my $partlength = $part-$prevpart;

			# Now split fasta and print partition
			foreach my $key (keys %fastaln){
				my $seq = $fastaln{$key};
				my $partseq = uc (substr ($seq, $prevpart, $partlength));
				my $len = length ($partseq);	
				my $lengthcheck = length($partseq)/3;
				if ($lengthcheck =~ /\.33/ || $lengthcheck =~ /\.66/ ){
					die "Partition $id\_parts_SPAN_$num\.cds.nonstop.fullpart.aln length $len is not multple of three! lengthcheck value is $lengthcheck\n";
				}

				my $countgap = () = $partseq =~ /\Q-/g;
				my $countn = () = $partseq =~ /\QN/g;
				my $diflen = $len - $countgap - $countn;	
				if ($diflen >= 24){
					print Resultspart ">$key\n$partseq\n";
					$totalseqs++;
				} else {
					print "$key has lower than 24 non-ambiguous nucleotides, discarded from alignment $outdir/$id\_parts_SPAN_$num\.cds.nonstop.aln\n";
				}	

			}	

			close Resultspart;

			if ($totalseqs < 17){ 
				print "Partition $outdir/$id\_parts_SPAN_$num\.cds.nonstop.aln contains lower than 17 sequences, discarding partition\n";
				system ("rm $outdir/$id\_parts_SPAN_$num\.cds.nonstop.aln");
			}

			$num++;
			$prevpart = $part;

		}

		# Print now here the final partition
		open (Resultspart, ">", "$outdir/$id\_parts_SPAN_$num\.cds.nonstop.aln");
		my $totalseqs = "0";

		# Check that the partition is multiple of three
		my $partcheck = ($prevpart)/3;
		if ($partcheck =~ /\.33/){
			$prevpart = $prevpart + 2;
		} elsif ($partcheck =~ /\.66/) {
			$prevpart = $prevpart +1;
		} else {
		    $prevpart = $prevpart;
		}

		# Now split fasta and print partition
		foreach my $key (keys %fastaln){
			my $seq = $fastaln{$key};
			my $partseq = uc (substr ($seq, $prevpart)); # It will retrieve everything through the end
			my $len = length ($partseq);	
			my $lengthcheck = length($partseq)/3;
			if ($lengthcheck =~ /\.33/ || $lengthcheck =~ /\.66/ ){
				die "Partition $id\_parts_SPAN_$num\.cds.nonstop.fullpart.aln length $len is not multple of three! lengthcheck value is $lengthcheck\n";
			}

			my $countgap = () = $partseq =~ /\Q-/g;
			my $countn = () = $partseq =~ /\QN/g;
			my $diflen = $len - $countgap - $countn;	
			if ($diflen >= 24){
				print Resultspart ">$key\n$partseq\n";
				$totalseqs++;
			} else {
				print "$key has lower than 24 non-ambiguous nucleotides, discarded from alignment $outdir/$id\_parts_SPAN_$num\.cds.nonstop.aln\n";
			}	

		}	

		close Resultspart;

		if ($totalseqs < 17){ 
			print "Partition $outdir/$id\_parts_SPAN_$num\.cds.nonstop.aln contains lower than 17 sequences, discarding partition\n";
			system ("rm $outdir/$id\_parts_SPAN_$num\.cds.nonstop.aln");
		}



=h
		# Old version
		# Multiple partition from GARD, generate the fasta files
		system ("hyphy /home/projects/ku_00039/people/joeviz/programs/hyphy-analyses/extract-partitions/extract-partitions.bf --msa $gardnex --output $outdir/$id\_parts  --extension cds.nonstop.fullpart.aln ENV=\"DATA_FILE_PRINT_FORMAT=9\"");

		# Add line to remove the last line in the partition fasta that contains a tree

		system ("ls $outdir/$id\_parts\*cds.nonstop.fullpart.aln > newpartitionfiles.txt");
		open (Partfiles, "<", "newpartitionfiles.txt");
		while (<Partfiles>){
			chomp;
			my $line3 = $_;
			next if ($line3 !~ /\S+/);

			my $partid = "";
			if ($line3 =~ /.*\/(\S+)\.cds/){
				$partid = $1;
			} else {
				die "Can't find OG id in $line3\n";
			}  

			system ("tail -1 $line3 > $outdir\/$partid\.treefile");
			system ("head -n -1 $line3 > temp.txt ; mv temp.txt $line3");


			### Exclude sequences with only gaps or missing data
			my %fasta;
			open(Filef, "<", $line3);
			while(<Filef>){
				chomp;
				my $line4 = $_;
				if ($line4 =~ />(\S+)/){
					$name = $1;
				} else {
					$fasta{$name} .= "$line4";
				}
			}
			close Filef;

			open (Resultsf, ">", "$outdir/$partid\.cds.nonstop.aln");
			foreach my $key (keys %fasta){
				my $seq = uc ($fasta{$key});
				my $len = length ($fasta{$key});
				my $countgap = () = $seq =~ /\Q-/g;
				my $countn = () = $seq =~ /\QN/g;
				my $diflen = $len - $countgap - $countn;	
				if ($diflen >= 12){
					print Resultsf ">$key\n$seq\n";
				} else {
					print "$key has lower than 12 non-ambiguous nucleotides, discarded from alignment $line3\n";
				}	
			}
			close Resultsf;

		}
		close Partfiles;
		system ("rm newpartitionfiles.txt");
=cut 

#		system ("python3 /home/projects/ku_00039/people/joeviz/Suz/ortholog_alignments/msa_mask_from_zorro.py $cdsdir\/$id\.cds.nonstop.aln $zorrofile $zorroscore $cdsdir\/$id\.cds.nonstop.zorromasked.aln ");
		#print ("python3 /home/projects/ku_00039/people/joeviz/Suz/ortholog_alignments/msa_mask_from_zorro.py $cdsdir\/$id\.cds.nonstop.aln $zorrofile 5 $cdsdir\/$id\.cds.nonstop.zorromasked.aln \n");

	}
#=cut

}
close File;
system ("rm tmp_protlist.txt");





