#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
my $dirname = dirname(__FILE__);

# Script to generate final files for a single genome, edit the path to specify each file for that genome

# usage: perl get_all_functional_annotation.pl (change lines where the input files are called)

my ($line, $name, $nameout);

#### Specify here the path to search for the input files, with GAGA-ID name as prefix

# Protein fasta containing the annotated genes that was used as input in the functional annotation, here one representative protein sequence per gene
my $annotdir = "/home/projects/ku_00039/people/joeviz/GAGA_genomes/Ectatomma/Final_gene_annotation/Ectatomma_ruidum_final_annotation_repfilt_addfunc.representative.pep.fasta";
# GO file generated in the previous step of the pipeline described in github
my $godir = "/home/projects/ku_00039/people/joeviz/GAGA_genomes/Ectatomma/functional_annotation/zijun/Ectatomma_ruidum_final_annotation_repfilt_addfunc.representative.pep.fasta.iprscan.gene.wego";
# Directory containing all functional annotation files from the previous step described in github (all_function_stat.pl)
my $funcannotdir = "/home/projects/ku_00039/people/joeviz/GAGA_genomes/Ectatomma/functional_annotation/zijun/"; # Directory with functional annotations from Zijun
# Annotation.xls file generated in the previous step (all_function_stat.pl)
my $funcannotdirupdated = "/home/projects/ku_00039/people/joeviz/GAGA_genomes/Ectatomma/functional_annotation/zijun/annotation.xls";
# Eggnog output file
my $eggnogdir = "/home/projects/ku_00039/people/joeviz/GAGA_genomes/Ectatomma/functional_annotation/Ectatomma_ruidum_final_annotation_repfilt_addfunc.representative/out_Ectatomma_ruidum_final_annotation_repfilt_addfunc.representative.emapper.annotations";
# Interpro scan output file
my $interprodir = "/home/projects/ku_00039/people/joeviz/GAGA_genomes/Ectatomma/functional_annotation/interpro/Ectatomma_ruidum_final_annotation_repfilt_addfunc.representative/Ectatomma_ruidum_final_annotation_repfilt_addfunc.representative.pep.fasta.tsv";
# Interpro scan output file parsed in the previous step of the pipeline (all_function_stat.pl)
my $interprozijun = "/home/projects/ku_00039/people/joeviz/GAGA_genomes/Ectatomma/functional_annotation/zijun/Ectatomma_ruidum_final_annotation_repfilt_addfunc.representative.pep.fasta.iprscan.xls";

#### OLD input to do all GAGA genome in the same script
#my $annotdir = "/home/projects/ku_00039/people/joeviz/GAGA_annotations/Final_GAGA_annotations/";
# File that will search (i.e: GAGA-0001): GAGA-0001/GAGA-0001_final_annotation_repfilt_addreannot_noparpse_representative.pep.fasta
#my $godir = "/home/projects/ku_00039/people/joeviz/GAGA_annotations/Functional_annotation/GO_annotations/";
# File: GAGA-0001_gene_merged.wego
#my $funcannotdir = "/home/projects/ku_00039/people/joeviz/functional_annotation/function_stat/";
#my $funcannotdirupdated = "/home/projects/ku_00039/people/joeviz/functional_annotation/function_stat_v2/function_stat/";
# Files: GAGA-0001/annotation.xls
# Unigene Swissprot       KEGG    COG     TrEMBL  Interpro
#my $eggnogdir = "/home/projects/ku_00039/people/joeviz/functional_annotation/eggnog-mapper/all_gaga_final_v2_parallel/";
#Files: GAGA-0001/out_GAGA-0001.emapper.annotations
#my $interprodir = "/home/projects/ku_00039/people/joeviz/functional_annotation/interpro_finalannot/";
#Files: GAGA-0001/GAGA-0001.fasta.tsv


# Open each proteome

	$line = $annotdir;
	my $gagaid = "";
	if ($line =~ /\S*\/(\S+?)\_final_annot/ || $line =~ /\S*\/(\S+?)\.fasta/){
		$gagaid = $1;
	} else {die "Can't find any id to use as name in $line\n";}

	print "Reading $gagaid\n\n";

	system ("mkdir -p $gagaid");
	system ("mkdir -p $gagaid\/all_functional_annotation_files");

	open (Results, ">", "$gagaid\/$gagaid\_final_annotation_repfilt_addreannot_noparpse_representative\_FunctionalAnnotation_summary.tsv");
	print Results "Protein id\tSwissprot\tEggNog\tKEGG\tCOG\tTrEMBL\tInterpro\tNumber of GOs\tNumber of GOs all combined\n";

#	open (ResultsGO, ">", "$gagaid\/$gagaid\_final_annotation_repfilt_addreannot_noparpse_representative\_FunctionalAnnotation_GO.annot");

	my %fannotsp; my %fannotegg; my %fannotkeg; my %fannotcog; my %fannottr; my %fannotipr; my %fannotgonum; my %fannotgonumcombined;

	my %kegannotcombined; my %goannotcombined;

	my $goinipr = "0"; my $goinegg = "0"; my $goinmerged = "0"; my $gocombined = "0";
	my $goiniprcounted = ""; my $goinmergedcounted = "";
	my $keginkeg = "0"; my $keginegg = "0"; my $kegcombined = "0";


	# Read all proteins	
	my %fasta;
	open (Infasta, "<", $line);
	while (<Infasta>){
		chomp;
		my $line2 = $_;
		next if ($line2 !~ /\S+/);
		if ($line2 =~ />(\S+)/){
			$name = $1;
		} else {
			$fasta{$name} .= $line2;
		}
	}
	close Infasta;


	# Read annotation summary
	
	my $annotfile = "$funcannotdirupdated";
	if (-f $annotfile) { 
#		print "the file exists\n"; 
		open (Anofile, "<", $annotfile);
		while (<Anofile>){
			chomp;
			my $line2 = $_;
			next if ($line2 !~ /\S+/);
			next if ($line2 =~ /Unigene\tSwissprot/);
			
			my @subl = split (/\t/, $line2);
			my $prot = $subl[0];
			$fannotsp{$prot} = $subl[1];
			$fannottr{$prot} = $subl[4];
			$fannotcog{$prot} = $subl[3];
			$fannotkeg{$prot} = $subl[2];

			if ($subl[2] =~ /(K\d\d\d\d\d)/){
				my $kegid = $1;
				my $kegpval = "";
				if ($subl[2] =~ /\/(\d[^\/]*\d)\//){ # Get the p-value
					$kegpval = $1;
					if ($kegpval <= /1e-20/){ # Only save KEGG hits with e-value lower than 1e-20
						$kegannotcombined{$prot} = $kegid;
						$keginkeg++;
					}
				} else {die "Cannot find p-value for KEGG hit $prot in $subl[2]\n";}
				#$kegannotcombined{$prot} = $1;
				#$keginkeg++;
			}
		}
		close Anofile;
	} else { 
		die "Can't find $gagaid annotation $annotfile file\n"; 
	}
	system ("cp $funcannotdir\/* $gagaid\/all_functional_annotation_files/");
#	system ("cp $funcannotdirupdated $gagaid\/all_functional_annotation_files/");	
#	system ("rm $gagaid\/all_functional_annotation_files/*fasta");


	# Read interpro file

	my %iprpval;

	my $iprfile = "$interprozijun";
	if (-f $iprfile) { 
#		print "the file exists\n"; 
		open (Anofile, "<", $iprfile);
		while (<Anofile>){
			chomp;
			my $line2 = $_;
			next if ($line2 !~ /\S+/);
			next if ($line2 =~ /Query_id\tSubject_id/);
			
			my @subl = split (/\t/, $line2);
			my $prot = $subl[0];
			if (exists $fannotipr{$prot}){
				if ($iprpval{$prot} !~ /\d/ && $subl[5] =~ /\d/){
					$fannotipr{$prot} = "$subl[1]\/$subl[5]\/$subl[6]";
					$iprpval{$prot} = $subl[5];
				} elsif ($subl[5] !~ /\d/){
					# Keep the previous
				} elsif ($subl[5] < $iprpval{$prot}){
					$fannotipr{$prot} = "$subl[1]\/$subl[5]\/$subl[6]";
					$iprpval{$prot} = $subl[5];
				}
			} else {
				$fannotipr{$prot} = "$subl[1]\/$subl[5]\/$subl[6]";
				$iprpval{$prot} = $subl[5];
			}	
		}
		close Anofile;
	} else { 
		die "Can't find $gagaid interpro annotation $iprfile file\n"; 
	}	


	# Copy original interpro

	my $iproriginal = "$interprodir";
	if (-f $iproriginal) { 
#		print "the file exists\n"; 
		system ("cp $iproriginal $gagaid\/all_functional_annotation_files/out_$gagaid\_interpro.tsv");
	} else { 
		die "Can't find $gagaid interpro tsv annotation $iproriginal file\n"; 
	}

	# Read GO from original interpro

	my $iproriginalfile = "$gagaid\/all_functional_annotation_files/out_$gagaid\_interpro.tsv";
	if (-f $iproriginalfile) { 
#		print "the file exists\n"; 
		open (Anofile, "<", $iproriginalfile);
		while (<Anofile>){
			chomp;
			my $line2 = $_;
			next if ($line2 !~ /\S+/);
			next if ($line2 =~ /Query_id\tSubject_id/);
			
			my @subl = split (/\t/, $line2);
			my $prot = $subl[0];
			next unless (exists $subl[13]); # Skip if last column with GO does not exist
			if ($subl[13] =~ /GO/){
				if ($goiniprcounted =~ / $prot /){
					#Already counter
				} else {
					$goinipr++;
					$goiniprcounted .= " $prot ";
				}

				my @subgo = split (/\|/, $subl[13]);
				foreach my $goipr (@subgo){
					if (exists $goannotcombined{$prot}){
						if ($goannotcombined{$prot} =~ /$goipr/ ){
							# already in
						}
						else {
							$goannotcombined{$prot} .= "\t$goipr";
						}
					} else {
						$goannotcombined{$prot} = "$goipr";
					}
				}
			}	
		}
		close Anofile;
	} else { 
		die "Can't find $gagaid interpro original annotation $iproriginalfile file\n"; 
	}	



	# Read eggnog table

#	my $eggnogfile = "$eggnogdir\/$gagaid\/out_$gagaid\.emapper.annotations";
	my $eggnogfile = "$eggnogdir";
	if (-f $eggnogfile) { 
#		print "the file exists\n"; 
		open (Anofile, "<", $eggnogfile);
		while (<Anofile>){
			chomp;
			my $line2 = $_;
			next if ($line2 !~ /\S+/);
			next if ($line2 =~ /^\#/);			
			
			my @subl = split (/\t/, $line2);
			my $prot = $subl[0];
			if ($subl[8] !~ /^\-$/){
				$fannotegg{$prot} = "$subl[8]\/$subl[7]";
			} else {
				$fannotegg{$prot} = "$subl[7]";
			}

			# GOs
			if ($subl[9] =~ /GO/){
				$goinegg++;

				my @subgo = split (/,/, $subl[9]);
				foreach my $goipr (@subgo){
					if (exists $goannotcombined{$prot}){
						if ($goannotcombined{$prot} =~ /$goipr/ ){
							# already in
						}
						else {
							$goannotcombined{$prot} .= "\t$goipr";
						}
					} else {
						$goannotcombined{$prot} = "$goipr";
					}
				}
			}	

			# KEGGs
			if ($subl[11] =~ /K/){
				$keginegg++;

				my @subkeg = split (/,/, $subl[11]);
				foreach my $kegfull (@subkeg){
					my $kogid = "";
					if ($kegfull =~ /(K\d\d\d\d\d)/ ){
						$kogid = $1;
					} else {
						die "Cannot find KEGG ID KNNNNN in $line2 in $eggnogfile\n";
					}
					if (exists $kegannotcombined{$prot}){
						if ($kegannotcombined{$prot} =~ /$kogid/ ){
							# already in
						}
						else {
							$kegannotcombined{$prot} .= "\t$kogid";
						}
					} else {
						$kegannotcombined{$prot} = "$kogid";
					}
				}
			}	

		}
		close Anofile;
	} else { 
		die "Can't find $gagaid EggNog annotation $eggnogfile file\n"; 
	}
	system ("cp $eggnogfile $gagaid\/all_functional_annotation_files/");

 
	# Read GO table

#	my $gofile = "$godir\/$gagaid\_gene_merged.wego";
	my $gofile = "$godir";
	if (-f $gofile) { 
#		print "the file exists\n"; 
		open (Anofile, "<", $gofile);
		while (<Anofile>){
			chomp;
			my $line2 = $_;
			next if ($line2 !~ /\S+/);
			next if ($line2 =~ /^\#/);			
			
			my @subl = split (/\t/, $line2);
			#my $prot = $subl[0];
			my $prot = shift @subl;
			$fannotgonum{$prot} = 0;
			foreach my $go (@subl){
				$fannotgonum{$prot}++;
#				print ResultsGO "$prot\t$go\n";

				if ($goinmergedcounted =~ / $prot /){
					#Already counted
				} else {
					$goinmerged++;
					$goinmergedcounted .= " $prot ";
				}
#=h
				if (exists $goannotcombined{$prot}){
					if ($goannotcombined{$prot} =~ /$go/ ){
						# already in
					}
					else {
						$goannotcombined{$prot} .= "\t$go";
					}
				} else {
					$goannotcombined{$prot} = "$go";
				}
#=cut
			}
		}
		close Anofile;
	} else { 
		die "Can't find $gagaid GO annotation $gofile file\n"; 
	}
#	system ("cp $gofile $gagaid\/$gagaid\_final_annotation_repfilt_addreannot_noparpse_representative\_FunctionalAnnotation_GO.tsv");
#	close ResultsGO;


	open (ResultsGOcomb, ">", "$gagaid\/$gagaid\_final_annotation_repfilt_addreannot_noparpse_representative\_FunctionalAnnotation_GO.annot");
	open (ResultsGOcombtab, ">", "$gagaid\/$gagaid\_final_annotation_repfilt_addreannot_noparpse_representative\_FunctionalAnnotation_GO.tsv");
# GO was to check, all good with the already merged
	open (ResultsKEGGcomb, ">", "$gagaid\/$gagaid\_final_annotation_repfilt_addreannot_noparpse_representative\_FunctionalAnnotation_KEGG.annot");
	open (ResultsKEGGcombtab, ">", "$gagaid\/$gagaid\_final_annotation_repfilt_addreannot_noparpse_representative\_FunctionalAnnotation_KEGG.tsv");

	foreach my $prot (sort keys %fasta){
#=h
		if (exists $goannotcombined{$prot}){
			$gocombined++;
			print ResultsGOcombtab "$prot\t$goannotcombined{$prot}\n";
			my @subgo = split (/\t/, $goannotcombined{$prot});
			foreach my $go (@subgo){
				$fannotgonumcombined{$prot}++;
				print ResultsGOcomb "$prot\t$go\n";
			}
		} else {
			#print ResultsGOcombtab "$prot\n"; # Print line but empty GO
		}
#=cut
		if (exists $kegannotcombined{$prot}){
			$kegcombined++;
			print ResultsKEGGcombtab "$prot\t$kegannotcombined{$prot}\n";
			my @subkeg = split (/\t/, $kegannotcombined{$prot});
			foreach my $keg (@subkeg){
				$kegannotcombined{$prot}++;
				print ResultsKEGGcomb "$prot\t$keg\n";
			}
		} else {
			#print ResultsKEGGcombtab "$prot\n"; # Print line but empty GO
		}
	}

#	close ResultsGOcomb;
#	close ResultsGOcombtab;
	close ResultsKEGGcomb;
	close ResultsKEGGcombtab;

	# Print summary table

	foreach my $prot (sort keys %fasta){
		my $swissprot = "NA";
		if (exists $fannotsp{$prot}){
			$swissprot = $fannotsp{$prot};
		}

		my $eggnog = "NA";
		if (exists $fannotegg{$prot}){
			$eggnog = $fannotegg{$prot};
			if ($eggnog =~ /^\-$/){
				$eggnog = "NA";
			}
		}

		my $kegg = "NA";
		if (exists $fannotkeg{$prot}){
			$kegg = $fannotkeg{$prot};
		}

		my $cog = "NA";
		if (exists $fannotcog{$prot}){
			$cog = $fannotcog{$prot};
		}

		my $trembl = "NA";
		if (exists $fannottr{$prot}){
			$trembl = $fannottr{$prot};
		}

		my $interpro = "NA";
		if (exists $fannotipr{$prot}){
			$interpro = $fannotipr{$prot};
		}

		my $gonumber = "0";
		if (exists $fannotgonum{$prot}){
			$gonumber = $fannotgonum{$prot};
		}	

		my $gonumbercombined = "0";
		if (exists $fannotgonumcombined{$prot}){
			$gonumber = $fannotgonumcombined{$prot};
		}	

		print Results "$prot\t$swissprot\t$eggnog\t$kegg\t$cog\t$trembl\t$interpro\t$gonumber\t$gonumbercombined\n";
	}

	close Results;

	print "Genes with GO annot in IPR: $goinipr\nIn EggNog: $goinegg\nIn merged file: $goinmerged\nCombined in all: $gocombined\n\n";
	print "Genes with KEGG annot in KEGG annot: $keginkeg\nIn EggNog: $keginegg\nCombined in all: $kegcombined\n\n";
	print "\n\n";


# Gzip at the end
#system ("gzip *\/*fasta");
#system ("gzip *\/*gff");














