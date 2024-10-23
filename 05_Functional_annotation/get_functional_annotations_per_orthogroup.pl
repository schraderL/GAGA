#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
my $dirname = dirname(__FILE__);

# usage: perl get_all_files.pl (change lines where the input files are called)

my ($line, $name, $nameout);


my $ogtable = "/home/projects/ku_00039/people/joeviz/orthofinder/run_allGAGA_final_annotations/toerda/Orthogroup_table_N0.tsv"; # Ortholog table from orthofinder
my $fannotdir = "/home/projects/ku_00039/people/joeviz/GAGA_annotations/Final_functional_annotation_wKEGG/"; # Functional annotation folder, it contains a folder for each ant genome (i.e.: GAGA-0001/) and it is generated in "get_all_functional_annotation.pl script"


# Load here functional annotations (tables and GOs)

system ("ls $fannotdir\/\*\/\*FunctionalAnnotation_summary.tsv > inputsums.txt");
#system ("ls $fannotdir\/GAGA-0001\/\*FunctionalAnnotation_summary.tsv > inputsums.txt"); # Only 1 genome to test

my %annotswiss; my %annotegg; my %annotipr; my %annotgos;  my %annotkegs;

open (File , "<", "inputsums.txt"); 
while (<File>) {
	chomp;
	$line = $_;
	next if ($line !~ /\S+/);
	my $gagaid = "";
	if ($line =~ /(GAGA\-\S\S\S\S)/ || $line =~ /(NCBI\-\S\S\S\S)/ || $line =~ /(OUT\-\S\S\S\S)/){
		$gagaid = $1;
	} else {die "Can't find GAGA id in $line\n";}

	# Read summary table
	open (Anofile, "<", $line);
	while (<Anofile>){
		chomp;
		my $line2 = $_;
		next if ($line2 !~ /\S+/);
		next if ($line2 =~ /Protein id\tSwissprot/);
			
		my @subl = split (/\t/, $line2);
		my $prot = $subl[0];

		if ($subl[1] =~ /\S+\/\S+\/(.*)/){
			$annotswiss{$gagaid}{$prot} = $1;
		}

		if ($subl[2] =~ /\S+\/(.*)/){
			$annotegg{$gagaid}{$prot} = $1;
		} elsif ($subl[2] !~ /^NA$/) {
			$annotegg{$gagaid}{$prot} = $subl[2];
		}

		if ($subl[6] =~ /(\S+)\/\S+\/(.*)/){
			$annotipr{$gagaid}{$prot} = "$1\/$2";
		}
	}
	close Anofile;

	# Read GO
	my $gofile = "$fannotdir\/\/$gagaid\/$gagaid\_final_annotation_repfilt_addreannot_noparpse_representative_FunctionalAnnotation_GO.tsv";
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
#			$annotgos{$gagaid}{$prot} = join('\t', @subl);
			@{$annotgos{$gagaid}{$prot}} = @subl;

		}
		close Anofile;
	} else { 
		die "Can't find $gagaid GO annotation $gofile file\n"; 
	}

	# Read KEGG
	my $kegfile = "$fannotdir\/\/$gagaid\/$gagaid\_final_annotation_repfilt_addreannot_noparpse_representative_FunctionalAnnotation_KEGG.tsv";
	if (-f $kegfile) { 
#		print "the file exists\n"; 
		open (Anofile, "<", $kegfile);
		while (<Anofile>){
			chomp;
			my $line2 = $_;
			next if ($line2 !~ /\S+/);
			next if ($line2 =~ /^\#/);			
			
			my @subl = split (/\t/, $line2);
			#my $prot = $subl[0];
			my $prot = shift @subl;
#			$annotgos{$gagaid}{$prot} = join('\t', @subl);
			@{$annotkegs{$gagaid}{$prot}} = @subl;

		}
		close Anofile;
	} else { 
		die "Can't find $gagaid KEGG annotation $kegfile file\n"; 
	}

}


# Open HOG table

open (Results, ">", "Orthogroups_functional_annotation_summary.tsv");
print Results "Orthogroup\tNumber of sequences\tNumber of species\tSwissprot\tNumber of seqs with that Swissprot\tEggNog\tNumber of seqs with that EggNog\tInterPro\tNumber of seqs with that InterPro\tNumber of GOs (at least 33\% sequences)\tTotal number of GOs\tNumber of KEGGs (at least 33\% sequences)\tTotal number of KEGGs\n";
open (ResultsGO, ">", "Orthogroups_functional_annotation_GO.annot");
open (ResultsGOtsv, ">", "Orthogroups_functional_annotation_GO.tsv");
open (ResultsGOall, ">", "Orthogroups_functional_annotation_GO_allnofilter.annot");
open (ResultsGOtsvall, ">", "Orthogroups_functional_annotation_GO_allnofilter.tsv");
open (ResultsKEGG, ">", "Orthogroups_functional_annotation_KEGG.annot");
open (ResultsKEGGtsv, ">", "Orthogroups_functional_annotation_KEGG.tsv");
open (ResultsKEGGall, ">", "Orthogroups_functional_annotation_KEGG_allnofilter.annot");
open (ResultsKEGGtsvall, ">", "Orthogroups_functional_annotation_KEGG_allnofilter.tsv");
my $header = "";
open (File , "<", $ogtable); 
while (<File>) {
	chomp;
	$line = $_;
	next if ($line !~ /\S+/);
	
	if ($line =~ /Parent Clade\t(.*)/){
		$header = $1;
		next;
	}
	my @gagaids = split (/\t/, $header);
	
	my @subl = split (/\t/, $line);
	my $ogid = "";
	if ($subl[0] =~ /\S+\.(HOG\S+)/){
		$ogid = $1;
	}

	my %hoggo; my %hogkeg; my %hogsprot; my %hogegg; my %hogipr;
	my $numsp = "0"; my $numseqs = "0";

	my $discard = shift @subl; $discard = shift @subl; $discard = shift @subl; 
	my $col = 0;
	foreach my $seqs (@subl){
		if ($seqs !~ /\S+/){
			$col++;
			next;
		}
		my @arrayseqs = split (/\,/, $seqs);
		foreach my $seq (@arrayseqs){
			next if ($seq !~ /\S+/);
			$seq =~ s/\s+//g;
			$numseqs++;
			my $gagaid = $gagaids[$col];
			# Get GOs
			if (exists $annotgos{$gagaid}{$seq}){
#				my $golist = $annotgos{$gagaid}{$seq};
#				my @goarray = split (/\t/, $golist);
				my @goarray = @{$annotgos{$gagaid}{$seq}};
				foreach my $eachgo (@goarray){
					$hoggo{$eachgo}++;
				}
			}
			# Get KEGGs
			if (exists $annotkegs{$gagaid}{$seq}){
#				my $keglist = $annotkegs{$gagaid}{$seq};
#				my @kegarray = split (/\t/, $keglist);
				my @kegarray = @{$annotkegs{$gagaid}{$seq}};
				foreach my $eachkeg (@kegarray){
					$hogkeg{$eachkeg}++;
				}
			}
			# Get other annotations
			if (exists $annotegg{$gagaid}{$seq}){
				my $annolist = $annotegg{$gagaid}{$seq};
				my @anoarray = split (/\t/, $annolist);
				foreach my $eachano (@anoarray){
					$hogegg{$eachano}++;
				}
			}
			if (exists $annotipr{$gagaid}{$seq}){
				my $annolist = $annotipr{$gagaid}{$seq};
				my @anoarray = split (/\t/, $annolist);
				foreach my $eachano (@anoarray){
					$hogipr{$eachano}++;
				}
			}
			if (exists $annotswiss{$gagaid}{$seq}){
				my $annolist = $annotswiss{$gagaid}{$seq};
				my @anoarray = split (/\t/, $annolist);
				foreach my $eachano (@anoarray){
					$hogsprot{$eachano}++;
				}
			}						
		}
		if ($numseqs > 0){
			$numsp++;
		}
		$col++;
	}

	# Print GOs
	my $gonum = 0;
	my $gonumall = 0;
	if ( !keys %hoggo){
		#Empty
	} else {
		print ResultsGOtsvall "$ogid\t";
		foreach my $key (sort keys %hoggo){
			my $times = $hoggo{$key};
			$gonumall++;
			print ResultsGOall "$ogid\t$key\n";
			print ResultsGOtsvall "$key\t";
			
			my $filter = $numseqs/3;
			if ($times >= $filter){
				if ($gonum == 0){
					print ResultsGOtsv "$ogid\t";
				}
				$gonum++;
				print ResultsGO "$ogid\t$key\n";
				print ResultsGOtsv "$key\t";				
			}
		}
	}
	if ($gonumall > 0){
		print ResultsGOtsvall "\n";
	}
	if ($gonum > 0){
		print ResultsGOtsv "\n";
	}
	# Print KEGGs
	my $kegnum = 0;
	my $kegnumall = 0;
	if ( !keys %hogkeg){
		#Empty
	} else {
		print ResultsKEGGtsvall "$ogid\t";
		foreach my $key (sort keys %hogkeg){
			my $times = $hogkeg{$key};
			$kegnumall++;
			print ResultsKEGGall "$ogid\t$key\n";
			print ResultsKEGGtsvall "$key\t";
			
			my $filter = $numseqs/3;
			if ($times >= $filter){
				if ($kegnum == 0){
					print ResultsKEGGtsv "$ogid\t";
				}
				$kegnum++;
				print ResultsKEGG "$ogid\t$key\n";
				print ResultsKEGGtsv "$key\t";				
			}
		}
	}
	if ($kegnumall > 0){
		print ResultsKEGGtsvall "\n";
	}
	if ($kegnum > 0){
		print ResultsKEGGtsv "\n";
	}

	my $sprotannot = "NA";
	my $sprotnum = 0;
	if ( !keys %hogsprot){
		#Empty
	} else {
		foreach my $key (sort keys %hogsprot){
			my $times = $hogsprot{$key};
			if ($times > $sprotnum){
				$sprotnum = $times;
				$sprotannot = $key;
			}
		}
	}

	my $iprannot = "NA";
	my $iprnum = 0;
	if ( !keys %hogipr){
		#Empty
	} else {
		foreach my $key (sort keys %hogipr){
			my $times = $hogipr{$key};
			if ($times > $iprnum){
				$iprnum = $times;
				$iprannot = $key;
			}
		}
	}

	my $eggannot = "NA";
	my $eggnum = 0;
	if ( !keys %hogegg){
		#Empty
	} else {
		foreach my $key (sort keys %hogegg){
			my $times = $hogegg{$key};
			if ($times > $eggnum){
				$eggnum = $times;
				$eggannot = $key;
			}
		}
	}

#print Results "Orthogroup\tNumber of sequences\tNumber of species\tSwissprot\tNumber of seqs with that Swissprot\tEggNog\tNumber of seqs with that EggNog\tInterPro\tNumber of seqs with that InterPro\tNumber of GOs\t";
	print Results "$ogid\t$numseqs\t$numsp\t$sprotannot\t$sprotnum\t$eggannot\t$eggnum\t$iprannot\t$iprnum\t$gonum\t$gonumall\t$kegnum\t$kegnumall\n";

}
close File;













