#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
my $dirname = dirname(__FILE__);

# Changed to rename ORs in scaffold order
# Changed to keep in the final GFF all OR > 100|50aa, name them as fragmentOR[number]

# usage: perl run_OR_classification.pl ABCENTH_clean.gff OUTNAME(NCBI-0001) GENOME
# perl run_OR_classification.pl testcleanfix.gff3 NCBI-0001 /home/projects/ku_00039/people/joeviz/GAGA_genomes/Genome_assemblies/GAGA_all_final_assemblies_oldv/NCBI-0001_Ooceraea_biroi_dupsrm_filt.fasta

# The script will use the following files that need to be previously generated: ABCENTH_clean.pep.fasta ABCENTH_clean.pep.fasta.tsv(interpro) ABCENTH_clean.pep.fasta.ORcoblast.txt ABCENTH_clean.pep.fasta.ORblast.txt  
# It requires to have the file with GAGA id and short name to rename the genes in $gagatable line 39.

my $orcompletelength = "350"; # Length to consider complete OR copies
my $minorlength = "200"; # Minimum length to retain partial ORs

my ($line, $name, $nameout);
my $gff = $ARGV[0];
my $prefix = "";
if ($gff =~ /(\S+)\.gff3/){
	$prefix = $1;
} else {die "Can't find prefix on gff3 file\n";}

my $gagaid = $ARGV[1];
my $genome = $ARGV[2];

my $outgffall = "$gagaid\_"."$prefix"."\_OR_renamed_all.gff3"; 
my $outprotall = "$gagaid\_"."$prefix"."\_OR_renamed_all"; 

my $outgff = "$gagaid\_"."$prefix"."\_OR_renamed_all_nofragment.gff3"; 
my $outprot = "$gagaid\_"."$prefix"."\_OR_renamed_all_nofragment"; 

my $outgffcomp = "$gagaid\_"."$prefix"."\_OR_renamed_complete.gff3"; 
my $outprotcomp = "$gagaid\_"."$prefix"."\_OR_renamed_complete"; 

# Read GAGA table
my $gagatable = "/home/projects/ku_00039/people/joeviz/GAGA_species_list.txt";
my %gagasp;
my %gagashort;
open(File, "<", $gagatable);
while(<File>){
	chomp;
	my $line = $_;
	$line =~ s/\r//g;
	next if ($line !~ /\S+/);
	my @subl = split (/\t/, $line);
	$gagasp{$subl[0]} = $subl[4];
	$gagashort{$subl[0]} = $subl[3];
}
close File;


# Reading Protein fasta

my %fasta;
my $protfile = "$prefix".".pep.fasta";
open(File, "<", $protfile);
while(<File>){
	chomp;
	$line = $_;
	next if ($line !~ /\S+/);	

	if ($line =~ />(\S+)/){
		$name = $1;
	} else {
		$fasta{$name} .= "$line";
	}
}
close File;



# Reading tsv

my %interpro;
my %pfam;
my $interprofile = "$prefix".".pep.fasta.tsv";
open (File, "<", $interprofile);
while(<File>){
	chomp;
	$line = $_;
	next if ($line !~ /\S+/);	
	my @subline = split (/\t/, $line);
	my $totalcol = scalar(@subline);
	next if ($totalcol <= 11); ## Avoid coil records

	$pfam{$subline[0]} .= "$subline[4] ";

	if (exists $interpro{$subline[0]}){
		if ($interpro{$subline[0]} =~ /$subline[11]\s/){
			next;
		} else {
			$interpro{$subline[0]} .= "$subline[11] ";
		}
	}	else {
		$interpro{$subline[0]} .= "$subline[11] ";
	}
}
close File;



# Reading blast results

my %orco;
my $blastorcofile = "$prefix".".pep.fasta.ORcoblast.txt";
open (File, "<", $blastorcofile);
while(<File>){
	chomp;
	$line = $_;
	next if ($line !~ /\S+/);	
	my @subline = split (/\t/, $line);

	if ($subline[10] < "1e-100"){
		$orco{$subline[0]} = "$subline[1]";
	}
}
close File;

my %orblast;
my $blastorfile = "$prefix".".pep.fasta.ORblast.txt";
open (File, "<", $blastorfile);
while(<File>){
	chomp;
	$line = $_;
	next if ($line !~ /\S+/);	
	my @subline = split (/\t/, $line);

	if ($subline[10] < "1e-5"){
		if (exists $orblast{$subline[0]}){
			if ($subline[10] < $orblast{$subline[0]}){
				$orblast{$subline[0]} = "$subline[10]";
			}
		} else {
			$orblast{$subline[0]} = "$subline[10]";
		}

	}
}
close File;

my %grblast;
my $blastgrfile = "$prefix".".pep.fasta.GRblast.txt";
open (File, "<", $blastgrfile);
while(<File>){
	chomp;
	$line = $_;
	next if ($line !~ /\S+/);	
	my @subline = split (/\t/, $line);

	if ($subline[10] < "1e-5"){
		if (exists $grblast{$subline[0]}){
			if ($subline[10] < $grblast{$subline[0]}){
				$grblast{$subline[0]} = "$subline[10]";
			}
		} else {
			$grblast{$subline[0]} = "$subline[10]";
		}

	}
}
close File;


# Reading GFF to obtain gene name from Sean pipeline

my %clusterorder;
#my %seannamerev;
my %seanname;
my %exonnumber;
open (GFFfile , "<", $gff); 
while (<GFFfile>) {
	chomp;
	$line = $_;
	next if ($line !~ /\S+/);
	next if ($line =~ /^#/);
	my @subline = split (/\t/, $line);


	if ($subline[2] =~ /mRNA/){

		my $genename = "";
		if ($subline[8] =~ /ID=([^;]+)/){
			$genename = $1;
		}
		else {print "ERROR in run_OR_classification: It fails detecting ID in GFF file $line\n";}

		my $genenamesean = "";
		if ($subline[8] =~ /transcript_id=([^;]+)/){
			$genenamesean = $1;
		}
		else {print "ERROR in run_OR_classification: It fails detecting gene ID Sean name in GFF file $line\n";}

		$seanname{$genename} = $genenamesean;

		my $seannameclust = "";
		my $seannameclustcoord = "";
		if ($genenamesean =~ /coords(\S+)\-(\d+)/ ){
			$seannameclust = $1;
			$seannameclustcoord = $2;				
		} else {die "Error in run_OR_clasffication: Can't find coords in $line\n";}

		$clusterorder{$seannameclust}{$seannameclustcoord}{$genenamesean} = $genename; 
#		$clusterorder{$seannameclust} = $genenamesean; 
#		$seannamerev{$genenamesean} = $genename;

	} elsif ($subline[2] =~ /CDS/){

		my $genename = "";
		if ($subline[8] =~ /Parent=([^;]+)/){
			$genename = $1;
		}
		else {print "ERROR in run_OR_classification: It fails detecting ID in GFF file $line\n";}

		$exonnumber{$genename}++;

	}
}
close GFFfile;


# Creating table and dividing proteins

my %newnames;
my $deleteproteins;
my $deleteproteinsl50;
my %hstatus;

open (Results, ">", "$gagaid\_ORs.txt");
print Results "GAGA ID\tSpecies name\tGene name\tOriginal name\tOR name\tLength\tExon number\tOR blast hit\tOR domain\tStatus Complete Partial Pseudogene\tFragment OR\n";

open (Resultssum, ">", "$gagaid\_ORs_summary.txt");
print Resultssum "GAGA ID\tSpecies name\tOR total number\tComplete\tComplete with domain\tPartial\tPartial with domain\tPartial longer 350aa\tPseudogene\tPseudogene with domain\tPseudogene longer 350aa\tLength lower than 200aa\tLength lower than 100aa\tLength lower than 50aa\tOrco\tOrco complete\tNine-exon ORs complete\tNine-exon ORs incomplete\tNine-exon ORs incomplete (length > 350aa)\tNine-exon ORs pseudogenes\tOR complete + partial with domain\tFinal retained ORs\tFinal retained complete ORs (complete + partial >350aa) \tFinal retained partial ORs\tFinal retained pseudogene ORs\tFinal retained complete 9-exon ORs\tFinal retained incomplete or pseudogene 9-exon ORs\tFragment ORs partial\tFragment ORs pseudogenes\n";

my $ornum = "1";
my $fragornum = "1";
my $totalor = "0"; my $com = "0"; my $comd = "0"; my $par = "0"; my $pard = "0"; my $parc = "0"; my $pse = "0"; my $psed = "0"; my $psec = "0"; my $l50 = "0"; my $l100 = "0"; my $l200 = "0"; my $orconum = "0"; my $orconumcom = "0"; my $ninec = "0"; my $ninei = "0"; my $nineic = "0"; my $ninep = "0";
my $compard = "0"; my $ftot = "0"; my $fcom = "0"; my $fpar = "0"; my $fparc = "0"; my $fcomp = "0"; my $fpse = "0"; my $fninec = "0"; my $fninei = "0"; my $frpartial = "0"; my $frpse = "0";

my %lsfasta;
my %lengthseq;
foreach my $prot (sort keys %fasta){
	my $length = length ($fasta{$prot});
	$lsfasta{$length}{$prot} = "$fasta{$prot}";
	$lengthseq{$prot} = $length;
}

#foreach my $len (sort { $b <=> $a } (keys %lsfasta)){
#	foreach my $prot (keys %{$lsfasta{$len}}){
foreach my $seanprotscaf (sort (keys %clusterorder)){
foreach my $seanprotcoord (sort { $a <=> $b } (keys %{$clusterorder{$seanprotscaf}})){
	foreach my $seanprot (keys %{$clusterorder{$seanprotscaf}{$seanprotcoord}}){

		my $prot = $clusterorder{$seanprotscaf}{$seanprotcoord}{$seanprot};
		my $len = $lengthseq{$prot};

		my $seq = $lsfasta{$len}{$prot};

		$totalor++;

		my $orcohit = "None";
		if (exists $orco{$prot}){
			if ($orco{$prot} =~ /Orco/){
				$orcohit = "Yes";
				$orconum++;
#				if ($status =~ /Complete/){
#					$orconumcom++;
#				}
			}
		}

		my $status = "";
		my $discard = "No"; # Label for fragment ORs

		if ($len >= $orcompletelength){
			$status = "Complete";
			$hstatus{$prot} = "Complete";
#				$com++;
		} else {
			$status = "Partial";
			$hstatus{$prot} = "Partial";
#			$par++;
			if ($len < 200){
				$l200++;
#				$deleteproteins .= "$prot ";	
				#$discard = "Yes";				
			}
			if ($len < 100){
				$l100++;
#				$deleteproteins .= "$prot ";	
				#$discard = "Yes";				
			}
			if ($len < 50){
				$l50++;
#				$deleteproteins .= "$prot ";   # Delete all protein shorter than 50aa
				$deleteproteinsl50 .= "$prot ";   # Delete all protein shorter than 50aa
#				$discard = "Yes";
			}
			if ($len < $minorlength){
				if ($orcohit =~ /None/){ # Avoid deleting ORCO even if it is partial
					$deleteproteins .= "$prot ";	# Delete all proteins shorter than the specified minimum length
					$discard = "Yes";		
				}		
			}				
		} 

		my $middleseq = substr ($seq, 5, -5);
		if ($middleseq =~ /X/){
			$status = "Pseudogene";
			$hstatus{$prot} = "Pseudogene";
#			$pse++;
		}


		my $sufix = "";
		my $sname = $seanname {$prot};
		if ($sname =~ /(\S+\d)\-RA/){
			$sufix = "";
#			$status = "Complete"; # Correct based on naming
			if ($status !~ /Complete/){ # Add sufix if not complete
				$sufix = "S"; # Use S "short" for short models without any tag from Sean pipline (P, I, C, N)
			}
		} elsif ($sname =~ /\S+\d+(\D+)\-RA/){
			$sufix = $1;	
			if ($sufix =~ /P/){
				$status = "Pseudogene";    # Correcting partial copies based on Sean output
				$hstatus{$prot} = "Pseudogene";
			} else {
				$status = "Partial";    # Correcting partial copies based on Sean output
				$hstatus{$prot} = "Partial";				
			}

#			if ($status =~ /Complete/){
#				$status = "Partial";    # Correcting partial copies based on Sean output
#				$hstatus{$prot} = "Partial";
#			}
		}

		my $orhit = "None";
		if (exists $orblast{$prot} && exists $grblast{$prot} ){
			if ($orblast{$prot} <= $grblast{$prot}){
				$orhit = "OR";
			} elsif ($orblast{$prot} > $grblast{$prot}){
				$orhit = "GR";
				$deleteproteinsl50 .= "$prot "; ### Delete protein identified as GR
				$discard = "Yes";
			}  
		} elsif (exists $orblast{$prot}){
			$orhit = "OR";
		} elsif (exists $grblast{$prot}){
			$orhit = "GR";
			$deleteproteinsl50 .= "$prot "; ### Delete protein identified as GR
			$discard = "Yes";
		}

		my $ordomain = "None";
		if (exists $interpro{$prot}){
			if ($interpro{$prot} =~ /IPR004117/){
				$ordomain = "OR";
				if ($status =~ /Pseudogene/){
					$psed++;
				} elsif ($status =~ /Complete/){
					$comd++;
				} elsif ($status =~ /Partial/){
					$pard++;
					$compard++;
				}

			} elsif ($interpro{$prot} =~ /IPR013604/){
				$ordomain = "GR";
				$deleteproteins .= "$prot "; ### Delete protein identified as GR
				$discard = "Yes";
			} else {
				print "Warning: $prot in $gagaid with $len length contains an unexpected domain $interpro{$prot}\n";
				$ordomain = "Other";
				if ($status =~ /Partial/ || $status =~ /Pseudogene/){
					$deleteproteins .= "$prot "; ### Delete protein identified with other domain, and are partial or pseudogenes
					$discard = "Yes";
				} elsif ($status =~ /Complete/ && $orhit !~ /OR/){
					$deleteproteins .= "$prot "; ### Delete complete protein identified with other domain, and does not hit any OR
					$discard = "Yes";
				}

			}
		}


		if ($status =~ /Complete/){
			$com++;
			$compard++;
		} elsif ($status =~ /Partial/){
			$par++;
			if ($len >= $orcompletelength){
				$status = "Complete"; ## Uncomment if we finally keep these also as complete copies
				$hstatus{$prot} = "Complete";
				$parc++;
			}

		} elsif ($status =~ /Pseudogene/){
			$pse++;
			if ($len >= $orcompletelength){
				#$status = "Complete"; ## Uncomment if we finally keep these also as complete copies (Likely not as putative pseudogenes)
				#$hstatus{$prot} = "Complete";
				$psec++;
			}
		}		

		# Final filter to discard proteins
		if ($status =~ /Partial/ || $status =~ /Pseudogene/){
			if ($ordomain =~ /None/ && $len < $minorlength){
				if ($orcohit =~ /None/){ # Avoid deleting ORCO even if it is partial and have no domain
					$deleteproteins .= "$prot "; ### Delete protein partial without domain, and length < specified in minorlength
					$discard = "Yes";		
				}	
			}
		}

		# Count 9 exon ORs
		if ($exonnumber{$prot} == 9){
			if ($status =~ /Complete/){
				$ninec++;
			} elsif ($status =~ /Partial/) {
				$ninei++;
				if ($len >= $orcompletelength){
					$nineic++;
				}				
			} elsif ($status =~ /Pseudogene/) {
				$ninep++;
			}
		}
		unless ($discard =~ /Yes/){
			if ($exonnumber{$prot} == 9){
				if ($status =~ /Complete/){
					$fninec++;
				} else {
					$fninei++;
				}	
			}
		}	

		# Get new names
		my $nname = "";
		if ($discard =~ /Yes/){
			#$nname = "Discarded";
			if ($len < 50){
				$nname = "Discarded";
			} else {
				$nname = "$gagashort{$gagaid}"."fragOR"."$fragornum"."$sufix";
				$fragornum++;
			}
		} else {
			if ($orcohit =~ /Yes/){
				$nname = "$gagashort{$gagaid}"."Orco"."$sufix";
				if ($status =~ /Complete/){
					$orconumcom++;
				}
				if ($orconum > 1){
					$nname = "$gagashort{$gagaid}"."Orco"."$orconum"."$sufix";
				}
			} else {
				$nname = "$gagashort{$gagaid}"."OR"."$ornum"."$sufix";
				$ornum++;
			}
		}
		$newnames{$prot} = "$nname";


		print Results "$gagaid\t$gagasp{$gagaid}\t$prot\t$sname\t$nname\t$len\t$exonnumber{$prot}\t$orhit\t$ordomain\t$status\t$discard\n";

	}

}}
close Results;


# Reading GFF and generating final file filtering short proteins and renaming ORs

open (Resultsgff, ">", "$outgff");
open (Resultsgffcomp, ">", "$outgffcomp");
open (Resultsgffall, ">", "$outgffall");

open (GFFfile , "<", $gff); 
while (<GFFfile>) {
	chomp;
	$line = $_;
	next if ($line !~ /\S+/);

	if ($line =~ /^#/){
		print Resultsgff "$line\n";
		next;
	}
	
	my @subline = split (/\t/, $line);

	if ($subline[2] =~ /CDS/){
		my $genename = "";
		my $rest = "";
		if ($subline[8] =~ /Parent=([^;]+)(\S+)/){
		#if ($subline[8] =~ /transcript_id \"(\S+)\"/){
			$genename = $1;
			$rest = $2;
		}
		else {die "ERROR in run_OR_classification.pl: It fails detecting Parent ID in $line\n";}

		my $nnamef = $newnames{$genename};

# 		Keeping now all genes as fragment ORs; and creating multiple GFFs
#		next if ($deleteproteins =~ /$genename /); ### Delete genes that do not pass the filter
		next if ($deleteproteinsl50 =~ /$genename /); ### Delete from GFF all genes shorted than 50aa

		print Resultsgffall "$subline[0]\t$subline[1]\t$subline[2]\t$subline[3]\t$subline[4]\t$subline[5]\t$subline[6]\t$subline[7]\tParent=$nnamef"."$rest\;Status=$hstatus{$genename}\n";

		unless ($deleteproteins =~ /$genename /){
			print Resultsgff "$subline[0]\t$subline[1]\t$subline[2]\t$subline[3]\t$subline[4]\t$subline[5]\t$subline[6]\t$subline[7]\tParent=$nnamef"."$rest\;Status=$hstatus{$genename}\n";

			if ($hstatus{$genename} =~ /Complete/){
				print Resultsgffcomp "$subline[0]\t$subline[1]\t$subline[2]\t$subline[3]\t$subline[4]\t$subline[5]\t$subline[6]\t$subline[7]\tParent=$nnamef"."$rest\;Status=$hstatus{$genename}\n";
			}
		}
		
	}
	elsif ($subline[2] =~ /mRNA/){
		my $genename = "";
		my $rest = "";
		if ($subline[8] =~ /ID=([^;]+)\;Parent=[^;]+(\;\S+)/){
		#if ($subline[8] =~ /transcript_id \"(\S+)\"/){
			$genename = $1;
			$rest = $2;
		}
		else {print "ERROR in run_OR_classification.pl: It fails detecting ID in $line\n";}

		my $nnamef = $newnames{$genename};

# 		Keeping now all genes as fragment ORs; and creating multiple GFFs
#		next if ($deleteproteins =~ /$genename /); ### Delete genes that do not pass the filter
		next if ($deleteproteinsl50 =~ /$genename /); ### Delete from GFF all genes shorted than 50aa

		print Resultsgffall "$subline[0]\t$subline[1]\tgene\t$subline[3]\t$subline[4]\t$subline[5]\t$subline[6]\t$subline[7]\tID=g"."$nnamef"."$rest;Status=$hstatus{$genename}\n";		
		print Resultsgffall "$subline[0]\t$subline[1]\t$subline[2]\t$subline[3]\t$subline[4]\t$subline[5]\t$subline[6]\t$subline[7]\tID=$nnamef\;Parent=g"."$nnamef"."$rest;Status=$hstatus{$genename}\n";		

		unless ($deleteproteins =~ /$genename /){
			print Resultsgff "$subline[0]\t$subline[1]\tgene\t$subline[3]\t$subline[4]\t$subline[5]\t$subline[6]\t$subline[7]\tID=g"."$nnamef"."$rest;Status=$hstatus{$genename}\n";		
			print Resultsgff "$subline[0]\t$subline[1]\t$subline[2]\t$subline[3]\t$subline[4]\t$subline[5]\t$subline[6]\t$subline[7]\tID=$nnamef\;Parent=g"."$nnamef"."$rest;Status=$hstatus{$genename}\n";		

			if ($hstatus{$genename} =~ /Complete/){
				print Resultsgffcomp "$subline[0]\t$subline[1]\tgene\t$subline[3]\t$subline[4]\t$subline[5]\t$subline[6]\t$subline[7]\tID=g"."$nnamef"."$rest;Status=$hstatus{$genename}\n";		
				print Resultsgffcomp "$subline[0]\t$subline[1]\t$subline[2]\t$subline[3]\t$subline[4]\t$subline[5]\t$subline[6]\t$subline[7]\tID=$nnamef\;Parent=g"."$nnamef"."$rest;Status=$hstatus{$genename}\n";		
			}
		}

		
		# Obtain final OR number statistics
		if ($deleteproteins =~ /$genename /){
			if ($hstatus{$genename} =~ /Complete/){
				# Nothing here
			} elsif ($hstatus{$genename} =~ /Partial/){
				$frpartial++;
			} elsif ($hstatus{$genename} =~ /Pseudogene/){
				$frpse++;
			}		

		} else {
			$ftot++;
			if ($hstatus{$genename} =~ /Complete/){
				$fcom++;
				$fcomp++;
			} elsif ($hstatus{$genename} =~ /Partial/){
				$fpar++;
				if ($lengthseq{$genename} >= $orcompletelength){ # It will not do nothing here as I relabelled those as complete in the status
					$fparc++;
					$fcomp++;
				}
			} elsif ($hstatus{$genename} =~ /Pseudogene/){
				$fpse++;
			}
		}

	}
}
close GFFfile;


#Encode proteins and CDS from the generated GFFs

system ("perl scripts/gff2fasta_v3.pl $genome $outgff $outprot ");
system ("sed \'s\/X\*\$\/\/\' $outprot\.pep.fasta > $outprot\.pep.fasta.tmp");
system ("mv $outprot\.pep.fasta.tmp $outprot\.pep.fasta");

system ("perl scripts/gff2fasta_v3.pl $genome $outgffall $outprotall ");
system ("sed \'s\/X\*\$\/\/\' $outprotall\.pep.fasta > $outprotall\.pep.fasta.tmp");
system ("mv $outprotall\.pep.fasta.tmp $outprotall\.pep.fasta");

system ("perl scripts/gff2fasta_v3.pl $genome $outgffcomp $outprotcomp ");
system ("sed \'s\/X\*\$\/\/\' $outprotcomp\.pep.fasta > $outprotcomp\.pep.fasta.tmp");
system ("mv $outprotcomp\.pep.fasta.tmp $outprotcomp\.pep.fasta");

#print Resultssum "$gagaid\t$gagasp{$gagaid}\t$totalor\t$com\t$comd\t$par\t$pard\t$parc\t$pse\t$psed\t$psec\t$l200\t$l100\t$l50\t$orconum\t$orconumcom\t$ninec\t$ninei\t$nineic\t$ninep\t$compard\t$ftot\t$fcom\t$fpar\t$fparc\t$fcomp\t$fpse\t$fninec\t$fninei\n";
print Resultssum "$gagaid\t$gagasp{$gagaid}\t$totalor\t$com\t$comd\t$par\t$pard\t$parc\t$pse\t$psed\t$psec\t$l200\t$l100\t$l50\t$orconum\t$orconumcom\t$ninec\t$ninei\t$nineic\t$ninep\t$compard\t$ftot\t$fcomp\t$fpar\t$fpse\t$fninec\t$fninei\t$frpartial\t$frpse\n";
close Resultssum;

