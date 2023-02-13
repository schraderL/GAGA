#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
my $dirname = dirname(__FILE__);

# Change to rename ORs in scaffold order
# Change to keep in the final GFF all OR > 100|50aa, name them as fragmentOR[number]

# usage: perl run_OR_classification_bitacora_GR.pl Bitacora.gff3 OUTNAME(NCBI-0001) GENOME
# perl run_OR_classification.pl testcleanfix.gff3 NCBI-0001 /home/projects/ku_00039/people/joeviz/GAGA_genomes/Genome_assemblies/GAGA_all_final_assemblies_oldv/NCBI-0001_Ooceraea_biroi_dupsrm_filt.fasta

# It will search for gff3 prefix the following files: ABCENTH_clean.pep.fasta ABCENTH_clean.pep.fasta.tsv(interpro) ABCENTH_clean.pep.fasta.ORcoblast.txt ABCENTH_clean.pep.fasta.ORblast.txt  

my $orcompletelength = "350"; # Length to consider complete OR copies
my $minorlength = "150"; # Minimum length to retain partial ORs
my $maxorlength = "1500"; # Maximum average length, it will print a warning with the sequences larger than this length

# Seqs shorter than "minorlength" but longer than 50aa will be labelled as fragment. 

my ($line, $name, $nameout);
my $gff = $ARGV[0];
my $prefix = "";
if ($gff =~ /(\S+)\.gff3/){
	$prefix = $1;
} else {die "Can't find prefix on gff3 file\n";}

my $gagaid = $ARGV[1];
my $genome = $ARGV[2];

my $outgffall = "$gagaid\_"."$prefix"."\_IR_iGluR_renamed_all.gff3"; 
my $outprotall = "$gagaid\_"."$prefix"."\_IR_iGluR_renamed_all"; 

my $outgffir = "$gagaid\_"."$prefix"."\_IR_renamed_all_complete.gff3"; 
my $outprotir = "$gagaid\_"."$prefix"."\_IR_renamed_all_complete"; 

my $outgffiglur = "$gagaid\_"."$prefix"."\_iGluR_renamed_complete.gff3"; 
my $outprotiglur = "$gagaid\_"."$prefix"."\_iGluR_renamed_complete"; 

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

	if ($subline[4] =~ /PF08395/){
		$pfam{$subline[0]}++;
	}
#	$pfam{$subline[0]} .= "$subline[4] ";

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

my %besthit;
my %dmelhit;
my %anthit;
my $blastirfile = "$prefix".".pep.fasta.IRblast.txt";
open (File, "<", $blastirfile);
while(<File>){
	chomp;
	$line = $_;
	next if ($line !~ /\S+/);	
	my @subline = split (/\t/, $line);

	next unless ($subline[10] <= "1e-5");

	if (exists $besthit{$subline[0]}){
		my @subl = split (/\t/, $besthit{$subline[0]});
		if ($subline[10] < $subl[1]){
			$besthit{$subline[0]} = "$subline[1]\t$subline[10]";
		}			
	} else {
		$besthit{$subline[0]} = "$subline[1]\t$subline[10]";
	}

	if ($subline[1] =~ /Dmel/){
#		next unless ($subline[10] <= "1e-20"); # High evalue to keep Dmel orthologs
		if (exists $dmelhit{$subline[0]}){
			my @subl = split (/\t/, $dmelhit{$subline[0]});
			if ($subline[10] < $subl[1]){
				$dmelhit{$subline[0]} = "$subline[1]\t$subline[10]";
			}			
		} else {
			$dmelhit{$subline[0]} = "$subline[1]\t$subline[10]";
		}
	}	

	if ($subline[1] =~ /Hs/ || $subline[1] =~ /Cf/){
		if (exists $anthit{$subline[0]}){
			my @subl = split (/\t/, $anthit{$subline[0]});
			if ($subline[10] < $subl[1]){
				$anthit{$subline[0]} = "$subline[1]\t$subline[10]";
			}			
		} else {
			$anthit{$subline[0]} = "$subline[1]\t$subline[10]";
		}
	}	

}
close File;


# Reading GFF to obtain gene order and exon number

my %clusterorder;
my %scaflengthsort; my %coordlengthsort;
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

#		$clusterorder{$seannameclust}{$seannameclustcoord}{$genenamesean} = $genename; 
		$clusterorder{$subline[0]}{$subline[3]}{$genename} = $genename; 

		$scaflengthsort{$genename} = $subline[0];
		$coordlengthsort{$genename} = $subline[3];


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
my $deleteproteins = "";
my $deleteproteinsl50 = "";
my %hstatus;

open (Results, ">", "$gagaid\_IR_iGluR.txt");
print Results "GAGA ID\tSpecies name\tGene name\tScaffold location\tIR/iGluR name\tIR/iGluR subfamily\tLength\tExon number\tBlast hit\tDomain\tDmel best hit\tAnt best hit\tStatus Complete Partial Pseudogene\tFragment GR\n";

open (Resultssum, ">", "$gagaid\_IR_iGluR_summary.txt");
print Resultssum "GAGA ID\tSpecies name\tIR/iGluR total number\tComplete\tComplete with domain\tPartial\tPartial with domain\tPartial longer 250aa\tPseudogene\tPseudogene with domain\tPseudogene longer 250aa\tLength lower than 150aa\tLength lower than 100aa\tLength lower than 50aa\t";
print Resultssum "NMDAR count\tKainate count\tAMPA count\tIR25a count\tIR8a count\tIR93a count\tIR68a count\tIR76b count\tIR21a count\tIR40a count\tIR75 count\t";
print Resultssum "IR complete + partial with domain\tFinal retained IRs/iGluRs\tFinal retained complete IRs/iGluRs (complete + partial >350aa)\tFinal retained complete de novo annot\tFinal retained partial IRs/iGluRs\tFinal retained partial de novo annot\tFinal retained pseudogene IRs/iGluRs\tFinal retained complete IRs\tFinal retained incomplete or pseudogene IRs\tFragment IR/iGluRs partial\tFragment IR/iGluRs pseudogenes\n";

my $irnum = "100"; # Starting to count IRs in 100
my $fragornum = "1";
my $totalor = "0"; my $com = "0"; my $comd = "0"; my $par = "0"; my $pard = "0"; my $parc = "0"; my $pse = "0"; my $psed = "0"; my $psec = "0"; my $l50 = "0"; my $l100 = "0"; my $l200 = "0"; 
my $compard = "0"; my $ftot = "0"; my $fcom = "0"; my $fpar = "0"; my $fparc = "0"; my $fcomp = "0"; my $fpse = "0"; my $fninec = "0"; my $fninei = "0"; my $frpartial = "0"; my $frpse = "0";
my $fircount = "0"; my $firpcount = "0"; 
my $nmda = "0"; my $kai = "0"; my $ampa = "0"; my $ir25a = "0"; my $ir8a = "0"; my $ir93a = "0"; my $ir68a = "0"; my $ir76b = "0"; my $ir21a = "0"; my $ir40a = "0";my $ir75 = "0";
my $denovoc = "0"; my $denovop = "0";

my %subfcount = (); my %subfcountfrag = ();

my %lsfasta;
my %lengthseq;
foreach my $prot (sort keys %fasta){
	my $length = length ($fasta{$prot});
	$lsfasta{$length}{$prot} = "$fasta{$prot}";
	$lengthseq{$prot} = $length;
}

foreach my $leng (sort { $b <=> $a } (keys %lsfasta)){
	foreach my $prot (keys %{$lsfasta{$leng}}){
#foreach my $scaf (sort (keys %clusterorder)){
#foreach my $coord (sort { $a <=> $b } (keys %{$clusterorder{$scaf}})){
#	foreach my $genename (keys %{$clusterorder{$scaf}{$coord}}){

#		my $prot = $clusterorder{$scaf}{$coord}{$genename};

		my $scaf = $scaflengthsort{$prot};
		my $coord = $coordlengthsort{$prot};

		my $len = "0";
		my $seq = "";
		if (exists $lengthseq{$prot}){
			$len = $lengthseq{$prot};
			$seq = $lsfasta{$len}{$prot};
		} else { # Avoid errors in empty genes from genewise
#			die "Error: Can't find length of $prot, located in $scaf $coord $genename\n";
			$seq = "NNNNNNNNNNNNNXNNNXNNNNNNNNN";
			$exonnumber{$prot} = 0;
		}

		$totalor++;



		my $status = "";
		my $sufix = "";
		my $discard = "No"; # Label for fragment ORs

		if ($len >= $orcompletelength){
			$status = "Complete";
			$hstatus{$prot} = "Complete";
#			$com++;
		} else {
			$status = "Partial";
			$hstatus{$prot} = "Partial";
			$sufix = "S";
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
#				if ($orcohit =~ /None/){ # Avoid deleting ORCO even if it is partial
					$deleteproteins .= "$prot ";	# Delete all proteins shorter than the specified minimum length
					$discard = "Yes";		
#				}		
			}				
		} 

		my $middleseq = substr ($seq, 5, -5);
		if ($middleseq =~ /X/){
			$status = "Pseudogene";
			$hstatus{$prot} = "Pseudogene";
			$sufix = "P";
#			$pse++;
		}


		my $irblasthit = "None";
		my $dmelblasthit = "None";
		my $antblasthit = "None";
		my $subfamily = "pIR";
		if (exists $besthit{$prot}){
			my @subl = split (/\t/, $besthit{$prot});
			$irblasthit = "$subl[0] $subl[1]"; 
		}
		if (exists $dmelhit{$prot}){
			my @subl = split (/\t/, $dmelhit{$prot});
			$dmelblasthit = "$subl[0] $subl[1]";

#my $nmda = "0"; my $kai = "0"; my $ampa = "0"; my $ir25a = "0"; my $ir8a = "0"; my $ir93a = "0"; my $ir68a = "0"; my $ir76b = "0"; my $ir21a = "0"; my $ir40a = "0";my $ir75 = "0";


			if ($subl[1] > "1e-30"){
				if (exists $anthit{$prot}){
					my @ahit = split (/\t/, $anthit{$prot});
					if ($ahit[1] < "1e-30"){
						if ($ahit[0] =~ /NMDAR1/ ){
							$subfamily = "NMDAR1";
							$nmda++;
						} elsif ($ahit[0] =~ /NMDAR2/ ){
							$subfamily = "NMDAR2";
							$nmda++;
						} elsif ($ahit[0] =~ /NMDAR3/ ){
							$subfamily = "NMDAR3";
							$nmda++;
						} elsif ($ahit[0] =~ /IR25a/ ){
							$subfamily = "IR25a";
							$ir25a++;
						} elsif ($ahit[0] =~ /IR8a/ ){
							$subfamily = "IR8a";
							$ir8a++;
						} elsif ($ahit[0] =~ /IR93a/ ){
							$subfamily = "IR93a";
							$ir93a++;
						} elsif ($ahit[0] =~ /IR68a/ ){
							$subfamily = "IR68a";
							$ir68a++;
						} elsif ($ahit[0] =~ /IR76b/ ){
							$subfamily = "IR76b";
							$ir76b++;
						} elsif ($ahit[0] =~ /IR21a/ ){
							$subfamily = "IR21a";
							$ir21a++;
						} elsif ($ahit[0] =~ /IR40a/ ){
							$subfamily = "IR40a";
							$ir40a++;
						} elsif ($ahit[0] =~ /IR75/ ){
							$subfamily = "IR75";
							$ir75++;
						} elsif ($ahit[0] =~ /IR/ ){
							$subfamily = "IR";
						}						

					} else {
						if ($ahit[1] > $subl[1]){ # Dmel shows better evalue
							if ($dmelblasthit =~ /NMDAR1/ ){
								$subfamily = "NMDAR1";
#								if ($subl[1] > "1e-100"){ $subfamily = "NMDAR3"; } # Assign as NMDAR3 those not hitting any Dmel NMDAR1 or 2 with eval < 10e-100
								$nmda++;
							} elsif ($dmelblasthit =~ /NMDAR2/ ){
								$subfamily = "NMDAR2";
#								if ($subl[1] > "1e-100"){ $subfamily = "NMDAR3"; } # Assign as NMDAR3 those not hitting any Dmel NMDAR1 or 2 with eval < 10e-100
								$nmda++;
							} elsif ($dmelblasthit =~ /AMPA/ ){
								$subfamily = "AMPA";
								$ampa++;
							} elsif ($dmelblasthit =~ /Kaina/ ){
								$subfamily = "KA";
								$kai++;							
							} else {
								$subfamily = "IR";
							}
						} else { # Ant hit shows better evalue
							if ($ahit[0] =~ /NMDAR1/ ){
								$subfamily = "NMDAR1";
								$nmda++;
							} elsif ($ahit[0] =~ /NMDAR2/ ){
								$subfamily = "NMDAR2";
								$nmda++;
							} elsif ($ahit[0] =~ /NMDAR3/ ){
								$subfamily = "NMDAR3";
								$nmda++;
							} elsif ($ahit[0] =~ /AMPA/ ){
								$subfamily = "AMPA";
								$ampa++;
							} elsif ($ahit[0] =~ /Kaina/ ){
								$subfamily = "KA";
								$kai++;							
							} else {
								$subfamily = "IR";
							}
						}
					}
				} else {
						if ($dmelblasthit =~ /NMDAR1/ ){
							$subfamily = "NMDAR1";
#							if ($subl[1] > "1e-100"){ $subfamily = "NMDAR3"; } # Assign as NMDAR3 those not hitting any Dmel NMDAR1 or 2 with eval < 10e-100
							$nmda++;
						} elsif ($dmelblasthit =~ /NMDAR2/ ){
							$subfamily = "NMDAR2";
#							if ($subl[1] > "1e-100"){ $subfamily = "NMDAR3"; } # Assign as NMDAR3 those not hitting any Dmel NMDAR1 or 2 with eval < 10e-100
							$nmda++;
						} elsif ($dmelblasthit =~ /AMPA/ ){
							$subfamily = "AMPA";
							$ampa++;
						} elsif ($dmelblasthit =~ /Kaina/ ){
							$subfamily = "KA";
							$kai++;							
						} else {
							$subfamily = "IR";
						}

				}
			} else {

				if ($dmelblasthit =~ /NMDAR/ ){
#					$subfamily = "NMDAR1";
					$nmda++;
					if (exists $anthit{$prot}){
						my @ahit = split (/\t/, $anthit{$prot});
						if ($ahit[1] > $subl[1]){
							if ($dmelblasthit =~ /NMDAR2/ ){
								$subfamily = "NMDAR2";
							} elsif ($dmelblasthit =~ /NMDAR1/ ){
								$subfamily = "NMDAR1";
							}							
						} else {
							if ($ahit[0] =~ /NMDAR2/ ){
								$subfamily = "NMDAR2";
							} elsif ($ahit[0] =~ /NMDAR1/ ){
								$subfamily = "NMDAR1";
							} elsif ($ahit[0] =~ /NMDAR3/ ){
								$subfamily = "NMDAR3";
							}								
						}
					} else {
						if ($dmelblasthit =~ /NMDAR2/ ){
							$subfamily = "NMDAR2";
						} elsif ($dmelblasthit =~ /NMDAR1/ ){
							$subfamily = "NMDAR1";
						}					
					}
#					if ($subl[1] > "1e-100"){ $subfamily = "NMDAR3"; } # Assign as NMDAR3 those not hitting any Dmel NMDAR1 or 2 with eval < 10e-100
#					$nmda++;
#				} elsif ($dmelblasthit =~ /NMDAR2/ ){
#					$subfamily = "NMDAR2";
#					if ($subl[1] > "1e-100"){ $subfamily = "NMDAR3"; } # Assign as NMDAR3 those not hitting any Dmel NMDAR1 or 2 with eval < 10e-100
#					$nmda++;
				} elsif ($dmelblasthit =~ /AMPA/ ){
					$subfamily = "AMPA";
					$ampa++;
				} elsif ($dmelblasthit =~ /Kaina/ ){
					$subfamily = "KA";
					$kai++;
				} elsif ($dmelblasthit =~ /IR25a/ ){
					$subfamily = "IR25a";
					$ir25a++;
				} elsif ($dmelblasthit =~ /IR8a/ ){
					$subfamily = "IR8a";
					$ir8a++;
				} elsif ($dmelblasthit =~ /IR93a/ ){
					if ($subl[1] > "1e-20"){
						$subfamily = "IR";
					} else {
						$subfamily = "IR93a";
						$ir93a++;
					}
				} elsif ($dmelblasthit =~ /IR68a/ ){
					if ($subl[1] > "1e-20"){
						$subfamily = "IR";
					} else {
						$subfamily = "IR68a";
						$ir68a++;
					}
				} elsif ($dmelblasthit =~ /IR76b/ ){
					if ($subl[1] > "1e-20"){
						$subfamily = "IR";
					} else {
						$subfamily = "IR76b";
						$ir76b++;
					}
				} elsif ($dmelblasthit =~ /IR21a/ ){
					if ($subl[1] > "1e-20"){
						$subfamily = "IR";
					} else {
						$subfamily = "IR21a";
						$ir21a++;
					}
				} elsif ($dmelblasthit =~ /IR40a/ ){
					if ($subl[1] > "1e-20"){
						$subfamily = "IR";
					} else {
						$subfamily = "IR40a";
						$ir40a++;
					}
				} elsif ($dmelblasthit =~ /IR75/ ){
					if ($subl[1] > "1e-20"){
						$subfamily = "IR";
					} else {
						$subfamily = "IR75";
						$ir75++;
					}
				} elsif ($dmelblasthit =~ /IR/ ){
					$subfamily = "IR";
				}
			}
		} 
		if (exists $anthit{$prot}){
			my @subl = split (/\t/, $anthit{$prot});
			$antblasthit = "$subl[0] $subl[1]";
			unless (exists $dmelhit{$prot}){
				if ($subl[1] <= "1e-30"){
					$subfamily = "IR";
				}
			}
		}



		my $ordomain = "None";
		if (exists $interpro{$prot}){

			my $ligchandomain = "No";
			if ($interpro{$prot} =~ /IPR001320/ && $interpro{$prot} =~ /IPR019594/ && $interpro{$prot} =~ /IPR001828/){
				$ligchandomain = "Yes";
				$ordomain = "ANF,Lig_chan,Lig_chan2";
				if ($subfamily =~ /pIR/){
					$subfamily = "piGluR"; # Re-assigning family with the domain info
				}
			} elsif ($interpro{$prot} =~ /IPR001320/ && $interpro{$prot} =~ /IPR001828/){
				$ligchandomain = "Yes";
				$ordomain = "ANF,Lig_chan";
				if ($subfamily =~ /pIR/){
					$subfamily = "piGluR"; # Re-assigning family with the domain info
				}
			} elsif ($interpro{$prot} =~ /IPR019594/ && $interpro{$prot} =~ /IPR001828/){
				$ligchandomain = "Yes";
				$ordomain = "ANF,Lig_chan2";
				if ($subfamily =~ /pIR/){
					$subfamily = "piGluR"; # Re-assigning family with the domain info
				}
			} elsif ($interpro{$prot} =~ /IPR001320/ && $interpro{$prot} =~ /IPR019594/){
				$ligchandomain = "Yes";
				$ordomain = "Lig_chan,Lig_chan2";
			} elsif ($interpro{$prot} =~ /IPR001320/){
				$ligchandomain = "Yes";
				$ordomain = "Lig_chan";
			} elsif ($interpro{$prot} =~ /IPR019594/){
				$ligchandomain = "Yes";
				$ordomain = "Lig_chan2";
			} elsif ($interpro{$prot} =~ /IPR001828/){
				$ligchandomain = "No";
				$ordomain = "ANF";
				if ($subfamily =~ /pIR/){
					$subfamily = "piGluR"; # Re-assigning family with the domain info
				}
			} else {
				print "Warning: $prot in $gagaid with $len length contains an unexpected domain $interpro{$prot}\n";
				$ordomain = "Other";
				if ($status =~ /Partial/ || $status =~ /Pseudogene/){
					$deleteproteins .= "$prot "; ### Delete protein identified with other domain, and are partial or pseudogenes
					$discard = "Yes";
				} elsif ($status =~ /Complete/ && $irblasthit =~ /None/){
					$deleteproteins .= "$prot "; ### Delete complete protein identified with other domain, and does not hit any IR
					$discard = "Yes";
				}

			}


			if ($ligchandomain =~ /Yes/){
				if ($status =~ /Pseudogene/){
					$psed++;
				} elsif ($status =~ /Complete/){
					$comd++;
				} elsif ($status =~ /Partial/){
					$pard++;
					$compard++;
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
#				if ($orcohit =~ /None/){ # Avoid deleting ORCO even if it is partial and have no domain
					$deleteproteins .= "$prot "; ### Delete protein partial without domain, and length < specified in minorlength
					$discard = "Yes";		
#				}	
			}
		}



		# Get new names
		my $nname = "";
		if ($discard =~ /Yes/){
			#$nname = "Discarded";
			#if ($len < 50){
			#	$nname = "Discarded";
			if ($deleteproteinsl50 =~ /$prot /){
				$nname = "Discarded";
			} else {
				$subfcountfrag{$subfamily}++;
				$subfcount{$subfamily}++;
				if ($subfamily =~ /IR$/){
					$nname = "$gagashort{$gagaid}"."frag$subfamily$irnum"."$sufix";
					$irnum++;
#				} elsif ($subfcountfrag{$subfamily} == 1){
#					$nname = "$gagashort{$gagaid}"."frag$subfamily"."$sufix";
				} else {
#					$nname = "$gagashort{$gagaid}"."frag$subfamily"."\.$subfcountfrag{$subfamily}"."$sufix";
					$nname = "$gagashort{$gagaid}"."frag$subfamily"."\.$subfcount{$subfamily}"."$sufix"; # Change to number it for all complete and partial
				}
#				$nname = "$gagashort{$gagaid}"."frag$subfamily"."$fragornum"."$sufix";
#				$fragornum++;
			}
		} else {
			$subfcount{$subfamily}++;

			if ($subfamily =~ /IR$/){
				$nname = "$gagashort{$gagaid}"."$subfamily$irnum"."$sufix";
				$irnum++;
			} elsif ($subfcount{$subfamily} == 1){
				$nname = "$gagashort{$gagaid}"."$subfamily"."$sufix";
			} else {
				$nname = "$gagashort{$gagaid}"."$subfamily"."\.$subfcount{$subfamily}"."$sufix";
			}
			if ($len >= $maxorlength){ # Print a warning in retained GRs longer than the max specified length
				print "Warning: $prot in $gagaid with length $len is longer than the max specified length\n";
			}
		}
		$newnames{$prot} = "$nname";


#print Results "GAGA ID\tSpecies name\tGene name\tScaffold location\tIR/iGluR name\tIR/iGluR subfamily\tLength\tExon number\tBlast hit\tDomain\tDmel best hit\tAnt best hit\tStatus Complete Partial Pseudogene\tFragment GR\n";
		print Results "$gagaid\t$gagasp{$gagaid}\t$prot\t$scaf\-$coord\t$nname\t$subfamily\t$len\t$exonnumber{$prot}\t$irblasthit\t$ordomain\t$dmelblasthit\t$antblasthit\t$status\t$discard\n";

	}

#}}
}
close Results;

# Reading GFF and generating final file filtering short proteins and renaming ORs

open (Resultsgffir, ">", "$outgffir");
open (Resultsgffiglur, ">", "$outgffiglur");
open (Resultsgffall, ">", "$outgffall");

open (GFFfile , "<", $gff); 
while (<GFFfile>){
	chomp;
	$line = $_;
	next if ($line !~ /\S+/);

	if ($line =~ /^#/){
		print Resultsgffall "$line\n";
		next;
	}
	
	my @subline = split (/\t/, $line);

	if ($subline[2] =~ /CDS/){
		my $genename = "";
		my $rest = "";
		if ($subline[8] =~ /Parent=([^;]+)(\S*)/){
			$genename = $1;
			$rest = $2;
		}
		else {die "ERROR in run_OR_classification.pl: It fails detecting Parent ID in $line\n";}

		my $nnamef = "";
		if (exists $newnames{$genename}){
			$nnamef = $newnames{$genename};
		} else {
			die "Can't find gene new name for $genename in $line\n";
		}

# 		Keeping now all genes as fragment ORs; and creating multiple GFFs
#		next if ($deleteproteins =~ /$genename /); ### Delete genes that do not pass the filter
		next if ($deleteproteinsl50 =~ /$genename /); ### Delete from GFF all genes shorter than 50aa or mis-identification such as GR

		print Resultsgffall "$subline[0]\t$subline[1]\t$subline[2]\t$subline[3]\t$subline[4]\t$subline[5]\t$subline[6]\t$subline[7]\tParent=$nnamef"."$rest\;Status=$hstatus{$genename}\n";

		unless ($deleteproteins =~ /$genename /){
#			print Resultsgff "$subline[0]\t$subline[1]\t$subline[2]\t$subline[3]\t$subline[4]\t$subline[5]\t$subline[6]\t$subline[7]\tParent=$nnamef"."$rest\;Status=$hstatus{$genename}\n";

			if ($hstatus{$genename} =~ /Complete/){
#				print Resultsgffcomp "$subline[0]\t$subline[1]\t$subline[2]\t$subline[3]\t$subline[4]\t$subline[5]\t$subline[6]\t$subline[7]\tParent=$nnamef"."$rest\;Status=$hstatus{$genename}\n";

				if ($nnamef =~ /$gagashort{$gagaid}\S*IR/){
					print Resultsgffir "$subline[0]\t$subline[1]\t$subline[2]\t$subline[3]\t$subline[4]\t$subline[5]\t$subline[6]\t$subline[7]\tParent=$nnamef"."$rest\;Status=$hstatus{$genename}\n";
				} else {
					print Resultsgffiglur "$subline[0]\t$subline[1]\t$subline[2]\t$subline[3]\t$subline[4]\t$subline[5]\t$subline[6]\t$subline[7]\tParent=$nnamef"."$rest\;Status=$hstatus{$genename}\n";				
				}

			}
		}
		
	}
	elsif ($subline[2] =~ /mRNA/){
		my $genename = "";
		my $rest = "";
		if ($subline[8] =~ /ID=([^;]+)\;Parent=[^;]+(\;\S*)/){
		#if ($subline[8] =~ /transcript_id \"(\S+)\"/){
			$genename = $1;
			$rest = $2;
		}
		else {print "ERROR in run_OR_classification.pl: It fails detecting ID in $line\n";}

		my $nnamef = "";
		if (exists $newnames{$genename}){
			$nnamef = $newnames{$genename};
		} else {
			die "Can't find gene new name for $genename in $line\n";
		}


# 		Keeping now all genes as fragment ORs; and creating multiple GFFs
#		next if ($deleteproteins =~ /$genename /); ### Delete genes that do not pass the filter
		next if ($deleteproteinsl50 =~ /$genename /); ### Delete from GFF all genes shorter than 50aa or mis-identification such as GR

		print Resultsgffall "$subline[0]\t$subline[1]\tgene\t$subline[3]\t$subline[4]\t$subline[5]\t$subline[6]\t$subline[7]\tID=g"."$nnamef"."$rest;Status=$hstatus{$genename}\n";		
		print Resultsgffall "$subline[0]\t$subline[1]\t$subline[2]\t$subline[3]\t$subline[4]\t$subline[5]\t$subline[6]\t$subline[7]\tID=$nnamef\;Parent=g"."$nnamef"."$rest;Status=$hstatus{$genename}\n";		

		unless ($deleteproteins =~ /$genename /){
#			print Resultsgff "$subline[0]\t$subline[1]\tgene\t$subline[3]\t$subline[4]\t$subline[5]\t$subline[6]\t$subline[7]\tID=g"."$nnamef"."$rest;Status=$hstatus{$genename}\n";		
#			print Resultsgff "$subline[0]\t$subline[1]\t$subline[2]\t$subline[3]\t$subline[4]\t$subline[5]\t$subline[6]\t$subline[7]\tID=$nnamef\;Parent=g"."$nnamef"."$rest;Status=$hstatus{$genename}\n";		

			if ($hstatus{$genename} =~ /Complete/){
#				print Resultsgffcomp "$subline[0]\t$subline[1]\tgene\t$subline[3]\t$subline[4]\t$subline[5]\t$subline[6]\t$subline[7]\tID=g"."$nnamef"."$rest;Status=$hstatus{$genename}\n";		
#				print Resultsgffcomp "$subline[0]\t$subline[1]\t$subline[2]\t$subline[3]\t$subline[4]\t$subline[5]\t$subline[6]\t$subline[7]\tID=$nnamef\;Parent=g"."$nnamef"."$rest;Status=$hstatus{$genename}\n";		

				if ($nnamef =~ /$gagashort{$gagaid}\S*IR/){
					print Resultsgffir "$subline[0]\t$subline[1]\tgene\t$subline[3]\t$subline[4]\t$subline[5]\t$subline[6]\t$subline[7]\tID=g"."$nnamef"."$rest;Status=$hstatus{$genename}\n";		
					print Resultsgffir "$subline[0]\t$subline[1]\t$subline[2]\t$subline[3]\t$subline[4]\t$subline[5]\t$subline[6]\t$subline[7]\tID=$nnamef\;Parent=g"."$nnamef"."$rest;Status=$hstatus{$genename}\n";		
				} else {
					print Resultsgffiglur "$subline[0]\t$subline[1]\tgene\t$subline[3]\t$subline[4]\t$subline[5]\t$subline[6]\t$subline[7]\tID=g"."$nnamef"."$rest;Status=$hstatus{$genename}\n";		
					print Resultsgffiglur "$subline[0]\t$subline[1]\t$subline[2]\t$subline[3]\t$subline[4]\t$subline[5]\t$subline[6]\t$subline[7]\tID=$nnamef\;Parent=g"."$nnamef"."$rest;Status=$hstatus{$genename}\n";							
				}

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
				if ($genename =~ /IR_iGluR.*\.t/){
					$denovoc++; # Counting new annotated genes in bitacora gemoma
				}
				if ($nnamef =~ /$gagashort{$gagaid}\S*IR/){
					$fircount++;
				}
			} elsif ($hstatus{$genename} =~ /Partial/){
				$fpar++;
				if ($genename =~ /IR_iGluR.*\.t/){
					$denovop++; # Counting new annotated genes in bitacora gemoma
				}
				if ($lengthseq{$genename} >= $orcompletelength){ # It will not do nothing here as I relabelled those as complete in the status
					$fparc++;
					$fcomp++;
				}
				if ($nnamef =~ /$gagashort{$gagaid}\S*IR/){
					$firpcount++;
				}
			} elsif ($hstatus{$genename} =~ /Pseudogene/){
				$fpse++;
				if ($nnamef = /$gagashort{$gagaid}\S*IR/){
					$firpcount++;
				}
			}
		}

	}
}
close GFFfile;


#Encode proteins and CDS from the generated GFFs

system ("perl /home/projects/ku_00039/people/joeviz/programs/bitacora/Scripts/Tools/gff2fasta_v3.pl $genome $outgffir $outprotir ");
system ("sed \'s\/X\*\$\/\/\' $outprotir\.pep.fasta > $outprotir\.pep.fasta.tmp");
system ("mv $outprotir\.pep.fasta.tmp $outprotir\.pep.fasta");

system ("perl /home/projects/ku_00039/people/joeviz/programs/bitacora/Scripts/Tools/gff2fasta_v3.pl $genome $outgffall $outprotall ");
system ("sed \'s\/X\*\$\/\/\' $outprotall\.pep.fasta > $outprotall\.pep.fasta.tmp");
system ("mv $outprotall\.pep.fasta.tmp $outprotall\.pep.fasta");

system ("perl /home/projects/ku_00039/people/joeviz/programs/bitacora/Scripts/Tools/gff2fasta_v3.pl $genome $outgffiglur $outprotiglur ");
system ("sed \'s\/X\*\$\/\/\' $outprotiglur\.pep.fasta > $outprotiglur\.pep.fasta.tmp");
system ("mv $outprotiglur\.pep.fasta.tmp $outprotiglur\.pep.fasta");

#print Resultssum "$gagaid\t$gagasp{$gagaid}\t$totalor\t$com\t$comd\t$par\t$pard\t$parc\t$pse\t$psed\t$psec\t$l200\t$l100\t$l50\t$orconum\t$orconumcom\t$ninec\t$ninei\t$nineic\t$ninep\t$compard\t$ftot\t$fcom\t$fpar\t$fparc\t$fcomp\t$fpse\t$fninec\t$fninei\n";
print Resultssum "$gagaid\t$gagasp{$gagaid}\t$totalor\t$com\t$comd\t$par\t$pard\t$parc\t$pse\t$psed\t$psec\t$l200\t$l100\t$l50\t";
print Resultssum "$nmda\t$kai\t$ampa\t$ir25a\t$ir8a\t$ir93a\t$ir68a\t$ir76b\t$ir21a\t$ir40a\t$ir75\t";
print Resultssum "$compard\t$ftot\t$fcomp\t$denovoc\t$fpar\t$denovop\t$fpse\t$fircount\t$firpcount\t$frpartial\t$frpse\n";
close Resultssum;


#print Resultssum "GAGA ID\tSpecies name\tIR/iGluR total number\tComplete\tComplete with domain\tPartial\tPartial with domain\tPartial longer 250aa\tPseudogene\tPseudogene with domain\tPseudogene longer 250aa\tLength lower than 150aa\tLength lower than 100aa\tLength lower than 50aa\t";
#print Resultssum "NMDAR count\tKainate count\tAMPA count\tIR25a count\tIR8a count\tIR93a count\tIR68a count\tIR76b count\tIR21a count\tIR40a count\tIR75 count\t";
#print Resultssum "IR complete + partial with domain\tFinal retained IRs/iGluRs\tFinal retained complete IRs/iGluRs (complete + partial >350aa)\tFinal retained complete de novo annot\tFinal retained partial IRs/iGluRs\tFinal retained partial de novo annot\tFinal retained pseudogene IRs/iGluRs\tFinal retained complete IRs\tFinal retained incomplete or pseudogene IRs\tFragment IR/iGluRs partial\tFragment IR/iGluRs pseudogenes\n";

