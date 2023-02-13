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

my $orcompletelength = "350"; # Length to consider complete GR copies
my $minorlength = "200"; # Minimum length to retain partial GRs
my $maxorlength = "535"; # Maximum average length, it will print a warning with the sequences larger than this length

# Seqs shorter than "minorlength" but longer than 50aa will be labelled as fragment. 

my ($line, $name, $nameout);
my $gff = $ARGV[0];
my $prefix = "";
if ($gff =~ /(\S+)\.gff3/){
	$prefix = $1;
} else {die "Can't find prefix on gff3 file\n";}

my $gagaid = $ARGV[1];
my $genome = $ARGV[2];

my $outgffall = "$gagaid\_"."$prefix"."\_GR_renamed_all.gff3"; 
my $outprotall = "$gagaid\_"."$prefix"."\_GR_renamed_all"; 

my $outgff = "$gagaid\_"."$prefix"."\_GR_renamed_all_nofragment.gff3"; 
my $outprot = "$gagaid\_"."$prefix"."\_GR_renamed_all_nofragment"; 

my $outgffcomp = "$gagaid\_"."$prefix"."\_GR_renamed_complete.gff3"; 
my $outprotcomp = "$gagaid\_"."$prefix"."\_GR_renamed_complete"; 

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

my %orco;
my $blastorcofile = "$prefix".".pep.fasta.ORcoblast.txt";
open (File, "<", $blastorcofile);
while(<File>){
	chomp;
	$line = $_;
	next if ($line !~ /\S+/);	
	my @subline = split (/\t/, $line);

	if ($subline[10] <= "1e-100"){
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

	if ($subline[10] <= "1e-5"){
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

my %dmelblast; my %dmelblasteval;
my $blastdmelfile = "$prefix".".pep.fasta.GRdmelblast.txt";
open (File, "<", $blastdmelfile);
while(<File>){
	chomp;
	$line = $_;
	next if ($line !~ /\S+/);	
	my @subline = split (/\t/, $line);

	if ($subline[10] <= "1e-40"){
		if (exists $dmelblast{$subline[0]}){
			if ($subline[10] < $dmelblasteval{$subline[0]}){
				$dmelblasteval{$subline[0]} = "$subline[10]";
				$dmelblast{$subline[0]} = "$subline[1]";
			}
		} else {
			$dmelblasteval{$subline[0]} = "$subline[10]";
			$dmelblast{$subline[0]} = "$subline[1]";
		}

	}
}
close File;

my %blast; my %length;
my %grblast;
my $blastgrfile = "$prefix".".pep.fasta.GRblast.txt";
open (File, "<", $blastgrfile);
while(<File>){
	chomp;
	$line = $_;
	next if ($line !~ /\S+/);	
	my @subline = split (/\t/, $line);

	if ($subline[10] <= "1e-5"){
		if (exists $grblast{$subline[0]}){
			if ($subline[10] < $grblast{$subline[0]}){
				$grblast{$subline[0]} = "$subline[10]";
			}
		} else {
			$grblast{$subline[0]} = "$subline[10]";
		}

		push (@{$blast{$subline[0]}}, join("\t",$subline[1], $subline[10], $subline[3], $subline[12], $subline[13], $subline[6], $subline[7], $subline[2])); # New line to save the whole blast results and explore for chimeric genes
		$length{$subline[0]} = $subline[12];
	}
}
close File;


### Finding chimeric/fused genes

my $minlengthcut = "30"; ## Minimum positions required to trim a protein (i.e. blast hits starting in position 10 will report the full sequence instead of trimming the first 10 aa)

foreach my $key (sort keys %blast) {
	my $blasthit2="";
	my $hitlvl = "0";
	my (@ini, @fin);
	$ini[0] = "99999999999999999999";
	$fin[0] = "1";
	foreach my $blastresult (@{$blast{$key}}) {
		my @subline = ();
		@subline = split (/\t/, $blastresult);
		my $hitblast = $key;
		my $filtro1 = ($subline[4]*2)/3;
		my $filtro2 = ($subline[3]*0.8);

#### NO FILTERING HERE		
#		if ($subline[2] < 51) { # If the alignment is lower than 50aa, and smaller than 2/3 QUERY length protein used, it should contain a similarity higher than 80%. If not it is removed as false positive (small domains hitting non-related proteins)
#			unless ($subline[2] >= $filtro1){
#				next if ($subline[7] < 80);
#			}
#		}

#		if ($subline[2] >= $filtro1 || $subline[2] >= $filtro2) { # BLAST filtering: alignment covering 2/3 of subject, or 80% of query
			$hitlvl = "2";
			my $n = 0;
			my $extrahit = 0;
			foreach my $i (@ini){
				my $f = $fin[$n];
				my $f2 = $f + 10;
				my $i2 = $i - 10;

				if ($i >= int($subline[5]) && $f <= int($subline[6])){
					$ini[$n] = int($subline[5]);
					$fin[$n] = int($subline[6]);
				}
				elsif ($i <= int($subline[5]) && $f < int($subline[6]) && $f2 >= int($subline[5])){
					$fin[$n] = int($subline[6]);
				}
				elsif ($i > int($subline[5]) && $f >= int($subline[6]) && $i2 <= int($subline[6])){
					$ini[$n] = int($subline[5]);
				}
				elsif ($f2 < int($subline[5]) || $i2 > int($subline[6])) {
					$extrahit++;
				}

				$n++;
			}

			if ($extrahit >= $n) {
					$ini[$n] = int($subline[5]);
					$fin[$n] = int($subline[6]);				
			}
#		}
	}

	if ($hitlvl == 2){
		my $hits = scalar(@ini);

		if ($hits == 1){

			# Length filter to avoid exluding a few ($minlengthcut) initial or end positions
			my $ipos = $ini[0];
			my $fpos = $fin[0];
			if ($ipos <= $minlengthcut){ # Initial position
				$ipos = 1;
			}
			my $filterend = $length{$key} - $minlengthcut;
			if ($fpos >= $filterend){ # Initial position
				$fpos = $length{$key};
			}

#			print Results "$key annot $ipos $fpos blastp\n";
		}
		else {
			my $n = 0;

			print "Warning $key in $gagaid with length $length{$key} should be split in $hits parts\n";
			
			foreach my $i (@ini){
				my $f = $fin[$n];

				# Length filter to avoid exluding a few ($minlengthcut) initial or end positions
				my $ipos = $ini[$n];
				my $fpos = $fin[$n];
				if ($ipos <= $minlengthcut){ # Initial position
					$ipos = 1;
				}
				my $filterend = $length{$key} - $minlengthcut;
				if ($fpos >= $filterend){ # Initial position
					$fpos = $length{$key};
				}

				my $nn= $n+1;
#				print Results "$key\_split$nn annot $ipos $fpos blastp\n";
				$n++;

			}
		}
	}

}

###


# Reading GFF to obtain gene order and exon number

my %clusterorder;
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

open (Results, ">", "$gagaid\_GRs.txt");
print Results "GAGA ID\tSpecies name\tGene name\tScaffold location\tGR name\tLength\tExon number\tBlast hit\tDomain\tDmel function\tStatus Complete Partial Pseudogene\tFragment GR\n";

open (Resultssum, ">", "$gagaid\_GRs_summary.txt");
print Resultssum "GAGA ID\tSpecies name\tGR total number\tComplete\tComplete with domain\tPartial\tPartial with domain\tPartial longer 350aa\tPseudogene\tPseudogene with domain\tPseudogene longer 350aa\tLength lower than 200aa\tLength lower than 100aa\tLength lower than 50aa\tOR count\tOrco\tOrco complete\tSingle-exon GRs complete\tSingle-exon GRs incomplete\tSingle-exon GRs incomplete (length > 350aa)\tSingle-exon GRs pseudogenes\tDmel sugar receptors\tDmel bitter receptors\tDmel CO2 receptors\tDmel fructose receptors\tGR complete + partial with domain\tFinal retained GRs\tFinal retained complete GRs (complete + partial >350aa)\tFinal retained complete de novo annot\tFinal retained partial GRs\tFinal retained partial de novo annot\tFinal retained pseudogene GRs\tFinal retained complete 1-exon GRs\tFinal retained incomplete or pseudogene 1-exon GRs\tFragment GRs partial\tFragment GRs pseudogenes\n";

my $ornum = "1";
my $fragornum = "1";
my $totalor = "0"; my $com = "0"; my $comd = "0"; my $par = "0"; my $pard = "0"; my $parc = "0"; my $pse = "0"; my $psed = "0"; my $psec = "0"; my $l50 = "0"; my $l100 = "0"; my $l200 = "0"; my $orconum = "0"; my $orconumcom = "0"; my $ninec = "0"; my $ninei = "0"; my $nineic = "0"; my $ninep = "0";
my $compard = "0"; my $ftot = "0"; my $fcom = "0"; my $fpar = "0"; my $fparc = "0"; my $fcomp = "0"; my $fpse = "0"; my $fninec = "0"; my $fninei = "0"; my $frpartial = "0"; my $frpse = "0";
my $grcount = "0";
my $bitter = "0"; my $sugar = "0"; my $co2 = "0"; my $gr43 = "0";
my $denovoc = "0"; my $denovop = "0";

my %lsfasta;
my %lengthseq;
foreach my $prot (sort keys %fasta){
	my $length = length ($fasta{$prot});
	$lsfasta{$length}{$prot} = "$fasta{$prot}";
	$lengthseq{$prot} = $length;
}

#foreach my $len (sort { $b <=> $a } (keys %lsfasta)){
#	foreach my $prot (keys %{$lsfasta{$len}}){
foreach my $scaf (sort (keys %clusterorder)){
foreach my $coord (sort { $a <=> $b } (keys %{$clusterorder{$scaf}})){
	foreach my $genename (keys %{$clusterorder{$scaf}{$coord}}){

		my $prot = $clusterorder{$scaf}{$coord}{$genename};

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
			$sufix = "P";
#			$pse++;
		}


		my $orhit = "None";
		if (exists $orblast{$prot} && exists $grblast{$prot} ){
			if ($orblast{$prot} <= $grblast{$prot}){
				$orhit = "OR";
			} elsif ($orblast{$prot} > $grblast{$prot}){
				$orhit = "GR";
#				$deleteproteinsl50 .= "$prot "; ### Delete protein identified as GR
#				$discard = "Yes";
#				$grcount++;
			}  
		} elsif (exists $orblast{$prot}){
			$orhit = "OR";
		} elsif (exists $grblast{$prot}){
			$orhit = "GR";
#			$deleteproteinsl50 .= "$prot "; ### Delete protein identified as GR
#			$discard = "Yes";
#			$grcount++;
		}

		my $ordomain = "None";
		if (exists $interpro{$prot}){
			if ($interpro{$prot} =~ /IPR013604/ || $interpro{$prot} =~ /IPR009318/){
				$ordomain = "GR";
				if ($status =~ /Pseudogene/){
					$psed++;
				} elsif ($status =~ /Complete/){
					$comd++;
				} elsif ($status =~ /Partial/){
					$pard++;
					$compard++;
				}

			} elsif ($interpro{$prot} =~ /IPR004117/){
				$ordomain = "OR";
#				if ($orhit !~ /OR/){ # There are OR that have the Gr domain but are still ORs, so using this to avoid deleting good proteins
#					$deleteproteinsl50 .= "$prot "; ### Delete protein identified as GR
#					$discard = "Yes";					
#				}
			} elsif ($interpro{$prot} =~ /IPR001320/ || $interpro{$prot} =~ /IPR019594/ || $interpro{$prot} =~ /IPR001828/){
				$ordomain = "IR";
				$deleteproteinsl50 .= "$prot "; ### Delete protein identified as IR-iGluR
				$discard = "Yes";
			} elsif ($interpro{$prot} =~ /IPR005055/){
				$ordomain = "CSP";
				$deleteproteinsl50 .= "$prot "; ### Delete protein identified as CSP
				$discard = "Yes";
			} elsif ($interpro{$prot} =~ /IPR006170/){
				$ordomain = "OBP";
				$deleteproteinsl50 .= "$prot "; ### Delete protein identified as OBP
				$discard = "Yes";
			} elsif ($interpro{$prot} =~ /IPR003172/){
				$ordomain = "NPC2";
				$deleteproteinsl50 .= "$prot "; ### Delete protein identified as NPC2
				$discard = "Yes";
			} elsif ($interpro{$prot} =~ /IPR002159/){
				$ordomain = "CD36";
				$deleteproteinsl50 .= "$prot "; ### Delete protein identified as CD36-SNMP
				$discard = "Yes";
			} elsif ($interpro{$prot} =~ /IPR001873/){
				$ordomain = "PPK";
				$deleteproteinsl50 .= "$prot "; ### Delete protein identified as PPK
				$discard = "Yes";
			} else {
				print "Warning: $prot in $gagaid with $len length contains an unexpected domain $interpro{$prot}\n";
				$ordomain = "Other";
				if ($status =~ /Partial/ || $status =~ /Pseudogene/){
					$deleteproteins .= "$prot "; ### Delete protein identified with other domain, and are partial or pseudogenes
					$discard = "Yes";
				} elsif ($status =~ /Complete/ && $orhit !~ /GR/){
					$deleteproteins .= "$prot "; ### Delete complete protein identified with other domain, and does not hit any OR
					$discard = "Yes";
				}

			}
		}
		if (exists $pfam{$prot}){ # Print warning in proteins with the same pfam domain repeated more than once: putative fused genes
			if ($pfam{$prot} > 1){
				print "Warning: $prot in $gagaid with length $len contains $pfam{$prot} GR pfam domains\n";
			}
		}

		if ($orhit =~ /OR/){
			$deleteproteinsl50 .= "$prot "; ### Delete protein identified as OR
			$discard = "Yes";
			$grcount++;
		} elsif ($ordomain =~ /OR/){
			if ($orhit !~ /GR/){
				$deleteproteinsl50 .= "$prot "; ### Delete protein identified as OR
				$discard = "Yes";
				$grcount++;
			}
			# If blast hits OR, it is ok as some domain are mis-identified
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

		# Count 9 exon ORs || Single exon GRs 
		if ($exonnumber{$prot} == 1){
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
			if ($exonnumber{$prot} == 1){
				if ($status =~ /Complete/){
					$fninec++;
				} else {
					$fninei++;
				}	
			}
		}	


		# Find Dmel function
		my $dmelhit = "None";
		if (exists $dmelblast{$prot}){
			if ($dmelblast{$prot} =~ /Gr64a/ || $dmelblast{$prot} =~ /Gr64b/ || $dmelblast{$prot} =~ /Gr64c/ || $dmelblast{$prot} =~ /Gr64d/ || $dmelblast{$prot} =~ /Gr64e/ || $dmelblast{$prot} =~ /Gr64f/ || $dmelblast{$prot} =~ /Gr5a/ || $dmelblast{$prot} =~ /Gr61a/){
				$dmelhit = "$dmelblast{$prot} Sugar";
				$sugar++;
			} elsif ($dmelblast{$prot} =~ /Gr21a/ || $dmelblast{$prot} =~ /Gr63a/){
				$dmelhit = "$dmelblast{$prot} CO2";
				$co2++;
			} elsif ($dmelblast{$prot} =~ /Gr43a/){
				$dmelhit = "$dmelblast{$prot} Fructose";
				$gr43++;
			} elsif ($dmelblast{$prot} =~ /Gr93a/ || $dmelblast{$prot} =~ /Gr66a/ || $dmelblast{$prot} =~ /Gr33a/){
				$dmelhit = "$dmelblast{$prot} Bitter";
				$bitter++;
			} else {
				$dmelhit = "$dmelblast{$prot}";
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
				$nname = "$gagashort{$gagaid}"."fragGR"."$fragornum"."$sufix";
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
				$deleteproteinsl50 .= "$prot "; # Delete ORco from final retained GRs
			} else {
				$nname = "$gagashort{$gagaid}"."GR"."$ornum"."$sufix";
				$ornum++;
			}
			if ($len >= $maxorlength){ # Print a warning in retained GRs longer than the max specified length
				print "Warning: $prot in $gagaid with length $len is longer than the max specified length\n";
			}
		}
		$newnames{$prot} = "$nname";


		print Results "$gagaid\t$gagasp{$gagaid}\t$prot\t$scaf\-$coord\t$nname\t$len\t$exonnumber{$prot}\t$orhit\t$ordomain\t$dmelhit\t$status\t$discard\n";

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
		if ($subline[8] =~ /Parent=([^;]+)(\S*)/){
		#if ($subline[8] =~ /transcript_id \"(\S+)\"/){
			$genename = $1;
			$rest = $2;
		}
		else {die "ERROR in run_OR_classification.pl: It fails detecting Parent ID in $line\n";}

#		my $nnamef = $newnames{$genename};
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
			print Resultsgff "$subline[0]\t$subline[1]\t$subline[2]\t$subline[3]\t$subline[4]\t$subline[5]\t$subline[6]\t$subline[7]\tParent=$nnamef"."$rest\;Status=$hstatus{$genename}\n";

			if ($hstatus{$genename} =~ /Complete/){
				print Resultsgffcomp "$subline[0]\t$subline[1]\t$subline[2]\t$subline[3]\t$subline[4]\t$subline[5]\t$subline[6]\t$subline[7]\tParent=$nnamef"."$rest\;Status=$hstatus{$genename}\n";
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
				if ($genename =~ /GR.*\.t/){
					$denovoc++; # Counting new annotated genes in bitacora gemoma
				}
			} elsif ($hstatus{$genename} =~ /Partial/){
				$fpar++;
				if ($genename =~ /GR.*\.t/){
					$denovop++; # Counting new annotated genes in bitacora gemoma
				}
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

system ("perl scripts/Tools/gff2fasta_v3.pl $genome $outgff $outprot ");
system ("sed \'s\/X\*\$\/\/\' $outprot\.pep.fasta > $outprot\.pep.fasta.tmp");
system ("mv $outprot\.pep.fasta.tmp $outprot\.pep.fasta");

system ("perl scripts/Tools/gff2fasta_v3.pl $genome $outgffall $outprotall ");
system ("sed \'s\/X\*\$\/\/\' $outprotall\.pep.fasta > $outprotall\.pep.fasta.tmp");
system ("mv $outprotall\.pep.fasta.tmp $outprotall\.pep.fasta");

system ("perl scripts/Tools/gff2fasta_v3.pl $genome $outgffcomp $outprotcomp ");
system ("sed \'s\/X\*\$\/\/\' $outprotcomp\.pep.fasta > $outprotcomp\.pep.fasta.tmp");
system ("mv $outprotcomp\.pep.fasta.tmp $outprotcomp\.pep.fasta");

#print Resultssum "$gagaid\t$gagasp{$gagaid}\t$totalor\t$com\t$comd\t$par\t$pard\t$parc\t$pse\t$psed\t$psec\t$l200\t$l100\t$l50\t$orconum\t$orconumcom\t$ninec\t$ninei\t$nineic\t$ninep\t$compard\t$ftot\t$fcom\t$fpar\t$fparc\t$fcomp\t$fpse\t$fninec\t$fninei\n";
print Resultssum "$gagaid\t$gagasp{$gagaid}\t$totalor\t$com\t$comd\t$par\t$pard\t$parc\t$pse\t$psed\t$psec\t$l200\t$l100\t$l50\t$grcount\t$orconum\t$orconumcom\t$ninec\t$ninei\t$nineic\t$ninep\t$sugar\t$bitter\t$co2\t$gr43\t$compard\t$ftot\t$fcomp\t$denovoc\t$fpar\t$denovop\t$fpse\t$fninec\t$fninei\t$frpartial\t$frpse\n";
close Resultssum;

