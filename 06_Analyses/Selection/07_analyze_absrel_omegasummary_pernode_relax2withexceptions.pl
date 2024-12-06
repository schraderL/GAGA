#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
my $dirname = dirname(__FILE__);

sub by_number {
    if ($a < $b){ -1 } elsif ($a > $b) { 1 } else { 0 }
}

# usage: perl analyze_absrel_omegasummary_pernode.pl hogtable outputname strictparameter hogqualfilter[None; Average; Percnine; Percsix]

## It goes through all the nodes, and reconcile the gene tree with the species tree. Strict mode: it only uses the nodes containing the expected number of species
# Relaxed mode: it uses all nodes from the gene tree, reconciling with the corresponding node in the species tree, even if only a fraction of the species are in the gene tree
# Additional option to relax reconciliation (named relax2), allow to reconcile a node if it contains some species that do not belong to that node (i.e.: leptanilla or any poneroid (only 1 or few) is in the formicoid clade, it reconciles in the formicoid clade rather than going to the root node)
	# Added 4 exceptions for relax2: node 164 (all would go down to 165); 194 (all would go to 198 with no dorylinae); 248 (all would go to 250 myrmicinae ignoring only two ectatomminae); See results in node 250, because it could also happen that results go to node 255

#my $srmode = "0.6"; # Define here the mode: Percentage of species in the "species tree node" required to be in the "gene tree node" to reconcile it. I.e.: 0.6 and 0.8 would be strict, lower than 0.5 would be more relaxed
my $srmode = $ARGV[2];
#my $relaxoutspecies = "0.02"; #| 0 = Strict | 0.02 = Relax2 ## Define here the percentage of species allowed not to be in that node in order to reconcile it. For instance, a node of 50 species, allows 1 species not in that node to reconcile that specific node. 
my $relaxoutspecies = $ARGV[4];

# Files with partitions and cleaned aln hmmclean 50% codon filt
my $cdsdir = "/home/projects/ku_00039/people/joeviz/orthology_trees/all_orthologs/Cleaned_alignments/run_all_withgard/hmmcleaner_50_aln/";
my $absreldir = "/home/projects/ku_00039/people/joeviz/orthology_trees/all_orthologs/Hyphy_absrel/absrel_withgardpartitions_hmmclean50aln/absrel_output_tables_dnds/"; # Directory containing Absrel output tables
my $alnext = "nonstop_hmmclean_gapfilt.aln"; # extension of the codon alignment files to find them

#my $pvalue = "0.001"; # Pvalue to filter significant positive selection, both in fdr and uncorrected pvalues
my $pvalue = $ARGV[5];

# File with orthogroups to use
my $hogfile = $ARGV[0];
#my $hogfile = "/home/projects/ku_00039/people/joeviz/orthology_trees/all_orthologs/Selection_tables/species_all_100percsp_orthogroups_singlecopy_counts.tsv";
#my $hogfile = "/home/projects/ku_00039/people/joeviz/orthology_trees/all_orthologs/Selection_tables/species_all_80percsp_orthogroups_singlecopy_counts.tsv";
#my $hogfile = "/home/projects/ku_00039/people/joeviz/orthology_trees/all_orthologs/Selection_tables/species_all_orthogroups.tsv";

my $filterminbranch = "150"; # Minimum number of branches required to include the dn/ds and positive selection results, otherwise the HOG is skipped

# File with the node and leafs in the tree
my $treetable = "/home/projects/ku_00039/people/joeviz/orthology_trees/all_orthologs/Selection_tables/Table_node_leaf_ids_withnames.txt";

# File with the alignment quality for each HOG, see lines to modify how to filter HOGs with low quality
my $zorrotable = "/home/projects/ku_00039/people/joeviz/orthology_alignments/all_orthogroups/codon_alignments/Zorro_average_scores.txt";
#my $qualfilter = "None"; # Options are: None; Average; Percnine; Percsix . Filter by average quality, by percentage of positions with score higher than 9, and the last is the sam but score of 6.
my $qualfilter = $ARGV[3];

# Output name
my $outputname = $ARGV[1];
#my $outputname = "Hyphy_absrel_omega_pernode_strict_100pc1to1.txt";
#my $outputname = "Hyphy_absrel_withgard_omega_perspecies_100pc1to1_fdr0001_strict60_relax2nof.txt";


## CODE

open (Scores, ">", $outputname);
print Scores "Species\tNodeID in species tree\tSpecies-Node name\tAverage dN\tAverage dN (excluding dN>2)\tAverage dS\tAverage dS (excluding dS>2)\tSum dN\tSum dN (excluding dN>2)\tSum dS\tSum dS (excluding dS>2)\tOmega (average (dN/dS)/3): calculated second/fourth column\tOmega (average (dN/dS)/3): calculated third/fifth column (excluding dN and dS >2)\tAverage omega ((dN/dS)/3): average value of the omega values in this node, calculated from the separate dn ds values\tAverage omega ((dN/dS)/3) (excluding omega>2)\tMedian omega ((dN/dS)/3): median value of the omega value in this node\tAverage omega (baseline omega): Average value of the baseline omega in this node\tAverage omega (baseline omega) (excluding omega>2)\tMedian omega (baseline omega)\tAverage omega (absrel omega): Average value of the absrel reported omega across the node\tAverage omega (absrel omega) (excluding omega>2)\tMedian omega (absrel omega)\tGenes under positive selection\tGenes under positive selection without FDR correction\tTotal number of genes tested\tTotal number of partitions tested\tProportion of genes under selection\tProportion of genes under selection across all partitions\tProportion of genes under selection without FDR correction across all partitions\tGene-HOG under selection\n";

open (Rtable, ">", "$outputname\_reconcilednodes.txt");
print Rtable "HOGID\tPercent of species covered in species tree node by gene tree node\tSpecies in the gene tree node\tSpecies from genetree node not in the species tree node\tSpecies in the species tree node\tGene tree node\tSpecies tree node\tGene tree leafs\tSpecies tree leafs\n";


my ($line, $name);

# Reading HOG table
my $hogtokeep = "";
open(File, "<", $hogfile);
while(<File>){
	chomp;
	$line = $_;
	next if ($line !~ /\S+/);
	next if ($line =~ /HOG\s/); # Skip header
	if ($line =~ /(HOG\S+)/){
		$hogtokeep .= " $1 ";
	}	
}
close File;

# Reading HOG alignment quality table
my $hogtokeepquality = "";
open(File, "<", $zorrotable);
while(<File>){
	chomp;
	$line = $_;
	next if ($line !~ /\S+/);
	next if ($line =~ /Average score/); # Skip header
	my @subl = split (/\t/, $line);

	if ($qualfilter =~ /None/){
		$hogtokeepquality .= " $subl[0] "; # No filter
	} elsif ($qualfilter =~ /Average/){
		if ($subl[1] >= 2){ # Filter by average quality, use 4
			$hogtokeepquality .= " $subl[0] ";
		}
	} elsif ($qualfilter =~ /Percnine/){
		my $perc = $subl[6]/$subl[3];
		if ($perc >= 0.10){ # Filter by number of well aligned positions >9
			$hogtokeepquality .= " $subl[0] ";
		}
	} elsif ($qualfilter =~ /Percsix/){
		my $perc = $subl[7]/$subl[3];
		if ($perc >= 0.50){ # Filter by number of well aligned positions >6
			$hogtokeepquality .= " $subl[0] ";
		}
	} else {
		die "Quality filter parameter is undefined\n";
	}
}
close File;

# Reading tree node table
my %treenodes; my %nodename;
open(File, "<", $treetable);
while(<File>){
	chomp;
	$line = $_;
	next if ($line !~ /\S+/);
	next if ($line =~ /Leafs/); # Skip header
	my @subl = split (/\t/, $line);
	$treenodes{$subl[0]} = $subl[1];
	$nodename{$subl[0]} = $subl[2];

}
close File;


# Initializing all variables

my %sumdn; my %sumds; my %sumdnfilt; my %sumdsfilt; 
my %numbranches; my %psel; my %pselfdr; my %signhogs;
my %numbranchesnopart; my %hogsinnode; # To count hogs, no partitions
my %totalspomega3; my %totalspomega3filt; # Filt will contain the omega values<2
my %totalfitomega; my %totalfitomegafilt;
my %totalabsrelomega; my %totalabsrelomegafilt;
my %arrayomega; my %arrayomegafit; my %arrayomegaabsrel;


# Reading Protein fasta

#my %protfasta;
system ("ls $cdsdir\/\*$alnext > tmp_protlist_$outputname\.txt");
open(File, "<", "tmp_protlist_$outputname\.txt");
while(<File>){
	chomp;
	$line = $_;
	next if ($line !~ /\S+/);	

	my $id = "";
	if ($line =~ /.*\/(\S+)\.cds/){
		$id = $1;
	} else {
		die "Can't find OG id in $line\n";
	}  

	my $hogid = "";
	if ($line =~ /(HOG\d\d\d\d\d\d\d)/){ # HOG0000492
		$hogid = $1;
	} else {
		die "Can't find HOG id in $line\n";
	}  

	#### Skipping here the HOG if it is not in the table (i.e. only single copy)
	next unless ($hogtokeep =~ /\s$hogid\s/); # Skip hogs not in hogtokeep list

	#### Skipping here the HOG if the alignment quality did not pass the filter stablished in this script
	next unless ($hogtokeepquality =~ /\s$hogid\s/); # Skip hogs not in hogtokeep list


	my $absrelout = "$absreldir\/$id\_absrel_wdnds.tsv";


	# Count the number of lines (i.e.: tested branches), and filter by the used $filterminbranch value.
	my $countlines = "0";
	open (Filejson, "<", $absrelout);
	while (<Filejson>){
		chomp;
		my $line2 = $_;
		next if ($line2 !~ /\S+/);	
		next if ($line2 =~ /Non-synonymous/);
		$countlines++;	
	}
	close Filejson;
	if ($countlines == 0){
		#die "Died because the number of branches is 0 for $absrelout\n";
		print "Warning: $id table is empty or it does not exist: $absrelout\n";
		next;
	}
	elsif ($countlines < $filterminbranch){
		next;
	}

	open (Filejson, "<", $absrelout);
	while (<Filejson>){
		chomp;
		my $line2 = $_;
		next if ($line2 !~ /\S+/);	
		next if ($line2 =~ /Non-synonymous/);

		$line2 =~ s/GAGA\_/GAGA\-/g; # Replace _ in GAGA-ID
		$line2 =~ s/OUT\_/OUT\-/g; # Replace _ in GAGA-ID
		$line2 =~ s/NCBI\_/NCBI\-/g; # Replace _ in GAGA-ID

		my @subl = split (/\t/, $line2);

		# Detect for which node or leaf the results are
		my $node = "NA";; 

		my $hognode = $subl[0];
		my $hogleafs = $subl[3];
		my @ahogleafs = split(/,/, $hogleafs);
		my $counthogleafs = @ahogleafs;

		if ($counthogleafs == 1){ # Leaf: GAGA protein id
			foreach my $treenode (keys %treenodes){
				my @atreenode = split (/\s/, $treenodes{$treenode});
				my $counttreeleafs = @atreenode;
				next if ($counttreeleafs > 1);
				if ($hognode =~ /$treenodes{$treenode}/){
					$node = $treenode;
				}
			}
		} else { # Node: Node id
			# Reconcile the node id with the species tree and get the corresponding node
			my $besthit = "nohit"; my $bhspcovered = "0"; my $bhspnotcovered = "0"; my $bhspnotcoveredlist = ""; my $bhsptotal = "1"; my $besthitleafs = "";

			foreach my $treenode (keys %treenodes){
				my $spingenetree = "0";
				my @atreenode = split (/\s/, $treenodes{$treenode});
				my $counttreeleafs = @atreenode;
				next if ($counttreeleafs == 1);

				my $genetreenodespid = " ";  my $counthogleafspecies = "0";
				my $spinspeciestree = "0"; my $spnotinspeciestree = "0";
				my $spnotinspeciestreelist = ""; # Variable to save species for debugging purposes
				foreach my $sp (@ahogleafs){ # Iterate through the species in the gene tree node
					next if ($sp !~ /\S+/); # Avoid spaces or other non-characters
					if ($sp =~ /(GAGA-\d\d\d\d)/ || $sp =~ /(NCBI-\d\d\d\d)/ || $sp =~ /(OUT-\d\d\d\d)/){
						my $spid = $1;
						if ($genetreenodespid =~ /$spid /){ # Species already counted, duplicate gene
							next;
						} else {
							$genetreenodespid .= "$spid ";
							$counthogleafspecies++;
							if ($treenodes{$treenode} =~ /$spid/){
								$spinspeciestree++; # Count if the species is present in the species tree node
							} else {
								$spnotinspeciestree++; # Count if the species is not present in the species tree node
								$spnotinspeciestreelist .= "$spid";
							}
						}
					}
				}

				# $spingenetree contains the number of species from the gene tree node, into the species tree node
				# $counthogleafs contain the total number of genes (! not species) in that gene tree node
				# $counthogleafspecies contain the total number of species in that gene tree node (each species can have 1 or more gene)

				# $spinspeciestree contains number of species from the gene tree covered in the species tree node
				# $spnotinspeciestree contains the number of species from the gene tree node that are NOT in the species tree node
				# $counttreeleafs contain the number of species in the species tree node

				#print Rtable "Doing $id\t$spingenetree\t$counthogleafspecies\t$spinspeciestree\t$spnotinspeciestree\t$counttreeleafs\t$hognode\t$hogleafs\t$treenode\t$treenodes{$treenode}\n"; # Debug

				my $percoutspeciesinclade = 100000;
				if ($spinspeciestree > 0){
					$percoutspeciesinclade = $spnotinspeciestree/$spinspeciestree; # Counting here the percentage based on the total species contained in this case
				}
				#$percoutspeciesinclade = $spnotinspeciestree/$counttreeleafs; # Counting here the percentage based on the total species contained in the species tree node
				#if ($spnotinspeciestree > 0){ # Strict reconciliation: Gene tree node contain a species that it is not in the species tree node. This node is not the corresponding
				if ($percoutspeciesinclade > $relaxoutspecies){ # Relax2 reconciliation: Gene tree node contain a species that it is not in the species tree node, but we check for a percentage to allow some bad assigned genes, or cases with leptanilla, so we use the defined percentage to allow it.
					next;
				}

				# Moved block to here to avoid iterations with species tree nodes that already do not cover the species of the gene tree (reduce computation time)
				foreach my $sp (@atreenode){ # Iterate through the species in species tree node
					next if ($sp !~ /\S+/); # Avoid spaces or other non-characters
					if ($hogleafs =~ /$sp/){
						$spingenetree++; # Count how many species are in the gene tree node
					}
				}


				if ($spinspeciestree == $spingenetree){ # Control, this should be the same count
					if ($besthit =~ /nohit/){ # First hit, save
						$besthit = $treenode;
						$besthitleafs = $treenodes{$treenode};
						$bhspcovered = $spingenetree;
						$bhspnotcovered = $spnotinspeciestree;
						$bhspnotcoveredlist = $spnotinspeciestreelist;
						$bhsptotal = $counttreeleafs;
					} elsif ($spingenetree < $bhspcovered){ # It can now happen with relax2 ($relaxoutspecies > 0), as the current hit is a younger node allowing one out species (based on the percentage paremeter) OLD: It should not happen
						#die "Number of species counted in the gene tree is lower than the previous hit in $line\tGene line is $line2\nCurrent species tree node is $treenode with $spingenetree species counted in gene tree, $spnotinspeciestree species that are in the gene tree node and not the species tree node: $spnotinspeciestreelist, and $counttreeleafs total species in the species tree node\nBest hit previously counted is species tree node $besthit with $bhspcovered species counted in gene tree, $bhspnotcovered species that are in the gene tree node and not the species tree node, and $bhsptotal total species in the species tree node\n\n";
						# the current is a younger node allowing one/two out species (based on the percentage parameter) 
						# Then save this one, with the following exceptions: nodes 164 (all would go down to 165); 194 (many would go to 198 with no dorylinae); 248 (all would go to 250 myrmicinae ignoring only two ectatomminae); 
						# EXCEPTIONS:
						if ($treenode =~ /^165$/ && $spnotinspeciestreelist =~ /GAGA-0392/){ 
							# Do not save this hit, as it should be the root. Node 165 + Leptanilla = Node 164.
						} elsif ($treenode =~ /^198$/ && $spnotinspeciestreelist =~ /GAGA-0534/){
							# Do not save this hit, as it should belong to Node 194. Node 198 + any dorylinae = Node 194.	
						} elsif ($treenode =~ /^198$/ && $spnotinspeciestreelist =~ /GAGA-0577/){
							# Do not save this hit, as it should belong to Node 194. Node 198 + any dorylinae = Node 194.	
						} elsif ($treenode =~ /^198$/ && $spnotinspeciestreelist =~ /GAGA-0379/){
							# Do not save this hit, as it should belong to Node 194. Node 198 + any dorylinae = Node 194.	
						} elsif ($treenode =~ /^198$/ && $spnotinspeciestreelist =~ /NCBI-0001/){
							# Do not save this hit, as it should belong to Node 194. Node 198 + any dorylinae = Node 194.	
						} elsif ($treenode =~ /^250$/ && $spnotinspeciestreelist =~ /GAGA-0306/){
							# Do not save this hit, as it should belong to Node 248. Node 250 + any Ectatomminae = Node 248.	
						} elsif ($treenode =~ /^250$/ && $spnotinspeciestreelist =~ /GAGA-0234/){
							# Do not save this hit, as it should belong to Node 248. Node 250 + any Ectatomminae = Node 248.	
						} elsif ($treenode =~ /^255$/ && $spnotinspeciestreelist =~ /NCBI-0012/){
							# Do not save this hit, as it should belong to Node 250. Node 255 + Pogonomyrmex-Manica-Myrmica = Node 250.	
						} elsif ($treenode =~ /^255$/ && $spnotinspeciestreelist =~ /NCBI-0114/){
							# Do not save this hit, as it should belong to Node 250. Node 255 + Pogonomyrmex-Manica-Myrmica = Node 250.	
						} else {
							$besthit = $treenode;
							$besthitleafs = $treenodes{$treenode};
							$bhspcovered = $spingenetree;
							$bhspnotcovered = $spnotinspeciestree;
							$bhspnotcoveredlist = $spnotinspeciestreelist;
							$bhsptotal = $counttreeleafs;							
						}
				

					} elsif ($spingenetree > $bhspcovered){ # It can now happen with relax2 ($relaxoutspecies > 0), as the previous best hit is a younger node allowing one out species (based on the percentage parameter) OLD: It should not happen
						#die "Number of species counted in the gene tree is higher than the previous hit in $line\tGene line is $line2\nCurrent species tree node is $treenode with $spingenetree species counted in gene tree, $spnotinspeciestree species that are in the gene tree node and not the species tree node, and $counttreeleafs total species in the species tree node\nBest hit previously counted is species tree node $besthit with $bhspcovered species counted in gene tree, $bhspnotcovered species that are in the gene tree node and not the species tree node, and $bhsptotal total species in the species tree node\n\n";
						# the previous best hit is a younger node allowing one/two out species (based on the percentage parameter) 
						# Don't save this one, with the following exceptions: nodes 164 (all would go down to 165); 194 (many would go to 198 with no dorylinae); 248 (all would go to 250 myrmicinae ignoring only two ectatomminae); 
						# EXCEPTIONS
						if ($besthit =~ /^165$/ && $bhspnotcoveredlist =~ /GAGA-0392/){ 
							# Save this hit, as it should be the root. Node 165 + Leptanilla = Node 164.
							$besthit = $treenode;
							$besthitleafs = $treenodes{$treenode};
							$bhspcovered = $spingenetree;
							$bhspnotcovered = $spnotinspeciestree;
							$bhspnotcoveredlist = $spnotinspeciestreelist;
							$bhsptotal = $counttreeleafs;								
						} elsif ($besthit =~ /^198$/ && $bhspnotcoveredlist =~ /GAGA-0534/){
							# Save this hit, as it should belong to Node 194. Node 198 + any dorylinae = Node 194.	
							$besthit = $treenode;
							$besthitleafs = $treenodes{$treenode};
							$bhspcovered = $spingenetree;
							$bhspnotcovered = $spnotinspeciestree;
							$bhspnotcoveredlist = $spnotinspeciestreelist;
							$bhsptotal = $counttreeleafs;	
						} elsif ($besthit =~ /^198$/ && $bhspnotcoveredlist =~ /GAGA-0577/){
							# Save this hit, as it should belong to Node 194. Node 198 + any dorylinae = Node 194.	
							$besthit = $treenode;
							$besthitleafs = $treenodes{$treenode};
							$bhspcovered = $spingenetree;
							$bhspnotcovered = $spnotinspeciestree;
							$bhspnotcoveredlist = $spnotinspeciestreelist;
							$bhsptotal = $counttreeleafs;	
						} elsif ($besthit =~ /^198$/ && $bhspnotcoveredlist =~ /GAGA-0379/){
							# Save this hit, as it should belong to Node 194. Node 198 + any dorylinae = Node 194.	
							$besthit = $treenode;
							$besthitleafs = $treenodes{$treenode};
							$bhspcovered = $spingenetree;
							$bhspnotcovered = $spnotinspeciestree;
							$bhspnotcoveredlist = $spnotinspeciestreelist;
							$bhsptotal = $counttreeleafs;	
						} elsif ($besthit =~ /^198$/ && $bhspnotcoveredlist =~ /NCBI-0001/){
							# Save this hit, as it should belong to Node 194. Node 198 + any dorylinae = Node 194.	
							$besthit = $treenode;
							$besthitleafs = $treenodes{$treenode};
							$bhspcovered = $spingenetree;
							$bhspnotcovered = $spnotinspeciestree;
							$bhspnotcoveredlist = $spnotinspeciestreelist;
							$bhsptotal = $counttreeleafs;	
						} elsif ($besthit =~ /^250$/ && $bhspnotcoveredlist =~ /GAGA-0306/){
							# Save this hit, as it should belong to Node 248. Node 250 + any Ectatomminae = Node 248.	
							$besthit = $treenode;
							$besthitleafs = $treenodes{$treenode};
							$bhspcovered = $spingenetree;
							$bhspnotcovered = $spnotinspeciestree;
							$bhspnotcoveredlist = $spnotinspeciestreelist;
							$bhsptotal = $counttreeleafs;	
						} elsif ($besthit =~ /^250$/ && $bhspnotcoveredlist =~ /GAGA-0234/){
							# Save this hit, as it should belong to Node 248. Node 250 + any Ectatomminae = Node 248.
							$besthit = $treenode;
							$besthitleafs = $treenodes{$treenode};
							$bhspcovered = $spingenetree;
							$bhspnotcovered = $spnotinspeciestree;
							$bhspnotcoveredlist = $spnotinspeciestreelist;
							$bhsptotal = $counttreeleafs;	
						} elsif ($besthit =~ /^255$/ && $bhspnotcoveredlist =~ /NCBI-0012/){
							# Save this hit, as it should belong to Node 250. Node 255 + Pogonomyrmex-Manica-Myrmica = Node 250.	
							$besthit = $treenode;
							$besthitleafs = $treenodes{$treenode};
							$bhspcovered = $spingenetree;
							$bhspnotcovered = $spnotinspeciestree;
							$bhspnotcoveredlist = $spnotinspeciestreelist;
							$bhsptotal = $counttreeleafs;	
						} elsif ($besthit =~ /^255$/ && $bhspnotcoveredlist =~ /NCBI-0114/){
							# Save this hit, as it should belong to Node 250. Node 255 + Pogonomyrmex-Manica-Myrmica = Node 250.	
							$besthit = $treenode;
							$besthitleafs = $treenodes{$treenode};
							$bhspcovered = $spingenetree;
							$bhspnotcovered = $spnotinspeciestree;
							$bhspnotcoveredlist = $spnotinspeciestreelist;
							$bhsptotal = $counttreeleafs;								
						} else {
							# Do not save remaining nodes here							
						}						

					} elsif ($bhsptotal eq $counttreeleafs){
						print "Number of leafs gene tree node is the same as the previous hit in $line\tGene line is $line2\nCurrent species tree node is $treenode with $spingenetree counts and $counttreeleafs leafs\nBest hit previously counted is species tree node $besthit with $bhspcovered count and $bhsptotal leafs\nCurrent species tree node is $treenode with $spingenetree species counted in gene tree, $spnotinspeciestree species that are in the gene tree node and not the species tree node, and $counttreeleafs total species in the species tree node\nBest hit previously counted is species tree node $besthit with $bhspcovered species counted in gene tree, $bhspnotcovered species that are in the gene tree node and not the species tree node, and $bhsptotal total species in the species tree node\n\n";
					} elsif ($counttreeleafs < $bhsptotal){ # Node contain less leafs and still adjust to the gene tree node, keep this one
						# Exceptions added because of the relaxed reconciliation added with parameter $relaxoutspecies; are nodes 164 (all would go down to 165); 194 (all would go to 198 with no dorylinae); 248 (all would go to 250 myrmicinae ignoring only two ectatomminae); See results in node 250, because it could also happen that results go to node 255
						if ($treenode =~ /^165$/ && $spnotinspeciestreelist =~ /GAGA-0392/){ 
							# Do not save this hit, as it should be the root. Node 165 + Leptanilla = Node 164.
						} elsif ($treenode =~ /^198$/ && $spnotinspeciestreelist =~ /GAGA-0534/){
							# Do not save this hit, as it should belong to Node 194. Node 198 + any dorylinae = Node 194.	
						} elsif ($treenode =~ /^198$/ && $spnotinspeciestreelist =~ /GAGA-0577/){
							# Do not save this hit, as it should belong to Node 194. Node 198 + any dorylinae = Node 194.	
						} elsif ($treenode =~ /^198$/ && $spnotinspeciestreelist =~ /GAGA-0379/){
							# Do not save this hit, as it should belong to Node 194. Node 198 + any dorylinae = Node 194.	
						} elsif ($treenode =~ /^198$/ && $spnotinspeciestreelist =~ /NCBI-0001/){
							# Do not save this hit, as it should belong to Node 194. Node 198 + any dorylinae = Node 194.	
						} elsif ($treenode =~ /^250$/ && $spnotinspeciestreelist =~ /GAGA-0306/){
							# Do not save this hit, as it should belong to Node 248. Node 250 + any Ectatomminae = Node 248.	
						} elsif ($treenode =~ /^250$/ && $spnotinspeciestreelist =~ /GAGA-0234/){
							# Do not save this hit, as it should belong to Node 248. Node 250 + any Ectatomminae = Node 248.	
						} elsif ($treenode =~ /^255$/ && $spnotinspeciestreelist =~ /NCBI-0012/){
							# Do not save this hit, as it should belong to Node 250. Node 255 + Pogonomyrmex-Manica-Myrmica = Node 250.	
						} elsif ($treenode =~ /^255$/ && $spnotinspeciestreelist =~ /NCBI-0114/){
							# Do not save this hit, as it should belong to Node 250. Node 255 + Pogonomyrmex-Manica-Myrmica = Node 250.	
						} else {
							$besthit = $treenode;
							$besthitleafs = $treenodes{$treenode};
							$bhspcovered = $spingenetree;
							$bhspnotcovered = $spnotinspeciestree;
							$bhspnotcoveredlist = $spnotinspeciestreelist;
							$bhsptotal = $counttreeleafs;							
						}

					} elsif ($counttreeleafs > $bhsptotal){ # Node contain more leafs, keep previous node
						# Don't save this node
					} else {
						print "Situation I didn't think in $line\tGene line is $line2\nCurrent species tree node is $treenode with $spingenetree species counted in gene tree, $spnotinspeciestree species that are in the gene tree node and not the species tree node, and $counttreeleafs total species in the species tree node\nBest hit previously counted is species tree node $besthit with $bhspcovered species counted in gene tree, $bhspnotcovered species that are in the gene tree node and not the species tree node, and $bhsptotal total species in the species tree node\n\n";
					}
				} else {
					die "Different number of counts in $line\n Species counted in species tree is $spinspeciestree for gene tree node: $hognode in $line2\n Species counted in gene tree is $spingenetree in $treenode\n";
				}

			}		
			# Strict mode, only corresponding a gene tree node with the species tree if it covers more than 60% or 80% of the species
			# Relax mode, corresponding a gene tree node with the species tree covering the species in the gene tree, with at least 20% or 40% of the species
			my $div = $bhspcovered/$bhsptotal;
			my $perchit = sprintf("%.2f", $div);
			print Rtable "$id\t$perchit\t$bhspcovered\t$bhspnotcovered\t$bhsptotal\t$hognode\t$besthit\t$hogleafs\t$besthitleafs\n";

			if ($perchit >= $srmode){ # Strict/relaxed variable
				$node = $besthit;
			} else {
				next;
			}
		}

		#

		# Step to filter omegas from the calculations that are too low or too high, due to impossibility to estimate dn or ds. Only if positive selection is not detected
		if ($subl[2] > $pvalue){
#			if ($subl[8] <= 1e-10 || $subl[9] <= 1e-10){ # Filter both dn and ds
#			if ($subl[8] <= 1e-10){ # Filter dn and ds
#			if ($subl[9] <= 1e-10){ # Filter ds
#				next;
#			}
		}

		$sumdn{$node} += $subl[8];
		if ($subl[8] <= 2){
			$sumdnfilt{$node} += $subl[8];
		}
		$sumds{$node} += $subl[9];
		if ($subl[9] <= 2){
			$sumdsfilt{$node} += $subl[9];
		}
	
		$totalfitomega{$node} += $subl[4];
		push (@{$arrayomegafit{$node}}, $subl[4]);
		$totalabsrelomega{$node} += $subl[6];
		push (@{$arrayomegaabsrel{$node}}, $subl[6]);
		$totalspomega3{$node} += $subl[10];
		push (@{$arrayomega{$node}}, $subl[10]);

		if ($subl[4] <= 2){
			$totalfitomegafilt{$node} += $subl[4];
		}
		if ($subl[6] <= 2){
			$totalabsrelomegafilt{$node} += $subl[6];
		}
		if ($subl[10] <= 2){
			$totalspomega3filt{$node} += $subl[10];
		}

		$numbranches{$node}++; # Count all branches, including partitions

		if (exists $hogsinnode{$node}){
				#if ($hogsinnode{$node} =~ /$hogid\-\-$subl[0]\, /){
				if ($hogsinnode{$node} =~ /$hogid\-\-/){ ### Using just the hog id to count the HOG, previous line would count additional branches in the tree for the same HOG
					my $test = "OK"; # Random variable, just to skip this one as the partition is already counted
				} else {
					$numbranchesnopart{$node}++; # Count only once each partition for the specific species/node
					$hogsinnode{$node} .= "$hogid\-\-$subl[0], ";
				}
		} else {
			$numbranchesnopart{$node}++;
			$hogsinnode{$node} .= "$hogid\-\-$subl[0], ";
		}

		if ($subl[1] <= $pvalue){ # Count positive selection, FDR
			#### ADD here if partitions exist, to only count it once
			if (exists $signhogs{$node}){
				if ($signhogs{$node} =~ /$hogid\-\-$subl[0]\, /){
					next;
				}
			}
			####
			$psel{$node}++;
			$pselfdr{$node}++;
			$signhogs{$node} .= "$hogid\-\-$subl[0], ";
		} elsif ($subl[2] <= $pvalue){
			$psel{$node}++;
		}
		 
	}
	close Filejson;

}
close File;


# Print scores for each species

foreach my $node (sort by_number keys(%treenodes)){

	my @leafs = split (/\s/, $treenodes{$node});
	my $countleafs = @leafs;
#	next if ($countleafs > 1); ### Only go through the leafs here, skip internal nodes

	# Variables
	my $avedn = "NA"; my $aveds = "NA"; my $avednfilt = "NA"; my $avedsfilt = "NA"; my $omegaave3 = "NA"; my $omegaave3filt = "NA"; my $spomega3 = "NA"; my $omegafit = "NA"; my $omegaabsrel = "NA"; my $spomega3filt = "NA"; my $omegafitfilt = "NA"; my $omegaabsrelfilt = "NA";
	my $medianspomega3 = "NA"; my $medianomegafit = "NA"; my $medianomegaabsrel = "NA";

	unless (exists $numbranches{$node}){ # Just in case that a node have no data
		print Scores "$treenodes{$node}\t$node\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n";
		next;
	} 
	if ($numbranches{$node} eq 0){ # Just in case that a node have no data
		print Scores "$treenodes{$node}\t$node\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n";
		next;
	} 

	$avedn = $sumdn{$node}/$numbranches{$node};
	$aveds = $sumds{$node}/$numbranches{$node};
	$avednfilt = $sumdnfilt{$node}/$numbranches{$node};
	$avedsfilt = $sumdsfilt{$node}/$numbranches{$node};	
	$omegaave3 = ($avedn/$aveds)/3;
	$omegaave3filt = ($avednfilt/$avedsfilt)/3;	
	$spomega3 = $totalspomega3{$node}/$numbranches{$node};	
	$omegafit = $totalfitomega{$node}/$numbranches{$node};
	$omegaabsrel = $totalabsrelomega{$node}/$numbranches{$node};
	$spomega3filt = $totalspomega3filt{$node}/$numbranches{$node};	
	$omegafitfilt = $totalfitomegafilt{$node}/$numbranches{$node};
	$omegaabsrelfilt = $totalabsrelomegafilt{$node}/$numbranches{$node};

	# Median
	my $mid = int @{$arrayomega{$node}}/2;
	my @sorted_values = sort by_number@{$arrayomega{$node}};
	if (@{$arrayomega{$node}} % 2) {
	    $medianspomega3 = $sorted_values[ $mid ];
	} else {
	    $medianspomega3 = ($sorted_values[$mid-1] + $sorted_values[$mid])/2;
	} 

	$mid = ""; @sorted_values = ();
	$mid = int @{$arrayomegafit{$node}}/2;
	@sorted_values = sort by_number@{$arrayomegafit{$node}};
	if (@{$arrayomegafit{$node}} % 2) {
	    $medianomegafit = $sorted_values[ $mid ];
	} else {
	    $medianomegafit = ($sorted_values[$mid-1] + $sorted_values[$mid])/2;
	} 

	$mid = ""; @sorted_values = ();
	$mid = int @{$arrayomegaabsrel{$node}}/2;
	@sorted_values = sort by_number@{$arrayomegaabsrel{$node}};
	if (@{$arrayomegaabsrel{$node}} % 2) {
	    $medianomegaabsrel = $sorted_values[ $mid ];
	} else {
	    $medianomegaabsrel = ($sorted_values[$mid-1] + $sorted_values[$mid])/2;
	} 

	my $pselfdrpercnopart = 0; my $pselfdrperc = 0; my $pselfdrnum = 0;
	if (exists $pselfdr{$node}){ # In case there is no positive selection in this node, it will give an error with uninitialized value
		$pselfdrpercnopart = $pselfdr{$node}/$numbranchesnopart{$node};
		$pselfdrperc = $pselfdr{$node}/$numbranches{$node};
		$pselfdrnum = $pselfdr{$node};
	}
	my $pselperc = 0; my $pselnum = 0;
	if (exists $psel{$node}){ # In case there is no positive selection in this node, it will give an error with uninitialized value
		$pselperc = $psel{$node}/$numbranches{$node};
		$pselnum = $psel{$node};
	}
	my $signhoglist = "";
	if (exists $signhogs{$node}){ # In case there is no positive selection in this node, it will give an error with uninitialized value
		$signhoglist = $signhogs{$node};
	}

	print Scores "$treenodes{$node}\t$node\t$nodename{$node}\t$avedn\t$avednfilt\t$aveds\t$avedsfilt\t$sumdn{$node}\t$sumdnfilt{$node}\t$sumds{$node}\t$sumdsfilt{$node}\t$omegaave3\t$omegaave3filt\t$spomega3\t$spomega3filt\t$medianspomega3\t$omegafit\t$omegafitfilt\t$medianomegafit\t$omegaabsrel\t$omegaabsrelfilt\t$medianomegaabsrel\t$pselfdrnum\t$pselnum\t$numbranchesnopart{$node}\t$numbranches{$node}\t$pselfdrpercnopart\t$pselfdrperc\t$pselperc\t$signhoglist\n";
	
}


system ("rm tmp_protlist_$outputname\.txt");
close Scores;
close Rtable;


# Change dots for commas to read it in excel
system("sed \'s/\\./,/g\' $outputname > $outputname\_comma.txt");


