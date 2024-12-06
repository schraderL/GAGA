#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
my $dirname = dirname(__FILE__);

### Create script to submit absrel in a subset of orthologs. Here using single copy, specified in the hogtokepp file

# Number of threads per alignment
my $nthreads=4;

my $alnfolder = "/home/projects/ku_00039/people/joeviz/orthology_trees/all_orthologs/Cleaned_alignments/run_all_withgard/hmmcleaner_50_aln/";
my $alnext = "nonstop_hmmclean_gapfilt.aln";

my $treefolder = "/home/projects/ku_00039/people/joeviz/orthology_trees/all_orthologs/Cleaned_alignments/run_all_withgard/hmmcleaner_50_aln/";
my $treeext = "hmmclean_gapfilt.rooted.treefile";

my $outdir="absrel_1to1_rooted";

### Hog file to filter only this hogs
my $hogfile = "/home/projects/ku_00039/people/joeviz/orthology_trees/all_orthologs/Selection_tables/species_all_80percsp_orthogroups_singlecopy_counts.tsv"; # List of orthogroups (HOG IDs) to subset
my $hogstokeep = "";
open(File, "<", $hogfile);
while(<File>){
    chomp;
    my $line = $_;
    next if ($line !~ /\S+/);
    next if ($line =~ /HOG\s/); # Skip header
    if ($line =~ /(HOG\S+)/){
        $hogstokeep .= " $1 ";
    }   
}
close File;
 

system ("mkdir -p $outdir");

my $counter = "1";

system ("ls $alnfolder\/\*$alnext > tmp_alns.txt");

open (Script, ">", "submit_all_hyphy.sh");
open (File, "<", "tmp_alns.txt");
while(<File>){
    chomp;  
    my $line = $_;
    next if ($line !~ /\S+/);
    my $filename = "";
    if ($line =~ /.*\/(\S+?)\.cds/){
        $filename = $1;
    } else {
        die "It cannot find filename in $line\n";
    }
    my $hogname = "";
    if ($line =~ /.*\/(HOG\d+)/){
        $hogname = $1;
    } else {
        die "It cannot find HOG name id in $line\n";
    }

    # Only process hogs in the orhtolog group of interest
    unless ($hogstokeep =~ / $hogname /){
#        print "Skipping no single copy $line\n";
        next;
    }

    my $treefile = "$treefolder\/$filename\.$treeext";
    unless (-s $treefile) { 
        print "$treefile not found or empty in tree folder for $filename in aln $line, skipping...\n";
        next; 
    }

    open (Scripthog, ">", "submit_hyphy_$filename\.sh");
    print Scripthog "#!/bin/sh
#PBS -W group_list=ku_00039 -A ku_00039
#PBS -m n
#PBS -l nodes=1:ppn=$nthreads:thinnode
#PBS -l mem=16gb
#PBS -l walltime=99:00:00

# Go to the directory from where the job was submitted (initial directory is \$HOME)
echo Working directory is \$PBS_O_WORKDIR
cd \$PBS_O_WORKDIR

### Here follows the user commands:
# Define number of processors
NPROCS=`wc -l < \$PBS_NODEFILE`
echo This job has allocated \$NPROCS nodes

# Load all required modules for the job
module load tools
module load ngs
module load fasttree/2.1.11
module load perl/5.24.0
module load openmpi/gcc
module load hyphy/2.5.38

# This is where the work is done
# Make sure that this script is not bigger than 64kb ~ 150 lines, otherwise put in seperat script and execute from here

hyphy CPU=$nthreads absrel --alignment $line --tree $treefile --output $outdir/$filename\_ABSREL.json > $outdir/$filename\_absrel.output

";

    close Scripthog;

    print Script "qsub submit_hyphy_$filename\.sh\n";

    $counter++;
    if ($counter > 100){
        $counter = 1;
        print Script "sleep 60s\n";
    }

}
close File;
close Script;

system ("rm tmp_alns.txt");


