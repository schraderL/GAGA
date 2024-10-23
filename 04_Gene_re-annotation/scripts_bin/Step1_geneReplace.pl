#!/usr/bin/perl -w
use strict;
use Data::Dumper;

if (@ARGV != 5) {
    warn "#Usge: perl $0 <gff1> <gff2> <mRNA.lst> <gene.lst> <prefix> \n";
    exit;
}

my $file_gff1 = shift;
my $file_gff2 = shift;
my $lst1 = shift;
my $lst2 = shift;
my $prefix = shift;

my %gff1;
my %gff2;
readGff($file_gff1,\%gff1,1);
readGff($file_gff2,\%gff2,0);

my %lst1;
my %lst2;
readList($lst1,\%lst1);
readList($lst2,\%lst2);

#print STDERR Dumper \%gff1;
#exit;

foreach my $gene_id (keys %gff1) {
    my $gene_line = $gff1{$gene_id}{gene};
    my @t1= split /\s+/,$gene_line;

    if(exists $gff2{$gene_id}{gene}) {
        my $gene_line2 = $gff2{$gene_id}{gene};
        my @t2= split /\s+/,$gene_line2;

        if($t2[3] < $t1[3]) {
            $t1[3] = $t2[3];
        }
        if($t2[4] > $t1[4]) {
            $t1[4] = $t2[4];
        }

        my $func = $lst2{$gene_id};
        $t1[8] .= "geneName=$func;";
    }
    print join("\t",@t1),"\n";
    foreach my $mRNA_id (keys %{$gff1{$gene_id}{mRNA}}) {
        if(exists $lst1{$mRNA_id}) {
            my $mRNA_id_new = $lst1{$mRNA_id};
            foreach my $line (@{$gff2{$gene_id}{mRNA}{$mRNA_id_new}}) {
                print "$line\n";
            }
	    delete $gff2{$gene_id}{mRNA}{$mRNA_id_new};
        }else {
            foreach my $line (@{$gff1{$gene_id}{mRNA}{$mRNA_id}}) {
                print "$line\n";
            }
	    
        }
    }
    foreach my $mRNA_id (keys %{$gff2{$gene_id}{mRNA}}) {
	    foreach my $line (@{$gff2{$gene_id}{mRNA}{$mRNA_id}}) {
		print "$line\n";
	    }
    }
    delete $gff2{$gene_id};
}


foreach my $gene_id (keys %gff2) {
    my $gene_line = $gff2{$gene_id}{gene};
    print "$gene_line\n";
    foreach my $mRNA_id (keys %{$gff2{$gene_id}{mRNA}}) {
        foreach my $line (@{$gff2{$gene_id}{mRNA}{$mRNA_id}}) {
            print "$line\n";
        }
    }
}

####################
##################
sub readGff {
    my ($file,$p,$tag) = @_;

    open IN, "<$file" or die "failed to open $file:$!\n";
    my $gene_id='';
    my $id='';
    while (<IN>) {
        next if(/#/);
        chomp;
        if(/\sgene\s/) {
            $gene_id = $1 if(/ID=([^;]+);/);
            $p->{$gene_id}{gene} = $_;
        }elsif(/\smRNA\s/) {
            $id = $1 if(/ID=([^;]+);/);
	    $gene_id = $1 if(/($prefix\_g\d+)/);
            push @{$p->{$gene_id}{mRNA}{$id}},$_;
        }elsif(/\sCDS\s/){
	    $gene_id = $1 if(/($prefix\_g\d+)/);
	    if($tag) {
	        $id = $1 if(/Parent=([^;]+);/);
            }
            push @{$p->{$gene_id}{mRNA}{$id}},$_;
        }
    }
    close (IN);
}

sub readList {
    my ($file,$p) = @_;

    open IN, "<$file" or die "failed to open $file:$!\n";
    while (<IN>) {
        chomp;
        my @t=split /\s+/;
        $p->{$t[0]} = $t[1];
    }
    close (IN);
}
