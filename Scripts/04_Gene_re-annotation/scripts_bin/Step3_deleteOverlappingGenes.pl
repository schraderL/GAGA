#!/usr/bin/perl -w
use strict;
use Data::Dumper;

if (@ARGV != 3) {
    warn "#Usge: perl $0 <gff1> <Overlap.lst.replace.filter3> <prefix> \n";
    exit;
}

my $file_gff1 = shift;
my $lst1 = shift;
my $prefix = shift;

my %gff1;
readGff($file_gff1,\%gff1);

my %lst1;
readList($lst1,\%lst1);

foreach my $gene_id (keys %gff1) {
    my $gene_line = $gff1{$gene_id}{gene};

    print "$gene_line\n";
    foreach my $mRNA_id (keys %{$gff1{$gene_id}{mRNA}}) {
        if(exists $lst1{$mRNA_id}) {
		;
        }else {
            foreach my $line (@{$gff1{$gene_id}{mRNA}{$mRNA_id}}) {
                print "$line\n";
            }
        }
    }
}

####################
##################
sub readGff {
    my ($file,$p) = @_;

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
            push @{$p->{$gene_id}{mRNA}{$id}},$_;
        }elsif(/\sCDS\s/){
            push @{$p->{$gene_id}{mRNA}{$id}},$_;
        }
    }
    close (IN);
}

sub readList {
    my ($file,$p) = @_;

    open IN, "<$file" or die "failed to open $file:$!\n";
    while (<IN>) {
	next if(/^$/);
        chomp;
        my @t=split /\s+/;
	next if(!defined $t[1]);
	if($t[0] =~ /$prefix\_g/) {
		$p->{$t[0]} = 1;
	}elsif($t[1] =~ /$prefix\_g/) {
		$p->{$t[1]} = 1;
	}
    }
    close (IN);
}
