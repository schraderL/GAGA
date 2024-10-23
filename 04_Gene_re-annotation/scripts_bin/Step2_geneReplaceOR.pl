#!/usr/bin/perl -w
use strict;
use Data::Dumper;

if (@ARGV != 3) {
    warn "#Usge: perl $0 <gff1> <gff2> <mRNA.lst.id> \n";
    exit;
}

my $file_gff1 = shift;
my $file_gff2 = shift;
my $lst1 = shift;

my %gff1;
my %gff2;
readGff($file_gff1,\%gff1);
readGff($file_gff2,\%gff2);

my %lst1;
readList($lst1,\%lst1);


foreach my $gene_id (keys %gff1) {
    my $gene_line = $gff1{$gene_id}{gene};
    my @t1= split /\s+/,$gene_line;

    print join("\t",@t1),"\n";
    my $gene_id_new = '';
    foreach my $mRNA_id (keys %{$gff1{$gene_id}{mRNA}}) {
        if(exists $lst1{$mRNA_id}) {
            my $mRNA_id_new = $lst1{$mRNA_id};
	    $gene_id_new = $mRNA_id_new;
	    $gene_id_new = "g$mRNA_id_new" ;

            foreach my $line (@{$gff2{$gene_id_new}{mRNA}{$mRNA_id_new}}) {
                $line =~ s/Parent=$gene_id_new/Parent=$gene_id/;
                print "$line\n";
            }
	    delete $gff2{$gene_id_new}{mRNA}{$mRNA_id_new};
	    delete $gff2{$gene_id_new};
        }else {
            foreach my $line (@{$gff1{$gene_id}{mRNA}{$mRNA_id}}) {
                print "$line\n";
            }
	    
        }
    }
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
        $p->{$t[1]} = $t[0];
    }
    close (IN);
}
