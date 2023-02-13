#!/usr/bin/perl 
use strict;
use File::Basename qw(basename);
use FindBin qw($Bin);
use Getopt::Long;

my $version=1.0; #
my $version=1.1; # add -cds_ns ,2009-4-13

die "Usage:perl $0 -perfect  -check_cds -mini_cds 150 -mini_seq 2000 -cds_ns 10 -verbose -outdir ./ file.gff file.fa\n" if (@ARGV<2);

my ($Perfect_gene,$Check_cds,$Mini_cds,$Max_Ns,$Verbose,$Outdir);
my ($Mini_seq);
my $CDS_ns;
GetOptions(
	"perfect"=>\$Perfect_gene,
	"check_cds"=>\$Check_cds,
	"mini_cds:i"=>\$Mini_cds,
	"mini_seq:i"=>\$Mini_seq,
	"max_ns:i"=>\$Max_Ns,
	"cds_ns:i"=>\$CDS_ns,
	"verbose"=>\$Verbose, #Globe para
	"outdir:s"=>\$Outdir,
);

die `pod2text $0` if (@ARGV<2);
$Outdir ||=".";
$Outdir=~s/\/$//;#Globe para
my $GFF_file=shift;
my $Genome_file=shift;

my $GFF_file_name=basename($GFF_file);
my $GFF_out_file=$Outdir."/".$GFF_file_name.".check.gff";

my %mRNA;
read_gff($GFF_file,\%mRNA);

check_redundance(\%mRNA);#test,ok

mini_cds(\%mRNA,$Mini_cds) if (defined $Mini_cds);#test ok

perfect_gene(\%mRNA) if (defined $Perfect_gene);#test ok

check_cds($Genome_file,\%mRNA) if (defined $Check_cds);

create_gff($GFF_out_file,\%mRNA);

sub create_gff{
	my ($file,$p)=@_;
	my $output;
	foreach my $chr( sort keys %$p){
		foreach my $id(sort keys %{$p->{$chr}}){
			$output.=join("\t",@{$p->{$chr}{$id}{mRNA}})."\n";
			@{$p->{$chr}{$id}{CDS}}= sort {$a->[3] <=>$b->[3]} @{$p->{$chr}{$id}{CDS}};
			foreach my $cds (@{$p->{$chr}{$id}{CDS}}){
				$output.=join("\t",@$cds)."\n";
			}
		}
	}
	open OUT,">$file" or die "Fail $file:$!";
	print OUT $output;
	close OUT;
	
}

sub read_gff{
	my ($file,$p)=@_;
	open IN,$file or die "Fail $file:$!";
	while(<IN>){
		next if /^#/;
		chomp;
		my @c=split(/\t/);
		@c[3,4]=@c[4,3] if ($c[3]>$c[4]);
		if ($c[2] eq 'mRNA' && $c[8]=~/ID=([^;\s]+)/){
			my $id=$1;
			@{$p->{$c[0]}{$id}{mRNA}}=@c;
			if ($Perfect_gene && $c[8]=~/Type=([^;\s]+)/){
				$p->{$c[0]}{$id}{type}=$1;#just for predict result.
			}elsif($Perfect_gene){
				die "You chose perfect_gene option ,but input gff file has no information of type. Column 9 of mRNA  in  gff file  should be in fellow format:\n\tID=id;Type=type\nType must be:multi-exon,single-exon,no-first,no-last or no-first-last\n";
			}
		}elsif($c[2] eq 'CDS' && $c[8]=~/Parent=([^;\s]+)/){
			push @{$p->{$c[0]}{$1}{CDS}},[@c];
		}else{
			#die "Format ERROR: Line $. of $file\n";
		}
	}
	close IN;
}

sub check_redundance{
	my ($p)=@_;
	my ($all_size,$pure_size,$redunt_size,$total_gene_len);
	foreach my $seq_name (sort keys %$p) {
                my $seq_p = $p->{$seq_name};
                my (@pos1,@pos2);
                foreach my $gene_name (sort keys %$seq_p) {
                        my $gene_p = $seq_p->{$gene_name};
                        my $strand = $gene_p->{mRNA}[6];
                        my $gene_start = $gene_p->{mRNA}[3];
                        my $gene_end = $gene_p->{mRNA}[4];
                        my $gene_len = abs($gene_end-$gene_start) + 1;
                        $total_gene_len += $gene_len; 
                        push @pos1, [$gene_start,$gene_end] if($strand eq '+');
                        push @pos2, [$gene_end,$gene_start] if($strand eq '-');
                        
                }
	                
                my ($num1,$num2,$num3) = Conjoin_fragment(\@pos1);
                $all_size += $num1;
                $pure_size += $num2;
                $redunt_size += $num3;

                my ($num1,$num2,$num3) = Conjoin_fragment(\@pos2);
                $all_size += $num1;
                $pure_size += $num2;
                $redunt_size += $num3;
	}
	print STDERR "\ncheck redundance on bp level:\n  $all_size (all) = $pure_size (pure) + $redunt_size (redunt)\n\n" if($Verbose);
	
}

##conjoin the overlapped fragments, and caculate the redundant size
##usage: conjoin_fragment(\@pos);
##               my ($all_size,$pure_size,$redunt_size) = conjoin_fragment(\@pos);
sub Conjoin_fragment{
        my $pos_p = shift; ##point to the two dimension input array
        my $new_p = [];         ##point to the two demension result array

        my ($all_size, $pure_size, $redunt_size) = (0,0,0);

        return (0,0,0) unless(@$pos_p);

        foreach my $p (@$pos_p) {
                        ($p->[0],$p->[1]) = ($p->[0] <= $p->[1]) ? ($p->[0],$p->[1]) : ($p->[1],$p->[0]);
                        $all_size += abs($p->[0] - $p->[1]) + 1;
        }

        @$pos_p = sort {$a->[0] <=>$b->[0]} @$pos_p;
        push @$new_p, (shift @$pos_p);

        foreach my $p (@$pos_p) {
                        if ( ($p->[0] - $new_p->[-1][1]) <= 0 ) { # conjoin
                                        if ($new_p->[-1][1] < $p->[1]) {
                                                        $new_p->[-1][1] = $p->[1];
                                        }

                        }else{  ## not conjoin
                                        push @$new_p, $p;
                        }
        }
        @$pos_p = @$new_p;

        foreach my $p (@$pos_p) {
                        $pure_size += abs($p->[0] - $p->[1]) + 1;
        }

        $redunt_size = $all_size - $pure_size;
        return ($all_size,$pure_size,$redunt_size);
}

sub perfect_gene{
	my ($p)=@_;
	my $partial_num=0;
	foreach my $chr(sort keys %$p){
		foreach my $id(sort keys %{$p->{$chr}}){
			if ($p->{$chr}{$id}{type}=~/^no-/){
				delete $p->{$chr}{$id};
				$partial_num++;
			}
		}
	}
	print STDERR "perfecting the gene set:\n  remove $partial_num incomplete genes\n\n" if($Verbose);
}

#start,stop,triple,Ns
sub check_cds{
	my ($file,$p)=@_;
	my ($ns_num,$partal_num,$cds_ns_num)=(0,0,0);
	open IN,$file or die "Fail $file:$!";
	$/=">"; <IN>; $/="\n";
	while(<IN>){
		my $chr=$1 if(/^(\S+)/);
                $/=">";
                my $seq=<IN>;
                chomp $seq;
                $seq=~s/\s//g;
                $/="\n";
		next if not exists $p->{$chr};
		my $seq_p = $p->{$chr};
		delete $p->{$chr} if ( defined $Mini_seq && length($seq)<$Mini_seq );
		foreach my $gene_name (sort keys %$seq_p) {
                        my $gene_p = $seq_p->{$gene_name};
                        my $strand = $gene_p->{mRNA}[6];
			if (defined $Max_Ns ){
				my $gene_str = substr($seq,$gene_p->{mRNA}[3],abs($gene_p->{mRNA}[4]-$gene_p->{mRNA}[3])+1);
				if( filter_ns($gene_str) ){
					$ns_num++;
					delete $seq_p->{$gene_name};
					next;
				}
			}
                        my @exon = sort {$a->[3]<=>$b->[3]} @{$gene_p->{CDS}};
                        my $cds_str;
                        foreach my $p (@exon) {
                                my $str_len = abs($p->[4] - $p->[3]) + 1;
                                my $str_start = ($p->[3] < $p->[4]) ?  $p->[3] : $p->[4];
                                my $str = substr($seq,$str_start-1,$str_len);
                                $cds_str .= $str;
                        }
			$cds_str=Complement_Reverse($cds_str)  if($strand eq '-');
                        if( ! check_CDS($cds_str) ){
				$partal_num++;
                                delete $seq_p->{$gene_name};
                        }
			if (defined $CDS_ns){
				if( filter_ns($cds_str) ){
					$cds_ns_num++;
					delete $seq_p->{$gene_name};
					next;
				}
			}
                }
	} 
	close IN;
	print STDERR "check cds model on the sequence:\n" if ($Verbose);
        print STDERR "  remove $Max_Ns Ns genes $ns_num\n" if ( defined $Max_Ns && $Verbose);
	print STDERR "  remove wrong genes $partal_num\n" if ($Verbose);
	print STDERR "  remove $CDS_ns Ns cds $cds_ns_num\n" if ( defined $CDS_ns && $Verbose);
}


sub check_CDS{
        my $seq=shift;
#	my $id=shift;
	my $len=length($seq);
        my ($start,$end,$mid,$triple)=(0,0,0,0);
        $mid=1;
        my $len=length($seq);
        $triple=1 if($len%3 == 0);
        $start=1 if($seq=~/^ATG/);
        $end=1 if($seq=~/TAA$|TAG$|TGA$/);
#	print ">$id\n$seq\n";
        for (my $i=3; $i<$len-3; $i+=3) {
                my $codon=substr($seq,$i,3);
                if( ($codon eq 'TGA') || ($codon eq 'TAG') || ($codon eq 'TAA' )){
			$mid=0;
	#		print $id."\t".$len."\t".($i+1)."-".($i+3)."\n"; 
	#		print $seq."\n";die;
		}
        }
        if ($start && $mid && $end && $triple ) {
                return 1;
        }else{
                return 0;
        }
}

sub filter_ns{
	my ($seq)=@_;
	return ($seq=~/[N]{$Max_Ns}/i);
}

sub mini_cds{
	my ($p,$mini_cds)=@_;
	my $small_cds_num=0;
	foreach my $seq_name (sort keys %$p) {
                my $seq_p = $p->{$seq_name};
                foreach my $gene_name (sort keys %$seq_p) {
                        my $gene_p = $seq_p->{$gene_name};
                        my $cds_size;
                        foreach my $p (@{$gene_p->{CDS}}) {
                                $cds_size += abs($p->[4] - $p->[3]) + 1;
                        }
                        if($cds_size < $mini_cds){
                                delete $seq_p->{$gene_name};
                                $small_cds_num++;
                        }
                        
                }
        }
        print STDERR "Remove cds less than $mini_cds bp\n remove $small_cds_num small cds genes\n\n" if($Verbose);
}

sub Complement_Reverse{
	my $seq=shift;
	$seq=~tr/AGCTagct/TCGAtcga/;
	$seq=reverse($seq);
	return $seq;

}

