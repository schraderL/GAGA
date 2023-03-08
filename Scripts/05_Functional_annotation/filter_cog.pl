#!/usr/bin/perl

=head1 Name

filter_cog.pl  --  filter out cog genes which do not have COG ID.

=head1 Description

The purpose of this program is to reduce computing time.

=head1 Version

  Author: Fan Wei, fanw@genomics.org.cn
  Version: 1.0,  Date: 2008-5-21
  Note:

=head1 Usage
 
 perl filter_cog.pl <whog_file> <myva_file>
  --verbose   output running progress information to screen  
  --help      output help information to screen  

=head1 Exmple

 filter_cog.pl whog myva > cog_clean.fa

=cut

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;


my ($verbose,$help);
my %cog;			##store which has cog number
GetOptions(
	"verbose"=>\$verbose,
	"help"=>\$help,
);
die `pod2text $0` if(@ARGV < 2 || $help);

my ($whog,$myva) = @ARGV;

read_whog(\$whog,\%cog);


filter(\$myva,\%cog);



###################################

sub read_whog{
	my($whog_file,$cog_ref)=@_;
	open IN,$$whog_file or die $!;
	while(<IN>){
		if (/^  \w+:/) {
			s/^\s+//;
			chomp;
			my @t = split /\s+/;
			shift @t;
			foreach  (@t) {
				$cog_ref->{$_} = 1;
			}
		}
	}
	close IN;
}

sub filter{
	my ($myva,$cog_ref)=@_;
	open IN,$$myva or die $!;
	while(<IN>){
		chomp;
		(my $title=$_)=~s/^>//;
		$/=">";
		chomp(my $seq=<IN>);
		$/="\n";
		next unless $cog_ref->{$title};
		print ">$title\n$seq";
	}
	close IN;
}
		