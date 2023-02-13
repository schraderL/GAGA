#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
my $dirname = dirname(__FILE__);

# Get a GFF with only conserved GRs and rename the GFF to then merge with genewise

# usage: perl get_conserved_GRs.pl ABCENTH_GRs.txt ABCENTH_clean_GR_renamed_all_nofragment.gff3 


my ($line, $name, $nameout);
my $gff = $ARGV[1];
my $intable = $ARGV[0];
my $prefix = "";
if ($intable =~ /(\S+)\_GRs.txt/){
	$prefix = $1;
} else {die "Can't find prefix in $intable file\n";}

my $outgffall = "$prefix"."\_GR_all_renamed.gff3"; 
my $outgff = "$prefix"."\_GR_Dmelconserved_renamed.gff3"; 

# Read input table

my $dmelseq = "";

open(File, "<", $intable);
while(<File>){
	chomp;
	my $line = $_;
	next if ($line !~ /\S+/);
	my @subl = split (/\t/, $line);

	unless ($subl[9] =~ /None/){
		$dmelseq .= " $subl[4] ";
	}

}
close File;


# Reading GFF and generating final files

open (Resultsgff, ">", "$outgff");
print Resultsgff "##gff-version 3\n";

open (Resultsgffall, ">", "$outgffall");
print Resultsgffall "##gff-version 3\n";

open (GFFfile , "<", $gff); 
while (<GFFfile>) {
	chomp;
	$line = $_;
	next if ($line !~ /\S+/);

	if ($line =~ /^#/){
#		print Resultsgff "$line\n";
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

		my $nnamef = "Abcenth\_$genename";


		print Resultsgffall "$subline[0]\t$subline[1]\t$subline[2]\t$subline[3]\t$subline[4]\t$subline[5]\t$subline[6]\t$subline[7]\tParent=$nnamef"."$rest\;\n";

		if ($dmelseq =~ / $genename /){
			print Resultsgff "$subline[0]\t$subline[1]\t$subline[2]\t$subline[3]\t$subline[4]\t$subline[5]\t$subline[6]\t$subline[7]\tParent=$nnamef"."$rest\;\n";

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

		my $nnamef = "Abcenth\_$genename";


		print Resultsgffall "$subline[0]\t$subline[1]\tgene\t$subline[3]\t$subline[4]\t$subline[5]\t$subline[6]\t$subline[7]\tID=g"."$nnamef"."$rest;\n";		
		print Resultsgffall "$subline[0]\t$subline[1]\t$subline[2]\t$subline[3]\t$subline[4]\t$subline[5]\t$subline[6]\t$subline[7]\tID=$nnamef\;Parent=g"."$nnamef"."$rest;\n";		

		if ($dmelseq =~ / $genename /){
			print Resultsgff "$subline[0]\t$subline[1]\tgene\t$subline[3]\t$subline[4]\t$subline[5]\t$subline[6]\t$subline[7]\tID=g"."$nnamef"."$rest;\n";		
			print Resultsgff "$subline[0]\t$subline[1]\t$subline[2]\t$subline[3]\t$subline[4]\t$subline[5]\t$subline[6]\t$subline[7]\tID=$nnamef\;Parent=g"."$nnamef"."$rest;\n";		

		}

	}
}
close GFFfile;

close Resultsgff;
close Resultsgffall;

=head

#Encode proteins and CDS from the generated GFFs

system ("perl /home/projects/ku_00039/people/joeviz/programs/bitacora/Scripts/Tools/gff2fasta_v3.pl $genome $outgff $outprot ");
system ("sed \'s\/X\*\$\/\/\' $outprot\.pep.fasta > $outprot\.pep.fasta.tmp");
system ("mv $outprot\.pep.fasta.tmp $outprot\.pep.fasta");

system ("perl /home/projects/ku_00039/people/joeviz/programs/bitacora/Scripts/Tools/gff2fasta_v3.pl $genome $outgffall $outprotall ");
system ("sed \'s\/X\*\$\/\/\' $outprotall\.pep.fasta > $outprotall\.pep.fasta.tmp");
system ("mv $outprotall\.pep.fasta.tmp $outprotall\.pep.fasta");

system ("perl /home/projects/ku_00039/people/joeviz/programs/bitacora/Scripts/Tools/gff2fasta_v3.pl $genome $outgffcomp $outprotcomp ");
system ("sed \'s\/X\*\$\/\/\' $outprotcomp\.pep.fasta > $outprotcomp\.pep.fasta.tmp");
system ("mv $outprotcomp\.pep.fasta.tmp $outprotcomp\.pep.fasta");

#print Resultssum "$gagaid\t$gagasp{$gagaid}\t$totalor\t$com\t$comd\t$par\t$pard\t$parc\t$pse\t$psed\t$psec\t$l200\t$l100\t$l50\t$orconum\t$orconumcom\t$ninec\t$ninei\t$nineic\t$ninep\t$compard\t$ftot\t$fcom\t$fpar\t$fparc\t$fcomp\t$fpse\t$fninec\t$fninei\n";
#print Resultssum "$gagaid\t$gagasp{$gagaid}\t$totalor\t$com\t$comd\t$par\t$pard\t$parc\t$pse\t$psed\t$psec\t$l200\t$l100\t$l50\t$orconum\t$orconumcom\t$ninec\t$ninei\t$nineic\t$ninep\t$compard\t$ftot\t$fcomp\t$fpar\t$fpse\t$fninec\t$fninei\t$frpartial\t$frpse\n";
print Resultssum "$gagaid\t$gagasp{$gagaid}\t$totalor\t$com\t$comd\t$par\t$pard\t$parc\t$pse\t$psed\t$psec\t$l200\t$l100\t$l50\t$grcount\t$orconum\t$orconumcom\t$ninec\t$ninei\t$nineic\t$ninep\t$sugar\t$bitter\t$co2\t$gr43\t$compard\t$ftot\t$fcomp\t$denovoc\t$fpar\t$denovop\t$fpse\t$fninec\t$fninei\t$frpartial\t$frpse\n";

close Resultssum;

