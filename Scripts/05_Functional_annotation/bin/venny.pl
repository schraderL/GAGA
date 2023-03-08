#!/usr/bin/perl

#Author: liaoxinhui@genomics.org.cn
#Date:2013-04-11 

use strict;
use warnings;
use Getopt::Long;
use FindBin '$Bin';

my ($infile, $header, $name, $font_family, $name_font_size, $num_font_size, $colour, $desc, $outdir, $imgname, $help);

GetOptions(
        "infile:s" => \$infile,
	"header!" =>\$header,
        "name:s" => \$name,
	"desc:s" => \$desc,
        "color!" => \$colour,
	"font_family:s" => \$font_family,
	"name_font_size:f" => \$name_font_size,
	"num_font_size:f" => \$num_font_size,
        "outdir:s" => \$outdir,
		"imgname:s" => \$imgname,
        "help|?" => \$help
);
if (!defined $infile || defined $help) {
	die <<USAGE

Description: Generate tables & graph for Venny (2D,3D,4D,5D)
Usage: perl $0 [options]
Options:
	* -infile			input files with IDs in column 1, separated by ","
	* -name				name of each category, separated by ","
	  -desc				description file for each ID, format: GeneID, description1, description2... separated by "\\t" (header needed)
	  -header			only active when files with header
	  -color			flag to filled categories with predefined colors
	  -font_family			font family, default is Verdana
	  -name_font_size		font size of name, default is 18
	  -num_font_size		font size of numbers, default is 13
	  -outdir			directory of output files, default is current directory
	  -imgname          output image name, default is "Venn-[input.number]D"
	  -help 			help information
E.g.:
	perl $0 -infile A.xls,B.xls.C.xls,D.xls -name A,B,C,D -color

USAGE
}
my @files = split /,+/,$infile;
my @names = split /,+/,$name;
die "\nWarning: Files' number not equal to Names' number!\nWarning: Please double check the parameters!\n\n" if (@files != @names);
$font_family ||= "Verdana";
$name_font_size ||= 18;
$num_font_size ||= 13;
$outdir ||= "./";
`mkdir -p $outdir` unless (-d $outdir);
$imgname ||= "Venn-".scalar(@names)."D";
my @colours = ("#0000ff","#ff0000","#00cc00","#cc9900","#0099ff");
#read ID files
my %all_ids;
foreach my $index (0..$#files) {
	open INFILE,$files[$index];
	<INFILE> if defined $header;
	my %check_dup = ();
	while (<INFILE>) {
		chomp;my @temp = split /\s+/;
		next if ($check_dup{$temp[0]});
		$check_dup{$temp[0]} = 1;
		$all_ids{$temp[0]} .= "$names[$index]_";
	}
	close INFILE;
}
#read description files
my %desc; my $desc_header;
if (defined $desc) {
	open DESC,$desc;
	my $head = <DESC>;
	my @head = split /\t+/,$head;shift @head; $desc_header = join "\t",@head;
	while (<DESC>) {
		chomp;my @temp = split /\t+/,$_;
		my $id = shift @temp;
		$desc{$id} = join "\t",@temp;
	}
	close DESC;
}
#grouping & print
my %all_groups;
foreach my $id (keys %all_ids) {
	chop $all_ids{$id};
	my $id_desc = $desc{$id} ? $desc{$id} : "NULL";
	$all_groups{$all_ids{$id}}{'ids'} .= (defined $desc) ? "$id\t$id_desc\n" : "$id\n";
	$all_groups{$all_ids{$id}}{'num'} += 1;
}
foreach my $group (keys %all_groups) {
	open OUTFILE,">$outdir/$group.List.xls";
	print OUTFILE (defined $desc) ? "GeneID\t$desc_header" : "GeneID\n";
	print OUTFILE $all_groups{$group}{'ids'};
	close OUTFILE;
}
#draw graph
&draw_venn2(\%all_groups,\@names) if @names == 2;
&draw_venn3(\%all_groups,\@names) if @names == 3;
&draw_venn4(\%all_groups,\@names) if @names == 4;
&draw_venn5(\%all_groups,\@names) if @names == 5;

exit 0;
#sub programs
sub draw_venn2 {
	my ($groups, $names) = @_;
	my ($g12,$g1,$g2) = (
		$$groups{"$$names[0]_$$names[1]"} ? $$groups{"$$names[0]_$$names[1]"}{'num'} : 0,
		$$groups{"$$names[0]"} ? $$groups{"$$names[0]"}{'num'} : 0,
		$$groups{"$$names[1]"} ? $$groups{"$$names[1]"}{'num'} : 0
	);
	my $stroke_width = defined $colour ? 0.1 : 2;
        my ($colour1,$colour2) = @colours;
        my $opacity = defined $colour ? 0.3 : 0;
	open SVG,">$outdir/$imgname.svg";
	print SVG <<SVG;
<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<svg version="1.1"
   xmlns="http://www.w3.org/2000/svg"
   xmlns:xlink="http://www.w3.org/1999/xlink"
   width="400" height="300"
>
<circle cx="150" cy="170" r="100" style="stroke:black;stroke-width:$stroke_width;fill:$colour1;fill-opacity:$opacity"/>
<circle cx="250" cy="170" r="100" style="stroke:black;stroke-width:$stroke_width;fill:$colour2;fill-opacity:$opacity"/>
<text x="150" y="50" font-size="$name_font_size" text-anchor="end" font-family="$font_family">$$names[0]</text>
<text x="250" y="50" font-size="$name_font_size" text-anchor="start" font-family="$font_family">$$names[1]</text>
<text x="100" y="170" font-size="$num_font_size" text-anchor="middle" font-family="$font_family">$g1</text>
<text x="200" y="170" font-size="$num_font_size" text-anchor="middle" font-family="$font_family">$g12</text>
<text x="300" y="170" font-size="$num_font_size" text-anchor="middle" font-family="$font_family">$g2</text>
</svg>
SVG
system("/usr/bin/convert  $outdir/$imgname.svg $outdir/$imgname.png");
}

sub draw_venn3 {
	my ($groups, $names) = @_;
	my ($g123,$g12,$g13,$g23,$g1,$g2,$g3) = (
		$$groups{"$$names[0]_$$names[1]_$$names[2]"} ? $$groups{"$$names[0]_$$names[1]_$$names[2]"}{'num'} : 0,
		$$groups{"$$names[0]_$$names[1]"} ? $$groups{"$$names[0]_$$names[1]"}{'num'} : 0,
		$$groups{"$$names[0]_$$names[2]"} ? $$groups{"$$names[0]_$$names[2]"}{'num'} : 0,
		$$groups{"$$names[1]_$$names[2]"} ? $$groups{"$$names[1]_$$names[2]"}{'num'} : 0,
		$$groups{"$$names[0]"} ? $$groups{"$$names[0]"}{'num'} : 0,
		$$groups{"$$names[1]"} ? $$groups{"$$names[1]"}{'num'} : 0,
		$$groups{"$$names[2]"} ? $$groups{"$$names[2]"}{'num'} : 0
	);
	my $stroke_width = defined $colour ? 0.1 : 2;
	my ($colour1,$colour2,$colour3) = @colours;
	my $opacity = defined $colour ? 0.3 : 0;
	open SVG,">$outdir/$imgname.svg";
print SVG <<SVG;
<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<svg version="1.1"
   xmlns="http://www.w3.org/2000/svg"
   xmlns:xlink="http://www.w3.org/1999/xlink"
   width="500" height="500"
>
<circle cx="200" cy="200" r="100" style="stroke:black;stroke-width:$stroke_width;fill:$colour1;fill-opacity:$opacity"/>
<circle cx="300" cy="200" r="100" style="stroke:black;stroke-width:$stroke_width;fill:$colour2;fill-opacity:$opacity"/>
<circle cx="250" cy="290" r="100" style="stroke:black;stroke-width:$stroke_width;fill:$colour3;fill-opacity:$opacity"/>
<text x="140" y="110" font-size="$name_font_size" text-anchor="end" font-family="$font_family">$$names[0]</text>
<text x="360" y="110" font-size="$name_font_size" text-anchor="start" font-family="$font_family">$$names[1]</text>
<text x="250" y="420" font-size="$name_font_size" text-anchor="middle" font-family="$font_family">$$names[2]</text>
<text x="150" y="190" font-size="$num_font_size" text-anchor="middle" font-family="$font_family">$g1</text>
<text x="250" y="160" font-size="$num_font_size" text-anchor="middle" font-family="$font_family">$g12</text>
<text x="350" y="190" font-size="$num_font_size" text-anchor="middle" font-family="$font_family">$g2</text>
<text x="190" y="270" font-size="$num_font_size" text-anchor="middle" font-family="$font_family">$g13</text>
<text x="250" y="240" font-size="$num_font_size" text-anchor="middle" font-family="$font_family">$g123</text>
<text x="310" y="270" font-size="$num_font_size" text-anchor="middle" font-family="$font_family">$g23</text>
<text x="250" y="340" font-size="$num_font_size" text-anchor="middle" font-family="$font_family">$g3</text>
</svg>
SVG
system("/usr/bin/convert $outdir/$imgname.svg $outdir/$imgname.png");
}

sub draw_venn4 {
	my ($groups, $names) = @_;
	my ($g1234,$g123,$g124,$g134,$g234,$g12,$g13,$g14,$g23,$g24,$g34,$g1,$g2,$g3,$g4) = (
		$$groups{"$$names[0]_$$names[1]_$$names[2]_$$names[3]"} ? $$groups{"$$names[0]_$$names[1]_$$names[2]_$$names[3]"}{'num'} : 0,
		$$groups{"$$names[0]_$$names[1]_$$names[2]"} ? $$groups{"$$names[0]_$$names[1]_$$names[2]"}{'num'} : 0,
		$$groups{"$$names[0]_$$names[1]_$$names[3]"} ? $$groups{"$$names[0]_$$names[1]_$$names[3]"}{'num'} : 0,
		$$groups{"$$names[0]_$$names[2]_$$names[3]"} ? $$groups{"$$names[0]_$$names[2]_$$names[3]"}{'num'} : 0,
		$$groups{"$$names[1]_$$names[2]_$$names[3]"} ? $$groups{"$$names[1]_$$names[2]_$$names[3]"}{'num'} : 0,
		$$groups{"$$names[0]_$$names[1]"} ? $$groups{"$$names[0]_$$names[1]"}{'num'} : 0,
		$$groups{"$$names[0]_$$names[2]"} ? $$groups{"$$names[0]_$$names[2]"}{'num'} : 0,
		$$groups{"$$names[0]_$$names[3]"} ? $$groups{"$$names[0]_$$names[3]"}{'num'} : 0,
		$$groups{"$$names[1]_$$names[2]"} ? $$groups{"$$names[1]_$$names[2]"}{'num'} : 0,
		$$groups{"$$names[1]_$$names[3]"} ? $$groups{"$$names[1]_$$names[3]"}{'num'} : 0,
		$$groups{"$$names[2]_$$names[3]"} ? $$groups{"$$names[2]_$$names[3]"}{'num'} : 0,
		$$groups{"$$names[0]"} ? $$groups{"$$names[0]"}{'num'} : 0,
		$$groups{"$$names[1]"} ? $$groups{"$$names[1]"}{'num'} : 0,
		$$groups{"$$names[2]"} ? $$groups{"$$names[2]"}{'num'} : 0,
		$$groups{"$$names[3]"} ? $$groups{"$$names[3]"}{'num'} : 0
	);
	my $stroke_width = defined $colour ? 0.1 : 2;
        my ($colour1,$colour2,$colour3,$colour4) = @colours;
        my $opacity = defined $colour ? 0.3 : 0;
        open SVG,">$outdir/$imgname.svg";
	print SVG <<SVG;
<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<svg version="1.1"
   xmlns="http://www.w3.org/2000/svg"
   xmlns:xlink="http://www.w3.org/1999/xlink"
   width="680" height="680"
>
<ellipse cx="250" cy="350" rx="190" ry="110" style="fill:$colour1;fill-opacity:$opacity;stroke:black;stroke-width:$stroke_width" transform="rotate(40,250,350)"/>
<ellipse cx="325" cy="295" rx="180" ry="85" style="fill:$colour2;fill-opacity:$opacity;stroke:black;stroke-width:$stroke_width" transform="rotate(40,325,295)"/>
<ellipse cx="345" cy="295" rx="180" ry="85" style="fill:$colour3;fill-opacity:$opacity;stroke:black;stroke-width:$stroke_width" transform="rotate(-40,345,295)"/>
<ellipse cx="425" cy="350" rx="190" ry="110" style="fill:$colour4;fill-opacity:$opacity;stroke:black;stroke-width:$stroke_width" transform="rotate(-40,425,350)"/>
<text x="120" y="190" font-size="$name_font_size" text-anchor="middle" font-family="$font_family">$$names[0]</text>
<text x="240" y="150" font-size="$name_font_size" text-anchor="middle" font-family="$font_family">$$names[1]</text>
<text x="420" y="150" font-size="$name_font_size" text-anchor="middle" font-family="$font_family">$$names[2]</text>
<text x="550" y="190" font-size="$name_font_size" text-anchor="middle" font-family="$font_family">$$names[3]</text>
<text x="240" y="190" font-size="$num_font_size" text-anchor="middle" font-family="$font_family">$g2</text>
<text x="420" y="190" font-size="$num_font_size" text-anchor="middle" font-family="$font_family">$g3</text>
<text x="160" y="300" font-size="$num_font_size" text-anchor="middle" font-family="$font_family">$g1</text>
<text x="220" y="250" font-size="$num_font_size" text-anchor="middle" font-family="$font_family">$g12</text>
<text x="335" y="240" font-size="$num_font_size" text-anchor="middle" font-family="$font_family">$g23</text>
<text x="450" y="250" font-size="$num_font_size" text-anchor="middle" font-family="$font_family">$g34</text>
<text x="510" y="300" font-size="$num_font_size" text-anchor="middle" font-family="$font_family">$g4</text>
<text x="275" y="300" font-size="$num_font_size" text-anchor="middle" font-family="$font_family">$g123</text>
<text x="335" y="350" font-size="$num_font_size" text-anchor="middle" font-family="$font_family">$g1234</text>
<text x="400" y="300" font-size="$num_font_size" text-anchor="middle" font-family="$font_family">$g234</text>
<text x="235" y="380" font-size="$num_font_size" text-anchor="middle" font-family="$font_family">$g13</text>
<text x="290" y="405" font-size="$num_font_size" text-anchor="middle" font-family="$font_family">$g134</text>
<text x="380" y="405" font-size="$num_font_size" text-anchor="middle" font-family="$font_family">$g124</text>
<text x="440" y="380" font-size="$num_font_size" text-anchor="middle" font-family="$font_family">$g24</text>
<text x="335" y="460" font-size="$num_font_size" text-anchor="middle" font-family="$font_family">$g14</text>
</svg>
SVG
system("/usr/bin/convert $outdir/$imgname.svg $outdir/$imgname.png");
}

sub draw_venn5 {
	my ($groups, $names) = @_;
	my ($g12345,$g1234,$g1235,$g1245,$g1345,$g2345,$g123,$g124,$g125,$g134,$g135,$g145,$g234,$g235,$g245,$g345,$g12,$g13,$g14,$g15,$g23,$g24,$g25,$g34,$g35,$g45,$g1,$g2,$g3,$g4,$g5) =(
		$$groups{"$$names[0]_$$names[1]_$$names[2]_$$names[3]_$$names[4]"} ? $$groups{"$$names[0]_$$names[1]_$$names[2]_$$names[3]_$$names[4]"}{'num'} : 0,
		$$groups{"$$names[0]_$$names[1]_$$names[2]_$$names[3]"} ? $$groups{"$$names[0]_$$names[1]_$$names[2]_$$names[3]"}{'num'} : 0,
		$$groups{"$$names[0]_$$names[1]_$$names[2]_$$names[4]"} ? $$groups{"$$names[0]_$$names[1]_$$names[2]_$$names[4]"}{'num'} : 0,
		$$groups{"$$names[0]_$$names[1]_$$names[3]_$$names[4]"} ? $$groups{"$$names[0]_$$names[1]_$$names[3]_$$names[4]"}{'num'} : 0,
		$$groups{"$$names[0]_$$names[2]_$$names[3]_$$names[4]"} ? $$groups{"$$names[0]_$$names[2]_$$names[3]_$$names[4]"}{'num'} : 0,
		$$groups{"$$names[1]_$$names[2]_$$names[3]_$$names[4]"} ? $$groups{"$$names[1]_$$names[2]_$$names[3]_$$names[4]"}{'num'} : 0,
		$$groups{"$$names[0]_$$names[1]_$$names[2]"} ? $$groups{"$$names[0]_$$names[1]_$$names[2]"}{'num'} : 0,
		$$groups{"$$names[0]_$$names[1]_$$names[3]"} ? $$groups{"$$names[0]_$$names[1]_$$names[3]"}{'num'} : 0,
		$$groups{"$$names[0]_$$names[1]_$$names[4]"} ? $$groups{"$$names[0]_$$names[1]_$$names[4]"}{'num'} : 0,
		$$groups{"$$names[0]_$$names[2]_$$names[3]"} ? $$groups{"$$names[0]_$$names[2]_$$names[3]"}{'num'} : 0,
		$$groups{"$$names[0]_$$names[2]_$$names[4]"} ? $$groups{"$$names[0]_$$names[2]_$$names[4]"}{'num'} : 0,
		$$groups{"$$names[0]_$$names[3]_$$names[4]"} ? $$groups{"$$names[0]_$$names[3]_$$names[4]"}{'num'} : 0,
		$$groups{"$$names[1]_$$names[2]_$$names[3]"} ? $$groups{"$$names[1]_$$names[2]_$$names[3]"}{'num'} : 0,
		$$groups{"$$names[1]_$$names[2]_$$names[4]"} ? $$groups{"$$names[1]_$$names[2]_$$names[4]"}{'num'} : 0,
		$$groups{"$$names[1]_$$names[3]_$$names[4]"} ? $$groups{"$$names[1]_$$names[3]_$$names[4]"}{'num'} : 0,
		$$groups{"$$names[2]_$$names[3]_$$names[4]"} ? $$groups{"$$names[2]_$$names[3]_$$names[4]"}{'num'} : 0,
		$$groups{"$$names[0]_$$names[1]"} ? $$groups{"$$names[0]_$$names[1]"}{'num'} : 0,
		$$groups{"$$names[0]_$$names[2]"} ? $$groups{"$$names[0]_$$names[2]"}{'num'} : 0,
		$$groups{"$$names[0]_$$names[3]"} ? $$groups{"$$names[0]_$$names[3]"}{'num'} : 0,
		$$groups{"$$names[0]_$$names[4]"} ? $$groups{"$$names[0]_$$names[4]"}{'num'} : 0,
		$$groups{"$$names[1]_$$names[2]"} ? $$groups{"$$names[1]_$$names[2]"}{'num'} : 0,
		$$groups{"$$names[1]_$$names[3]"} ? $$groups{"$$names[1]_$$names[3]"}{'num'} : 0,
		$$groups{"$$names[1]_$$names[4]"} ? $$groups{"$$names[1]_$$names[4]"}{'num'} : 0,
		$$groups{"$$names[2]_$$names[3]"} ? $$groups{"$$names[2]_$$names[3]"}{'num'} : 0,
		$$groups{"$$names[2]_$$names[4]"} ? $$groups{"$$names[2]_$$names[4]"}{'num'} : 0,
		$$groups{"$$names[3]_$$names[4]"} ? $$groups{"$$names[3]_$$names[4]"}{'num'} : 0,
		$$groups{"$$names[0]"} ? $$groups{"$$names[0]"}{'num'} : 0,
		$$groups{"$$names[1]"} ? $$groups{"$$names[1]"}{'num'} : 0,
		$$groups{"$$names[2]"} ? $$groups{"$$names[2]"}{'num'} : 0,
		$$groups{"$$names[3]"} ? $$groups{"$$names[3]"}{'num'} : 0,
		$$groups{"$$names[4]"} ? $$groups{"$$names[4]"}{'num'} : 0,
	);
	my $stroke_width = defined $colour ? 0.1 : 2;
        my ($colour1,$colour2,$colour3,$colour4,$colour5) = @colours;
        my $opacity = defined $colour ? 0.3 : 0;
        open SVG,">$outdir/$imgname.svg";
        print SVG <<SVG;
<?xml version="1.0" standalone="no"?>
<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN"
  "http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd">
<svg version="1.1" width="600" height="600" viewBox="-362 -388 746 742" xmlns:xlink="http://www.w3.org/1999/xlink" xmlns="http://www.w3.org/2000/svg" preserveAspectRatio="xMidYMid meet" zoomAndPan="magnify">
 <title>Radially-symmetrical Five-set Venn Diagram</title>
 <desc>Devised by Branko Gruenbaum and rendered by CMG Lee.</desc>
 <defs>
  <ellipse id="ellipse" cx="36" cy="-56" rx="160" ry="320"/>
  <g id="ellipses">
   <use xlink:href="#ellipse" fill="$colour1"/>
   <use xlink:href="#ellipse" fill="$colour2" transform="rotate(72)"/>
   <use xlink:href="#ellipse" fill="$colour3" transform="rotate(144)"/>
   <use xlink:href="#ellipse" fill="$colour4" transform="rotate(216)"/>
   <use xlink:href="#ellipse" fill="$colour5" transform="rotate(288)"/>
  </g>
 </defs>
 <use xlink:href="#ellipses" fill-opacity="$opacity" stroke="black" stroke-width="$stroke_width"/>
 <g text-anchor="middle" font-family="$font_family" font-size="$num_font_size">
  <text x="30" y="-300" dy="0.7ex" font-size="$name_font_size" startOffset="0">$$names[0]</text>
  <text x="300" y="-60" dy="0.7ex" font-size="$name_font_size" startOffset="0">$$names[1]</text>
  <text x="160" y="280" dy="0.7ex" font-size="$name_font_size" startOffset="0">$$names[2]</text>
  <text x="-220" y="220" dy="0.7ex" font-size="$name_font_size" startOffset="0">$$names[3]</text>
  <text x="-280" y="-130" dy="0.7ex" font-size="$name_font_size" startOffset="0">$$names[4]</text>
  <text x="25" y="-270" dy="0.7ex" startOffset="0">$g1</text>
  <text x="280" y="-100" dy="0.7ex" startOffset="0">$g2</text>
  <text x="170" y="230" dy="0.7ex" startOffset="0">$g3</text>
  <text x="-160" y="240" dy="0.7ex" startOffset="0">$g4</text>
  <text x="-280" y="-90" dy="0.7ex" startOffset="0">$g5</text>
  <text x="180" y="-130" dy="0.7ex" startOffset="0">$g12</text>
  <text x="40" y="230" dy="0.7ex" startOffset="0">$g13</text>
  <text x="100" y="-200" dy="0.7ex" startOffset="0">$g14</text>
  <text x="-80" y="-215" dy="0.7ex" startOffset="0">$g15</text>
  <text x="190" y="125" dy="0.7ex" startOffset="0">$g23</text>
  <text x="-190" y="120" dy="0.7ex" startOffset="0">$g24</text>
  <text x="230" y="40" dy="0.7ex" startOffset="0">$g25</text>
  <text x="-60" y="220" dy="0.7ex" startOffset="0">$g34</text>
  <text x="-170" y="-150" dy="0.7ex" startOffset="0">$g35</text>
  <text x="-222" y="0" dy="0.7ex" startOffset="0">$g45</text>
  <text x="90" y="150" dy="0.7ex" startOffset="0">$g123</text>
  <text x="148" y="-153" dy="0.7ex" startOffset="0">$g124</text>
  <text x="170" y="-20" dy="0.7ex" startOffset="0">$g125</text>
  <text x="-33" y="208" dy="0.7ex" startOffset="0">$g134</text>
  <text x="-93" y="-193" dy="0.7ex" startOffset="0">$g135</text>
  <text x="20" y="-180" dy="0.7ex" startOffset="0">$g145</text>
  <text x="-120" y="120" dy="0.7ex" startOffset="0">$g234</text>
  <text x="190" y="100" dy="0.7ex" startOffset="0">$g235</text>
  <text x="-211" y="32" dy="0.7ex" startOffset="0">$g245</text>
  <text x="-150" y="-80" dy="0.7ex" startOffset="0">$g345</text>
  <text x="-30" y="160" dy="0.7ex" startOffset="0">$g1234</text>
  <text x="140" y="80" dy="0.7ex" startOffset="0">$g1235</text>
  <text x="120" y="-100" dy="0.7ex" startOffset="0">$g1245</text>
  <text x="-60" y="-140" dy="0.7ex" startOffset="0">$g1345</text>
  <text x="-160" y="20" dy="0.7ex" startOffset="0">$g2345</text>
  <text x="0" y="0" dy="0.7ex" startOffset="0">$g12345</text>
 </g>
</svg>
SVG
system("/usr/bin/convert $outdir/$imgname.svg $outdir/$imgname.png");
}
