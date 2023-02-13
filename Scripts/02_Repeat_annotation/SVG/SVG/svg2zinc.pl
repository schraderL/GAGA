#!/usr/bin/perl -w
#
#      svg2zinc.pl, a Perl script to display SVG file, or to generate
#                   scripts, modules or images from SVG files
#                   The result depends on the selected Backend
# 
#      Copyright (C) 2002-2004
#      Centre d'Études de la Navigation Aérienne
#      IntuiLab 2004
#
#      Authors: Christophe Mertz <mertz at intuilab at com>
#
# $Id: svg2zinc.pl,v 1.22 2004/05/01 09:19:33 mertz Exp $
#-----------------------------------------------------------------------------------

use vars qw( $VERSION);
($VERSION) = sprintf("%d.%02d", q$Revision: 1.22 $ =~ /(\d+)\.(\d+)/);

use strict;
use Carp;
use XML::Parser;
use SVG::SVG2zinc;
use Getopt::Long;

use File::Basename;

################ commande line options treatment
my ($out, $displayResult);
my $verbose = 0;
my $render = 1;
my $namespace = 0;
my $backend;

if (defined $ARGV[0]) {
    $backend = shift @ARGV;
    
    my $module = $backend;
    $module =~ s!::!\/!g ;
    if (!&findINC($backend.".pm") and !&findINC( "SVG/SVG2zinc/Backend/".$module.".pm")) {
	&usage ("$backend not found as a module or as SVG::SVG2zinc::Backend::$backend in perl path: @INC");
    }
} else {
    &usage ("required first parameter: Backend name") unless (defined $ARGV[0]);
}

Getopt::Long::Configure("pass_through");
GetOptions("help" => \&usage,
	   "out:s" => \$out,
	   "verbose" => \$verbose,
	   "render:i" => \$render,
	   "namespace" => \$namespace,
           );

my @otherParams;
while (@ARGV) {
    last if $ARGV[0] !~ /^-\w+/ or $#ARGV eq 0;
    push @otherParams, (shift @ARGV => shift @ARGV);
}

#print "otherparam = @otherParams\n";
    
my $svgfile;
if ($#ARGV >= 0) {
    $svgfile = shift @ARGV;
} else {
    &usage("missing svgfile param");
}

my @outParam;
if ($#ARGV >= 0) {
    @outParam = ( -out => shift @ARGV);
}

if (@outParam eq 0) {
    if ( $backend eq "PerlScript" or
	 $backend eq "TclScript" or
	 $backend eq "PerlClass" or
	 $backend eq "Image" ) {
	&usage ("The backend $backend requires and outfile as last parameter of $0");
    } elsif ($backend ne "Display" and $backend ne "Print") {
	warn "You did not specify an out file for the user defined Backend '$backend'. This may be an error\n";
    }
}

# print "out = @outParam\n";

sub usage {
    my ($error) = @_;
    print "$error\n\n" if defined $error;
    print basename($0), " v$VERSION\tSVG::SVG2zinc.pm v$SVG::SVG2zinc::VERSION\n";
    print "Usage:\n";
    print basename($0), " Backend options svgfile [outfile]\n";
    print " where options are :\n";
    print "  -help           : to get this little help\n";
    print "  -verbose        : to get more warning or prints\n";
    print "  -render (0|1|2) : to select Zinc rendering mode (default to 1)\n";
    print "  -namespace      : to treat XML namespace\n";
    print "  Backend may be Display, PerlScript, PerlClass, Image, Print\n";
    print "  or any user defined backend\n";
    exit;
}

&usage ("too much parameters: @ARGV")  unless $#ARGV < 0;
&usage ("Bad value ($render) for -render option. Must be 0, 1 or 2") unless ($render == 0 or $render == 1 or $render == 2);



my $file;
my $zinc;
my $mw;
my ($WIDTH,$HEIGHT) = (600,600);
my $top_group;

&SVG::SVG2zinc::parsefile($svgfile, $backend,
			  -verbose => $verbose,
			  -render => $render,
			  -namespace => $namespace,
			  @otherParams, @outParam);


__END__

=head1 NAME

svg2zinc.pl - displays a svg file or generates perl script, perl module or jpg,gif,png images

=head1 SYNOPSIS

B<svg2zinc.pl> -help

B<svg2zinc.pl> B<Display> [options] [-render 0|1|2] svgfile

B<svg2zinc.pl> B<PerlScript> [options] [-render 0|1|2] svgfile outfile

B<svg2zinc.pl> B<PerlClass> [options] [-prefix a_prefix] svgfile outfile

B<svg2zinc.pl> B<Image> [options] [-w pixels] [-h pixels] [-ratio %] svgfile outfile

B<svg2zinc.pl> B<UserBackend> [options] [userBackendOptions] svgfile [outfile]


=head1 DESCRIPTION

svg2zinc.pl is a perl script which examplifie the use of SVG::SVG2zinc to read
and interprete SVG file. Depending on the B<Backend> svg2zinc.pl can be used to:

=over

=item B<Display>

Displays a SVG file in a TkZinc windows.

=item B<PerlScript>

Generates a perl script able to display the SVG file in a TkZinc windows

=item B<PerlClass>

Generates a perl class module usable by other TkZinc based applications

=item B<TclSscript>

Generates a tcl script able to display the SVG file in a TkZinc windows  module usable by other TkZinc based applications. BEWARE: it is still higly experimental

=item B<Image>

Generates an image file of usual bitmap format such as gif, jpeg or png, depending on
the outfile suffix. Accepted suffixes are those of the imagemagic package.
B<-w> option defines the image width, B<-h> option defines the image
height and B<-ratio> option defines the image ratio towards the svg defined size.
This requires the imagemagic package

=item B<Print>

Just print on stdout some incomplete lines of code similar of the code produced (or evaluated)
with other backends. This backend is not very usefull except for SVG2zinc debug or as an exemple
of oversimplified backend.

=item B<UserBackend>

Generate data or file according to the provided B<UserBackend>. This backend can be either
a F<UserBackend.pm> file or a sub class of SVG::SVG2zinc::Backend. In both cases, they must be
reachable in the perl pathes. The options given to this backend
should be documented in this backend. To write such a backend you should read the
SVG::SVG2zinc::Backend(3pm) man pages and look to other backends as examples.

In this case, the [userBackendOptions] options must be given in the following form:
B<-option1 val1 -option2 val2> with fullnames for options and a value for each option.
These options are the same than those accepted by the F<UserBackend.pm> module

=back


=head1 OPTIONS

The generic options available for all backends are :

=over

=item B<-verbose>

to get more information on stdout

=item B<-namespace>

To treat xml file where SVG is a namespace

=back


=head1 BUGS and LIMITATIONS

Mainly the same bugs and limitations of SVG::SVG2zinc(3pm)

=head1 SEE ALSO

SVG::SVG2zinc(3pm) SVG::SVG2zinc::Backend(3pm) Tk::Zinc(3pm) 

TkZinc is available at www.openatc.org/zinc

=head1 AUTHORS

Christophe Mertz <mertz at intuilab dot com> with some help from Daniel Etienne <etienne at cena dot fr>

=head1 COPYRIGHT
    
CENA (C) 2002-2004 IntuiLab 2004

This program is free software; you can redistribute it and/or modify it under the term of the LGPL licence.

=cut

