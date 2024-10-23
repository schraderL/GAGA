package SVG::SVG2zinc::Backend::Tcl;

#	Backend Class for SVG2zinc
# 
#	Copyright 2003-2004
#	Centre d'Études de la Navigation Aérienne
#
#	Author: Christophe Mertz <mertz at intuilab dot com>
#
#       A module for code translation from perl to tcl generation
#
# $Id: Tcl.pm,v 1.2 2004/05/01 09:19:34 mertz Exp $
#############################################################################


use vars qw( $VERSION);
($VERSION) = sprintf("%d.%02d", q$Revision: 1.2 $ =~ /(\d+)\.(\d+)/);

@ISA = qw( Exporter );
@EXPORT = qw( perl2tcl );

use strict;
use Carp;


sub perl2tcl {
    my (@lines) = @_;
    my @res;
    foreach my $l (@lines) {

	$l =~ s/->(\w*)\((.*)\)/\$w\.zinc $1 $2/g;  # ->add(....) => $w.zinc add ...

	$l =~ s/\s*,\s*/ /g;  # replacing commas by spaces
	$l =~ s/\s*=>\s*/ /g;  # replacing => by spaces

 	$l =~ s/\s*\'([^\s;]+)\'\s*/ $1 /g ;  # removing single-quotes around string without spaces
	$l =~ s/\s*\"([^\s;]+)\"\s*/ $1 /g ;  # removing double-quotes around string without spaces
	$l =~ s/([\"\s])\#/$1\\\#/g ;  # prefixing # by a slash

	$l =~ s/\[/\{/g;  # replacing [ by }
	$l =~ s/\]/\}/g;  # replacing ] by }
	$l =~ s/\{\s+/\{/g;  # removing spaces after {
	$l =~ s/\s+\}/\}/g;  # removing spaces before }

	$l =~ s/-tags \{(\S+)\}/-tags $1/g;  # -tags {toto}  ==>> -tags toto
	$l =~ s/\'/\"/g;  # replacing all single quotes by double quotes

	$l = &hack($l);
    
	$l =~ s/\s+/ /g;  # dangerous: removing multiple occurences of blanks

	$l =~ s/^\s+//;   # removing blanks at the beginning
	$l =~ s/\s+$//;   # removing trailing blanks
	$l =~ s/\s*;$//;  # removing trailing ;
	push @res, $l;
    }
    return (@res);
}


# this routine is used to do some special code transformation,
# due to soem discrepancies between tcl/tk and perl/tk
# the following code is more or less dependant from the generated
# code by SVG2zinc.pm
#
# We assume is code has already been tcl-ised
sub hack {
    my ($l) = @_;

    if ($l =~ /^\$w\.zinc fontCreate/) {
	# this works because I know how fontCreate is used in SVG2zinc
	$l =~ s/\$w\.zinc fontCreate/font create/;
	$l =~ s/-weight medium/-weight normal/;
    }

    return $l;
}

1;

