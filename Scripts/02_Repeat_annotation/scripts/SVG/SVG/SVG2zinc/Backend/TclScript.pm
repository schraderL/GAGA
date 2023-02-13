package SVG::SVG2zinc::Backend::TclScript;

#	Backend Class for SVG2zinc
# 
#	Copyright 2003
#	Centre d'Études de la Navigation Aérienne
#
#	Author: Christophe Mertz <mertz at intuilab dot com>
#
#       A concrete class for code generation for Tcl Scripts
#
# $Id: TclScript.pm,v 1.3 2004/05/01 09:19:34 mertz Exp $
#############################################################################

use strict;
use Carp;

use SVG::SVG2zinc::Backend;
use SVG::SVG2zinc::Backend::Tcl;
use File::Basename;

use vars qw( $VERSION @ISA  );
@ISA = qw( SVG::SVG2zinc::Backend );
($VERSION) = sprintf("%d.%02d", q$Revision: 1.3 $ =~ /(\d+)\.(\d+)/);

sub new {
    my ($class, %passed_options) = @_;
    my $self = {};
    bless $self, $class;
    $self->{-render} = defined $passed_options{-render} ? delete $passed_options{-render} : 1;
    $self->_initialize(%passed_options);
    return $self;
}


sub treatLines {
    my ($self,@lines) = @_;
    foreach my $l (@lines) {
	$self->printLines( &perl2tcl($l) );
    }
}

sub fileHeader {
    my ($self) = @_;
    my $svgfile = $self->{-in};
    my $svgfilename = basename($svgfile);
    $svgfilename =~ s/\./_/g;
    my ($svg2zincPackage) = caller;
    my $VERSION = eval ( "\$".$svg2zincPackage."::VERSION" );
    $self->printLines('#!/bin/sh
# the next line restarts using wish \
    exec wish "$0" "$@"
');

    $self->printLines("

####### This Tcl script file has been generated
####### from $svgfile
####### by SVG::SVG2zinc.pm Version: $VERSION

");

    $self->printLines('
#
# Locate the zinc top level directory.
#
set zincRoot [file join [file dirname [info script]] ..]

#
# And adjust the paths accordingly.
#
lappend auto_path $zincRoot
set zinc_library $zincRoot

package require Tkzinc 3.2

## here we should import img for reading jpeg, png, gif files

');

    my $render = $self->{-render};
    $self->printLines(
<<HEADER
set w .$svgfilename
## catch {destroy \$w}
toplevel \$w
wm title \$w $svgfilename
wm iconname \$w $svgfilename

###########################################
# Zinc
##########################################
zinc \$w.zinc -width 600 -height 600 -font 9x15 -borderwidth 0 -backcolor grey90 -render $render

pack \$w.zinc

set topGroup [\$w.zinc add group 1]


HEADER
    );
}


sub fileTail {
    my ($self) = @_;
    $self->printLines(
<<'TAIL'
### translating ojects for making them all visibles

#set bbox [$w.zinc bbox $topGroup]

$w.zinc translate $topGroup 200 150


##### bindings for moving rotating scaling the items

bind $w.zinc <ButtonPress-1>  "press motion %x %y"
bind $w.zinc <ButtonRelease-1>  release
bind $w.zinc <ButtonPress-2>  "press zoom %x %y"
bind $w.zinc <ButtonRelease-2>  release
bind $w.zinc <ButtonPress-3>  "press mouseRotate %x %y"
bind $w.zinc <ButtonRelease-3>  release


set curX 0
set curY 0
set curAngle 0

proc press {action x y} {
    global w curAngle curX curY

    set curX $x
    set curY $y
    set curAngle [expr atan2($y, $x)]
    bind $w.zinc <Motion> "$action %x %y"
}

proc motion {x y} {
    global w topGroup curX curY

    foreach {x1 y1 x2 y2} [$w.zinc transform $topGroup \
			       [list $x $y $curX $curY]] break
    $w.zinc translate $topGroup [expr $x1 - $x2] [expr $y1 - $y2]
    set curX $x
    set curY $y
}

proc zoom {x y} {
    global w curX curY

    if {$x > $curX} {
	set maxX $x
    } else {
	set maxX $curX
    }
    if {$y > $curY} {
	set maxY $y
    } else {
	set maxY $curY
    }
    if {($maxX == 0) || ($maxY == 0)} {
	return;
    }
    set sx [expr 1.0 + (double($x - $curX) / $maxX)]
    set sy [expr 1.0 + (double($y - $curY) / $maxY)]
    $w.zinc scale __svg__1 $sx $sx
    set curX $x
    set curY $y
}

proc mouseRotate {x y} {
    global w curAngle

    set lAngle [expr atan2($y, $x)]
    $w.zinc rotate __svg__1 [expr $lAngle - $curAngle]
    set curAngle  $lAngle
}

proc release {} {
    global w

    bind $w.zinc <Motion> {}
}
TAIL
);

    $self->close;
}


1;


__END__

=head1 NAME

SVG:SVG2zinc::Backend::TclScript - a backend class for generating Tcl script

=head1 SYNOPSIS

 use SVG:SVG2zinc::Backend::TclScript;

 $backend = SVG:SVG2zinc::Backend::TclScript->new(
	       -out => filename_or_handle,
               -in => svgfilename,
	       -verbose => 0|1,
	       -render => 0|1|2,
	       );

 $backend->fileHeader();

 $backend->treatLines("lineOfCode1", "lineOfCode2",...);

 $backend->comment("comment1", "comment2", ...);

 $backend->printLines("comment1", "comment2", ...);

 $backend->fileTail();

=head1 DESCRIPTION

SVG:SVG2zinc::Backend::TclScript is a class for generating Tcl script to display SVG files. The generated script is based on TkZinc.

For more information, you should look at SVG:SVG2zinc::Backend(3pm).

The new method accepts parameters described in the SVG:SVG2zinc::Backend class and the following additionnal parameters:

=over

=item B<-render>

The render value of the TkZinc widget. 0 means no openGL, 1 or 2 for openGL. Defaults to 1.

=back

=head1 SEE ALSO

SVG::SVG2zinc::Backend and SVG::SVG2zinc(3pm)

=head1 BUGS and LIMITATIONS

This is higly experimental. Only few tests... The author is not a Tcl coder!

The Tk::Zinc::SVGExtension perl module provided with SVG::SVG2zinc should be converted in Tcl and imported by (or included in) the generated Tcl script.

=head1 AUTHORS

Christophe Mertz <mertz at intuilab dot com>

=head1 COPYRIGHT
    
CENA (C) 2003-2004

This program is free software; you can redistribute it and/or modify it under the term of the LGPL licence.

=cut

