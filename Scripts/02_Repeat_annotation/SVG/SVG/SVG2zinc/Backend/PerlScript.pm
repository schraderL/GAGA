package SVG::SVG2zinc::Backend::PerlScript;

#	Backend Class for SVG2zinc
# 
#	Copyright 2003-2004
#	Centre d'Études de la Navigation Aérienne
#
#	Author: Christophe Mertz <mertz at intuilab dot com>
#
#       A concrete class for code generation for Perl Scripts
#
# $Id: PerlScript.pm,v 1.16 2004/05/01 09:19:34 mertz Exp $
#############################################################################

use strict;
use Carp;

use SVG::SVG2zinc::Backend;
use File::Basename;

use vars qw( $VERSION @ISA  );
@ISA = qw( SVG::SVG2zinc::Backend );
($VERSION) = sprintf("%d.%02d", q$Revision: 1.16 $ =~ /(\d+)\.(\d+)/);


sub treatLines {
    my ($self,@lines) = @_;
    foreach my $l (@lines) {
	$l =~ s/->/\$_zinc->/g;
	$self->printLines($l);
    }
}

sub fileHeader {
    my ($self) = @_;
    my $svgfile = $self->{-in};
    my ($svg2zincPackage) = caller;
    my $VERSION = eval ( "\$".$svg2zincPackage."::VERSION" );
    $self->printLines("#!/usr/bin/perl -w

####### This file has been generated from $svgfile by SVG::SVG2zinc.pm Version: $VERSION
");


    $self->printLines(
<<'HEADER'
use Tk::Zinc 3.295;
use Tk::Zinc::Debug;
use Tk::PNG;  # only usefull if loading png file
use Tk::JPEG; # only usefull if loading png file

use Tk::Zinc::SVGExtension;

my $mw = MainWindow->new();
HEADER
		      );
    my $svgfilename = basename($svgfile);
    $self->printLines("
\$mw->title('$svgfile');
my (\$WIDTH, \$HEIGHT) = (800, 600);
" );
    my $render = $self->{-render};
    $self->printLines("
my \$zinc = \$mw->Zinc(-width => \$WIDTH, -height => \$HEIGHT,
		     -borderwidth => 0,
                     -backcolor => 'white', # why white?
		     -render => $render,
		      )->pack(qw/-expand yes -fill both/);;
");

    $self->printLines(
<<'HEADER'
if (Tk::Zinc::Debug->can('init')) {
    # for TkZinc >= 3.2.96
    &Tk::Zinc::Debug::init($zinc, -optionsToDisplay => "-tags", -optionsFormat => "row");
} else {
    # for TkZinc <= 3.2.95
    &Tk::Zinc::Debug::finditems($zinc);
    &Tk::Zinc::Debug::tree($zinc, -optionsToDisplay => "-tags", -optionsFormat => "row");
}

my $top_group = 1; ###$zinc->add('group', 1);

my $_zinc=$zinc;

{ ###

HEADER
    );
}


sub fileTail {
    my ($self) = @_;
    $self->printLines(
<<'TAIL'
	   }

### on va retailler et translater les objets créés!

my @bbox = $_zinc->bbox($top_group);
$_zinc->translate($top_group, -$bbox[0], -$bbox[1]);
@bbox = $_zinc->bbox($top_group);
my $ratio = 1;
$ratio = $WIDTH / $bbox[2] if ($bbox[2] > $WIDTH);
$ratio = $HEIGHT/$bbox[3] if ($HEIGHT/$bbox[3] lt $ratio);
$zinc->scale($top_group, $ratio, $ratio);

### on ajoute quelques binding bien pratiques pour la mise au point

$_zinc->Tk::bind('<ButtonPress-1>', [\&press, \&motion]);
$_zinc->Tk::bind('<ButtonRelease-1>', [\&release]);
$_zinc->Tk::bind('<ButtonPress-2>', [\&press, \&zoom]);
$_zinc->Tk::bind('<ButtonRelease-2>', [\&release]);

# $_zinc->Tk::bind('<ButtonPress-3>', [\&press, \&mouseRotate]);
# $_zinc->Tk::bind('<ButtonRelease-3>', [\&release]);
$_zinc->bind('all', '<Enter>',
	[ sub { my ($z)=@_; my $i=$z->find('withtag', 'current');
			my @tags = $z->gettags($i);
			pop @tags; # pour enlever 'current'
			print "$i (", $z->type($i), ") [@tags]\n";}] );

&Tk::MainLoop;


##### bindings for moving, rotating, scaling the items
my ($cur_x, $cur_y, $cur_angle);
sub press {
    my ($zinc, $action) = @_;
    my $ev = $zinc->XEvent();
    $cur_x = $ev->x;
    $cur_y = $ev->y;
    $cur_angle = atan2($cur_y, $cur_x);
    $zinc->Tk::bind('<Motion>', [$action]);
}

sub motion {
    my ($zinc) = @_;
    my $ev = $zinc->XEvent();
    my $lx = $ev->x;
    my $ly = $ev->y;
    my @res = $zinc->transform($top_group, [$lx, $ly, $cur_x, $cur_y]);
    $zinc->translate($top_group, $res[0] - $res[2], $res[1] - $res[3]);
    $cur_x = $lx;
    $cur_y = $ly;
}

sub zoom {
    my ($zinc, $self) = @_;
    my $ev = $zinc->XEvent();
    my $lx = $ev->x;
    my $ly = $ev->y;
    my $maxx;
    my $maxy;
    my $sx;
    my $sy;
    
    if ($lx > $cur_x) {
	$maxx = $lx;
    } else {
	$maxx = $cur_x;
    }
    if ($ly > $cur_y) {
	$maxy = $ly
    } else {
	$maxy = $cur_y;
    }
    return if ($maxx == 0 || $maxy == 0);
    $sx = 1.0 + ($lx - $cur_x)/$maxx;
    $sy = 1.0 + ($ly - $cur_y)/$maxy;
    $cur_x = $lx;
    $cur_y = $ly;
    $zinc->scale($top_group, $sx, $sx); #$sy);
}

sub mouseRotate {
    my ($zinc) = @_;
    my $ev = $zinc->XEvent();
    my $lx = $ev->x;
    my $ly = $ev->y;
    my $langle = atan2($ly, $lx);
    $zinc->rotate($top_group, -($langle - $cur_angle));
    $cur_angle = $langle;
}

sub release {
    my ($zinc) = @_;
    $zinc->Tk::bind('<Motion>', '');
}
TAIL
);
    $self->close;
}


1;

__END__

=head1 NAME

SVG:SVG2zinc::Backend::PerlScript - a backend class generating Perl script displaying the content of a SVG file

=head1 SYNOPSIS

 use SVG:SVG2zinc::Backend::PerlScript;

 $backend = SVG:SVG2zinc::Backend::PerlScript->new(
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

SVG:SVG2zinc::Backend::PerlScript is a class for generating perl script which displays the content of a SVG file. The generated script requires Tk::Zinc.

For more information, you should look at SVG::SVG2zinc::Backend(3pm).

The generated perl script uses the Tk::Zinc::Debug tool, so it is easy to inspect items created in Tk::Zinc. Use the <ESC> key to get some help when the cursor is in the Tk::Zinc window. 

The B<new> method accepts parameters described in the SVG:SVG2zinc::Backend class and the following additionnal parameter:

=over

=item B<-render>

The render option of the Tk::Zinc widget. A value of 0 means no openGL, 1 or 2 for openGL. Defaults to 1.

=back

=head1 SEE ALSO

SVG::SVG2zinc::Backend(3pm) and SVG::SVG2zinc(3pm)

=head1 AUTHORS

Christophe Mertz <mertz at intuilab dot com>

=head1 COPYRIGHT
    
CENA (C) 2003-2004 IntuiLab (C) 2004

This program is free software; you can redistribute it and/or modify it under the term of the LGPL licence.

=cut

