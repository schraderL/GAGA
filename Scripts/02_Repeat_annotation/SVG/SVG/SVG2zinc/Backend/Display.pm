package SVG::SVG2zinc::Backend::Display;

#	Backend Class for SVG2zinc to display a svg file in a Tk::Zinc canvas
# 
#	Copyright 2003
#	Centre d'Études de la Navigation Aérienne
#
#	Author: Christophe Mertz <mertz at intuilab dot com>
#
# $Id: Display.pm,v 1.10 2004/05/01 09:19:33 mertz Exp $
#############################################################################

use SVG::SVG2zinc::Backend;

@ISA = qw( SVG::SVG2zinc::Backend );

use vars qw( $VERSION);
($VERSION) = sprintf("%d.%02d", q$Revision: 1.10 $ =~ /(\d+)\.(\d+)/);

use strict;
use Carp;
use Tk::Zinc::SVGExtension;

eval (require Tk::Zinc);
if ($@) {
    die "$@\nSVG::SVG2zinc::Backend requires Tk::Zinc to be installed\n";
} elsif (eval ('$Tk::Zinc::VERSION !~ /^\d\.\d+$/ or $Tk::Zinc::VERSION < 3.295') ) {
    die "Tk::Zinc must be at least 3.295";
}


sub new {
    my ($class, %passed_options) = @_;
    my $self = {};
    bless $self, $class;
    $self->_initialize(%passed_options);
    return $self;
}

my $zinc;
my ($WIDTH, $HEIGHT);
my $top_group;
sub _initialize {
    my ($self, %passed_options) = @_;
    $WIDTH = delete $passed_options{-width};
    $WIDTH = 600 unless defined $WIDTH;
    $HEIGHT = delete $passed_options{-height};
    $HEIGHT = 600 unless defined $HEIGHT;
    
    $self->SUPER::_initialize(%passed_options);

    require Tk::Zinc::Debug; # usefull for browsing items herarchy
    my $mw = MainWindow->new();
    my $svgfile = $self->{-in};
    $mw->title($svgfile);
    $zinc = $mw->Zinc(-width => $WIDTH, -height => $HEIGHT,
		      -borderwidth => 0,
		      -render => $self->{-render},
		      -backcolor => "white", ## why white?
		      )->pack(qw/-expand yes -fill both/);

    if (Tk::Zinc::Debug->can('init')) {
	# for TkZinc >= 3.2.96
	&Tk::Zinc::Debug::init($zinc, -optionsToDisplay => "-tags", -optionsFormat => "row");
    } else {
	# for TkZinc <= 3.2.95
	&Tk::Zinc::Debug::finditems($zinc);
	&Tk::Zinc::Debug::tree($zinc, -optionsToDisplay => "-tags", -optionsFormat => "row");
    }
}


sub treatLines {
    my ($self,@lines) = @_;
    my $verbose = $self->{-verbose};
    foreach my $l (@lines) {
	my $expr = $l;
	$expr =~ s/->/\$zinc->/g;
	my $unused = $zinc; ## due to a perl bug, this is needed so that $zinc will be known in $expr
	my $r = eval ($expr);
	if ($@) {
#	    &myWarn ("While evaluationg:\n$expr\nAn Error occured: $@\n");
	    print ("While evaluationg:\n$expr\nAn Error occured: $@\n");
	} elsif ($verbose) {
	    if ($l =~ /^->add/) {
		print "$r == $expr\n" if $verbose;
	    } else {
		print "$expr\n" if $verbose;
	    }
	}
    }
}


sub fileHeader {
#    my ($self) = @_;
}


my $zoom;
sub fileTail {
    # resizing to make them all visible
    $top_group = $zinc->find ('withtag', ".1");
    my @bbox = $zinc->bbox($top_group);
    $zinc->translate($top_group, -$bbox[0], -$bbox[1]) if defined $bbox[0] and $bbox[1];
    @bbox = $zinc->bbox($top_group);
    my $ratio = 1;
    $ratio = $WIDTH / $bbox[2] if ($bbox[2] and $bbox[2] > $WIDTH);
    $ratio = $HEIGHT/ $bbox[3] if ($bbox[3] and $HEIGHT/$bbox[3] lt $ratio);

    $zoom=1;
    $zinc->scale($top_group, $ratio, $ratio);

    # adding some usefull callbacks
    $zinc->Tk::bind('<ButtonPress-1>', [\&press, \&motion]);
    $zinc->Tk::bind('<ButtonRelease-1>', [\&release]);
    
    $zinc->Tk::bind('<ButtonPress-2>', [\&press, \&zoom]);
    $zinc->Tk::bind('<ButtonRelease-2>', [\&release]);

    $zinc->Tk::bind('<Control-ButtonPress-1>', [\&press, \&mouseRotate]);
    $zinc->Tk::bind('<Control-ButtonRelease-1>', [\&release]);
    $zinc->bind('all', '<Enter>',
                [ sub { my ($z)=@_; my $i=$z->find('withtag', 'current');
                        my @tags = $z->gettags($i);
                        pop @tags; # to remove the tag 'current'
                        print "$i (", $z->type($i), ") [@tags]\n";}] );

  Tk::MainLoop;
}

##### bindings for moving, rotating, scaling the displayed items
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
    $zinc->translate($top_group, ($res[0] - $res[2])*$zoom, ($res[1] - $res[3])*$zoom);
    $cur_x = $lx;
    $cur_y = $ly;
}

sub zoom {
    my ($zinc, $self) = @_;
    my $ev = $zinc->XEvent();
    my $lx = $ev->x;
    my $ly = $ev->y;
    my ($maxx, $maxy);
    
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
    my $sx = 1.0 + ($lx - $cur_x)/$maxx;
    my $sy = 1.0 + ($ly - $cur_y)/$maxy;
    $cur_x = $lx;
    $cur_y = $ly;
    $zoom = $zoom * $sx;
    $zinc->scale($top_group, $sx, $sx); #$sy);
}

sub mouseRotate {
    my ($zinc) = @_;
    my $ev = $zinc->XEvent();
    my $langle = atan2($ev->y, $ev->x);
    $zinc->rotate($top_group, -($langle - $cur_angle), $cur_x, $cur_y);
    $cur_angle = $langle;
}

sub release {
    my ($zinc) = @_;
    $zinc->Tk::bind('<Motion>', '');
}


sub displayVersion {
    print $0, " : Version $VERSION\n\tSVG::SVG2zinc.pm Version : $SVG::SVG2zinc::VERSION\n";
    exit;
}


1;


__END__

=head1 NAME

SVG:SVG2zinc::Backend::Display - a backend class for displaying SVG file

=head1 DESCRIPTION

SVG:SVG2zinc::Backend::Display is a class for displaying SVG files.

For more information, you should look at SVG:SVG2zinc::Backend(3pm).

The new method accepts parameters described in the SVG:SVG2zinc::Backend class and the following additionnal parameters:

=over

=item B<-render>

The render value of the Tk::Zinc widget. 0 means no openGL, 1 or 2 for openGL. Defaults to 1.

=back

=head1 SEE ALSO

SVG::SVG2zinc::Backend(3pm) and SVG::SVG2zinc(3pm)

=head1 AUTHORS

Christophe Mertz <mertz at intuilab dot com>

=head1 COPYRIGHT
    
CENA (C) 2003-2004 IntuiLab 2004

This program is free software; you can redistribute it and/or modify it under the term of the LGPL licence.

=cut

