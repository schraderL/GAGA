#!/usr/bin/perl -w
#
#      this sample script is used to test perl classes 
#      generated by the PerlClass backend from SVG::SVG2zinc
# 
#      Copyright (C) 2003
#      Centre d'�tudes de la Navigation A�rienne
#
#      Authors: Christophe Mertz <mertz at intuilab dot com>
#
# $Id: testGeneratedPerlClasses.pl,v 1.3 2004/05/01 09:19:33 mertz Exp $
#-----------------------------------------------------------------------------------



use strict;
use Tk::Zinc;
use Tk::Zinc::Debug;
use Getopt::Long;
use Carp;

use vars qw( $VERSION);
($VERSION) = sprintf("%d.%02d", q$Revision: 1.3 $ =~ /(\d+)\.(\d+)/);

################ le traitement des options de la ligne de commande
my ($out, $displayResult);
my $render = 1;
GetOptions("help" => \&usage,
	   "render:i" => \$render,
	   "version" => \&displayVersion,
           );

sub usage {
    my ($error) = @_;
    print "$0 [-help] [-version] [-render 0|1|2] class1 class2 ...\n";
    print "  to display svg code transformed in perl Classes\n";
    print "$error\n" if defined $error;
    exit;
}

#&usage unless ($#ARGV < 0);
&usage ("Bad value ($render) for -render option. Must be 0, 1 or 2") unless ($render == 0 or $render == 1 or $render == 2);

my ($WIDTH,$HEIGHT) = (600,600);

my $mw = MainWindow->new();
$mw->title("testing @ARGV");
my $zinc = $mw->Zinc(-width => $WIDTH, -height => $HEIGHT,
		     -borderwidth => 0,
		     -render => $render,
		     -backcolor => "white", ## Pourquoi blanc?
		     )->pack(qw/-expand yes -fill both/);;

if (Tk::Zinc::Debug->can('init')) {
    # for TkZinc >= 3.2.96
    &Tk::Zinc::Debug::init($zinc, -optionsToDisplay => "-tags", -optionsFormat => "row");
} else {
    # for TkZinc <= 3.2.95
    &Tk::Zinc::Debug::finditems($zinc);
    &Tk::Zinc::Debug::tree($zinc, -optionsToDisplay => "-tags", -optionsFormat => "row");
}

my $top_group = $zinc->add('group', 1, -tags => ['topsvggroup']);

foreach my $fileOrClass (@ARGV) {
    my $file;
    my $class;
    if ($fileOrClass =~ /(.*)\.pm$/) {
	$class = $1;
	$file = findINC($fileOrClass);
	$file or die "unable to find $fileOrClass in your perl path: @INC\n";
    } else {
	$class = $fileOrClass;
	$file = findINC("$fileOrClass.pm");
	$file or die "unable to find $fileOrClass.pm in your perl path: @INC\n";
    }
    require $file or carp ("unable to locate $class");

    $zinc->add('group', $top_group, -tags => [$class]);
    my $expr = "$class->new(-zinc => \$zinc, -topgroup => \$class);";
    # print "expr=$expr\n";
    eval($expr);
    carp ("$@\nwhen loading $class") if $@;
}


my $zoom=1;
$zinc->Tk::bind('<ButtonPress-1>', [\&press, \&motion]);
$zinc->Tk::bind('<ButtonRelease-1>', [\&release]);

$zinc->Tk::bind('<ButtonPress-2>', [\&press, \&zoom]);
$zinc->Tk::bind('<ButtonRelease-2>', [\&release]);

$zinc->Tk::bind('<Shift-ButtonPress-2>', [\&press, \&mouseRotate]);
$zinc->Tk::bind('<Shift-ButtonRelease-2>', [\&release]);
$zinc->bind('all', '<Enter>',
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
    $zinc->rotate($top_group, -($langle - $cur_angle));
    $cur_angle = $langle;
}

sub release {
    my ($zinc) = @_;
    $zinc->Tk::bind('<Motion>', '');
}


sub displayVersion {
    print $0, " : Version $VERSION\n";
    exit;
}

sub findINC {
 my $file = join('/',@_);
 my $dir;
 $file  =~ s,::,/,g;
 foreach $dir (@INC)
  {
      my $path; 
      return $path if (-e ($path = "$dir/$file"));
  }
 return undef;
}


__END__

=head1 NAME

testGeneratedPerlClasses.pl - sample script to test generated perl classes (with the PerlClass backend)

=head1 SYNOPSIS

B<testGeneratedPerlClasses.pl> -help

B<testGeneratedPerlClasses.pl> [-render 0|1|2] Class1 Class2 ...


=head1 DESCRIPTION

testGeneratedPerlClasses.pl is a perl script which loads perl one or more classes generated 
with the PerlClass backend of the SVG::SVG2zinc module. It is usefull for testing that 
the generated classes can displays properly in a Tk::Zinc widget. 

=head1 SEE ALSO

SVG::SVG2zinc(3pm) SVG::SVG2zinc::Backend(3pm) Tk::Zinc(3pm) 

TkZinc is available at www.openatc.org/zinc

=head1 AUTHORS

Christophe Mertz <mertz at intuilab dot com>

=head1 COPYRIGHT
    
CENA (C) 2003-2004

This program is free software; you can redistribute it and/or modify it under the term of the LGPL licence.

=cut
