package SVG::SVG2zinc::Backend::Image;

#	Backend Class for SVG2zinc to generate image files
# 
#	Copyright 2003
#	Centre d'Études de la Navigation Aérienne
#
#	Author: Christophe Mertz <mertz at intuilab dot com>
#
# $Id: Image.pm,v 1.5 2004/05/01 09:19:33 mertz Exp $
#############################################################################

use SVG::SVG2zinc::Backend;

@ISA = qw( SVG::SVG2zinc::Backend );

use vars qw( $VERSION);
($VERSION) = sprintf("%d.%02d", q$Revision: 1.5 $ =~ /(\d+)\.(\d+)/);

use strict;
use Carp;
use Tk::Zinc::SVGExtension;

eval (require Tk::Zinc);
if ($@) {
    print "$@\nSVG::SVG2zinc::Backend requires Tk::Zinc to be installed\n";
}

sub new {
    # testing that 'import' is available
    # is this test portable?
    my $ret = `which import`;
    croak ("## You need the 'import' command from 'imagemagic' package to use the Image backend.\n")
	if !$ret;
    
    my ($class, %passed_options) = @_;
    my $self = {};
    bless $self, $class;
    $self->_initialize(%passed_options);
    return $self;
}

my $zinc;
sub _initialize {
    my ($self, %passed_options) = @_;
    if (defined $passed_options{-ratio}) {
	if ($passed_options{-ratio} !~ /^\d+%$/) {
	    croak ("## -ratio should look like nn%");
	} else {
	    $self->{-ratio} = delete $passed_options{-ratio};
	}
    }
    if (defined $passed_options{-width}) {
	$self->{-width} = delete $passed_options{-width};
    }
    if (defined $passed_options{-height}) {
	$self->{-height} = delete $passed_options{-height};
    }

    $self->SUPER::_initialize(%passed_options);

    my $mw = MainWindow->new();
    my $svgfile = $self->{-in};
    $mw->title($svgfile);
    my $render = (defined $self->{-render}) ? $self->{-render} : 1;
    $zinc = $mw->Zinc(-borderwidth => 0,
		      -render => $render,
		      -backcolor => "white", ## why white?
		      )->pack(qw/-expand yes -fill both/);
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


sub fileTail {
    my ($self) = @_;
    my $outfile = $self->{-out};

    # to find the top group containing width and height
    my $svgGroup = $zinc->find('withtag', 'svg_top') ; 

    my $tags = join " ", $zinc->gettags($svgGroup);
#    print "svgGroup=$svgGroup  => $tags\n";
    my ($width) = $tags =~ /width=(\d+)/ ;
    my ($height) = $tags =~ /height=(\d+)/ ;
#    print "height => $height width => $width\n";
    
    $zinc->configure (-width => $width, -height => $height);
    $zinc->update;

    my $requiredWidth = $self->{-width};
    my $requiredHeigth = $self->{-height};
    my $importParams="";
    if (defined $requiredWidth and defined $requiredHeigth) {
	$importParams=" -resize $requiredWidth"."x$requiredHeigth";
    } elsif (defined $self->{-ratio}) {
	$importParams=" -resize ".$self->{-ratio};
    }
#    print "importParams=$importParams\n";

    ## following are for comments:
    my ($svg2zincPackage) = caller;
    my $VERSION = eval ( "\$".$svg2zincPackage."::VERSION" );
    my $svgfile = $self->{-in};

    my $command = "import -window " . $zinc->id . $importParams ." -comment 'created with SVG::SVG2zinc from $svgfile v$VERSION (c) CENA 2003 C.Mertz.' $outfile";
#    print "command=$command\n";
    my $return = system ($command);

    if ($return) {
	## -1 when import is not available
	print "## To use the Image Backend you need the 'import' command\n";
	print "## from the 'imagemagick' package on your system\n";
    }
}


1;


__END__

=head1 NAME

SVG:SVG2zinc::Backend::Image - a backend class for generating image file from SVG file

=head1 DESCRIPTION

SVG:SVG2zinc::Backend::Image is a backend class for generating image file from SVG files. It uses the 'import' command included in the ImageMagick package.

For more information, you should look at SVG:SVG2zinc::Backend(3pm).

The new method accepts parameters described in the SVG:SVG2zinc::Backend class and the following additionnal parameters:

=over

=item B<none>

=back

=head1 SEE ALSO

SVG::SVG2zinc::Backend(3pm) and SVG::SVG2zinc(3pm)

=head1 BUGS and LIMITATIONS

This backend generates images files from the content of a displayed Tk::Zinc window. The size (in pixels) of the generated image is thus limited to the maximal size of a window on your system. 

=head1 AUTHORS

Christophe Mertz <mertz at intuilab dot com>

=head1 COPYRIGHT
    
CENA (C) 2003-2004 IntuiLab 2004

This program is free software; you can redistribute it and/or modify it under the term of the LGPL licence.

=cut

