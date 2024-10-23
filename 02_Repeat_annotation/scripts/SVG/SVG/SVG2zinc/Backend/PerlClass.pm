package SVG::SVG2zinc::Backend::PerlClass;

#	Backend Class for SVG2zinc
# 
#	Copyright 2003-2004
#	Centre d'Études de la Navigation Aérienne
#
#	Author: Christophe Mertz <mertz at intuilab dot com>
#
#       An concrete class for code generation for Perl Class
#
# $Id: PerlClass.pm,v 1.5 2004/05/01 09:19:33 mertz Exp $
#############################################################################

use SVG::SVG2zinc::Backend;
use File::Basename;

@ISA = qw( SVG::SVG2zinc::Backend );

use vars qw( $VERSION);
($VERSION) = sprintf("%d.%02d", q$Revision: 1.5 $ =~ /(\d+)\.(\d+)/);

use strict;
use Carp;

sub _initialize {
    my ($self, %passed_options) = @_;
    $self->{-topgroup} = '$self->{-topgroup}'; # this code is used by SVG2zinc
    $self->SUPER::_initialize(%passed_options);
    return $self;
}

sub treatLines {
    my ($self,@lines) = @_;
    foreach my $l (@lines) {
	$l =~ s/^(\s*)->/$1\$_zinc->/g;
	$l =~ s/(\W)->/$1\$_zinc->/g;
	$self->printLines($l);
    }
}

sub fileHeader {
    my ($self) = @_;
    my $file = $self->{-in}; # print "file=$file\n";
    my ($svg2zincPackage) = caller;
    my $VERSION = eval ( "\$".$svg2zincPackage."::VERSION" );
    my ($package_name) = basename ($self->{-out}) =~ /(.*)\.pm$/ ;
    
    $self->printLines("package $package_name;

####### This file has been generated from $file by SVG::SVG2zinc.pm Version: $VERSION

");
    $self->printLines(
<<'HEADER'
use Tk;
use Tk::Zinc 3.295;
use Tk::PNG;  # only usefull if loading png file
use Tk::JPEG; # only usefull if loading png file
use Tk::Zinc::SVGExtension;
use strict;

use Carp;
		      

sub new {
    my ($class, %passed_options) = @_;
    my $self = {};
    bless $self, $class;

    my $_zinc = $passed_options{-zinc};
    croak ("-zinc option is mandatory at instanciation") unless defined $_zinc;

    if (defined $passed_options{-topgroup}) {
	$self->{-topgroup} = $passed_options{-topgroup};
    } else {
	$self->{-topgroup} = 1;
    }
    

# on now items creation!
HEADER
);
}


sub fileTail {
    my ($self) = @_;
    $self->comment ("", "Tail of SVG2zinc::Backend::PerlScript", "");
    $self->printLines(
<<'TAIL'
return $self;
}

1;
TAIL
);
    $self->close;
}


1;

