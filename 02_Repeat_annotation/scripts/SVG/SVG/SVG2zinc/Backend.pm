package SVG::SVG2zinc::Backend;

#	Backend for SVG2zinc
#
#	Copyright 2003-2004
#	Centre d'Études de la Navigation Aérienne
#
#	Author: Christophe Mertz <mertz at intuilab dot com>
#
#       An abstract class for code generation
#       Concrete sub-classes can generate code for perl (script / module), tcl,
#       printing, or direct execution
#
# $Id: Backend.pm,v 1.12 2004/05/01 09:19:33 mertz Exp $
#############################################################################

use strict;
use Carp;
use FileHandle;

use vars qw( $VERSION);
($VERSION) = sprintf("%d.%02d", q$Revision: 1.12 $ =~ /(\d+)\.(\d+)/);

sub new {
    my ($class, %passed_options) = @_;
    my $self = {};
    bless $self, $class;
    $self->_initialize(%passed_options);
    return $self;
}

my %new_options = (
		   -render => 1,
		   -out => 1,
		   -in => 1,
		   -verbose => 1,
);

sub _initialize {
    my ($self, %passed_options) = @_;
    foreach my $opt (keys (%passed_options)) {
	if (defined ($new_options{$opt})) {
	    $self->{$opt} = $passed_options{$opt};
	} else {
	    carp ("Warning: option $opt unknown for a ".ref($self)."\n");
	}
    }
    $self->{-render} = 1 unless defined $self->{-render};
    croak("undefined mandatory -in options") unless defined $self->{-in};
    if (defined $self->{-out} and $self->{-out}) {
	my $out = $self->{-out};
	if (ref($out) eq 'GLOB') {
	    ## nothing to do, the $out is supposed to be open!?
	} else {
	    my $fh = FileHandle->new("> " . $out);
	    if ($fh) {
		$self->{-filehandle} = $fh;
	    } else {
		carp ("unable to open " . $out);
	    }
	}
    }
    return $self;
}

# used by SVG2zinc to know the zinc group in which the svg topgroup  
# by default: 1
# currently default need be overriden by PerlClass only, as far as I now!?
sub _topgroup {
    my ($self) = @_;
    if ($self->{-topgroup}) {
	return $self->{-topgroup};
    } else {
	return 1;
    }
}

# returns true if code is put in a file
sub inFile {
    my ($self) = @_;
    return (defined $self->{-filehandle});
}

sub printLines {
    my ($self, @lines) = @_;
    if ($self->inFile) {
	my $fh = $self->{-filehandle};
	foreach my $l (@lines) {
	    print $fh "$l\n";
	}
    } else {
	carp "printLines cannot print if no outfile has been given\n";
    }
}

sub treatLines {
    my ($self, @lines) = @_;
    if ($self->inFile) {
	$self->printLines(@lines);
    }
}


## in case of file generation, should print a comment
## the default is to print comment starting with #
sub comment {
    my ($self, @lines) = @_;
    if ($self->inFile) {
	foreach my $l (@lines) {
	    $self->printLines("##  $l");
	}
    }
}

sub close {
    my ($self) = @_;
    if ($self->inFile) {
	$self->{-filehandle}->close;
    }
}

sub fileHeader {
    my ($self) = @_;
    $self->comment ("", "default Header of SVG::SVG2zinc::Backend", "");
}


sub fileTail {
    my ($self) = @_;
    $self->comment ("", "default Tail of SVG::SVG2zinc::Backend", "");
    $self->close;
}


	
1;


__END__

=head1 NAME

SVG:SVG2zinc::Backend - a virtual class SVG:SVG2zinc svg reader. Sub-class are specialized for different type of generation

=head1 SYNOPSIS

package SVG:SVG2zinc::Backend::SubClass

use SVG:SVG2zinc::Backend;

## some methods definition

....

 ## when using a specialized backend:

 use SVG:SVG2zinc::Backend::SubClass;

 $backend = SVG:SVG2zinc::Backend::SubClass->new(
	       -out => filename_or_handle,
               -in => svgfilename,
	       -verbose => 0|1,
	       [otheroptions],
	       );

 $backend->fileHeader();

 $backend->treatLines("lineOfCode1", "lineOfCode2",...);

 $backend->comment("comment1", "comment2", ...);

 $backend->printLines("comment1", "comment2", ...);

 $backend->fileTail();

=head1 DESCRIPTION

SVG:SVG2zinc::Backend is a perl virtual class which should be specialized in sub-classes. It defines
a common interface ot classes which can for example generate perl code with Tk::Zinc, display
SVG file in a Tk::Zinc widget, convert svg file in image files (e.g. png) or generate tcl code
to be used with TkZinc etc...

A backend should provide the following methods:

=over

=item B<new>

This creation class method should accept pairs of (-option => value) as well as the following arguments:

=over

=item B<-out>

A filename or a filehandle ready for writing the output. In same rare cases
(e.g. the Display backend which only displays the SVG file on the screen,
this option will not be used)

=item B<-in>

The svg filename. It should be used in comments only in the generated file

=item B<-verbose>

It will be used for letting the backend being verbose

=back

=item B<fileHeader>

Generates the header in the out file, if needed. This method should be called just after creating a backend and prior any treatLines or comment method call.

=item B<treatLines>

Processes the given arguments as lines of code. The arguments are very close to Tk::Zinc perl code. When creating a new backend, using the B<Print> backend can help understanding what are exactly these arguments.

=item B<comment>

Processes the given arguments as comments. Depending on the backend, this method must be redefined so that arguments are treated as comments, or just skipped.

=item B<printLines>

Print in an outfile the given arguments as lines of text. This method should not be re-defined, but used by any Backend which generates code.

=item B<fileTail>

Generate the tail in the out file if needed and closes the out file. This must be the last call.

=back

A backend can use the printLines method to print lines in the generated file.

=head1 SEE ALSO

SVG::SVG2zinc::Backend::Display(3pm), SVG::SVG2zinc::Backend::PerlScript(3pm),
SVG::SVG2zinc::Backend::TclScript(3pm), SVG::SVG2zinc::Backend::PerlClass(3pm) code
as examples of SVG::SVG2zinc::Backend subclasses.

SVG::SVG2zinc(3pm)

=head1 AUTHORS

Christophe Mertz <mertz at intuilab dot com> with some help from Daniel Etienne <etienne at cena dot fr>

=head1 COPYRIGHT
    
CENA (C) 2003-2004

This program is free software; you can redistribute it and/or modify it under the term of the LGPL licence.

=cut

