package SVG::Graph::Glyph::pictogram;

use base SVG::Graph::Glyph;
use strict;
use Data::Dumper;

our @lifts;

=head2 draw

 Title   : draw
 Usage   :
 Function:
 Example :
 Returns :
 Args    :


=cut

sub draw{
  my ($self,@args) = @_;

  my($x1,$x2,$y1,$y2) = (undef,undef,undef,undef);
  #this will be the width of a bar, we need to calculate the
  #minimum distance between any two values in the data,
  #and use that as the step size.
  my $step = undef;
  my %xcoords = map {$_->x => 1} $self->group->data;

  ($x1,$x2) = (undef,undef);
  foreach my $xcoord (sort {$a <=> $b} keys %xcoords){

	$x1 = $xcoord and next unless defined($x1);
	$x2 = $xcoord;
	
	if(!defined($step)){
	  $step = $x2 - $x1;
	} else {
	  $step = $x2 - $x1 < $step ? $x2 - $x1 : $step;
	}

	$x1 = $x2;
  }

  foreach my $data ($self->group->data_chunks){
	my @newlifts = $self->_glyphs($step,\@lifts,[$data->data]);
	my $i = 0;
	foreach my $newlift (@newlifts){
	  $lifts[$i] += $newlift;
	  $i++;
	}
  }
}

sub _glyphs {
  my $self = shift;
  my($step,$lifts,$data) = @_;

  my @ls = ();
  my $l;

  #FIXME this logic is fubar, but it keeps things in proper scale.

  my $yscale = $self->ysize / ($self->group->ymax);
  my $xscale = $self->xsize / (scalar(@$data));

  my $fontsize = int($self->xsize / (scalar(@$data) * 0.96));
  foreach my $datum (sort {$a->x <=> $b->x} @$data){


	my $xsc = 1;
	my $ysc = $datum->y * $yscale;
warn $ysc;
	my $x1 = (($datum->x - 1) * $fontsize) + $self->xoffset;
	my $y1 = $self->ysize + $self->yoffset - ($self->yoffset / 8);

	if($datum->y > 0){
	  $self->canvas->text(
						  x             => 0,
						  y             => 0,
#						  'transform'   => "scale($xsc,$ysc), translate(". $x1/$xsc .",". $y1/$ysc .")",
						  'transform'   => "scale($xsc,$ysc), translate(". $x1/$xsc .",". 0 .")",
						  'font-size'   => $fontsize,
						  'font-family' => 'Helvetica',
						  'text-anchor' => 'right',
						  'fill'        => $self->color->{$datum->label},
						 )->cdata($datum->label);
	}
  }
}

=head2 color

 Title   : color
 Usage   : $obj->color($newval)
 Function: 
 Example : 
 Returns : value of color (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub color{
    my $self = shift;

    return $self->{'color'} = shift if @_;
    return $self->{'color'};
}

1;
