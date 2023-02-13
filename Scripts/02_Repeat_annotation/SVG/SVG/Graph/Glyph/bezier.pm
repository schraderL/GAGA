package SVG::Graph::Glyph::bezier;

use base SVG::Graph::Glyph;
use strict;
use Math::Spline;

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

   my $id = 'n'.sprintf("%07d",int(rand(9999999)));
   my $group = $self->svg->group(id=>"bezier$id");

   my $xscale = $self->xsize / $self->group->xrange;
   my $yscale = $self->ysize / $self->group->yrange;


   my @p = map {[$_->x,$_->y]} sort {$a->x <=> $b->x} $self->group->data;
   my $c = spline_generate(@p);

   my @c;
   foreach my $i (1..$#p){
	 my($x1,$y1) = @{$p[$i-1]};
	 my($x2,$y2) = @{$p[$i]  };

	 my $e = 0;
	 my($cx,$cy,$cx_sum,$cy_sum);
	 for(my $d = 2; $d <= 100 ; $d += 2){
	   next if $d == 0;
	   my $delta = ($x2 - $x1) / $d;
	   my $y1d = spline_evaluate($x1+$delta,$c,@p);
	   my $y2d = spline_evaluate($x2-$delta,$c,@p);

	   my($int,$t1,$t2) = Intersection(Segment([$x1,$y1],[$x1+$delta,$y1d]),Segment([$x2-$delta,$y2d],[$x2,$y2]));
	   next unless ref($int);
	   next unless $x1 <= $int->[0] and $x2 >= $int->[0];

	   $e++;
	   $cx_sum += $int->[0];
	   $cy_sum += $int->[1];
	 }

	 $x1 = (($x1 - $self->group->xmin) * $xscale) + $self->xoffset;
	 $y1 = ($self->ysize - ($y1 - $self->group->ymin) * $yscale) + $self->yoffset;
	 $x2 = (($x2 - $self->group->xmin) * $xscale) + $self->xoffset;
	 $y2 = ($self->ysize - ($y2 - $self->group->ymin) * $yscale) + $self->yoffset;

#warn $cx_sum,"\t",$e;
#next unless $e;
	 if($e){
	   $cx = $cx_sum / $e;
	   $cy = $cy_sum / $e;

	   $cx = (($cx - $self->group->xmin) * $xscale) + $self->xoffset;
	   $cy = ($self->ysize - ($cy - $self->group->ymin) * $yscale) + $self->yoffset;

	   $group->path(d=>"M$x1,$y1 Q$cx,$cy $x2,$y2",
					style=>{$self->_style},
				   );
	 } else {
	   $group->line(x1=>$x1,y1=>$y1,x2=>$x2,y2=>$y2,style=>{$self->_style});
	 }
   }


}

=head2 spline_generate

 Title   : spline_generate
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub spline_generate{
   my (@points) = @_;

   my ($i, $delta, $temp, @factors, @coeffs);
   $coeffs[0] = $factors[0] = 0;

   # Decomposition phase of the tridiagonal system of equations
   for ($i = 1; $i < @points - 1; $i++) {

#ad
	 next unless ($points[$i+1][0] - $points[$i-1][0]);

	 $delta = ($points[$i][0] - $points[$i-1][0]) / (($points[$i+1][0] - $points[$i-1][0]));


	 $temp = $delta * ($coeffs[$i-1] || 0) + 2;


	 $coeffs[$i] = ($delta - 1) / @points;

#ad
	 next unless ($points[$i+1][0] - $points[$i][0]);
	 next unless ($points[$i][0] - $points[$i-1][0]);

	 $factors[$i] = ($points[$i+1][1] - $points[$i][1]) / (($points[$i+1][0] - $points[$i][0]))
	                -
	                ($points[$i][1] - $points[$i-1][1]) / (($points[$i][0] - $points[$i-1][0]));


	 $factors[$i] = ( 6 * $factors[$i] / (($points[$i+1][0] - $points[$i-1][0]))
					  -
					  $delta * $factors[$i-1] ) / $temp;
   }

   # Backsubstitution phase of the tridiagonal system
   #
   $coeffs[$#points] = 0;
   for ($i = @points - 2; $i >= 0; $i--) {
	 $coeffs[$i] = ($coeffs[$i] || 0) * ($coeffs[$i+1] || 0) + ($factors[$i] || 0);
   }
   return \@coeffs;
}

=head2 spline_evaluate

 Title   : spline_evaluate
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub spline_evaluate{
  my ($x, $coeffs, @points) = @_;
  my ($i, $delta, $mult);

  # Which section of the spline are we in?
  #
  for ($i = @points - 2; $i >= 1; $i--) {
	last if $x >= $points[$i][0];
  }

  $delta = $points[$i+1][0] - $points[$i][0];

#ad
  return 0 unless $delta;

  $mult = ( $coeffs->[$i]/2 ) + ($x - $points[$i][0]) * ($coeffs->[$i+1] - $coeffs->[$i]) / (6 * $delta);
  $mult *= $x - $points[$i][0];
  $mult += ($points[$i+1][1] - $points[$i][1]) / ($delta);
  $mult -= ($coeffs->[$i+1] + 2 * $coeffs->[$i]) * $delta / 6;
  return $points[$i][1] + $mult * ($x - $points[$i][0]);
}

=head2 Segment

 Title   : Segment
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub Segment {
  [ $_[0], $_[1] ];
}

=head2 VectorSum

 Title   : VectorSum
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub VectorSum {
  my ($x, $y) = (0, 0);
  for (@_) {
    $x += $_->[0];
    $y += $_->[1];
  }
  [ $x, $y ];
}

=head2 ScalarProd

 Title   : ScalarProd
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub ScalarProd {
  my $s = shift;
  my $v = shift;
  [ $v->[0] * $s, $v->[1] * $s ];
}

=head2 Minus

 Title   : Minus
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub Minus {
  my $v = shift;
  [ - $v->[0], - $v->[1] ];
}

=head2 Intersection

 Title   : Intersection
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub Intersection {
  my $l1 = shift;
  my $l2 = shift;
  my ($bp1, $op1) = @$l1;
  my ($bp2, $op2) = @$l2;
  my $v1 = VectorSum($op1, Minus($bp1));
  my $v2 = VectorSum($op2, Minus($bp2));
  my ($x1, $y1) = @$v1;
  my ($x2, $y2) = @$v2;

  my $DEN = $x1 * $y2 - $x2 * $y1;

  # Lines are parallel.
  return undef if $DEN == 0;

  my ($bx1, $by1) = @$bp1;
  my ($bx2, $by2) = @$bp2;
  
  my $t1 = (($bx2 - $bx1) * $y2 - ($by2 - $by1) * $x2) / $DEN;
  my $t2 = (($bx2 - $bx1) * $y1 - ($by2 - $by1) * $x1) / $DEN;

  my $RESULT = VectorSum($bp1, ScalarProd($t1, $v1));
  ($RESULT, $t1, $t2);
}

1;
