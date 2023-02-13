package SVG::Graph::Glyph::heatmap;

use base SVG::Graph::Glyph;
use strict;

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
  # my $frame_transform = $self->frame_transform;
   my $group = $self->svg->group(id=>"heatmap$id");

   my $xscale = $self->xsize / $self->group->xrange;
   my $yscale = $self->ysize / $self->group->yrange;

   my($x1,$x2,$y1,$y2);
   my @data = sort {$a->x <=> $b->x} $self->group->data;

   my %datum = ();
   foreach my $datum ($self->group->data){
	 $datum{$datum->x}{$datum->y} = $datum->z;
   }

   my %x = map {$_->x => 1} $self->group->data;
   my %y = map {$_->y => 1} $self->group->data;
   my @x = sort {$a<=>$b} keys %x;
   my @y = sort {$a<=>$b} keys %y;

#   my $zmean = $self->group->zmean;
   my $zmean = $self->group->zmean;
   my $zmax = $self->group->zmax;
   my $zmin = $self->group->zmin;

   #$zmean = ($zmean - $zmin) / ($zmax - $zmin);
   $zmean = 0.5;

   my($xdim,$ydim);
   for(my $i = 1 ; $i <= $#x + 1; $i++){
	 $x1 = $x[$i-1];
	 $x2 = $x[$i];

	 $x1 = (($x1 - $self->group->xmin) * $xscale) + $self->xoffset;
	 $x2 = (($x2 - $self->group->xmin) * $xscale) + $self->xoffset;
	 $xdim ||= abs($x2-$x1);

	 for(my $j = 1 ; $j <= $#y + 1; $j++){

	   $y1 = $y[$j-1];
	   $y2 = $y[$j];

	   $y1 = ($self->ysize - ($y1 - $self->group->ymin) * $yscale) + $self->yoffset;
	   $y2 = ($self->ysize - ($y2 - $self->group->ymin) * $yscale) + $self->yoffset;
	   $ydim ||= abs($y2-$y1);

	   my $z = $datum{$x[$i-1]}{$y[$j-1]};
	   
	   #$z = $zmean;
	   
	   if(defined $z){

		 #my $zval = ($z - $zmean) / $zmean; #($zmax - $zmean);
		 
		 my $zval = ($z - $zmin)/ ($zmax-$zmin);
		 my $zfp = ($zval - 0.5) / 0.5;
		 my $zfm = (0.5 - $zval) / 0.5;
		 
		 #my $zval = $z;
		 #warn "$zval $zfp $zfm $zmean $zmin $zmax";
		 
		
		 my($rm,$gm,$bm)    = ($self->rgb_m->[0] , $self->rgb_m->[1] , $self->rgb_m->[2]);
		 #my($rm, $gm,$bm)  = (0,0,0);
		 # 255, -255, 0
		 my($rhd,$ghd,$bhd) = ($self->rgb_h->[0] - $self->rgb_m->[0] , $self->rgb_h->[1] - $self->rgb_m->[1] , $self->rgb_h->[2] - $self->rgb_m->[2]);
		 # 0, -255, 255
		 my($rld,$gld,$bld) = ($self->rgb_l->[0] - $self->rgb_m->[0] , $self->rgb_l->[1] - $self->rgb_m->[1] , $self->rgb_l->[2] - $self->rgb_m->[2]);
		 my($r, $g, $b) = (0, 0, 0);

		 #warn "$rld $gld $bld";

		 #if value above mean and red high brighter than red mean
		 if($zval > $zmean && $rhd > 0){
		   $r = int($rm + ($rhd * $zfp));
		 } elsif($zval > $zmean && $rhd < 0){
		   $r = int($rm + ($rhd * $zfp));
		 }
		 if($zval < $zmean && $rld > 0){
		   $r = int($rm + ($rld * $zfm));
		 } elsif($zval < $zmean && $rld < 0) {
		   $r = int($rm + ($rld * $zfm));
		 }

		 #if value above mean and green high brighter than green mean
		 if($zval > $zmean && $ghd > 0){
		   $g = int($gm + ($ghd * $zfp));
		 } elsif($zval > $zmean && $ghd < 0){
		   $g = int($gm + ($ghd * $zfp));
		 } 
		 if($zval < $zmean && $gld > 0){
		   $g = int($gm + ($gld * $zfm));
		 } elsif($zval < $zmean && $gld < 0) {
		   $g = int($gm + ($gld * $zfm));
		 }

		 #if value above mean and blue high brighter than blue mean
		 if($zval > $zmean && $bhd > 0){
		   $b = int($bm + ($bhd * $zfp));
		 } elsif($zval > $zmean && $bhd < 0){
		   $b = int($bm + ($bhd * $zfp));
		 } 
		 if($zval < $zmean && $bld > 0){
		   $b = int($bm + ($bld * $zfm));
		 } elsif($zval < $zmean && $bld < 0) {
		   $b = int($bm + ($bld * $zfm));
		 }
		 
		 #if($zval == $zmean) {
		 #  $r = int($rm);
		 #  $g = int($gm);
		 #  $b = int($bm);
		 #}

		 #warn "red $r, gr $g, blu $b";

#		 $r = int($z > $zmean ? $self->rgb_h->[0] * ($z-$zmean)/($zmax-$zmean) :
#				  $self->rgb_l->[0] * ($zmean-$z)/($zmean-$zmin));
#		 $g = int($z > $zmean ? $self->rgb_h->[1] * ($z-$zmean)/($zmax-$zmean) :
#				  $self->rgb_l->[1] * ($zmean-$z)/($zmean-$zmin));
#		 $b = int($z > $zmean ? $self->rgb_h->[2] * ($z-$zmean)/($zmax-$zmean) :
#				  $self->rgb_l->[2] * ($zmean-$z)/($zmean-$zmin));

#warn "$z [$zmin,$zmean,$zmax] rgb($r,$g,$b)";
#		 $group->rect(x=>$x1,y=>$y1-abs($y2-$y1),width=>abs($x2-$x1),height=>abs($y2-$y1),style=>{'stroke'=>"rgb($r,$g,$b)",'fill'=>"rgb($r,$g,$b)"});
		 
$group->rect(x=>$x1,y=>$y1-$ydim,width=>$xdim,height=>$ydim,style=>{'stroke'=>"rgb($r,$g,$b)",'fill'=>"rgb($r,$g,$b)"});
	   }
	 }
   }
}

=head2 rgb_l

 Title   : rgb_l
 Usage   : $obj->rgb_l($newval)
 Function: 
 Example : 
 Returns : value of rgb_l (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub rgb_l{
    my $self = shift;

    return $self->{'rgb_l'} = shift if @_;
    return $self->{'rgb_l'} || [0,0,0];
}

=head2 rgb_m

 Title   : rgb_m
 Usage   : $obj->rgb_m($newval)
 Function: 
 Example : 
 Returns : value of rgb_m (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub rgb_m{
    my $self = shift;

    return $self->{'rgb_m'} = shift if @_;
    return $self->{'rgb_m'} || [128,128,128];
}

=head2 rgb_h

 Title   : rgb_h
 Usage   : $obj->rgb_h($newval)
 Function: 
 Example : 
 Returns : value of rgb_h (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub rgb_h{
    my $self = shift;

    return $self->{'rgb_h'} = shift if @_;
    return $self->{'rgb_h'} || [255,255,255];
}


1;
