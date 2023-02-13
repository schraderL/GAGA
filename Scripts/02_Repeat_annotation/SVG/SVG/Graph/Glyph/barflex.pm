package SVG::Graph::Glyph::barflex;

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
   my $group = $self->svg->group(id=>"bar$id");

   my $xscale = $self->xsize / $self->group->xrange;
   my $yscale = $self->ysize / $self->group->yrange;

   my($x1,$x2,$y1,$y2) = (undef,undef,undef,undef);
   foreach my $datum (sort {$a->x <=> $b->x} $self->group->data){
	 if(!defined($x1) and !defined($y1)){
	   $x1 = (($datum->x - $self->group->xmin) * $xscale) + $self->xoffset;
	   $y1 = ($self->ysize - ($datum->y - $self->group->ymin) * $yscale) + $self->yoffset;
	   next;
	 }

	 $x2 = (($datum->x - $self->group->xmin) * $xscale) + $self->xoffset;
	 $y2 = ($self->ysize - ($datum->y - $self->group->ymin) * $yscale) + $self->yoffset;

#	 $group->line(x1=>$x1,y1=>$y1,x2=>$x2,y2=>$y2,style=>{$self->_style});
	 $group->rectangle(x=>$x1,
					   y=>$y1,
					   width=>$x2-$x1,
					   height=>$self->ysize - $y1 + $self->yoffset,
					   style=>{$self->_style});

	 $x1 = $x2;
	 $y1 = $y2;
   }
}

1;
