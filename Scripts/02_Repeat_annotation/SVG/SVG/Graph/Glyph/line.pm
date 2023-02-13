package SVG::Graph::Glyph::line;

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
   my $group = $self->svg->group(id=>"line$id");

   my $xscale = $self->xsize / $self->group->xrange;
   my $yscale = $self->ysize / $self->group->yrange;

   my($x1,$x2,$y1,$y2);
   foreach my $datum (sort {$a->x <=> $b->x} $self->group->data){
	 if(!defined($x1) and !defined($y1)){
	   $x1 = (($datum->x - $self->group->xmin) * $xscale) + $self->xoffset;
	   $y1 = ($self->ysize - ($datum->y - $self->group->ymin) * $yscale) + $self->yoffset;
	   next;
	 }

	 $x2 = (($datum->x - $self->group->xmin) * $xscale) + $self->xoffset;
	 $y2 = ($self->ysize - ($datum->y - $self->group->ymin) * $yscale) + $self->yoffset;

	 $group->line(x1=>$x1,y1=>$y1,x2=>$x2,y2=>$y2,style=>{$self->_style});

	 $x1 = $x2;
	 $y1 = $y2;
   }
}

1;
