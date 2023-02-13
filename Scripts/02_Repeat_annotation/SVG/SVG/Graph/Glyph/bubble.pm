package SVG::Graph::Glyph::bubble;

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
   my $group = $self->svg->group(id=>"bubble$id");

   my $xscale = $self->xsize / $self->group->xrange;
   my $yscale = $self->ysize / $self->group->yrange;
   my $zscale = 30           / $self->group->zrange;

   foreach my $datum($self->group->data){
	 my $cx = (($datum->x - $self->group->xmin) * $xscale) + $self->xoffset;
	 my $cy = (($self->xsize - ($datum->y - $self->group->ymin) * $yscale)) + $self->yoffset;

	 my $r  = $zscale * ($datum->z + 0.001);

	 $group->circle(cx=>$cx,cy=>$cy,r=>$r,style=>{$self->_style});
   }
}

1;
