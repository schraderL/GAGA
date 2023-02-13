package SVG::Graph::Glyph::wedge;

use base SVG::Graph::Glyph;
use strict;
use constant PI => 3.14159;

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

   my %fill = (1=>'red',
			   2=>'orange',
			   3=>'yellow',
			   4=>'green',
			   5=>'blue',
			   6=>'indigo',
			   7=>'violet',
			  );

   my $total = 0;
   my $wedge_count = 0;
   foreach my $datum ($self->group->data){
	 $wedge_count++;
	 $total += $datum->x;
	 die __PACKAGE__." can't take negative values" if $datum->x < 0;
   }

   my $cx = ($self->xsize / 2) + $self->xoffset;
   my $cy = ($self->ysize / 2) + $self->yoffset;
   my $r  = ($self->xsize) / 2;
   my $s  = $r;

   my $offset = 0;
   my $old_offset;
   my $wedge = 0;
   foreach my $datum ($self->group->data){
	 $wedge++;
	 $old_offset = $offset;
	 $offset += ($datum->x / $total);

	 my $v = $old_offset * 2 * PI;
	 my $w = $offset     * 2 * PI;

	 my $x1 = $cx + cos($v)*$r;
	 my $y1 = $cy + sin($v)*$s;
	 my $x2 = $cx + cos($w)*$r;
	 my $y2 = $cy + sin($w)*$s;

	 my $large = $datum->x < $total / 2 ? 0 : 1;

	 my %extra = ();

	 my $id = 'n'.sprintf("%07d",int(rand(9999999)));
	 my $group = $self->svg->group(id=>"wedge$id",%extra);

#	 $group->line(x1=>$cx,y1=>$cy , x2=>$x2 , y2=> $y2,style=>{'stroke-width'=>1,'stroke'=>'black'});
#	 $group->ellipse(cx=>$cx,cy=>$cy,rx=>$r,ry=>$s,style=>{'fill-opacity'=>0,'stroke'=>'black'});

	 $group->path(d=>"M$cx,$cy L$x1,$y1 A$r,$s 0 $large 1 $x2,$y2 L$cx,$cy",
					  style=>{'fill-opacity'=>0.4,'fill'=>$fill{$wedge},stroke=>$fill{$wedge},'stroke-width'=>1},
					 );
   }
}

1;
