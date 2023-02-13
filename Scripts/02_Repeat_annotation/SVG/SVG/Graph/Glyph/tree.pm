package SVG::Graph::Glyph::tree;

use base SVG::Graph::Glyph;
use Data::Dumper;
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
   my $group = $self->svg->group(id=>"tree$id");

   my($data) = $self->group->data;
   my $root = $data->root;
   my $xsize   = $self->xsize;
   my $ysize   = $self->ysize;
   my $yscale  = $ysize / $root->depth;
   my $xscale  = $xsize / (scalar($root->leaves_under));
   my $xoffset = $self->xoffset;
   my $yoffset = $self->yoffset;


   $group->line(x1=>$xoffset+($xsize/2),x2=>$xoffset+($xsize/2),y1=>$yoffset,y2=>$yoffset+($yscale*$root->branch_length),style=>{$self->_style});

   my $maxy = $self->getmaxy(node=>$root, yscale=>$yscale,  yoffset=>$yoffset+($yscale*$root->branch_length));
   
   #root branch
   $self->_draw(group=>$group,node=>$root,
				xsize=>$xsize,
				yscale=>$yscale,
				xoffset=>$xoffset,yoffset=>$yoffset+($yscale*$root->branch_length),
				maxy=>$maxy
			   );
}

sub getmaxy {
  my($self,%a) = @_;
  my $maxy = $a{yoffset};
  foreach my $d ($a{node}->daughters){
	my $temp = $self->getmaxy(node=>$d, yscale=>$a{yscale}, yoffset => $a{yoffset} + ($d->branch_length * $a{yscale}));
	if ($temp > $maxy) {
	  $maxy = $temp;
	}
  }
  
  return $maxy;

}

=head2 _draw

 Title   : _draw
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub _draw {
  my($self,%a) = @_;

  my $xoffset = $a{xoffset};
  my $xoffset_old;
  my $maxy = $a{maxy};
  my $c = 1;

  foreach my $d ($a{node}->daughters){

	my %style = $a{node}->_style ? $a{node}->_style : $self->_style;

	my $frac = $d->leaves_under / $a{node}->leaves_under;
	$xoffset_old = $xoffset;

	#horizontal bar spanning the domain of a child subtree
	if(!$d->left_sisters and !$d->right_sisters){
	} elsif(!$d->left_sisters){
	  $a{group}->line(x1=>$xoffset+($frac*$a{xsize}/2),
					  y1=>$a{yoffset},
					  x2=>$xoffset+($frac*$a{xsize}),
					  y2=>$a{yoffset},
					  style=>{%style});

#warn join ' ', $a{node}->_style;

	} elsif(!$d->right_sisters){
	  $a{group}->line(x1=>$xoffset,
					  y1=>$a{yoffset},
					  x2=>$xoffset+($frac*$a{xsize}/2),
					  y2=>$a{yoffset},
					  style=>{%style});
	} else {
	  $a{group}->line(x1=>$xoffset,
					  y1=>$a{yoffset},
					  x2=>$xoffset+($frac*$a{xsize}),
					  y2=>$a{yoffset},
					  style=>{%style});
	}



	if($d->depth_under <= 1) {
	  $a{group}->line(x1=>(($frac*$a{xsize}+$xoffset)+($xoffset))/2,x2=>(($frac*$a{xsize}+$xoffset)+($xoffset))/2,
					  y1=>$a{yoffset},
					  y2=>$maxy,
					  style=>{%style}
					 );
	}
	else {
	  #root branch of the child subtree
	  $a{group}->line(x1=>(($frac*$a{xsize}+$xoffset)+($xoffset))/2,x2=>(($frac*$a{xsize}+$xoffset)+($xoffset))/2,
					  y1=>$a{yoffset},
					  y2=>$a{yoffset}+($a{yscale}*$d->branch_length),
					  style=>{%style}
					 );
	}
	
	#root branch label
	my($cx,$cy) = ((($frac*$a{xsize}+$xoffset)+($xoffset))/2 , $maxy);
	#my($cx,$cy) = ((($frac*$a{xsize}+$xoffset)+($xoffset))/2 , $a{yoffset}+($a{yscale}*$d->branch_length));
	$a{group}->text(x=>$cx,
					y=>$cy,
#					style=>{'font-size'=>'15px'},
					style=>{$self->_style,'stroke-width'=>0.5},
					transform=>"rotate(90,$cx,$cy)"
				 )->cdata($d->branch_label());
#->cdata($d->name);


	$xoffset += $frac*$a{xsize};
	$c++;

	$self->_draw(group=>$a{group},node=>$d,
				 xsize=>$xoffset - $xoffset_old,
				 yscale=>$a{yscale},
				 xoffset => $xoffset_old,
				 yoffset => $a{yoffset} + ($d->branch_length * $a{yscale}),
				 maxy => $maxy
				);
  }


}

1;
