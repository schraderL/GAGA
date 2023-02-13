package SVG::Graph::Glyph::axis;

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
   #my $frame_transform = $self->frame_transform;
   my $group = $self->svg->group(id=>"axis$id");

   my $xscale = $self->xsize / $self->group->xrange;
   my $yscale = $self->ysize / $self->group->yrange;

   my $xintercept = $self->x_intercept || 0;
   my $yintercept = $self->y_intercept || 0;

   #draw x axis line
   $group->line(x1=>$self->xoffset,
				y1=>$self->yoffset + $self->ysize - ($yintercept * $yscale),
				x2=>$self->xoffset + $self->xsize,
				y2=>$self->yoffset + $self->ysize - ($yintercept * $yscale),
				style=>{$self->_style}
			   );

   #draw y axis line
   $group->line(x1=>$self->xoffset + ($xintercept * $xscale),
				y1=>$self->yoffset,
				x2=>$self->xoffset + ($xintercept * $xscale),
				y2=>$self->yoffset + $self->ysize,
				style=>{$self->_style}
			   );

   my @x_tick_labels = ();
   my @x_intertick_labels = ();
   if($self->x_tick_labels){
	@x_tick_labels = @{ $self->x_tick_labels };
   }
   if($self->x_intertick_labels){
	@x_intertick_labels = @{ $self->x_intertick_labels };
   }

   my @y_tick_labels = ();
   my @y_intertick_labels = ();
   if($self->y_tick_labels){
	@y_tick_labels = @{ $self->y_tick_labels };
   }
   if($self->y_intertick_labels){
	@y_intertick_labels = @{ $self->y_intertick_labels };
   }

   #x ticks
   if($self->x_absolute_ticks and $self->x_fractional_ticks){
	 die "x axis can't have both absolute and fractional ticks";
   } elsif($self->x_absolute_ticks){

	 for(my $tick = $self->group->xmin ; $tick <= $self->group->xmax; $tick += $self->x_absolute_ticks){
	   my $tickpos = ($tick - $self->group->xmin) * $xscale;

	   $group->line(x1=>$self->xoffset+$tickpos,
					y1=>$self->yoffset+$self->ysize - ($self->y_intercept * $yscale),
					x2=>$self->xoffset+$tickpos,
					y2=>$self->yoffset+$self->ysize - ($self->y_intercept * $yscale)+($self->group->_parent_svg->margin/8),
					style=>{$self->_style}
				   );


	   #tick label
	   my $x_tick_label = shift @x_tick_labels;
	   my $x = $self->xoffset + $tickpos;
	   my $y = $self->yoffset+$self->ysize - ($self->y_intercept * $yscale)+($self->group->_parent_svg->margin/4);
	   $group->text(
					x=>$x,
					y=>$y,
					transform=>"rotate(90,$x,$y)",
				   )->cdata($x_tick_label);

	   #intertick label
	   my $x_intertick_label = shift @x_intertick_labels;
	   $tickpos = (((2 * $tick) - 1) / 2) * $xscale;
	   $x = $self->xoffset + $tickpos;
	   $group->text(
					x=>$x,
					y=>$y,
					transform=>"rotate(90,$x,$y)",
				   )->cdata($x_intertick_label);

	 }
   } elsif($self->x_fractional_ticks){
	 my $inc = $self->group->xrange / $self->x_fractional_ticks;
	 for(my $tick = $self->group->xmin ; $tick <= $self->group->xmax; $tick += $inc){
	   my $tickpos = ($tick - $self->group->xmin) * $xscale;

	   $group->line(x1=>$self->xoffset+$tickpos,
					y1=>$self->yoffset+$self->ysize - ($self->y_intercept * $yscale),
					x2=>$self->xoffset+$tickpos,
					y2=>$self->yoffset+$self->ysize - ($self->y_intercept * $yscale)+($self->group->_parent_svg->margin/8),
					style=>{$self->_style}
				   );

	 }
   }

   #y ticks
   if($self->y_absolute_ticks and $self->y_fractional_ticks){
	 die "y axis can't have both absolute and fractional ticks";
   } elsif($self->y_absolute_ticks){
	 for(my $tick = $self->group->ymin ; $tick <= $self->group->ymax; $tick += $self->y_absolute_ticks){
	   my $tickpos = ($tick - $self->group->ymin) * $yscale;

	   $group->line(x1=>$self->xoffset + ($self->x_intercept * $xscale),
					y1=>$self->yoffset+$self->ysize-$tickpos,
					x2=>$self->xoffset + ($self->x_intercept * $xscale)-($self->group->_parent_svg->margin/8),
					y2=>$self->yoffset+$self->ysize-$tickpos,
					style=>{$self->_style}
				   );


	   #tick label
	   my $y_tick_label = shift @y_tick_labels;
#	   my $x = $self->xoffset + $tickpos;
	   my $x = $self->xoffset + ($self->x_intercept * $xscale) - ($self->group->_parent_svg->margin);
#	   my $y = $self->yoffset+$self->ysize+($self->group->_parent_svg->margin/4);
	   my $y = $self->yoffset + $self->ysize - $tickpos;
	   $group->text(
					x=>$x,
					y=>$y,
#					transform=>"rotate(90,$x,$y)",
				   )->cdata($y_tick_label);

	   #intertick label
	   my $y_intertick_label = shift @y_intertick_labels;
	   $tickpos = (((2 * $tick) - 1) / 2) * $xscale;
	   $x = $self->xoffset + $tickpos;
	   $group->text(
					x=>$x,
					y=>$y,
#					transform=>"rotate(90,$x,$y)",
				   )->cdata($y_intertick_label);

	 }
   } elsif($self->y_fractional_ticks){
	 my $inc = $self->group->yrange / $self->y_fractional_ticks;
	 for(my $tick = $self->group->ymin ; $tick <= $self->group->ymax ; $tick += $inc){
	   my $tickpos = ($tick - $self->group->ymin) * $xscale;

	   $group->line(x1=>$self->xoffset + ($self->x_intercept * $xscale),
					y1=>$self->yoffset+$self->ysize-$tickpos,
					x2=>$self->xoffset + ($self->x_intercept * $xscale)-($self->group->_parent_svg->margin/8),
					y2=>$self->yoffset+$self->ysize-$tickpos,
					style=>{$self->_style}
				   );
	 }
   }

}

=head2 x_intercept

 Title   : x_intercept
 Usage   : $obj->x_intercept($newval)
 Function: 
 Example : 
 Returns : value of x_intercept (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub x_intercept{
    my $self = shift;

    return $self->{'x_intercept'} = shift if @_;
    return $self->{'x_intercept'};
}

=head2 y_intercept

 Title   : y_intercept
 Usage   : $obj->y_intercept($newval)
 Function: 
 Example : 
 Returns : value of y_intercept (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub y_intercept{
    my $self = shift;

    return $self->{'y_intercept'} = shift if @_;
    return $self->{'y_intercept'};
}

=head2 x_fractional_ticks

 Title   : x_fractional_ticks
 Usage   : $obj->x_fractional_ticks($newval)
 Function: 
 Example : 
 Returns : value of x_fractional_ticks (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub x_fractional_ticks{
    my $self = shift;

    return $self->{'x_fractional_ticks'} = shift if @_;
    return $self->{'x_fractional_ticks'};
}

=head2 y_fractional_ticks

 Title   : y_fractional_ticks
 Usage   : $obj->y_fractional_ticks($newval)
 Function: 
 Example : 
 Returns : value of y_fractional_ticks (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub y_fractional_ticks{
    my $self = shift;

    return $self->{'y_fractional_ticks'} = shift if @_;
    return $self->{'y_fractional_ticks'};
}

=head2 x_absolute_ticks

 Title   : x_absolute_ticks
 Usage   : $obj->x_absolute_ticks($newval)
 Function: 
 Example : 
 Returns : value of x_absolute_ticks (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub x_absolute_ticks{
    my $self = shift;

    return $self->{'x_absolute_ticks'} = shift if @_;
    return $self->{'x_absolute_ticks'};
}

=head2 y_fractional_ticks

 Title   : y_fractional_ticks
 Usage   : $obj->y_fractional_ticks($newval)
 Function: 
 Example : 
 Returns : value of y_fractional_ticks (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub y_absolute_ticks{
    my $self = shift;

    return $self->{'y_absolute_ticks'} = shift if @_;
    return $self->{'y_absolute_ticks'};
}

=head2 x_intertick_labels

 Title   : x_intertick_labels
 Usage   : $obj->x_intertick_labels($newval)
 Function: 
 Example : 
 Returns : value of x_intertick_labels (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub x_intertick_labels{
    my $self = shift;

    return $self->{'x_intertick_labels'} = shift if @_;
    return $self->{'x_intertick_labels'};
}

=head2 x_tick_labels

 Title   : x_tick_labels
 Usage   : $obj->x_tick_labels($newval)
 Function: 
 Example : 
 Returns : value of x_tick_labels (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub x_tick_labels{
    my $self = shift;

    return $self->{'x_tick_labels'} = shift if @_;
    return $self->{'x_tick_labels'};
}

=head2 y_intertick_labels

 Title   : y_intertick_labels
 Usage   : $obj->y_intertick_labels($newval)
 Function: 
 Example : 
 Returns : value of y_intertick_labels (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub y_intertick_labels{
    my $self = shift;

    return $self->{'y_intertick_labels'} = shift if @_;
    return $self->{'y_intertick_labels'};
}

=head2 y_tick_labels

 Title   : y_tick_labels
 Usage   : $obj->y_tick_labels($newval)
 Function: 
 Example : 
 Returns : value of y_tick_labels (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub y_tick_labels{
    my $self = shift;

    return $self->{'y_tick_labels'} = shift if @_;
    return $self->{'y_tick_labels'};
}

1;
