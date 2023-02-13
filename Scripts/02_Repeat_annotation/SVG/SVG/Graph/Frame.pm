package SVG::Graph::Frame;

use base SVG::Graph::Data;
use strict;
use Data::Dumper;

=head2 new

 Title   : new
 Usage   : you should not be calling new, see add_frame
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub new {
  my($class, %args) = @_;
  my $self = bless {}, $class;
  $self->init(%args);
  return $self;
}

=head2 init

 Title   : init
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub init {
  my($self, %args) = @_;

#  die "you must provide a 'data' arg to new()" unless $args{data};
  $self->_parent_svg($args{svg});

  foreach my $arg (keys %args){
	my $meth = 'add_'.$arg;
        $self->$meth($args{$arg});
  }

  my $id = 'n'.sprintf("%07d",int(rand(9999999)));
  my $group;

  if($self->frame_transform eq "top") {
	$group=$self->_parent_svg->svg->group(id=> $id);
  }
  elsif($self->frame_transform eq "left") {

	my $scaley = $self->xsize/$self->ysize;
	my $scalex = $self->ysize/$self->xsize;
   
	my $translateoffx = -$self->xoffset;
	my $translateoffy = -$self->yoffset;

	my $translatey = $self->ysize + $self->yoffset;
	my $translatex = $self->xoffset;
	
	$group = $self->_parent_svg->svg->group(id => $id, transform =>"translate($translatex, $translatey) scale($scaley, $scalex)  rotate(-90) translate($translateoffx, $translateoffy)");
  }
  elsif($self->frame_transform eq "right") {
	
	my $scalex = $self->xsize/$self->ysize;
	my $scaley = $self->ysize/$self->xsize;

	my $translateoffx = -$self->xoffset;
	my $translateoffy = -$self->yoffset;

	my $translatey = $self->yoffset;
	my $translatex = $self->xsize + $self->xoffset;

	$group=$self->_parent_svg->svg->group(id => $id, transform => "translate($translatex, $translatey) scale($scalex, $scaley) rotate(90) translate($translateoffx, $translateoffy)");

  }
  elsif($self->frame_transform eq "bottom") {

	my $translateoffx = -$self->xoffset;
	my $translateoffy = -$self->yoffset;

	my $translatex = $self->xsize + $self->xoffset;
 	my $translatey = $self->ysize + $self->yoffset;
	$group=$self->_parent_svg->svg->group(id => $id, transform => "translate($translatex, $translatey) rotate(180) translate($translateoffx, $translateoffy)");

  }
  else {
	$group=$self->_parent_svg->svg->group(id=> $id);
  }
  
  $self->svg($group);
  $self->is_changed(1);
}

=head2 add_glyph

 Title   : add_glyph
 Usage   : $frame->add_glyph( 'glyph_name', glyph_args=>arg)
 Function: adds a glyph to the Frame object
 Returns : a SVG::Graph::Glyph::glyph_name object
 Args    : glyph dependent


=cut

sub add_glyph {
  my($self, $glyphtype, %args) = @_;

  my $class = 'SVG::Graph::Glyph::'.$glyphtype || 'generic';
  eval "require $class"; if($@){ die "couldn't load $class: $@" };

  my $glyph = $class->new(%args, svg => $self->svg, group => $self,
						  xsize=>$self->xsize,
						  ysize=>$self->ysize,
						  xoffset=>$self->xoffset,
						  yoffset=>$self->yoffset,
						 );
  push @{$self->{glyphs}}, $glyph;
  return $glyph;
}

=head2 add_frame

 Title   : add_frame
 Usage   : my $frame = $graph->add_frame
 Function: adds a Frame to the current Frame
 Returns : a SVG::Graph::Frame object
 Args    : (optional) frame_transform => 'top' default orientation
                                         'bottom' rotates graph 180 deg (about the center)
                                         'right' points top position towards right
                                         'left' points top position towards left
 


=cut

sub add_frame {
  my($self,$frames) = @_;

  my $epitaph = "only SVG::Graph::Frame objects accepted";

  #die $epitaph unless ref $frames->{frame} eq 'SVG::Graph::Frame';
  
  my $frame_arg = $frames->{frame};

  if (ref($frame_arg) eq 'ARRAY') {
	foreach my $frame (@$frame_arg) {
	  die $epitaph unless ref $frame eq 'SVG::Graph::Frame';
	  push @{$self->{frames}}, $frame; 
	}
	  
  }
  elsif(ref($frame_arg) eq __PACKAGE__) {
	push @{$self->{frames}}, $frame_arg;
  }
  else {
	my $frame = SVG::Graph::Frame->new(svg=>$self->_parent_svg,
									   _parent_frame=>$self,
									   xoffset=>$self->_parent_svg->margin,
									   yoffset=>$self->_parent_svg->margin,
									   xsize=>$self->_parent_svg->width  - (2 * $self->_parent_svg->margin),
									   ysize=>$self->_parent_svg->height - (2 * $self->_parent_svg->margin),
									   frame_transform=>$frames->{frame_transform}
									  );	
	
	$frame->stack($self->stack) if $self->stack;
	$frame->ystat($self->ystat) if $self->stack;
	
	push @{$self->{frames}}, $frame;
	return $frame;
  }

  $self->is_changed(1);
}

=head2 frames

 Title   : frames
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub frames {
  my $self = shift;
  return $self->{frames} ? @{$self->{frames}} : ();
}

=head2 add_data

 Title   : add_data
 Usage   : $frame->add_data($data)
 Function: adds a SVG::Graph::Data object to the current Frame
 Returns : none
 Args    : SVG::Graph::Data object


=cut

sub add_data {
  my($self,@datas) = @_;

  my $epitaph = "only SVG::Graph::Data objects accepted";

  foreach my $data (@datas){
	if(ref $data eq 'ARRAY'){
	  foreach my $d (@$data){
		die $epitaph unless ref $d eq 'SVG::Graph::Data' || ref $data eq 'SVG::Graph::Data::Tree';
		push @{$self->{data}}, $d;
	  }
	} else {
	  die $epitaph unless ref $data eq 'SVG::Graph::Data' || ref $data eq 'SVG::Graph::Data::Tree';
	  push @{$self->{data}}, $data;
	}
  }

  $self->is_changed(1);
}

=head2 all_data

 Title   : all_data
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub all_data {
  my $self = shift;
  my $flag = shift;

  if(($self->_parent_frame && $flag) || !$self->_parent_frame){
	my @data = $self->data;

	#recurse down into subframes...
	foreach my $subframe ($self->frames){
	  push @data, $subframe->all_data(1);
	}

	return map {$_->can('data') ? $_->all_data(1) : $_ } @data;
  } elsif($self->_parent_frame) {
	return $self->_parent_frame->all_data;
  }
}

=head2 data

 Title   : data
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub data {
  my $self = shift;

  #these are SVG::Graph::Data objects
  my @data = $self->{data} ? @{$self->{data}} : ();

  #recurse down into subframes...
  foreach my $subframe ($self->frames){
	push @data, $subframe->data;
  }

  return map {$_->can('data') ? $_->data : $_ } @data;
}

=head2 glyphs

 Title   : glyphs
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub glyphs {
  my $self = shift;
  return $self->{glyphs} ? @{$self->{glyphs}} : ();
}

=head2 data_chunks

 Title   : data_chunks
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub data_chunks {
  my $self = shift;

  my @data = $self->{data} ? @{$self->{data}} : ();

  #recurse down into subframes...
  foreach my $subframe ($self->frames){
	push @data, $subframe->data_chunks;
  }

  return @data;
}

=head2 draw

 Title   : draw
 Usage   : should not directly call this method, see SVG::Graph->draw
 Function: depends on child glyph implementations
 Example :
 Returns : 
 Args    :


=cut

sub draw {
  my($self, $svg) = @_;

  foreach my $frame ($self->frames){
#warn $frame;
	$frame->draw($self);
  }

  foreach my $glyph ($self->glyphs){
#warn $glyph;
	$glyph->draw($self);
  }
}

=head2 _recalculate_stats

 Title   : _recalculate_stats
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub _recalculate_stats{
   my ($self,@args) = @_;
   return undef unless $self->is_changed;

   my $xstat = Statistics::Descriptive::Full->new();
   my $ystat = Statistics::Descriptive::Full->new();
   my $zstat = Statistics::Descriptive::Full->new();

   #right now we only support y-stacking.  this may need to be extended in the future
   if($self->stack){
	 my @ystack;
	 foreach my $data ($self->data_chunks){
	   my $i = 0;
	   foreach my $datum ($data->data){
		 $ystack[$i] += $datum->y;
		 $i++;
	   }
	 }

	 $ystat->add_data($_) foreach @ystack;
   } else {
	 $ystat->add_data(map {ref($_) && $_->can('y') ? $_->y : $_} map {$_->can('data') ? $_->data : $_->y} $self->all_data);

   }

   $xstat->add_data(map {ref($_) && $_->can('x') ? $_->x : $_} map {$_->can('data') ? $_->data : $_->x} $self->all_data);
   $zstat->add_data(map {ref($_) && $_->can('z') ? $_->z : $_} map {$_->can('data') ? $_->data : $_->z} $self->all_data);

   $self->xstat($xstat);
   $self->ystat($ystat);
   $self->zstat($zstat);

   $self->is_changed(0);
}

=head2 _parent_svg

 Title   : _parent_svg
 Usage   : $obj->_parent_svg($newval)
 Function: 
 Example : 
 Returns : value of _parent_svg (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub _parent_svg{
    my $self = shift;

    return $self->{'_parent_svg'} = shift if @_;
    return $self->{'_parent_svg'};
}

=head2 _parent_frame

 Title   : _parent_frame
 Usage   : $obj->_parent_frame($newval)
 Function: 
 Example : 
 Returns : value of _parent_frame (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub add__parent_frame{return shift->_parent_frame(@_)}
sub _parent_frame{
    my $self = shift;

    return $self->{'_parent_frame'} = shift if @_;
    return $self->{'_parent_frame'};
}

=head2 svg

 Title   : svg
 Usage   : $obj->svg($newval)
 Function: 
 Example : 
 Returns : value of svg (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub svg{
    my $self = shift;

    return $self->{'svg'} = shift if @_;
    return $self->{'svg'};
}

=head2 xsize

 Title   : xsize
 Usage   : $obj->xsize($newval)
 Function: 
 Example : 
 Returns : value of xsize (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub add_xsize {return shift->xsize(@_)}
sub xsize{
    my $self = shift;

    return $self->{'xsize'} = shift if @_;
    return $self->{'xsize'};
}

=head2 ysize

 Title   : ysize
 Usage   : $obj->ysize($newval)
 Function: 
 Example : 
 Returns : value of ysize (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub add_ysize {return shift->ysize(@_)}
sub ysize{
    my $self = shift;

    return $self->{'ysize'} = shift if @_;
    return $self->{'ysize'};
}

=head2 xoffset

 Title   : xoffset
 Usage   : $obj->xoffset($newval)
 Function: 
 Example : 
 Returns : value of xoffset (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub add_xoffset {return shift->xoffset(@_)}
sub xoffset{
    my $self = shift;

    return $self->{'xoffset'} = shift if @_;
    return $self->{'xoffset'};
}

=head2 yoffset

 Title   : yoffset
 Usage   : $obj->yoffset($newval)
 Function: 
 Example : 
 Returns : value of yoffset (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub add_yoffset {return shift->yoffset(@_)}
sub yoffset{
    my $self = shift;

    return $self->{'yoffset'} = shift if @_;
    return $self->{'yoffset'};
}

=head2 xmin

 Title   : xmin
 Usage   : $obj->xmin($newval)
 Function: 
 Example : 
 Returns : value of xmin (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub xmin{
    my $self = shift;

    return $self->{'xmin'} = shift if @_;
    return $self->{'xmin'} if defined $self->{'xmin'};
    $self->_recalculate_stats();
    return $self->xstat->min;

}

=head2 xmax

 Title   : xmax
 Usage   : $obj->xmax($newval)
 Function: 
 Example : 
 Returns : value of xmax (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub xmax{
    my $self = shift;

    return $self->{'xmax'} = shift if @_;
    return $self->{'xmax'} if defined $self->{'xmax'};
    $self->_recalculate_stats();
    return $self->xstat->max;

}

=head2 ymin

 Title   : ymin
 Usage   : $obj->ymin($newval)
 Function: 
 Example : 
 Returns : value of ymin (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub ymin{
    my $self = shift;

    return $self->{'ymin'} = shift if @_;
    return $self->{'ymin'} if defined $self->{'ymin'};
    $self->_recalculate_stats();
    return $self->ystat->min;

}

=head2 ymax

 Title   : ymax
 Usage   : $obj->ymax($newval)
 Function: 
 Example : 
 Returns : value of ymax (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub ymax{
    my $self = shift;

    return $self->{'ymax'} = shift if @_;
    return $self->{'ymax'} if defined $self->{'ymax'};
    $self->_recalculate_stats();
    return $self->ystat->max;

}

=head2 xrange

 Title   : xrange
 Usage   : $obj->xrange($newval)
 Function: 
 Example : 
 Returns : value of xrange (a scalar)


=cut

sub xrange{
    my $self = shift;

    return $self->xmax - $self->xmin;
}

=head2 yrange

 Title   : yrange
 Usage   : $obj->yrange($newval)
 Function: 
 Example : 
 Returns : value of yrange (a scalar)


=cut

sub yrange{
    my $self = shift;

    return $self->ymax - $self->ymin;
}

=head2 stack

 Title   : stack
 Usage   : $obj->stack($newval)
 Function: 
 Example : 
 Returns : value of stack (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub stack{
    my $self = shift;

    return $self->{'stack'} = shift if @_;
    return $self->{'stack'};
}

sub add_frame_transform{return shift->frame_transform(@_)}
sub frame_transform {
    my $self = shift;

    return $self->{'frame_transform'} = shift if @_;
    return $self->{'frame_transform'} if defined $self->{'frame_transform'};
}

1;
