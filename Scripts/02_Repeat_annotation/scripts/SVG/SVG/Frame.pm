package SVG::Graph::Group;

use base SVG::Graph::Data;
use strict;
use Data::Dumper;

=head2 new

 Title   : new
 Usage   :
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

  my $id = 'n'.sprintf("%07d",int(rand(9999999)));
  my $group = $self->_parent_svg->svg->group(id=>$id);

  foreach my $arg (keys %args){
	my $meth = 'add_'.$arg;
        $self->$meth($args{$arg});
  }
  $self->svg($group);
  $self->is_changed(1);
}

=head2 add_glyph

 Title   : add_glyph
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


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

=head2 add_group

 Title   : add_group
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub add_group {
  my($self,@groups) = @_;

  my $epitaph = "only SVG::Graph::Group objects accepted";

  if(scalar(@groups)){
	foreach my $group (@groups){
	  if(ref $group eq 'ARRAY'){
		foreach my $g (@$group){
		  die $epitaph unless ref $g eq 'SVG::Graph::Group';
		  push @{$self->{groups}}, $g;
		}
	  } else {
		die $epitaph unless ref $group eq 'SVG::Graph::Group';
		push @{$self->{groups}}, $group;
	  }
	}
  } else {
   my $group = SVG::Graph::Group->new(svg=>$self->_parent_svg,
									  _parent_group=>$self,
									  xoffset=>$self->_parent_svg->margin,
									  yoffset=>$self->_parent_svg->margin,
									  xsize=>$self->_parent_svg->width  - (2 * $self->_parent_svg->margin),
									  ysize=>$self->_parent_svg->height - (2 * $self->_parent_svg->margin),
									 );	

   $group->stack($self->stack) if $self->stack;
   $group->ystat($self->ystat) if $self->stack;

   push @{$self->{groups}}, $group;
	return $group;
  }

  $self->is_changed(1);
}

=head2 groups

 Title   : groups
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub groups {
  my $self = shift;
  return $self->{groups} ? @{$self->{groups}} : ();
}

=head2 add_data

 Title   : add_data
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


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

  if(($self->_parent_group && $flag) || !$self->_parent_group){
	my @data = $self->data;

	#recurse down into subgroups...
	foreach my $subgroup ($self->groups){
	  push @data, $subgroup->all_data(1);
	}

	return map {$_->can('data') ? $_->all_data(1) : $_ } @data;
  } elsif($self->_parent_group) {
	return $self->_parent_group->all_data;
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

  #recurse down into subgroups...
  foreach my $subgroup ($self->groups){
	push @data, $subgroup->data;
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

  #recurse down into subgroups...
  foreach my $subgroup ($self->groups){
	push @data, $subgroup->data_chunks;
  }

  return @data;
}

=head2 draw

 Title   : draw
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub draw {
  my($self, $svg) = @_;

  foreach my $group ($self->groups){
#warn $group;
	$group->draw($self);
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

=head2 _parent_group

 Title   : _parent_group
 Usage   : $obj->_parent_group($newval)
 Function: 
 Example : 
 Returns : value of _parent_group (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub add__parent_group{return shift->_parent_group(@_)}
sub _parent_group{
    my $self = shift;

    return $self->{'_parent_group'} = shift if @_;
    return $self->{'_parent_group'};
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

1;
