package SVG::Graph::Glyph;

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

  my $id = sprintf("%07d",int(rand(9999999)));
  my($glyphname) = ref($self) =~ /::([^:]+)$/;


  foreach my $arg (keys %args){
	my $meth = $arg;
	#caveat, all constructor args will now pass
	if($self->can($meth)){
	  $self->$meth($args{$arg});
	} else {
	  $self->_style($arg => $args{$arg});
	}
  }

  $self->canvas($self->svg->group(id=>"$glyphname$id")) unless $self->canvas;
}

=head2 _style

 Title   : _style
 Usage   : $obj->_style($newval)
 Function: 
 Example : 
 Returns : 
 Args    : 


=cut

sub _style{
    my $self = shift;
	my($key,$val) = @_;

	if(defined($key) and not defined($val)){
	  return $self->{'_style'}{$key};
	} elsif(defined($key) and defined($val)){
	  $self->{'_style'}{$key} = $val;
	  return $val;
	} else {
	  return $self->{'_style'} ? %{$self->{'_style'}} : ();
	}
}

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

   die "method undefined by ".__PACKAGE__;
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

=head2 group

 Title   : group
 Usage   : $obj->group($newval)
 Function: 
 Example : 
 Returns : value of group (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub group{
    my $self = shift;

    return $self->{'group'} = shift if @_;
    return $self->{'group'};
}

=head2 xsize

 Title   : xsize
 Usage   : $obj->xsize($newval)
 Function: 
 Example : 
 Returns : value of xsize (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

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

sub yoffset{
    my $self = shift;

    return $self->{'yoffset'} = shift if @_;
    return $self->{'yoffset'};
}

=head2 xscale

 Title   : xscale
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub xscale{
   my ($self,@args) = @_;

   return $self->xsize / $self->group->xrange;
}

=head2 yscale

 Title   : yscale
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub yscale{
   my ($self,@args) = @_;

   return $self->ysize / $self->group->yrange;
}

=head2 canvas

 Title   : canvas
 Usage   : $obj->canvas($newval)
 Function: 
 Example : 
 Returns : value of canvas (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub canvas{
    my $self = shift;

    return $self->{'canvas'} = shift if @_;
    return $self->{'canvas'};
}


1;
