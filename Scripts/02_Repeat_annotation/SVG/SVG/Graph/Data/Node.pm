package SVG::Graph::Data::Node;

use strict;
use base qw(Tree::DAG_Node);

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

  $self->SUPER::_init;

  foreach my $arg (keys %args) {
	my $meth = $arg;
	if($self->can($meth)){
	  $self->$meth($args{$arg});
	} else {
	  $self->_style($arg => $args{$arg});
	}
  }

}

=head2 depth

 Title   : depth
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub depth{
   my ($self,@args) = @_;

   my $depth = $self->branch_length;

   my $maxdepth = 0;
   foreach my $daughter ($self->daughters){
	 my $ddepth = $daughter->depth;
	 $maxdepth = $ddepth > $maxdepth ? $ddepth : $maxdepth;
   }

   return $depth + $maxdepth;
}

=head2 branch_length

 Title   : branch_length
 Usage   : $obj->branch_length($newval)
 Function: 
 Example : 
 Returns : value of branch_length (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub branch_length{
    my $self = shift;

    return $self->{'branch_length'} = shift if @_;
    return $self->{'branch_length'} || 1;
}

=head2 branch_type

 Title   : branch_type
 Usage   : $obj->branch_type($newval)
 Function: 
 Example : 
 Returns : value of branch_type (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub branch_type{
    my $self = shift;

    return $self->{'branch_type'} = shift if @_;
    return $self->{'branch_type'};
}

=head2 branch_label

 Title   : branch_label
 Usage   : $obj->branch_label($newval)
 Function: 
 Example : 
 Returns : value of branch_label (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub branch_label{
    my $self = shift;

    return $self->{'branch_label'} = shift if @_;
    return $self->{'branch_label'};
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

1;
