package SVG::Graph::Data::Tree;

use strict;
use SVG::Graph::Data::Node;

=head2 new

 Title   : new
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub new {
  my($class, @args) = @_;
  my $self = bless {}, $class;
  $self->init(@args);
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

  $self->root(SVG::Graph::Data::Node->new);
  $self->root->name('root');

  foreach my $arg (keys %args) {
	$self->$arg($args{$arg});
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

   return $self->root->depth;
}

=head2 root

 Title   : root
 Usage   : $obj->root($newval)
 Function: 
 Example : 
 Returns : value of root (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub root{
    my $self = shift;

    return $self->{'root'} = shift if @_;
    return $self->{'root'};
}

=head2 new_node

 Title   : new_node
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub new_node{
   my ($self,@args) = @_;

   return SVG::Graph::Data::Node->new(@args);
}

1;
