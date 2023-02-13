package SVG::Graph::File;

use strict;

=head2 new

 Title   : new
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub new
{
   my ($class,%args) = @_;
   my $self = bless {}, $class;
   $self->init(%args);
   return($self);
}

=head2 init

 Title   : init
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub init
{
  my ($self, %args) = @_;
  foreach my $arg (keys %args){
	my $meth = $arg;
	if($self->can($meth)){
	  $self->$meth($args{$arg});
	} else {
	  $self->_style($arg => $args{$arg});
	}
  }
}


=head2 read_data

 Title   : read_data
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub read_data{
   my ($self,@args) = @_;

   die "method undefined by ".__PACKAGE__;
}

=head2 write_data

 Title   : write_data
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub write_data{
   my ($self,@args) = @_;

   die "method undefined by ".__PACKAGE__;
}


=head2 _style

 Title   : _style
 Usage   : $obj->_style($newval)
 Function: 
 Example : 
 Returns : 
 Args    : 


=cut

sub _style
{
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
