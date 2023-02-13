package SVG::Graph::Data::Datum;

use strict;
#use overload
#  '""'  => \&label,
#  '<=>' => sub { my($x,$y) = &check; $x <=> $y },
#  '+'   => sub { my($x,$y) = &check; $x+$y     },
#  '-'   => sub { my($x,$y) = &check; $x-$y     },
#  '*'   => sub { my($x,$y) = &check; $x*$y     },
#  '/'   => sub { my($x,$y) = &check; $x/$y     },
#  ;

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
  foreach my $arg (keys %args) {
	$self->$arg($args{$arg});
  }
}

=head2 x

 Title   : x
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub x {
  my($self,$arg) = @_;
  $self->{x} = $arg if defined $arg;
  return $self->{x};
}

=head2 y

 Title   : y
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub y {
  my($self,$arg) = @_;
  $self->{y} = $arg if defined $arg;
  return $self->{y};
}

=head2 z

 Title   : z
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub z {
  my($self,$arg) = @_;
  $self->{z} = $arg if defined $arg;
  return $self->{z};
}

=head2 label

 Title   : label
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub label {
  my($self,$arg) = @_;
  $self->{label} = $arg if defined $arg;
  return $self->{label};
}

=head2 check

 Title   : check
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub check {
  my($x,$y) = @_;

  $x = $x->x if ref $x eq __PACKAGE__;
  $y = $y->x if ref $y eq __PACKAGE__;

  return($x,$y);
}

1;
