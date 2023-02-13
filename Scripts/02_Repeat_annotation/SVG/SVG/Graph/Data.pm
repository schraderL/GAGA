package SVG::Graph::Data;

use strict;
use Statistics::Descriptive;
use Data::Dumper;

=head2 new

 Title   : new
 Usage   : my $data = SVG::Graph::Data->new
 Function: creates a new SVG::Graph::Data object
 Returns : a SVG::Graph::Data object
 Args    : (optional) array of SVG::Graph::Data::Datum objects


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
  foreach my $arg (keys %args){
	my $meth = 'add_'.$arg;
        $self->$meth($args{$arg});
  }
  $self->is_changed(1);
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

  return $self->{data} ? @{$self->{data}} : ();
}

=head2 add_data

 Title   : add_data
 Usage   : $data->add_data($datum)
 Function: adds a Datum object to the current Data object
 Returns : none
 Args    : SVG::Graph::Data::Datum object


=cut

sub add_data {
  my($self,@data) = @_;

  my $epitaph = "only SVG::Graph::Data::Datum objects accepted";

  foreach my $data (@data){
	if(ref $data eq 'ARRAY'){
	  foreach my $d (@$data){
		die $epitaph unless ref $d eq 'SVG::Graph::Data::Datum';
		push @{$self->{data}}, $d;
	  }
	} else {
	  die $epitaph unless ref $data eq 'SVG::Graph::Data::Datum';
	  push @{$self->{data}}, $data;
	}
  }

  $self->is_changed(1);
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

   #x
   my $xstat = Statistics::Descriptive::Full->new();
   $xstat->add_data(map {$_->x} $self->data);
   $self->xstat($xstat);

   #y
   my $ystat = Statistics::Descriptive::Full->new();
   $ystat->add_data(map {$_->y} $self->data);
   $self->ystat($ystat);

   #z
   my $zstat = Statistics::Descriptive::Full->new();
   $zstat->add_data(map {$_->z} $self->data);
   $self->zstat($zstat);

   $self->is_changed(0);
}

=head2 xstat

 Title   : xstat
 Usage   : $obj->xstat($newval)
 Function: 
 Example : 
 Returns : value of xstat (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub xstat{
    my $self = shift;

    return $self->{'xstat'} = shift if @_;
    return $self->{'xstat'};
}

=head2 ystat

 Title   : ystat
 Usage   : $obj->ystat($newval)
 Function: 
 Example : 
 Returns : value of ystat (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub ystat{
    my $self = shift;

    return $self->{'ystat'} = shift if @_;
    return $self->{'ystat'};
}

=head2 zstat

 Title   : zstat
 Usage   : $obj->zstat($newval)
 Function: 
 Example : 
 Returns : value of zstat (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub zstat{
    my $self = shift;

    return $self->{'zstat'} = shift if @_;
    return $self->{'zstat'};
}


sub xmean       {$_[0]->_recalculate_stats; return $_[0]->xstat->mean}
sub xmode       {$_[0]->_recalculate_stats; return $_[0]->xstat->mode}
sub xmedian     {$_[0]->_recalculate_stats; return $_[0]->xstat->median}
sub xmin        {$_[0]->_recalculate_stats; return $_[0]->xstat->min}
sub xmax        {$_[0]->_recalculate_stats; return $_[0]->xstat->max}
sub xrange      {$_[0]->_recalculate_stats; return $_[0]->xstat->sample_range}
sub xstdv       {$_[0]->_recalculate_stats; return $_[0]->xstat->standard_deviation}
sub xpercentile {$_[0]->_recalculate_stats; return $_[0]->xstat->percentile($_[1])}

sub ymean       {$_[0]->_recalculate_stats; return $_[0]->ystat->mean}
sub ymode       {$_[0]->_recalculate_stats; return $_[0]->ystat->mode}
sub ymedian     {$_[0]->_recalculate_stats; return $_[0]->ystat->median}
sub ymin        {$_[0]->_recalculate_stats; return $_[0]->ystat->min}
sub ymax        {$_[0]->_recalculate_stats; return $_[0]->ystat->max}
sub yrange      {$_[0]->_recalculate_stats; return $_[0]->ystat->sample_range}
sub ystdv       {$_[0]->_recalculate_stats; return $_[0]->ystat->standard_deviation}
sub ypercentile {$_[0]->_recalculate_stats; return $_[0]->ystat->percentile($_[1])}

sub zmean       {$_[0]->_recalculate_stats; return $_[0]->zstat->mean}
sub zmode       {$_[0]->_recalculate_stats; return $_[0]->zstat->mode}
sub zmedian     {$_[0]->_recalculate_stats; return $_[0]->zstat->median}
sub zmin        {$_[0]->_recalculate_stats; return $_[0]->zstat->min}
sub zmax        {$_[0]->_recalculate_stats; return $_[0]->zstat->max}
sub zrange      {$_[0]->_recalculate_stats; return $_[0]->zstat->sample_range}
sub zstdv       {$_[0]->_recalculate_stats; return $_[0]->zstat->standard_deviation}
sub zpercentile {$_[0]->_recalculate_stats; return $_[0]->zstat->percentile($_[1])}

=head2 is_changed

 Title   : is_changed
 Usage   : $obj->is_changed($newval)
 Function: 
 Example : 
 Returns : value of is_changed (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub is_changed{
    my $self = shift;

    return $self->{'is_changed'} = shift if @_;
    return $self->{'is_changed'};
}

=head2 svg

 Title   : svg
 Usage   : $obj->svg($newval)
 Function: 
 Example : 
 Returns : value of svg (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub add_svg {return shift->svg(@_)};
sub svg{
    my $self = shift;

    return $self->{'svg'} = shift if @_;
    return $self->{'svg'};
}


1;
