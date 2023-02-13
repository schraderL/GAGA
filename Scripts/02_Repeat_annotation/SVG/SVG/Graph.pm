package SVG::Graph;

=head1 NAME

SVG::Graph - Visualize your data in Scalable Vector Graphics (SVG) format.

=head1 SYNOPSIS

  use SVG::Graph;
  use SVG::Graph::Data;
  use SVG::Graph::Data::Datum;

  #create a new SVG document to plot in...
  my $graph = SVG::Graph->new(width=>600,height=>600,margin=>30);

  #and create a frame to hold the data/glyphs
  my $frame = $graph->add_frame;

  #let's plot y = x^2
  my @data = map {SVG::Graph::Data::Datum->new(x=>$_,y=>$_^2)}
                 (1,2,3,4,5);
  my $data = SVG::Graph::Data->new(data => \@data);

  #put the xy data into the frame
  $frame->add_data($data);

  #add some glyphs to apply to the data in the frame
  $frame->add_glyph('axis',        #add an axis glyph
    'x_absolute_ticks' => 1,       #with ticks every one
                                   #unit on the x axis
    'y_absolute_ticks' => 1,       #and ticks every one
                                   #unit on the y axis

    'stroke'           => 'black', #draw the axis black
    'stroke-width'     => 2,       #and 2px thick
  );

  $frame->add_glyph('scatter',     #add a scatterplot glyph
    'stroke' => 'red',             #the dots will be outlined
                                   #in red,
    'fill'   => 'red',             #filled red,
    'fill-opacity' => 0.5,         #and 50% opaque
  );

  #print the graphic
  print $graph->draw;

=head1 DESCRIPTION

SVG::Graph is a suite of perl modules for plotting data.  SVG::Graph
currently supports plots of one-, two- and three-dimensional data, as
well as N-ary rooted trees.  Data may be represented as:

 Glyph Name	Dimensionality supported
          	1d	2d	3d	tree
 --------------------------------------------------------
 Axis				x
 Bar Graph		x
 Bubble Plot			x
 Heatmap Graph			x
 Line Graph		x
 Pie Graph	x
 Scatter Plot		x
 Spline Graph		x
 Tree					x

SVG::Graph 0.02 is a pre-alpha release. Keep in mind that many of the
glyphs are not very robust. 

=head1 PLOTTING

You need to create a SVG::Graph::Frame instance from the parent
SVG::Graph instance for each set of data to be plotted.  Datasets
can be hierarchical, and to represent this, SVG::Graph::Frame
instances can themselves contain subframes.  SVG::Graph::Frame can
contain:

 - multiple subframes as instances of SVG::Graph::Frame
 - a single SVG::Graph::Data instance
 - multiple SVG::Graph::Glyph instances with which to render
   the attached SVG::Graph::Data instance, and all SVG::Graph::Data
   instances attached to SVG::Graph::Frame subinstances

See L<SVG::Graph::Frame> and L<SVG::Graph::Glyph> for details.

=head2 ONE DATA SET

 1. create an SVG::Graph instance
 2. create an SVG::Graph::Frame instance by calling
    SVG::Graph::add_frame();
 3. create an SVG::Graph::Data instance, containing
    an SVG::Graph::Data::Datum instance for each data point.
 4. Attach the SVG::Graph::Data instance to your SVG::Graph::Frame
    using SVG::Graph::Frame::add_data();
 5. Attach glyphs to the SVG::Graph::Frame instance using
    SVG::Graph::Frame::add_glyph();
 6. Call SVG::Graph::draw();

=head2 MULTIPLE DATA SETS

 1. create an SVG::Graph instance
 2. create an SVG::Graph::Frame instance by calling
    SVG::Graph::add_frame();
 3. create an SVG::Graph::Data instance, containing
    an SVG::Graph::Data::Datum instance for each data point.
 4. Attach the SVG::Graph::Data instance to your SVG::Graph::Frame
    using SVG::Graph::Frame::add_data();
 5. Attach glyphs to the SVG::Graph::Frame instance using
    SVG::Graph::Frame::add_glyph();
 6. repeat [2-5] for each additional data set to be added.
    add_frame() can be called on SVG::Graph to add top-level data
    sets, or SVG::Graph::Frame to add hierarchical data sets.
 7. Call SVG::Graph::draw();

=head1 FEEDBACK

Send an email to the svg-graph-developers list.  For more info,
visit the project page at http://www.sf.net/projects/svg-graph

=head1 AUTHORS

 Allen Day,      <allenday@ucla.edu>
 Chris To,       <crsto@ucla.edu>

=head1 CONTRIBUTORS

 James Chen,     <chenj@seas.ucla.edu>
 Brian O'Connor, <boconnor@ucla.edu>

=head1 SEE ALSO

L<SVG>

=cut

use SVG;
use SVG::Graph::Frame;

use Data::Dumper;
use strict;
our $VERSION = '0.02';

=head2 new

 Title   : new
 Usage   : my $graph = SVG::Graph->new(width=>600,
                                       height=>600,
                                       margin=>20);
 Function: creates a new SVG::Graph object
 Returns : a SVG::Graph object
 Args    : width => the width of the SVG 
           height => the height of the SVG
           margin => margin for the root frame


=cut

sub new{
   my ($class,@args) = @_;

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

sub init{
  my($self, %args) = @_;

  foreach my $arg (keys %args){
	my $meth = $arg;
        $self->$meth($args{$arg});
  }

  #allow passing of an existing SVG
  if(!$self->svg){
	$self->svg(SVG->new(xmlns=>"http://www.w3.org/2000/svg",width=>$self->width,height=>$self->height));
  }
}

=head2 width

 Title   : width
 Usage   : $obj->width($newval)
 Function: 
 Example : 
 Returns : value of width (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub width{
    my $self = shift;

    return $self->{'width'} = shift if @_;
    return $self->{'width'};
}

=head2 height

 Title   : height
 Usage   : $obj->height($newval)
 Function: 
 Example : 
 Returns : value of height (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub height{
    my $self = shift;

    return $self->{'height'} = shift if @_;
    return $self->{'height'};
}

=head2 margin

 Title   : margin
 Usage   : $obj->margin($newval)
 Function: 
 Example : 
 Returns : value of margin (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub margin{
    my $self = shift;

    return $self->{'margin'} = shift if @_;
    return $self->{'margin'};
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

=head2 add_frame

 Title   : add_frame
 Usage   : my $frame = $graph->add_frame
 Function: adds a Frame to the current Graph
 Returns : a SVG::Graph::Frame object
 Args    : a hash.  usable keys:
             frame_transform (optional)
               'top' default orientation
               'bottom' rotates graph 180 deg (about the center)
               'right' points top position towards right
               'left' points top position towards left

=cut

sub add_frame{
   my ($self,%args) = @_;

   my $margin = $self->margin || 0;
   my $height = $self->height || 0;
   my $width  = $self->width  || 0;
   my $xoffset = $self->xoffset || 0;
   my $yoffset = $self->yoffset || 0;

   my $frame = SVG::Graph::Frame->new(svg=>$self,
									  xoffset=>$xoffset + $margin,
									  yoffset=>$yoffset + $margin,
									  xsize=>$width  - (2 * $margin),
									  ysize=>$height - (2 * $margin),
									  frame_transform=>$args{frame_transform}
									 );

   #print STDERR Dumper($frame);

   push @{$self->{frames}}, $frame;
   return $frame;
}

=head2 frames

 Title   : frames
 Usage   : get/set
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub frames{
   my ($self,@args) = @_;

   return $self->{frames} ? @{$self->{frames}} : ();
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

=head2 draw

 Title   : draw
 Usage   : $graph=>draw
 Function: depends on child glyph implementations 
 Returns : xmlifyied SVG object
 Args    : none


=cut

sub draw{
   my ($self,@args) = @_;

   foreach my $frame ($self->frames){
	 $frame->draw;
   }

   return $self->svg->xmlify;
}


1;
__END__
