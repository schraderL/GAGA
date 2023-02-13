package SVG::SVG2zinc::Conversions;

use Math::Trig;
use Math::Bezier::Convert;
use strict;
use Carp;

use vars qw( $VERSION @ISA @EXPORT );

($VERSION) = sprintf("%d.%02d", q$Revision: 1.10 $ =~ /(\d+)\.(\d+)/);

@ISA = qw( Exporter );

@EXPORT = qw( InitConv
	      removeComment convertOpacity
	      createNamedFont
	      defineNamedGradient namedGradient namedGradientDef existsGradient
	      extractGradientTypeAndStops addTransparencyToGradient 
	      colorConvert
	      pathPoints points
	      cleanName
	      float2int sizesConvert sizeConvert
	      transform
	      );

# some variables to be initialized at the beginning

my ($warnProc, $lineNumProc); # two proc
my %fonts; # a hashtable to identify all used fonts 
my %gradients;

sub InitConv {
    ($warnProc, $lineNumProc) = @_;
    %fonts = ();
    %gradients = ();
    return 1;
}

sub myWarn{
    &{$warnProc}(@_);
}

### remove SVG comments in the form  /* */  in $str
### returns the string without these comments
sub removeComment {
    my ($str) = @_;
#    my $strOrig = $str;
    return "" unless defined $str;

    while ($str =~ s|(.*)(/\*.*\*/){1}?|$1|) {
#	print "begin='$str'\n";
    }
#    print "'$strOrig' => '$str'\n";
    $str =~ s/^\s*// ;
    return $str;
}

## returns an opacity value between 0 and 1
## returns 1 if the argument is undefined
sub convertOpacity {
    my ($opacity) = @_;
    $opacity = 1 unless defined $opacity;
    $opacity = 0 if $opacity<0;
    $opacity = 1 if $opacity>1;
    return $opacity;
}


######################################################################################
#   fontes management
######################################################################################

# the following hashtable is used to maps SVG font names to X font names
# BUG: obvioulsy this hashtable should be defined in the system or at
# least as a configuration file or in the SVG2zinc parser parameters
my %fontsMapping =
    ( 'comicsansms' => "comic sans ms",
#      'helvetica'   => "arial", # "verdana",
      'arialmt' => "arial",
      );

sub createNamedFont {
    my ($fullFamily, $size, $weight) = @_;
    $fullFamily = "verdana" if $fullFamily eq "";
    my $family = lc($fullFamily);

    $weight = "normal" unless $weight;  ## valeur par défaut

    if ( $size =~ /(.*)pt/ ) {
	## size in points
	$size = $1;
    } elsif ( $size =~ /(\d*(.\d*)?)\s*$/ ) {
	## size in pixel
	## BUG: generates a bug in TkZinc when render != 0 (TBC)
	$size = -$1;
    }

    $size = &float2int($size); # I round the font size, at least until we have vectorial font in Tk::Zinc
    
    if ( $family =~ /(\w*)-bold/ ) {
	$family = $1;
	$weight = "bold";  # this might be in contradiction with the wieght defined in SVG (??)
    } else {
	$weight = "medium";
    }
    $family = $fontsMapping{$family} if defined $fontsMapping{$family};
#    print "FontFamily: '$fullFamily' => '$family'\n";
    
    my $fontKey = join "_", ($family, $size, $weight);
    if (!defined $fonts{$fontKey}) {
	$fonts{$fontKey} = $fontKey;
	print "In createNamedFont, a new font: $fontKey\n";
	return ($fontKey, "->fontCreate('$fontKey', -family => \"$family\", -size => $size, -weight => \"$weight\");");
    } else {
	return ($fontKey,"");
    }
    
} # end of createNamedFont

######################################################################################
#   gradients management
######################################################################################
# my %gradients;

## Check if the new gradient does not already exists (with another name)
## In this case, the hash is extended with an "auto-reference"
##   $gradients{newName} = "oldName"
## and the function returns 0
## Otherwise, add an entry in the hastable
##   $gradients{newName} = "newDefinition"
## and returns 1
sub defineNamedGradient {
    my ($newGname, $newGradDef) = @_;
    my $prevEqGrad;
    $newGradDef =~ s/^\s*(.*\S)\s*$/$1/ ; # removing trailing/leading blank
    $newGradDef =~ s/\s*\|\s*/ \| /g ;    # inserting blanks around the |
    $newGradDef =~ s/\s\s+/ /g;  # removing multiple occurence of blanks
#    print "CLEANED grad='$newGradDef'\n";
    foreach my $gname (keys %gradients) {
	if ($gradients{$gname} eq $newGradDef) {
	    ## such a gradient already exist with another name
	    $gradients{$newGname} = $gname;
#	    print "GRADIENT: $newGname == $gname\n";

#	    $res .= "\n###### $newGname => $gname"; ### 

	    return 0;
	}
    }
    ## there is no identical gradient with another name
    ## we add the definition in the hashtable
    $gradients{$newGname} = $newGradDef;
    return $newGradDef;
}

## returns the name of a gradient, by following if necessary
## "auto-references" in the hashtable
sub namedGradient {
    my ($gname) = @_;
    my $def = $gradients{$gname};
    return $gname unless defined $def; 
    ## to avoid looping if the hashtable is buggy:
    return $gname if !defined $gradients{$def} or $def eq $gradients{$def}; 
    return &namedGradient($gradients{$gname});
}

## returns the definition associated to a named gradient, following if necessary
## "auto-references" in the hashtable
sub namedGradientDef {
    my ($gname) = @_;
    my $def = $gradients{$gname};
    return "" unless defined $def; 
    ## to avoid looping if the hashtable is buggy:
    return $def if !defined $gradients{$def} or $def eq $gradients{$def}; 
    return $gradients{&namedGradient($gradients{$gname})};
}

# returns 1 if the named has an associated gradient
sub existsGradient {
    my ($gname) = @_;
    if (defined $gradients{$gname}) {return 1} else {return 0}; 
}

## this function returns both the radial type with its parameters AND
## a list of stops characteristics as defined in TkZinc
## usage: ($radialType, @stops) = &extractGradientTypeAndStops(<namedGradient>);
## this func assumes that <namedGradient> DOES exist
sub extractGradientTypeAndStops {
    my ($namedGradient) = @_;
    my $gradDef = &namedGradientDef($namedGradient);
    my @defElements = split (/\s*\|\s*/ , $gradDef);
    my $gradientType;
    $gradientType = shift @defElements;
    return ($gradientType, @defElements); 
}

## combines the opacity to every parts of a named gradient
## if some parts of the gradients are themselves partly transparent, they are combined
## if $opacity is 1, returns directly $gname
## else returns a new definition of a gradient
sub addTransparencyToGradient {
    my ($gname,$opacity) = @_;
    return $gname if $opacity == 100;
    &myWarn ("ATTG: ERROR $gname\n"), return $gname if !&namedGradientDef($gname); ## this cas is certainly an error in the SVG source file!
    my ($gradientType, @stops) = &extractGradientTypeAndStops($gname);

    my @newStops;
    foreach my $stop (@stops) {
	my $newStop="";
	if ($stop =~ /^([^\s;]+)\s*;\s*(\d+)\s*(\d*)\s*$/   #   red;45 50   or red;45
	    ) {
	    my ($color,$trans,$pos) = ($1,$2,$3);
#	    print "$stop => '$color','$trans','$pos'\n";
	    my $newtransp = &float2int($trans*$opacity/100);
	    if ($pos) {
		$newStop="$color;$newtransp $pos";
	    } else {
		$newStop="$color;$newtransp";
	    }
	} elsif ($stop =~ /^(\S+)\s+(\d+)$/) {             # red 50
	    my ($color,$pos) = ($1,$2);
#	    print "$stop => '$color','$pos'\n";
	    my $newtransp = &float2int($opacity);
	    $newStop="$color;$newtransp $pos";
	} elsif ($stop =~ /^(\S+)$/) {
	    my ($color) = ($1);
#	        print "$stop => '$color'\n";
	    my $newtransp = &float2int($opacity);
	    $newStop="$color;$newtransp";
	} else {
	    &myWarn ("In addTransparencyToGradient: bad gradient Elements: '$stop'\n");
	}
	push @newStops, $newStop;
    }
    return ( $gradientType . " | " . join (" | ", @newStops));
} # end of addTransparencyToGradient


######################################################################################
# color conversion
######################################################################################
# a hash table to define non-X SVG colors
# THX to Lemort for bug report and correction!
my %color2color = ('lime'    => 'green',
		   'Lime'    => 'green',
		   'crimson' => '#DC143C',
		   'Crimson' => '#DC143C',
		   'aqua'    => '#00ffff',
		   'Aqua'    => '#00ffff',
		   'fuschia' => '#ff00ff',
		   'Fuschia' => '#ff00ff',
		   'fuchsia' => '#ff00ff',
		   'Fuchsia' => '#ff00ff',
		   'indigo'  => '#4b0082',
		   'Indigo'  => '#4b0082',
		   'olive'   => '#808000',
		   'Olive'   => '#808000',
		   'silver'  => '#c0c0c0',
		   'Silver'  => '#c0c0c0',
		   'teal'    => '#008080',
		   'Teal'    => '#008080',
		   'green'   => '#008000',
		   'Green'   => '#008000',
		   'grey'    => '#808080',
		   'Grey'    => '#808080',
		   'gray'    => '#808080',
		   'Gray'    => '#808080',
		   'maroon'  => '#800000',
		   'Maroon'  => '#800000',
		   'purple'  => '#800080',
		   'Purple'  => '#800080',
		   );

#### BUG: this is certainly only a partial implementation!!
sub colorConvert {
    my ($color) = @_;
    if ($color =~ /^\s*none/m) {
	return 'none';
    } elsif ($color =~ /rgb\(\s*(.+)\s*\)/ ) {
	## color like "rgb(...)"
	my $rgbs = $1;
	if ($rgbs =~ /([\d.]*)%\s*,\s*([\d.]*)%\s*,\s*([\d.]*)%/ ) {
  	    ## color like "rgb(1.2% , 45%,67.%)"
	    my ($r,$g,$b) = ($1,$2,$3);
	    $color = sprintf ("#%02x%02x%02x",
			      sprintf ("%.0f",2.55*$r),
			      sprintf ("%.0f",2.55*$g),
			      sprintf ("%.0f",2.55*$b));
	    return $color;
	} elsif ($rgbs =~ /(\d*)\s*,\s*(\d*)\s*,\s*(\d*)/ ) {
	    ## color like "rgb(255, 45,67)"
	    my ($r,$g,$b) = ($1,$2,$3);
	    $color = sprintf "#%02x%02x%02x", $r,$g,$b;
	    return $color;
	    } else {
	    &myWarn ("Unknown rgb color coding: $color\n");
	}
    } elsif ($color =~ /^url\(\#(.+)\)/ ) {
	## color like "url(#monGradient)"
	$color = $1;
	my $res = &namedGradient($color);
	return $res; #&namedGradient($1);
    } elsif ( $color =~ /\#([0-9a-fA-F]{3}?)$/ ) {
	## color like #fc1 => #ffcc11
	$color =~ s/([0-9a-fA-F])/$1$1/g ;
        # on doubling the digiys, because Tk does not do it properly
	return $color;
    } else {
	## named colors!
	## except those in the %color2color, all other should be defined in the
	## standard rgb.txt file
	my $converted = $color2color{lc($color)}; # THX to Lemort for bug report!
	if (defined $converted) {
	    return $converted;
	} else {
	    return $color;
	}
    }
} # end of colorConvert

######################################################################################
# path points commands conversion
######################################################################################


# &pathPoints (\%attrs)
# returns a boolean and a list of table references
# - the boolean is true is the path has more than one contour or if it must be closed
# - every table referecne pints to a table of strings, each string describing coordinates
# possible BUG: in Tk::Zinc when a curve has more than one contour, they are all closed
#  how is it in SVG?
sub pathPoints {
    my ($ref_attrs) = @_;
    my $str = $ref_attrs->{d};
#    print "#### In PathPoints : $str\n";
    my ($x,$y) = (0,0); # current values
    my $closed = 1;
    my $atLeastOneZ=0; # true if at least one z/Z command. The curve must then be closed
    my @fullRes;
    my @res ;
    my ($firstX, $firstY); # for memorizing the first point for a 'm' command after a 'z'!
    my ($prevContrlx,$prevContrly); # useful for the s/S commande

    # I use now a repetitive search on the same string, without allocating
    # a $last string for the string end; with very long list of points, such
    # as iceland.svg, we can gain 30% in this function and about 3s over 30s
    while ( $str =~ m/\s*([aAmMzZvVhHlLcCsSqQtT])\s*([^aAmMzZvVhHlLcCsSqQtT]*)\s*/g ) {
	my ($command, $args)=($1,$2);
	&myWarn ("!!!! Ill-formed path command: '", substr($str,pos($str), 40), "...'\n") unless defined $command ;
#	print "Command=$command args=$args x=$x y=$y\n";
	if ($command eq "M") {   ## moveto absolute
	    if (!$closed) {
		## creating a new contour
		push @fullRes, [ @res ];
		$atLeastOneZ = 1;
		@res = ();
	    }
	    my @points = &splitPoints($args);
	    ($prevContrlx,$prevContrly) = (undef,undef);
	    $firstX = $points[0];
	    $firstY = $points[1];
	    while (@points) {
		$x = shift @points;
		$y = shift @points;
		push @res , "[$x, $y]";
	    }
	    next;
	} elsif ($command eq "m") {  ## moveto relative
	    if (!$closed) {
		## creating a new contour
		push @fullRes, [ @res ];
		$atLeastOneZ = 1;
		@res = ();
	    }
	    my @dxy = &splitPoints($args);
	    $firstX = $x+$dxy[0];
	    $firstY = $y+$dxy[1];
# 	    print "m command: $args  => @dxy ,$x,$y\n";
	    while (@dxy) {
		## trying to minimize the number of operation
		## to speed a bit this loop
		$x += shift @dxy;
		$y += shift @dxy;
		push @res, "[$x, $y]";
	    }
	    next;
	} elsif ($command eq 'z' or $command eq 'Z') {
	    push @fullRes, [ @res ];
	    $closed = 1;
	    $atLeastOneZ = 1;
	    @res = ();
	    $x=$firstX;
	    $y=$firstY;
	    next;
	}
	# as a command will/should follow, the curve is no more closed
	$closed = 0;
	if ($command eq "V") {  ## vertival lineto absolute
	    ($y) = $args =~ /(\S+)/m ;  ## XXXX what about multiple y !?
	    push @res , "[$x, $y]";
	} elsif ($command eq "v") {  ## vertical lineto relative
	    my ($dy) = $args =~ /(\S+)/m ;  ## XXXX what about multiple dy !?
	    $y += $dy;
	    push @res , "[$x, $y]";
	} elsif ($command eq "H") {  ## horizontal lineto absolute
	    ($x) = $args =~ /(\S+)/m ;  ## XXXX what about multiple x !?
	    push @res , "[$x, $y]";
	} elsif ($command eq "h") {  ## horizontal lineto relative
	    my ($dx) = $args =~ /(\S+)/m ;  ## XXXX what about multiple dx !?
	    $x += $dx;
	    push @res , "[$x, $y]";
	} elsif ($command eq "L") {  ## lineto absolute
	    my @points = &splitPoints($args);
	    while (@points) {
		$x = shift @points;
		$y = shift @points;
		push @res , "[$x, $y]";
	    }
	} elsif ($command eq "l") {  ## lineto relative
	    ### thioscommand can have more than one point as arguments
	    my @points = &splitPoints($args);
	    # for (my $i = 0; $i < $#points; $i+=2)
	    # is not quicker than the following while
  	    while (@points) {
		## trying to minimize the number of operation
		## to speed a bit this loop
  		$x += shift @points;
  		$y += shift @points;
  		push @res , "[$x, $y]";
  	    }
	} elsif ($command eq "C" or $command eq "c") {  ## cubic bezier
	    &myWarn ("$command command in a path must not be the first one") ,last
		if (scalar @res < 1);
	    my @points = &splitPoints($args);
	    while (@points) {
		&myWarn ("$command command must have 6 coordinates x N times") ,last
		    if (scalar @points < 6);
		my $x1 = shift @points;
		my $y1 = shift @points;
		$prevContrlx = shift @points;
		$prevContrly = shift @points;
		my $xf = shift @points;
		my $yf = shift @points;
		if ($command eq "c") { $x1+=$x; $y1+=$y; $prevContrlx+=$x; $prevContrly+=$y; $xf+=$x; $yf+=$y}
		push @res, ( "[$x1, $y1, 'c'], [$prevContrlx, $prevContrly, 'c'], [$xf, $yf]");
		$x=$xf;
		$y=$yf;
	    }
	} elsif ($command eq "S" or $command eq "s") {  ## cubic bezier with opposite last control point
	    &myWarn ("$command command in a path must not be the first one") ,last
		if (scalar @res < 1);
#	    print "$command command : $args\n";
	    my @points = &splitPoints($args);
	    if ($command eq "s") {
		for (my $i=0; $i <= $#points; $i += 2) {
		    $points[$i] += $x;
		}
		for (my $i=1; $i <= $#points; $i += 2) {
		    $points[$i] += $y;
		}
	    }
	    while (@points) {
		&myWarn ("$command command must have 4 coordinates x N times; skipping @points") ,last
		    if (scalar @points < 4);
		my $x1 = (defined $prevContrlx) ? $prevContrlx : $x;
		$x1 = 2*$x-$x1;
		my $y1 = (defined $prevContrly) ? $prevContrly : $y;
		$y1 = 2*$y-$y1;
		$prevContrlx = shift @points;
		$prevContrly = shift @points;
		$x = shift @points;
		$y = shift @points;
		push @res, ( "[$x1, $y1, 'c'], [$prevContrlx, $prevContrly, 'c'], [$x, $y]");
	    }


	} elsif ($command eq "Q" or $command eq "q") {  ## quadratic bezier
	    &myWarn ("$command command in a path must not be the first one") ,last
		if (scalar @res < 1);
	    my @points = &splitPoints($args);
	    if ($command eq "q") {
		for (my $i=0; $i <= $#points; $i += 2) {
		    $points[$i] += $x;
		}
		for (my $i=1; $i <= $#points; $i += 2) {
		    $points[$i] += $y;
		}
	    }
	    while (@points) {
		&myWarn ("$command command must have 4 coordinates x N times") ,last
		    if (scalar @points < 4);
		$prevContrlx = shift @points;
		$prevContrly = shift @points;
		
		my $last_x = $x;
		my $last_y = $y;

		$x = shift @points;
		$y = shift @points;

		# the following code has been provided by Lemort@intuilab.com
		my @coordsToConvert = ($last_x,$last_y, $prevContrlx, $prevContrly,$x,$y);
		my @convertCoords = Math::Bezier::Convert::quadratic_to_cubic(@coordsToConvert);
		# removing the first point, already present
		splice(@convertCoords, 0, 2);
		
		while (@convertCoords) {
		    my ($ctrl1_x, $ctrl1_y) = splice(@convertCoords, 0, 2);
		    my ($ctrl2_x, $ctrl2_y) = splice(@convertCoords, 0, 2);
		    my ($pt_x, $pt_y) = splice(@convertCoords, 0, 2);

		    push @res, ("[$ctrl1_x, $ctrl1_y, 'c'], [$ctrl2_x, $ctrl2_y, 'c'], [$pt_x, $pt_y]");
		}

	    }

	} elsif ($command eq "T" or $command eq "t") {  ## quadratic bezier with opposite last control point?!
	    &myWarn ("$command command in a path must not be the first one") ,last
		if (scalar @res < 1);
	    my @points = &splitPoints($args);

	    if ($command eq "t") {
		for (my $i=0; $i <= $#points; $i += 2) {
		    $points[$i] += $x;
		}
		for (my $i=1; $i <= $#points; $i += 2) {
		    $points[$i] += $y;
		}
	    }
	    while (@points) {
		&myWarn ("$command command must have 2 coordinates x N times") ,last
		    if (scalar @points < 2);
		my $x1 = (defined $prevContrlx) ? $prevContrlx : $x;
		$prevContrlx = 2*$x-$x1;
		my $y1 = (defined $prevContrly) ? $prevContrly : $y;
		$prevContrly = 2*$y-$y1;
		
		my $last_x = $x;
		my $last_y = $y;

		$x = shift @points;
		$y = shift @points;

		# the following code has been provided by Lemort@intuilab.com
		my @coordsToConvert = ($last_x,$last_y, $prevContrlx, $prevContrly,$x,$y);
		my @convertCoords = Math::Bezier::Convert::quadratic_to_cubic(@coordsToConvert);
		# removing the first point, already present
		splice(@convertCoords, 0, 2);
		
		while (@convertCoords) {
		    my ($ctrl1_x, $ctrl1_y) = splice(@convertCoords, 0, 2);
		    my ($ctrl2_x, $ctrl2_y) = splice(@convertCoords, 0, 2);
		    my ($pt_x, $pt_y) = splice(@convertCoords, 0, 2);

		    push @res, ("[$ctrl1_x, $ctrl1_y, 'c'], [$ctrl2_x, $ctrl2_y, 'c'], [$pt_x, $pt_y]");
		}

	    }
	} elsif ($command eq 'a' or $command eq 'A') {
	    my @points = &splitPoints($args);
	    while (@points) {
		&myWarn ("bad $command command parameters: @points\n") if (scalar @points < 7);
#	    print "($x,$y) $command command: @points\n";
		if ($command eq 'a') {
		    $points[5] += $x;
		    $points[6] += $y;
		}
#	    print "($x,$y) $command command: @points\n";
		my @coords = &arcPathCommand ( $x,$y, @points[0..6] );
		push @res, @coords;
		$x = $points[5];
		$y = $points[6];
		last if (scalar @points == 7);
		@points = @points[7..$#points]; ### XXX à tester!
	    }
	} else {
	    &myWarn ("!!! bad path command: $command\n");
	}
    }
    if (@res) {
	return ( $atLeastOneZ, [@res], @fullRes);
    } else { return ( $atLeastOneZ, @fullRes) }
} # end of pathPoints




# this function can be called many many times; so it has been "optimized"
# even if a bit less readable
sub splitPoints {
    $_ = shift;
    ### adding a space before every dash (-) when the dash preceeds by a digit
    s/(\d)-/$1 -/g;
    ### adding a space before à dot (.) when more than one real are not separated;
    ### e.g.:  '2.3.45.6.' becomes '2.3 .45 .5'
    while ( scalar s/\.(\d+)\.(\d+)/\.$1 \.$2/) {
    }
    return  split ( /[\s,]+/ );
}



sub arcPathCommand {
    my ($x1,$y1, $rx,$ry, $x_rot, $large_arc_flag,$sweep_flag, $x2,$y2) = @_;
    return ($x2,$y2) if ($rx == 0 and $ry == 0);
    $rx = -$rx if $rx < 0;
    $ry = -$ry if $ry < 0;

    # computing the center
    my $phi = deg2rad($x_rot);

    # compute x1' and y1' (formula F.6.5.1)
    my $deltaX = ($x1-$x2)/2;
    my $deltaY = ($y1-$y2)/2;
    my $xp1 =  cos($phi)*$deltaX + sin($phi)*$deltaY;
    my $yp1 = -sin($phi)*$deltaX + cos($phi)*$deltaY;
#    print "xp1,yp1= $xp1 , $yp1\n";

    # the radius_check has been suggested by lemort@intuilab.com
    # checking that radius are correct
    my $radius_check = ($xp1/$rx)**2 + ($yp1/$ry)**2;
    
    if ($radius_check > 1) {
        $rx *= sqrt($radius_check);
        $ry *= sqrt($radius_check);
    }

    # compute the sign:  (formula F.6.5.2)
    my $sign = 1;
    $sign = -1 if $large_arc_flag eq $sweep_flag;
    # compute the big square root  (formula F.6.5.2)
#    print "denominator: ", ( ($rx*$ry)**2 - ($rx*$yp1)**2 - ($ry*$xp1)**2 ),"\n";
    my $bigsqroot = (
		     abs( ($rx*$ry)**2 - ($rx*$yp1)**2 - ($ry*$xp1)**2 ) ### ABS ?!?!
		     /
		     ( ($rx*$yp1)**2 + ($ry*$xp1)**2 )
		     );
    # computing c'x and c'y  (formula F.6.5.2) 
    $bigsqroot = $sign * sqrt ($bigsqroot);
    my $cpx = $bigsqroot * ($rx*$yp1/$ry);
    my $cpy = $bigsqroot * (- $ry*$xp1/$rx);

    # compute cx and cy  (formula F.6.5.3)
    my $middleX = ($x1+$x2)/2;
    my $middleY = ($y1+$y2)/2;
    my $cx = cos($phi)*$cpx - sin($phi)*$cpy + $middleX;
    my $cy = sin($phi)*$cpx + cos($phi)*$cpy + $middleY;

    # computing theta1   (formula F.6.5.5)
    my $XX = ($xp1-$cpx)/$rx;
    my $YY = ($yp1-$cpy)/$ry;
    my $theta1 = rad2deg (&vectorProduct ( 1,0,
					   $XX,$YY));
    # computing dTheta (formula F.6.5.6)
    my $dTheta = rad2deg (&vectorProduct ( $XX,$YY,
					   (-$xp1-$cpx)/$rx,(-$yp1-$cpy)/$ry ));
    # Next To be implemented!!
#    printf "cx,cy=%d,%d\ttheta1,dtheta=%d,%d\trx,ry=%d,%d\n",$cx,$cy,$theta1,$dTheta,$rx,$ry;
    if (!$sweep_flag and $dTheta>0) {
	$dTheta-=360;
    }
    if ($sweep_flag and $dTheta<0) {
	$dTheta+=360;
    }
    return  join (",", &computeArcPoints($cx,$cy,$rx,$ry,
				       $phi,deg2rad($theta1),deg2rad($dTheta))), "\n";
}

sub computeArcPoints {
    my ($cx,$cy,$rx,$ry,$phi,$theta1,$dTheta) = @_;
    my $Nrad = 3.14/18;
    my $N = &float2int(abs($dTheta/$Nrad));
     my $cosPhi = cos($phi);
    my $sinPhi = sin($phi);
#    print "N,dTheta: $N,$dTheta\n";
    my $dd = $dTheta/$N;
    my @res;
    for (my $i=0; $i<=$N; $i++) {
	my $a = $theta1 + $dd*$i;
	my $xp = $rx*cos($a);
	my $yp = $ry*sin($a);
	my $x1 = $cosPhi*$xp - $sinPhi*$yp + $cx;
	my $y1 = $sinPhi*$xp + $cosPhi*$yp + $cy;
	push @res, "[$x1, $y1]";
    }
    return @res;
}

## vectorial product
sub vectorProduct {
    my ($x1,$y1, $x2,$y2) = @_;
    my $sign = 1;
    $sign = -1 if ($x1*$y2 - $y1*$x2) < 0;

    return $sign * acos ( ($x1*$x2 + $y1*$y2)
			  /
			  sqrt ( ($x1**2 + $y1**2) * ($x2**2 + $y2**2) )
			  );
}

######################################################################################
# points conversions for polygone / polyline
######################################################################################

# &points (\%attrs)
# converts the string, value of an attribute points
# to a string of coordinate list for Tk::Zinc
sub points {
    my ($ref_attrs) = @_;
    my $str = $ref_attrs->{points};
    # suppressing leading and trailing blanks:
    ($str) = $str =~ /^\s*   # leading blanks
                      (.*\S) # 
		      \s*$   # trailing blanks
		     /x;

    $str =~ s/([^,])[\s]+([^,])/$1,$2/g ;  # replacing blanks separators by a comma
    return $str;
}

######################################################################################
# cleaning an id to make it usable as a TkZinc Tag
######################################################################################

## the following function cleans an id, ie modifies it so that it
## follows the TkZinc tag conventions.
## BUG: the cleanning is far from being complete
sub cleanName {
    my $id = shift;
     # to avoid numeric ids
    if ($id =~ /^\d+$/) {
#	&myWarn ("id: $id start with digits\n");
	$id = "id_".$id;
    }
    # to avoid any dots in a tag
    if ($id =~ /\./) {
#	&myWarn ("id: $id contains dots\n");
	$id =~ s/\./_/g ;
    }
    return $id;
}

################################################################################
# size conversions
################################################################################

## get a list of "size" attributes as listed in @attrs (e.g.: x y width height...)
## - convert all in pixel
## - return 0 for attributes listed in @attrs and not available in %{$ref_attrs}
sub sizesConvert {
    my ($ref_attrs,@attrs) = @_;
    my %attrs = %{$ref_attrs};
    my @res;
    foreach my $attr (@attrs) {
	my $value;
	if (!defined ($value = $attrs{$attr}) ) {
	    push @res,0;
#	    print "!!!! undefined attr: $attr\n";
	} else {
	    push @res,&sizeConvert ($value);
	}
    }
    return @res;
} # end of sizesConvert

# currently, to simplify this code, I suppose the screen is 100dpi!
# at least the generated code is currently independant from the host
# where is is supposed to run
# maybe this should be enhanced
sub sizeConvert {
    my ($value) = @_;
    if ($value =~ /(.*)cm/) {
	return $1 * 40;   ## approximative pixel / cm
    } elsif ($value =~ /(.*)mm/) {
	return $1 * 4;   ## approximative pixel / mm
    } elsif ($value =~ /(\d+)px/) {
	return $1;   ## exact! pixel / pixel
    } elsif ($value =~ /(.*)in/) {
	return &float2int($1 * 100);   ## approximative pixel / inch
    } elsif ($value =~ /(.*)pt/) {
	return &float2int($1 * 100 / 72);   ## approximative pixel / pt  (a pt = 1in/72)
    } elsif ($value =~ /(.*)pc/) {
	return &float2int($1 * 100 / 6);   ##   (a pica = 1in/6)
    } elsif ($value =~ /(.*)%/) {
	return $1/100;   ## useful for coordinates using % 
	                 ## in lienar gradient (x1,x2,y2,y2)
    } elsif ($value =~ /(.*)em/) { # not yet implemented
	&myWarn ("em unit not yet implemented in sizes");
	return $value;
    } elsif ($value =~ /(.*)ex/) { # not yet implemented
	&myWarn ("ex unit not yet implemented in sizes");
	return $value;
    } else {
	return $value;
    }
} # end of sizeConvert


sub float2int {
    return sprintf ("%.0f",$_[0]);
}


# process a string describing transformations
# returns a list of string describing transformations
# to be applied to Tk::Zinc item Id
sub transform {
    my ($id, $str) = @_;
    return () if !defined $str;
    &myWarn ("!!! Need an Id for applying a transformation\n"), return () if !defined $id;
    my @fullTrans;
    while ($str  =~ m/\s*(\w+)\s*\(([^\)]*)\)\s*/g) {
	my ($trans, $params) = ($1,$2);
	my @params = split (/[\s,]+/, $params);
	if ($trans eq 'translate') {
	    $params[1] = 0 if scalar @params == 1; ## the 2nd paramter defaults to 0
	    my $translation = "->translate($id," . join (",",@params) . ");"  ;
	    push @fullTrans, $translation;
	} elsif ($trans eq 'rotate') {
	    $params[0] = deg2rad($params[0]);
	    my $rotation = "->rotate($id," . join (",",@params) . ");";
	    push @fullTrans, $rotation;
	} elsif ($trans eq 'scale') {
	    $params[1] = $params[0] if scalar @params == 1; ## the 2nd scale parameter defaults to the 1st
	    my $scale = "->scale($id," . join (",",@params) . ");";
	    push @fullTrans,$scale;
	} elsif ($trans eq 'matrix') {
	    my $matrixParams = join ',',@params;
	    my $matrix = "->tset($id, $matrixParams);";
	    push @fullTrans, $matrix;
	} elsif ($trans eq 'skewX'){
	    my $skewX = "->skew($id, " . deg2rad($params[0]) . ",0);";
#	    print "skewX=$skewX\n";
	    push @fullTrans, $skewX;
	} elsif ($trans eq 'skewY'){
	    my $skewY = "->skew($id, 0," . deg2rad($params[0]) . ");";
#	    print "skewY=$skewY\n";
	    push @fullTrans, $skewY;
	} else {
	    &myWarn ("!!! Unknown transformation '$trans'\n");
	}
#	$str = $rest;
    }
    return reverse @fullTrans;
} # end of transform

1;
