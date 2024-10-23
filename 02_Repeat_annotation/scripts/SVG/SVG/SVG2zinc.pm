package SVG::SVG2zinc;

#
#	convertisseur SVG->TkZinc
# 
#	Copyright 2002-2003
#	Centre d'Études de la Navigation Aérienne
#
#	Author: Christophe Mertz <mertz at intuilab dot com> 
#                                previously <mertz at cena dot fr>
#          with many helps from
#          Alexandre Lemort <lemort at intuilab dot com>
#          Celine Schlienger <celine at intuilab dot com>
#          Stéphane Chatty <chatty at intuilab dot com>
#
# $Id: SVG2zinc.pm,v 1.42 2004/05/01 09:19:32 mertz Exp $
#############################################################################
#
# this is the main module of the a converter from SVG file
# to either perl script/module (an eventually other scripting language)
# It is also usable to display SVG graphic file in Tk::Zinc

#############################################################################
# limitations are now listed in the POD at the end of this file
#############################################################################

use strict;
use XML::Parser;
use Carp;
use Math::Trig;
use Tk::PNG;
use Tk::JPEG;
use English;
use File::Basename;

use SVG::SVG2zinc::Conversions;

use vars qw($VERSION $REVISION @ISA @EXPORT);
@EXPORT = qw( parsefile findINC );

$REVISION = q$Revision: 1.42 $ ;
$VERSION = "0.10";

# to suppress some stupid warning usefull for debugging only
my $warn=0;

my $verbose;

my $current_group;
my @prev_groups = ();
my %current_context;
my @prev_contexts = ();

my $itemCount = 0;
my $effectiveItemCount = 0;  # to know if some groups are empty (cf &defs et &defs_)
my $prefix;                  # prefix used in tags associated to generated items
my $colorSep = ";";

sub InitVars {
    @prev_groups = ();
    %current_context = ();
    @prev_contexts = ();
    
    $itemCount = 0;
    $effectiveItemCount = 0;
    $colorSep = ";";
}

# This hash table indicates all non-implemented extensions
# Normaly, the href extension is the only implemented extension listed in the SVG entity
# The hash-value associated to a not implemented etension is 0
# The hash-value is then set to 1 when an warning message has been printed once
my %notImplementedExtensionPrefix; 


# events on "graphics and container elements"
my @EVENT_ON_GRAPHICS = qw/onfocusin onfocusout onactivate onclick
    onmousedown onmouseup onmouseover onmousemove onmouseout onload/ ;
# events on "Document-level event attributes"
my @EVENT_ON_DOC = qw /onunload onabort onerror onresize onscroll onzoom/;
# events "Animation event attributes"
my @EVENT_ON_ANIM = qw /onbegin onend onrepeat/ ;

my %EVENT_ON_GRAPHICS = map { $_ => 1 } @EVENT_ON_GRAPHICS;
my %EVENT_ON_DOC      = map { $_ => 1 } @EVENT_ON_DOC;
my %EVENT_ON_ANIM     = map { $_ => 1 } @EVENT_ON_ANIM;


### @STYLE_ATTRS and %STYLE_ATTRS are "constant" array and hash used in
#### &analyze_style , &analyze_text_style , &groupContext , &attrs_implemented
my @STYLE_ATTRS = qw(class style display fill fill-opacity fill-rule stroke
		     stroke-width stroke-opacity opacity font-size font-family
		     font-weight stroke-linejoin stroke-linecap
		     ) ;
my %STYLE_ATTRS = map { $_ => 1 } @STYLE_ATTRS;

#### not implemented / not implementable attributes
#### these attributes will generate only limited warning
#### used in &attrs_implemented
my @STYLE_ATTRS_NYI = qw (stroke-miterlimit stroke-dasharray
			  gradientUnits gradientTransform spreadMethod
			  clip-rule clip-path
			  name
			  ) ;   # what is the foolish name?
my %STYLE_ATTRS_NYI = map { $_ => 1 } @STYLE_ATTRS_NYI;

#### not yet implemented tags (to avoid many many error messages)
#### this list could be used to clearly distinguishe TAGS
#### not yet implemented or not implementable.
#### This list is curently not used! consider it as a piece of comment!
my @NO_YET_IMPLEMENTED_TAGS = qw ( midPointStop
				   filter feColorMatrix feComponentTransfer feFuncA
				   );

my $fileDir; ## in fact this could be a part of an url
             ## but we currently only get files in the some directories

my $backend; ## the backend used to produce/interpret  perl or tcl or whatever...

my $expat;
sub parsefile {
    my ($svgfile, $backendName, %args) = @_;

    # some init
    &InitVars;
    $fileDir = dirname($svgfile)."/";
    
    # the group where to create items, defaulted to 1
#    $current_group = defined $args{-group} ? $args{-group} : 1;

    # verbosity, defaulted to 0
    $verbose = defined $args{-verbose} ? $args{-verbose}: 0;

    # init of the prefix used to prefix tags. defaulted to the empty string
    $prefix =  defined $args{-prefix} ? $args{-prefix} : "";
    delete $args{-prefix}; # this option is not propagated to Backend

#    print "The prefix is '$prefix'\n";
    # should we treat XML namespace?
    my $namespace = defined $args{-namespace} ? $args{-namespace} : 0;
    delete $args{-namespace}; # this option is not propagated to Backend
    
    ## init of some global variables used by Conversions.pm
    &SVG::SVG2zinc::Conversions::InitConv(\&myWarn, \&current_line);

    
    my $filename;
    if ($filename = &findINC($backendName.".pm")) {
	# print " loading $filename\n";
	eval {require "$filename"};
    } elsif ($filename = &findINC("SVG/SVG2zinc/Backend",$backendName.".pm")) {
	# print " loading $filename\n";
	eval {require "$filename"};
	$backendName = "SVG::SVG2zinc::Backend::$backendName";
    } else {
	die "unable to find Backend $backendName in perl path @INC";
    }
    if ($@) {
	die "while loading Backend $backendName:\n$@\n";
    }

    $backend=$backendName->new(-in => $svgfile, %args);

    $current_group = $backend->_topgroup;
    $backend->fileHeader;
    my $parser = new XML::Parser(
				 Style => 'SVG2zinc',
				 Namespaces => $namespace,  # well this works for dia shape dtd!
				 Pkg => 'SVG::SVG2zinc',
				 ErrorContext => 3,
				 );
    $parser->setHandlers(Char    => \&Char,
			 Init    => \&Init,
			 Final   => \&Final,
			 XMLDecl => \&XMLDecl,
			 );
    my $svg=$parser->parsefile($svgfile);
    $backend->fileTail;
    &print_warning_for_not_implemented_attr;
} # end of parsefile

## as it seems that some svg files are using differencies between dtd 1.0 and 1.1
## we need to know which version of the dtd we are using (defaulted to 1.0)
my $dtdVersion;
sub XMLDecl {
    my ($parser, $Version, $Encoding, $Standalone) = @_;
#    $Standalone = '_undef' unless defined $Standalone;
#    print "XMLDecl: $parser, $Version, $Encoding, $Standalone\n";
    if (defined $Version) {
	$dtdVersion = $Version;
    } else {
	$dtdVersion = 1.0;
    }
} # end of XMLDecl



# the svg tags are translated in group items.
# If the SVG tag contains both width and height properties
# they will be reported in the generated group as tags :
# 'height=xxx' 'width=xxx'
sub svg {
    my ($parser, $elementname, %attrs) = @_;
    %attrs = &expandAttributes ($elementname,%attrs);
    my ($name,$natural) = &name ($elementname, $attrs{id});
    delete $attrs{xmlns}; # this attribute is mandatory, but useless for SVG2zinc

    my ($width,$height)=&sizesConvert( \%attrs , qw (width height)); #! this defines the Zinc size!
    # case when the width or height is defined in %
    # the % refers to the size of an including document
    undef $width if defined $attrs{width} and $attrs{width} =~ /%/ ;
    undef $height if defined $attrs{height} and $attrs{height}=~ /%/ ;
#    print "WIDTH,HEIGHT = $width    $height\n";
    my $widthHeightTags="";
    if (defined $width and defined $height) {
	$widthHeightTags = ", 'width=" . &float2int($width) . 
	    "', 'height=" . &float2int($height) . "'";
    }
    if (!@prev_contexts) { # we are in the very top svg group!
	$widthHeightTags .= ", 'svg_top'";
    }
    my $res = "->add('group',$current_group, -tags => [$name$widthHeightTags], -priority => 10";
    unshift @prev_contexts, \%current_context;
    my $prop;
    ($prop, %current_context) = &groupContext ($name, %attrs);
    $res .= $prop . ");";
    
    unshift @prev_groups, $current_group;
    $current_group = $name;
    
    foreach my $attr (keys %attrs) {
	if ($attr =~ /^xmlns:(.+)/ ) {
	    my $extensionPrefix = $1;
	    # this xlink extension is only partly implemented
	    # (ie. when the url refers an image file in the same directory than the SVG file)
	    next if ($extensionPrefix eq 'xlink');  
	    print "$extensionPrefix is not implemented\n";
	    $notImplementedExtensionPrefix{$extensionPrefix} = 0;
	}
    }
    
    &attrs_implemented ( 'svg', $name, [qw ( id width height viewBox preserveAspectRatio
					     xmlns),
					# the following attributes are not currently implementable
					qw ( enable-background overflow )
					], %attrs );
    &stackPort($name, $width,$height, $attrs{viewBox}, $attrs{preserveAspectRatio});
    &display ($res);
} # end of svg

my @portStack;
sub stackPort {
#    my ($name, $width,$height,$viewbox,$aspectRatio)=@_;
    unshift @portStack, [ @_ ];
}

## to treat the viewbox, preserveAspectRatio attributes
## of the svg, symbol, image, foreignObject... entities
sub viewPortTransforms {
    my $portRef = shift @portStack;
    my ($name, $width,$height,$viewbox,$aspectRatio)=@{$portRef};
    $viewbox = "" unless defined $viewbox;
    $aspectRatio = "" unless defined $aspectRatio;
    $width = "" unless defined $width;
    $height = "" unless defined $height;
#    print "In $name: width=$width height=$height viewbox=$viewbox aspectRatio=$aspectRatio\n";
    if ($viewbox and $width and $height ) {
	my $expr = "->adaptViewport($name, $width,$height, '$viewbox', '$aspectRatio');";
#	print "Expr = $expr\n";
	&display($expr);
#	if (!$aspectRatio or $aspectRatio eq "none") {
#	    my $translateX = $minx;
#	    my $translateY = $miny;
#	    my $scaleX= $width /  ($portWidth - $minx);
#	    my $scaleY= $height / ($portHeight - $miny);
#	    @transfs = ("->translate($name, $translateX, $translateY);",
#			"->scale($name, $scaleX, $scaleY);");
#           &display(@transfs);
    }
}


sub svg_ {
    my ($parser, $elementname) = @_;
    print "############ End of $elementname:\n" if $verbose;
    &viewPortTransforms;
    $current_group = shift @prev_groups;
    %current_context = %{shift @prev_contexts};
}

# just to avoid useless warning messages
sub desc {}
sub desc_ { }

# just to avoid useless warning messages
sub title {}
sub title_ { }

# just to avoid useless warning messages in svg tests suites
sub Paragraph {}
sub Paragraph_ { }

## return either the id of the object or a name of the form '__<elementtype>__<$counter>'
## the returned named includes single quotes!
## it also increments two counters:
##  - the itemCount used for naming any item
##  - the effectiveItemCount for counting graphic items only
##    This counter is used at the end of a defs to see if a group
##    must be saved, or if the group is just empty
sub name {
    my ($type, $id) = @_;
    print "############ In $type:\n" if $verbose;
    $itemCount++;
    $effectiveItemCount++ if (defined $id and
			      $type ne 'defs' and
			      $type ne 'switch' and
			      $type ne 'g' and
			      $type ne 'svg' and
			      $type !~ /Gradient/
			      );
    if (defined $id) {
	$id = &cleanName ($id);
	return ("'$id'", 1);
    } else {
	return ("'" . $prefix . "__$type"."__$itemCount'",0);
    }
} # end of name

sub g {
    my ($parser, $elementname, %attrs) = @_;
    %attrs = &expandAttributes ($elementname,%attrs);
    my ($name,$natural) = &name ($elementname, $attrs{id});
#    print "Group: $name\n";
    my $res = "->add('group',$current_group, -tags => [$name], -priority => 10";
    unshift @prev_groups, $current_group;
    $current_group = $name;
    unshift @prev_contexts, \%current_context;
    my $prop;
#    print "GROUP1 context: ", (%attrs),"\n";
    ($prop, %current_context) = &groupContext ($name, %attrs);
#    print "g1:",join(" ",( %current_context )), "\n";
    $res .= $prop . ");";
    &attrs_implemented ( 'g', $name, [qw ( id transform ) , @EVENT_ON_GRAPHICS ], %attrs ); ### les styles attrs sont à traiter à part!
    &display ($res,
	      &transform($name, $attrs{transform}));
    &treatGroupEvent ($name, %attrs);
} # end of g

## returns true if the parameter is an EVENT_ON_GRAPHICS (ie. applies only to group-like tags)
sub isGroupEvent {
    my ($attr) = @_;
    return $EVENT_ON_GRAPHICS{$attr} or 0;
}

## should bing callbacks to group, depending on events and scripts...
## not yet implemented
sub treatGroupEvent {
    my ($objname, %attr) = (@_);
    foreach my $event (@EVENT_ON_GRAPHICS) {
	my $value = $attr{$event};
	next unless defined $value;
#	print "## $objname HAS EVENT $event = $value\n";
	# XXX what should I do here?
    }
}

sub groupContext {
    my ($name, %attrs) = @_;
    my %childrenContext;
    my $prop = "";
    foreach my $attr (keys %attrs) {
	my $value = $attrs{$attr};
#	print "IN $name : $attr := $value\n";
	if (!defined $value) {
	    &myWarn ("!! Undefined value for attribute $attr in group $name !?");
	    next;
	} elsif (&isGroupEvent ($attr)) {
	    next;
	}
	$value = &removeComment($value);
	if ($attr eq 'opacity') { # attributes to apply directly to the group
	    $value = &convertOpacity ($value);
	    $prop = sprintf ", -alpha => %i", &float2int($value * 100); 
	} elsif ($attr eq 'id' or $attr eq 'transform') { # attributes treated before!
	    next;
	} elsif ($attr eq 'display' and $value eq 'none') {
	    ## beware: the visibility attribut is inheritated but can be modidied by a child
	    ## I put it in the %childrenContext and children will manage it
	    $prop .= ", -visible => 0, -sensitive => 0";
	    &myWarn ("!! The following group is not visible: $name  !?\n");
	} elsif (&isAnExtensionAttr($attr)) {
	    next;
	} elsif ($attr eq 'viewBox' or $attr eq 'preserveAspectRatio'  or $attr eq 'height' or $attr eq 'width') {
	    ### hack which works fine for managing the viewport!!
	} elsif (!defined $STYLE_ATTRS{$attr}) { # this attribute is not implemented!
	    if (defined $STYLE_ATTRS_NYI{$attr}) {
		&not_implemented_attr($attr);
	    } else {
		&myWarn ("!!! Unimplemented attribute '$attr' (='$value') in group $name\n");
	    }
	    next;
	} else { # all other attributes will be applied to children
	    $childrenContext{$attr} = $value;
	}
    }
    print "children context: ", join (", ", (%childrenContext)) , "\n" if $verbose;
    return ($prop, %childrenContext);
} # end of groupContext
	    

sub g_ {
    my ($parser, $elementname) = @_;
    print "############ End of $elementname:\n" if $verbose;
    $current_group = shift @prev_groups;
    %current_context = %{shift @prev_contexts};
}

## A switch is implemented as a group.
## BUG: In fact, we should select either the first if the tag is implemented
## or the secund sub-tag if not.
## In practice, the first sub-tag is not implemented in standard SVG, so we
## we forget it and take the second one.
## A problem will appear if the first tag is implemented, because, in this case
## we will instanciantes both the first and second
sub switch {
    my ($parser, $elementname, %attrs) = @_;
    %attrs = &expandAttributes ($elementname,%attrs);
    my ($name) = &name ($elementname, $attrs{id});
    $name =~ s/\'//g;
    $attrs{id} = $name;
    &g($parser, $elementname, %attrs);
} # end of switch

sub switch_ {
#    my ($parser, $elementname) = @_;
    &g_;
}

# a clipath is a not-visible groupe whose items define a clipping area
# usable with the clip-path attribute
# BUG: currently, the clipping is not implemented, but at least clipping
## items are put in a invisible sub-group and are not displayed
sub clipPath {
    my ($parser, $elementname, %attrs) = @_;
    %attrs = &expandAttributes ($elementname,%attrs);
    my ($name,$natural) = &name ($elementname, $attrs{id});
    print "In clippath $name NYI\n";
    ## should we verify that the clippath has an Id?
    my $res = "->add('group',$current_group, -tags => [$name, '__clipPath'], -priority => 10, -atomic => 1, -visible => 0";
    unshift @prev_groups, $current_group;
    $current_group = $name;
    unshift @prev_contexts, \%current_context;
    my $prop;
    ($prop, %current_context) = &groupContext ($name, %attrs);
    $res .= $prop . ");";
#    &attrs_implemented ( 'g', $name, [qw ( id transform ) , @EVENT_ON_GRAPHICS ], %attrs ); ### les styles attrs sont à traiter à part!
    &display ($res, &transform($name, $attrs{transform}));
#    &treatGroupEvent ($name, %attrs);
} # end of clippath

sub clipPath_ {
    my ($parser, $elementname) = @_;
    print "############ End of $elementname:\n" if $verbose;
    $current_group = shift @prev_groups;
    %current_context = %{shift @prev_contexts};
} # end of clippath_

# a symbol is a non-visible group which will be instancianted (cloned) 
# latter in a <use> tag
sub symbol {
    my ($parser, $elementname, %attrs) = @_;
    %attrs = &expandAttributes ($elementname,%attrs);
    my ($name,$natural) = &name ($elementname, $attrs{id});
    ## should we verify that the clippath has an Id?
    ## should we verify that <symbols> is defined inside a <defs> tag?
    my $res = "->add('group',$current_group, -tags => [$name], -priority => 10, -atomic => 1, -visible => 0";
    unshift @prev_groups, $current_group;
    $current_group = $name;
    unshift @prev_contexts, \%current_context;
    my $prop;
    ($prop, %current_context) = &groupContext ($name, %attrs);
    $res .= $prop . ");";
#    &attrs_implemented ( 'g', $name, [qw ( id transform ) , @EVENT_ON_GRAPHICS ], %attrs ); ### les styles attrs sont à traiter à part!
    &display ($res, &transform($name, $attrs{transform}));
#    &treatGroupEvent ($name, %attrs);
} # end of symbol

sub symbol_ {
    my ($parser, $elementname) = @_;
    print "############ End of $elementname:\n" if $verbose;
    $current_group = shift @prev_groups;
    %current_context = %{shift @prev_contexts};
} # end of symbol_

# this will clone and make visible either symbols or other items based on the Id refered by the xlink:href attribute
sub use {
    my ($parser, $elementname, %attrs) = @_;
    %attrs = &expandAttributes ($elementname,%attrs);
    my ($name,$natural) = &name ($elementname, $attrs{id});
#    my @attrs = %attrs; print "############ Start of $elementname: @attrs\n" if $verbose;
    my $ref = $attrs{'xlink:href'};
    if (!defined $ref) {
	&myWarn ("!! $elementname must have a xlink:href attribute\n");
	return;
    }
    $ref =~ s/\#//;
    my $cleanedId = &cleanName($ref); # to make the name zinc compliant
    my $res = "->clone('$cleanedId', -visible => 1, -tags => [$name, 'cloned_$cleanedId']";
    $res .= &analyze_style (\%attrs);
    $res .=");";
    my ($x,$y,$width,$height) = ($attrs{x},$attrs{y},$attrs{width},$attrs{height});
    my @transforms = "->chggroup($name, $current_group);";
    if (defined $x) {
	push @transforms, "->translate($name, $x,$y);";
    }
    &display ($res,@transforms);
}

{   ## start of defs block to share $res and other variables between many functions

    ## XXX: BUG this code DOES NOT allow recursive defs! (this is also probably a bug in the SVG file)
    my $defsCounter = 0;
    my $insideGradient = 0;  ## should never exceed 1!
    my $res; # the current gradient/object being defined
    my $defsId; # the group id containing items to be cloned
    # this group will be deleted later if it is empty
    
    my $effectiveItem;
    ## a <defs> will generate the creation of an invisible group in Tk::Zinc
    ## to be cloned latter in a <use> tag
    ## This group can be potentialy empty and in this cas it would be better
    ## not to create it, or at least delete it latter if it is empty
    sub defs {
	my ($parser, $elementname, %attrs) = @_;
	%attrs = &expandAttributes ($elementname,%attrs);
	($defsId) = &name ($elementname, $attrs{id});
	$defsId =~ s/\'//g;
	$attrs{id} = $defsId;
	&g($parser, $elementname, %attrs);
	&display("->itemconfigure('$defsId', -visible => 0);");
	$defsCounter++;
	$effectiveItem = $effectiveItemCount;
	print "############ $elementname: $defsId\n" if $verbose;
    }

sub defs_ {
    my ($parser, $elementname) = @_;
    $defsCounter++;
#    print "end of defs $defsId:", $effectiveItemCount , $effectiveItem, "\n";
    &g_;
    if ($effectiveItemCount == $effectiveItem) {
	&display ("->remove('$defsId');");
    }
}


######################################################################
### CSS : Cascading Style Sheet
######################################################################
{ ### css
    my @styles;
    my %classes;
    my %elementClasses;
    my %idClasses;
    my $in_css=0;
sub nextStyle {
    my $text = shift;
    push @styles,$text;
#    print "Style: $text\n";
}

## returns a reference to a hash-table defining pair of (attribut value) describing
## a CSS style associated to a Class
## returns undef if such class is not defined
sub getClass {
    my $class = shift;
    my $ref_styles = $classes{$class};
#    print "in getClass: $class  ",%classes, "\n";
#    my %styles = %{$ref_styles}; print "in getClass: $class ", (%styles), "\n";
    return ($ref_styles);
}

## returns a reference to a hash-table defining pair of (attribut value) describing
## a CSS style associated to an element type
## returns undef if such element type is not defined
sub getElementClass {
    my $element = shift;
    my $ref_styles = $elementClasses{$element};
#    my %styles = %{$ref_styles};
#    print "in getElementClass: $element ", (%styles), "\n";
    return ($ref_styles);
}

## returns a reference to a hash-table defining pair of (attribut value) describing
## a CSS style associated to an Id
## returns undef if such class is not defined
sub getIdClass {
    my $id = shift;
    my $ref_styles = $idClasses{$id};
#    my %styles = %{$ref_styles};
#    print "in getIdClass: $id ", (%styles), "\n";
    return ($ref_styles);
}

sub style {
    my ($parser, $elementname, %attrs) = @_;
    if ($attrs{type} eq "text/css") {
	$in_css=1;
    }
} # end of style

sub style_ {
    my ($parser, $elementname) = @_;
    my $str = "";
    foreach my $s (@styles) {
	$s = &removeComment($s);
	$s =~ s/\s(\s+)//g ; # removing blocks of blanks
	$str .= " " . $s;
    }
#    print "in style_: $str\n";
    while ($str) {
#	print "remaning str in style_: $str\n";
	if ($str =~ /^\s*\.(\S+)\s*\{\s*([^\}]*)\}\s*(.*)/ ) {
	    # class styling
	    my ($name,$value) = ($1,$2);
	    $str = $3;
#	    $value =~ s/\s+$//;
	    print "STYLE of class: '$name' => '$value'\n";
	    ## and now do something!
	    my %style = &expandStyle($value);
	    $classes{$1} = \%style;
	} elsif ( $str =~ /^\s*\#([^\{]+)\s*\{\s*([^\}]*)\}\s*(.*)/ ) {
	    my ($ids,$value) = ($1,$2);
	    $str = $3;
	    print "STYLE of ids: '$ids' => '$value'\n";
	    ## and now do something!
	} elsif ( $str =~ /^\s*\[([^\{]+)\]\s*\{\s*([^\}]*)\}\s*(.*)/ ) {
	    my ($attr_val,$value) = ($1,$2);
	    $str = $3;
	    print "STYLE of attr_values: '$attr_val' => '$value'\n";
	    ## and now do something!
	} elsif  ( $str =~ /^\s*\@font-face\s*\{\s*[^\}]*\}\s*(.*)/ ) {
	    print "STYLE of font-face", substr($str, 0, 100),"....\n";
	    $str = $1;
	} elsif  ( $str =~ /^\s*([^\s\{]+)\s*\{\s*([^\}]*)\}\s*(.*)/ ) {
	    my ($name,$value) = ($1,$2);
	    $str = $3;
	    print "STYLE of tags: '$name' => '$value'\n";
	    ## and now do something... NYI
	} else {
	    &myWarn ("unknown style : $str\nskipping this style");
	    return;
	}
    }
    $in_css=0;
    @styles=();
} # end of style_

} ### end of css

######################################################################
### gradients
######################################################################

my $gname;
my @stops;
my @inheritedStops;
my $angle;
my $center;

sub radialGradient {
    my ($parser, $elementname, %attrs) = @_;
    %attrs = &expandAttributes ($elementname,%attrs);
    my ($name,$natural) = &name ($elementname, $attrs{id});
    &myWarn ("!! $elementname must have an id\n") unless $natural;
    $gname = substr ($name,1,-1); # remove quote (') at the very beginning and end of $name
    $insideGradient ++;
    &myWarn ("Gradient '$gname' definition inside a previous gradient definition. This is bug in svg source\n")
	    unless $insideGradient == 1;
    $res="->gname(";
    @stops = ();
    @inheritedStops = ();
    if (defined $attrs{'xlink:href'}) {
	my $unused;
	my $link = delete $attrs{'xlink:href'};
	if ($link =~ /^\#(.+)$/) {
	    $link = $1;
	    ($unused, @inheritedStops) = &extractGradientTypeAndStops ($link);
	} else {
	    # BUG??: we only treat internal links like #gradientName
	    carp "bad link towards a gradient: $link";
	}
    }
    my ($fx,$fy,$cx,$cy, $r) = &sizesConvert( \%attrs , qw (fx fy cx cy r));
    # BUG: a serious limitation is that TkZinc (3.2.6i) does not support
    #      the cx, cy and r parameters

    if (defined $cx and $cx == $fx) { delete $attrs{cx}; }  # to avoid needless warning of &attrs_implemented
    if (defined $cy and $cy == $fy) { delete $attrs{cy}; }  # to avoid needless warning of &attrs_implemented
    &attrs_implemented ( 'radialGradient', $name, [qw ( id fx fy)], %attrs );

    $fx = &float2int(($fx -0.5) * 100);
    $fy = &float2int(($fy -0.5) * 100);
    $center = "$fx $fy";
}

sub radialGradient_ {
    $insideGradient --;
    if (!@stops) {
	if (@inheritedStops) {
	    @stops = @inheritedStops;
	} else {
	    carp ("Bad gradient def: nor stops, neither xlink;href");
	}
    }
    my $gradientDef = "=radial $center | " . join (" | ", @stops);
    $res .= "\"" . $gradientDef . "\", \"$gname\");";  ### BUG: limits: x y!
    # si il faut appliquer une transparence sur un gradient on est très embêté!
    &defineNamedGradient($gname, $gradientDef) ;
#    print "RADIAL='$res'\n";
    &display($res) ;
    @stops = ();
}

sub linearGradient {
    my ($parser, $elementname, %attrs) = @_;
    %attrs = &expandAttributes ($elementname,%attrs);
    my ($name,$natural) = &name ($elementname, $attrs{id});
    &myWarn ("!! $elementname must have an id\n") unless $natural;
    $gname = substr ($name,1,-1); # remove quote (') at the very beginning and end of $name
    $insideGradient ++;
    &myWarn ("Gradient '$gname' definition inside a previous gradient Definition. This will bug\n")
	    unless $insideGradient == 1;
    $res="->gname("; 
    @stops = ();
    @inheritedStops = ();
    if (defined $attrs{'xlink:href'}) {
	my $unused;
	my $link = delete $attrs{'xlink:href'};
	if ($link =~ /^\#(.+)$/) {
	    $link = $1;
	    ($unused, @inheritedStops) = &extractGradientTypeAndStops ($link);
	} else {
	    # BUG??: we only treat internal links like #gradientName
	    carp "bad link towards a gradient: $link";
	}
    }
    &attrs_implemented ( 'linearGradient', $name, [qw ( x1 x2 y1 y2 id )], %attrs );
    my ($x1,$x2,$y1,$y2) = &sizesConvert( \%attrs , qw (x1 x2 y1 y2));
    if ( ($y2 - $y1) or ($x2 - $x1) ) {
	my $atan = - rad2deg (atan2 ($y2-$y1,$x2-$x1));
	$angle = &float2int($atan);
    } else {
	$angle = 0;
    };
#    print "ANGLE = $angle\n";
}

sub linearGradient_ {
    $insideGradient --;
    if (!@stops) {
	if (@inheritedStops) {
	    @stops = @inheritedStops;
	} else {
	    carp ("Bad gradient def: nor stops, neither xlink;href");
	}
    }
    my $gradientDef = "=axial $angle | " . join (" | ", @stops);
    $res .=  "\"" . $gradientDef . "\", \"$gname\");";
    # si il faut appliquer une transparence sur un gradient on est très embêté!
    &defineNamedGradient($gname, $gradientDef) ;
    &display($res) ;
}

sub stop {
    my ($parser, $elementname, %attrs) = @_;
    %attrs = &expandAttributes ($elementname,%attrs);
#    my ($name,$natural) = &name ($elementname, $attrs{id});  # no name is needed!
    &myWarn ("$elementname should be defined inside <linearGradient> or <radialGradiant>\n") unless $insideGradient;

    my $style = delete $attrs{'style'};
    if (defined $style) {
	my %keyvalues = &expandStyle($style);
	%attrs = (%attrs , %keyvalues);
    }
    my $offset = $attrs{'offset'};
    my $color = $attrs{'stop-color'};
    if (!defined $color) {
	&myWarn ("!! Undefined  stop-color in a <stop>\n");
    } elsif (!defined $offset) {
	&myWarn ("!! Undefined offset in a <stop>\n");
    } else {
	if ($offset =~ /([\.\d]+)%/){  
	    $offset = &float2int($1);
#	} elsif ($offset =~ /^([.\d]+)$/) {
#	    $offset = &float2int($1*100);
	} else {
	    $offset = &float2int($offset*100);
	}
	my $color=&colorConvert($color);
	if (defined (my $stopOpacity = $attrs{'stop-opacity'})) {
	    $stopOpacity = &float2int($stopOpacity*100);
	    push @stops, "$color$colorSep$stopOpacity $offset";
	} else {
	    push @stops, "$color $offset";
	}
    }
} # end of stop

} # end of gradient closure


my %convertFormat = (
    'jpg' => 'jpeg',
    'jpeg' => 'jpeg',
    'png' => 'png',
);

sub image {
    my ($parser, $elementname, %attrs) = @_;
    %attrs = &expandAttributes ($elementname,%attrs);
    my ($name,$natural) = &name ($elementname, $attrs{id});
#    &myWarn ("!! $elementname must have an id\n") unless $natural;

    my $group = $current_group;
    my @RES;
    if (my $opacity = $attrs{'opacity'}) {
        # creating an intermediate group for managing the transparency
	# BUG: we could used the attribute -color := white:$opacity
	$opacity = &convertOpacity ($opacity);
	if ($opacity != 1) {
	    ## on crée un groupe pour gérer la transparence
	    my $opacity = &float2int(100 * $opacity);
	    my $newgroup = substr ($name,0,-1) . "transparency'";
	    push @RES , "->add('group', $current_group, -alpha => $opacity, -tags => [ $newgroup ], -priority => 10);\n";
	    $group = $newgroup;
	}
    }
    my $res = "";
    my $ref = "";
    if ($ref = $attrs{'xlink:href'}) {

	if ($ref =~ /^data:image\/(\w+);base64,(.+)/) {

	    # this code has been provided by A. Lemort from Intuilab
	    # for uuencoded inline image
	    my $format = $1;
	    my $data = $2;
#	    print ("data:image: '", substr($ref,0,30), "....' format=$format\n");
	    $ref = "data:image/$format;base64"; # $ref is used later in a tag of the icon
	    $format = $convertFormat{lc($format)};
	    $res .= "->add('icon',$group, -image => ->Photo(-data => '$data', -format => '$format')";
	} elsif ($ref =~ /^data:;base64,(.+)/) {
	    ## the following piece of code works more or less ?!
 	    ## BUG: there is a pb with scaling (ex: data-svg/vero_data/propal_crea_boutons.svg) 
	    my $data = $1;
#	    print ("data:; '", substr($ref,0,30), "....' NO format!\n");
	    $ref = "data:;base64"; # $ref is used later in a tag of the icon
	    $res .= "->add('icon',$group, -image => ->Photo(-data => '$data')";
	} else {
	    # It's a file
	    # print "Including image : $fileDir$ref\n";
	    if (open REF, "$fileDir$ref") {
		close REF;
#		print "group='$group' ref='$ref' filedir=$fileDir\n";
		$res .= "->add('icon',$group, -image => ->Photo('$ref', -file => '$fileDir$ref')";
	    } else {
		&myWarn ("When parsing the image '$name': no such file: '" . substr ("$fileDir$ref", 0,50) . "'\n") ;
		return;
	    }
	}
    } else {
	&myWarn ("Unable to parse the image '$name'") ;
	return;
    }

    $res .= ", -tags => [$name, '$ref'], -composescale => 1, -composerotation => 1, -priority => 10);";
    push @RES, $res ;

    my ($x, $y, $width, $height) = &sizesConvert ( \%attrs , qw (x y width height));
    if ($width == 0 or $height == 0) {
	&myWarn ("Skipping a 0 sized image: '$name' size is $width x $height\n");
    } elsif ($width < 0 or $height < 0) {
	&myWarn ("Error in the size of the image '$name' : $width x $height\n");
    } else {
	push @RES, "->adaptViewport($name, $width,$height);";
    }
    if ($x or $y) {
	push @RES, "->translate($name, $x,$y);";
    }
    
    &attrs_implemented ( 'image', $name, [qw ( x y width height id )], %attrs );
    &display (@RES,
	      &transform($name, $attrs{transform}) );
} # end of image


sub line {
    my ($parser, $elementname, %attrs) = @_;
    %attrs = &expandAttributes ($elementname,%attrs);
    my ($name,$natural) = &name ($elementname, $attrs{id});
    my $res = "->add('curve',$current_group,[$attrs{x1},$attrs{y1},$attrs{x2},$attrs{y2}], -priority => 10";
    $res .= ", -tags => ['line'";
    $res .= ", $name" if ($natural or $attrs{transform});
    $res .= "]";
    $res .= &analyze_style (\%attrs);
    $res .=");";
    &attrs_implemented ( 'line', $name, [qw (x1 y1 x2 y2 style id transform )], %attrs );
    &display ($res,
	      &transform($name, $attrs{transform}) );
} # end of line


sub Char {
    my ($expat, $text) = @_;
    return if !defined $text;
    my $type = ($expat->context)[-1];
    return if !defined $type;
    chomp $text;
    return if (!$text && ($text ne "0")); # empty text!
    # $text =~ s/([\x80-\xff])/sprintf "#x%X;", ord $1/eg;
    # $text =~ s/([\t\n])/sprintf "#%d;", ord $1/eg;
#    print "$type: $text\n";
    if ($type eq 'tspan' or $type eq 'text') {
#	print "[$text]\n";
	&nextText ($text);
    } elsif ($type eq 'style') {
	&nextStyle ($text);
    }
} # end of char


## this lexical block allows &text, &nextTetx, &tspan, and &text_ to share common variables
{
    my $res;
    my @transforms;
    my @texts;
    sub text {
	my ($parser, $elementname, %attrs) = @_;
	%attrs = &expandAttributes ($elementname,%attrs);
	my ($name,$natural) = &name ($elementname, $attrs{id});
	my ($x,$y)=&sizesConvert( \%attrs , qw (x y));
	$res = "->add('text',$current_group, -position => [$x,$y], -tags => ['text'";
	$res .= ", $name" if ($natural or $attrs{transform});
	$res .= "], -anchor => 'sw'"; ## XXX est-ce bien la bonne ancre?
	## XXX BUG? je ne suis pas sur que la ligne suivante soit indispensable?!
	$res .= &analyze_text_style (\%attrs);
	@texts = ();
	@transforms = reverse &transform($name, $attrs{transform});
	&attrs_implemented ( 'text', $name, [qw (x y id transform text-anchor font-family font-size)], %attrs );
    } # end of text

sub nextText {
    my $txt = shift;
    push @texts,$txt;
}


## BUG: <tspan> tags can be used to modiofy many graphics attributs of the part of the text
## such as colors, fonte, size and position...
## this is currently hard to implement as in Tk::Zinc a text item can only have one color, one size...
sub tspan {
    my ($expat, $elementname, %attrs) = @_;
#    my @attrs = %attrs; print "TSPAN: $elementname @attrs\n";
    $res .= &analyze_text_style (\%attrs);
} # end of tspan

sub text_ {
    my ($parser, $elementname, @rest) = @_;
#    my $text = join ('+++',@texts); print "TEXT_ : $text\n";
    for (my $i=0 ; $i <= $#texts ; $i++) {
	$texts[$i] =~ s/\'/\\'/g ;  #'
    }
    my $theText = join ('', @texts);
    $res .= ", -text => '$theText', -priority => 10);";
    &display ($res, @transforms);
} # end of test_

} ## end of text lexical block

sub polyline {
    my ($parser, $elementname, %attrs) = @_;
    %attrs = &expandAttributes ($elementname,%attrs);
    my ($name,$natural) = &name ($elementname, $attrs{id});
    my $res = "->add('curve',$current_group,[" . &points(\%attrs);
    $res .= "], -tags => ['polyline'";
    $res .= ", $name" if ($natural or $attrs{transform});
    $res .= "], -priority => 10";
    $res .= &analyze_style (\%attrs);
    $res .=");";
    &attrs_implemented ( 'polyline', $name, [qw (points style transform id )], %attrs );
    &display ($res,
	      &transform($name, $attrs{transform}) );
}

sub rect {
    my ($parser, $elementname, %attrs) = @_;
    %attrs = &expandAttributes ($elementname,%attrs);
    my ($name,$natural) = &name ($elementname, $attrs{id});
    my ($x,$y,$width,$height)=&sizesConvert( \%attrs , qw (x y width height));
    my $res = "->add('rectangle',$current_group,[$x,$y,"
	.($x+$width).",".($y+$height)."], -tags => ['rect'";
    $res .= ", $name" if ($natural or $attrs{transform});
    $res .= "], -priority => 10";
    # by default, rectangles are filled (cf example svg/use02_p87.svg
    # from svg specifs). The value is set here, and can be overidden later
    # in the  &analyze_style
    $res .= ", -filled => 1" unless defined $attrs{fill} and $attrs{fill} eq 'none';
    delete $attrs{'stroke-linejoin'}; ## BUG: due to TkZinc limitation: no -joinstyle for rectangle 
    $res .= &analyze_style (\%attrs);
    $res .=");";
    &attrs_implemented ( 'rect', $name, [qw (id x y width height style transform )], %attrs );
    &display ($res,
	      &transform($name, $attrs{transform}) );
}


sub ellipse {
    my ($parser, $elementname, %attrs) = @_;
    %attrs = &expandAttributes ($elementname,%attrs);
    my ($name,$natural) = &name ($elementname, $attrs{id});
    my ($cx,$cy,$rx,$ry)=&sizesConvert( \%attrs , qw (cx cy rx ry));
    my $res = "->add('arc',$current_group,[". ($cx-$rx) . ", ". ($cy-$ry) .
	", " . ($cx+$rx) . ", ". ($cy+$ry) . "], -tags => ['ellipse',";
    $res .= ", $name" if ($natural or $attrs{transform});
    $res .= "], -priority => 10";
    # by default, ellipses are filled
    # from svg specifs). The value is set here, and can be overidden later
    # in the  &analyze_style
    $res .= ", -filled => 1" unless defined $attrs{fill} and $attrs{fill} eq 'none';
    delete $attrs{'stroke-linejoin'}; ## BUG: due to TkZinc limitation: no -joinstyle for arc 
    $res .= &analyze_style (\%attrs);
    $res .=");";
    &attrs_implemented ( 'ellipse', $name, [qw (cx cy rx ry style transform id )], %attrs );
    &display ($res,
	      &transform($name, $attrs{transform}) );
}

sub circle {
    my ($parser, $elementname, %attrs) = @_;
    %attrs = &expandAttributes ($elementname,%attrs);
    my ($name,$natural) = &name ($elementname, $attrs{id});
    my ($cx,$cy,$r)=&sizesConvert( \%attrs , qw (cx cy r));
    my $res = "->add('arc',$current_group,[". ($cx-$r) . ", ". ($cy-$r) .
	", " . ($cx+$r) . ", ". ($cy+$r) . "], -tags => ['circle'";
    $res .= ", $name" if ($natural or $attrs{transform});
    $res .= "], -priority => 10";
    # by default, circles are filled
    # from svg specifs). The value is set here, and can be overidden later
    # in the  &analyze_style
    $res .= ", -filled => 1" unless defined $attrs{fill} and $attrs{fill} eq 'none';
    $res .= &analyze_style (\%attrs);
    $res .=");";
    delete $attrs{'stroke-linejoin'}; ## BUG: due to TkZinc limitation: no -joinstyle for arc 
    &attrs_implemented ( 'circle', $name, [qw ( cx cy r transform id )], %attrs );
    &display ($res,
	      &transform($name, $attrs{transform}) );
}


sub polygon {
    my ($parser, $elementname, %attrs) = @_;
    %attrs = &expandAttributes ($elementname,%attrs);
    my ($name,$natural) = &name ($elementname, $attrs{id});
    my $res = "->add('curve',$current_group,[" . &points(\%attrs);
    $res .= "], -closed => 1, -tags => ['polygon'";
    $res .= ", $name" if ($natural or $attrs{transform});
    $res .= "], -priority => 10";
    # by default, polygones are filled
    # from svg specifs). The value is set here, and can be overidden later
    # in the  &analyze_style
    $res .= ", -filled => 1" unless defined $attrs{fill} and $attrs{fill} eq 'none';
    $res .= &analyze_style (\%attrs);
    $res .= ");";
    &attrs_implemented ( 'polygone', $name, [qw ( points style transform id )], %attrs );
    &display ($res,
	      &transform($name, $attrs{transform}) );
}


sub path {
    my ($parser, $elementname, %attrs) = @_;
    %attrs = &expandAttributes ($elementname,%attrs);
#    my @attrs=%attrs; print "PATH attr=@attrs\n";
    my ($name,$natural) = &name ($elementname, $attrs{id});
    my $res = "->add('curve',$current_group,[";
    my ($closed, @listOfListpoints) = &pathPoints(\%attrs);
    my $refPoints = shift @listOfListpoints;
    $res .= join (", ", @{$refPoints});
    $res .= "], -tags => [$name], -priority => 10";
    # by default, paths are filled (cf exemple other-samples/logo_intuilab_illustrator.svg)
    # The value is set here, and can be overidden later
    # in the  &analyze_style
    $res .= ", -filled => 1" unless defined $attrs{fill} and $attrs{fill} eq 'none';
    if ( defined $attrs{'fill-rule'} ) {
	$res .= ", -fillrule => 'nonzero'" unless $attrs{'fill-rule'} eq 'evenodd';
	delete $attrs{'fill-rule'};
    }
    $res .= ", -closed => $closed";
    $res .= &analyze_style (\%attrs);
    $res .= ");";
    # and process other contours
    my @contours = ();
    foreach my $refPoints (@listOfListpoints) {
	my @points = @{$refPoints};
#	print "AN OTHER CONTOUR:  @points\n";
	my $contour = "->contour($name, 'add', 0, [";
	$contour .= join (", ", @points);
	$contour .= "]);";
	push @contours, $contour;
    }
    &attrs_implemented ( 'path', $name, [qw ( d style stroke-linejoin stroke-linecap transform id)], %attrs );
    &display ($res, @contours,
	      &transform($name, $attrs{transform}) );
} # end of path




sub expandAttributes {
    my ($elementName, %rawAttrs) = @_;
    my (%styleKeyValues, %classKeyValues, %elementKeyValues, %idKeyValues);
    my $style = delete $rawAttrs{'style'};
    if (defined $style) {
	%styleKeyValues = &expandStyle($style);
    }
    my $class = delete $rawAttrs{'class'};
    if (defined $class) {  ## for the css
	my $ref_styles = &getClass($class);
	if (defined $ref_styles) {
	    %classKeyValues = %{$ref_styles};
	} else {
	    &myWarn ("class attribute refers an illegal style: '$class'\n");
	}
    }
    my $ref_styles = &getElementClass($elementName);
    if (defined $ref_styles) {
	%elementKeyValues = %{$ref_styles};
    }
    my $id = $rawAttrs{id};
    if (defined $id) {
	my $ref_style = &getIdClass($id);
	if (defined $ref_style) {
	    %idKeyValues = %{$ref_styles};
	}
    } 
   return (%rawAttrs, %elementKeyValues, %classKeyValues, %styleKeyValues, %idKeyValues); ## the order is important!
}

### CM 19/1/03: This function could be really simplified (CM 09/09/3 why??? I do not remember!)
## analyze styles attached to an item (non text item) and on any of its groups
sub analyze_style {
    my ($ref_attr) = @_;
    my %ref_attr = %{$ref_attr};
    my %attrs = ( %current_context , %ref_attr );
#    print "analyze_style:",join(" ",( %attrs )), "\n";
    my %directkeyvalues;
    foreach my $attr (@STYLE_ATTRS) {
	my $value = $attrs{$attr};
	if (defined $value) {  
#	    print ("$attr := ", $value, "\n");
	    $directkeyvalues{$attr} = &removeComment($value);
	}
    }
    return &analyze_style_hash (\%directkeyvalues);
} # end of analyze_style;


## analyze styles attached to a text item  and on any of its groups
sub analyze_text_style {
    my ($ref_attr) = @_;
    my %attrs = ( %current_context , %{$ref_attr} );
    my $res = "";
    my $style = delete $attrs{'style'};
    if (defined $style) {
	my %keyvalues = &expandStyle($style);
	$res = &analyze_text_style_hash (\%keyvalues);
    }
    my %directkeyvalues;
    foreach my $attr (@STYLE_ATTRS) {
	my $value = $attrs{$attr};
	if (defined $value) {  
#	    print ("Analyzetext_style: $attr := ", $value, "\n");
	    $directkeyvalues{$attr} = &removeComment($value);
	}
    }
    $res .= &analyze_text_style_hash (\%directkeyvalues);
    return $res;
} # end of analyze_text_style;


## expanses the attribute = "prop:val;prop2:val2" in a hashtable like {prop => val, prop2 => val2, ...}
## and return this hash (BUG: may be it should return a reference!)
sub expandStyle {
    my ($style) = @_;
    return () unless defined $style;
    my %keyvalues;
    $style = &removeComment ($style);
    foreach my $keyvalue ( split ( /\s*;\s*/ , $style) ) {
	my ($key, $value) = $keyvalue =~ /(.*)\s*:\s*(.*)/ ;
#	print "Style: key = $key   value=$value\n";
	$keyvalues{$key} = $value;
    }
    return %keyvalues;
} # end of expandStyle


## Analyze attributes contained in the hashtable given as ref
## This hashtable {attribut =>value...} must contain all
## attributes to analyze
## returns a string containing the TkZinc attributes
sub analyze_style_hash {
    my ($ref_keyvalues) = @_;
    my %keyvalues = %{$ref_keyvalues};
    my $res = "";
    my $opacity = &convertOpacity(delete $keyvalues{'opacity'});

    ## we are treating now the stroke color and its transparency
    my $stroke = delete $keyvalues{'stroke'};
    my $strokeOpacity = delete $keyvalues{'stroke-opacity'};
    $strokeOpacity = 1 if !defined $strokeOpacity;
    $strokeOpacity = &float2int( &convertOpacity($strokeOpacity)*$opacity*100);
    if (defined $stroke) {
	my $color = &colorConvert($stroke);
#	print "stroke=$stroke  <=> '$color'\n";
	if ($color eq "none") {
	    $res .= ", -linewidth => 0";
	    delete $keyvalues{'stroke-width'};
	} elsif ( $strokeOpacity != 100 ) {
	    if ( &existsGradient($color) ) {
		# so, apply a transparency to a Tk::Zinc named gradient
		my $newColor = &addTransparencyToGradient($color,$strokeOpacity);
		$res .= ", -linecolor => \"$newColor\", -filled => 1";
	    } else {
		$res .= ", -linecolor => \"$color$colorSep$strokeOpacity\"";
	    }
	} else {
	    $res .= ", -linecolor => \"$color\"";
	}
    } elsif ( $strokeOpacity != 1 ) { # no stroke color, but opacity
	## what should I do?!
    }
	    
    ## we are treating now the fill color and its transparency
    my $fill = delete $keyvalues{'fill'};
    my $fillOpacity = delete $keyvalues{'fill-opacity'};
    $fillOpacity = 1 if !defined $fillOpacity;
    $fillOpacity = &float2int( &convertOpacity($fillOpacity)*$opacity*100);
    delete $keyvalues{'fill-opacity'};
    if (defined $fill) {  
	my $color = &colorConvert($fill);
	if ($color eq "none") {
	    $res .= ", -filled => 0";
	    delete $keyvalues{'fill-opacity'};
	} elsif ( $fillOpacity != 100 ) {
#	    print "fillOpacity=$fillOpacity\n";
	    if ( &existsGradient($color) ) {
		# so, apply a transparency to a Tk::Zinc named gradient
		my $newColor = &addTransparencyToGradient($color,$fillOpacity);
		$res .= ", -fillcolor => \"$newColor\", -filled => 1";
		## we must define the contour color, else it will be defaulted to black
		$res .= ", -linecolor => \"$newColor\"," unless defined $stroke;
	    } else {
		$res .= ", -fillcolor => \"$color$colorSep$fillOpacity\", -filled => 1";
		## we must define the contour color, else it will be defaulted to black
		$res .= ", -linecolor => \"$color$colorSep$fillOpacity\"," unless defined $stroke;
	    }
	} else {
	    $res .= ", -fillcolor => \"$color\", -filled =>1";
	    ## we must define the contour color, else it will be defaulted to black
	    $res .= ", -linecolor => \"$color\"" unless defined $stroke;
	}
    }

    # all other attributes now
    foreach my $key (sort keys %keyvalues) {
	my $value = $keyvalues{$key};
	next if (!defined $value);
#	print "KEY=$key VALUE=$value\n";
	if ($key eq 'stroke-width') {
	    if ( defined $keyvalues{stroke} and $keyvalues{stroke} eq 'none' ) {
		delete $keyvalues{stroke};
		next;
	    }
	    $value = &sizeConvert($value);
	    if ($value == 0 and $dtdVersion eq "1.0") {
		$value = 0.1;      # BUG? a widht of 0 is the smallest possible width in SVG 1.0 [true or false?]
	    }
	    $res .= ", -linewidth => $value";
	} elsif ($key eq 'display') {
	    if ($value eq 'none') {
		$res .= ", -visible => 0, -sensitive => 0";
	    }
	    ## We do not treat the other possible values for display as defined in CSS2?!
	} elsif ($key eq 'visibility') {
            ## BUG? if a "not-visible" <g> group contains a visible graphic element
	    ## this element WILL NOT be visible in TkZinc , but should be visible in SVG!!
	    ## Cf specif svg p. 284
	    if ($value eq 'hidden' or $value eq 'collapse') {
		$res .= ", -visible => 0";
	    }
	    ## We do not treat the other possible values for display as defined in CSS2?!
	} elsif ($key eq 'stroke-linecap') {
	    if ($value eq 'butt' or $value eq 'round') {
		$res .= ", -capstyle => \"$value\"";
	    } elsif ($value eq 'square') {
		$res .= ", -capstyle => \"projecting\"";
	    } else {
		&myWarn ("!! bad value for $key style : $value\n");
	    }
	} elsif ($key eq 'stroke-linejoin') {
	    ($value) = $value =~ /(\w+)/ ;  ## pour enlever d'eventuel blancs 
	    $res .= ", -joinstyle => \"$value\"";
	} elsif ($key eq 'fill-rule') {
	    ### this attributes is for shape only and is analyzed in &path
	} elsif ($key eq 'font-size') {
	    ### this attributes is for text only and is analyzed in &analyze_text_style_hash
	} else {
	    &myWarn ("Unknown Style (in analyze_style_hash): $key (value is $value)\n") if $warn;
	}
    }
    return $res;
} # end of analyze_style_hash


## We do not treat yet relative size of text e.g. : font-size = %120
sub analyze_text_style_hash {
    my ($ref_keyvalues) = @_;
    my %keyvalues = %{$ref_keyvalues};
#    print "analyze_text_style_hash: ", %keyvalues,"\n";
    my $res = "";
    my $opacity = &convertOpacity($keyvalues{opacity});
    delete $keyvalues{'opacity'};

    my $fontFamily="";
    my $fontSize ="";
    my $fontWeight ="";
    foreach my $key (sort keys %keyvalues) {
	my $value = $keyvalues{$key};
#	print "$key  ==>>  $value\n";
	next if (!defined $value);  # in this case, the SVG code is invalide (TBC)
	if ($key eq 'text-anchor') {
	    if ($value eq 'start') {
		$res .= ", -alignment => 'left'";
	    } elsif ($value eq 'end') {
		$res .= ", -alignment => 'right'";
	    } elsif ($value eq 'middle') {
		$res .= ", -alignment => 'center'"}
	} elsif ($key eq 'display') {
	    if ($value eq 'none') {
		$res .= ", -visible => 0, -sensitive => 0";
	    }
	    ## We do not treat the other possible values for display as defined in CSS2?!
	} elsif ($key eq 'visibility') {
            ## BUG? if a "not-visible" <g> group contains a visible graphic element
	    ## this element WILL NOT be visible in TkZinc , but should be visible in SVG!!
	    ## Cf specif svg p. 284
	    if ($value eq 'hidden' or $value eq 'collapse') {
		$res .= ", -visible => 0";
	    }
	    ## We do not treat the other possible values for display as defined in CSS2?!
	} elsif ($key eq 'font-family') {
	    $value =~ s/\'//g;  # on removing quotes around the fonte name
	    $fontFamily = $value;
#	    print "font-family  ==>>  $fontFamily\n";
	} elsif ($key eq 'font-size') {
	    $fontSize = $value;
	} elsif ($key eq 'font-weight') {
	    $fontWeight = $value;
#	    print "font-weight  ==>>  $fontWeight\n";
	} elsif ($key eq 'fill') {
	    my $fillOpacity;
	    my $color = &colorConvert($value);
	    if ($color eq 'none') {
		# $res .= ", -filled => 0"; # this is the default value in Tk::Zinc
	    } elsif ( ($fillOpacity = $keyvalues{'fill-opacity'} or $opacity != 1) ) {
		$fillOpacity = &convertOpacity($fillOpacity) * $opacity;
		delete $keyvalues{'fill-opacity'};
		if ( &existsGradient($color) ) {
		    #  so, apply a transparency to a Tk::Zinc named gradient
		    my $newColor = &addTransparencyToGradient($color,$fillOpacity);
		    $res .= ", -color => \"$newColor\"";
		} else {
		    $res .= ", -color => \"$color$colorSep$fillOpacity\"";
		}
	    } else {
		$res .= ", -color => \"$color\"";
	    }
	} else {
	    &myWarn ("Unknown Style of text: $key (value is $value)\n") if $warn;
	}
    }
    if ($fontFamily or $fontSize or $fontWeight) {
        ## to be extended to all other fonts definition parameters
	## NB: fontWeight is not used yet!
	my ($fontKey,$code) = &createNamedFont ($fontFamily, $fontSize, "");
	&display($code) if $code;
	$res .= ", -font => \"$fontKey\"";
    }
    return $res;
} # end of analyze_text_style_hash




## print warnings for all used attributes unkonwn or not implemented
sub attrs_implemented {
    my ($type, $name, $ref_attrs_implemented, %attrs) = @_;
    my %attrs_implemented;
    foreach my $attr (@{$ref_attrs_implemented}) {
	$attrs_implemented{$attr}=1;
    }
    my %expandStyle = &expandStyle ($attrs{style});
    my %attributes = ( %expandStyle, %attrs);
    foreach my $attr ( keys %attributes ) {
#	print "attr: $attr  $attributes{$attr}\n";
	if (!&isAnExtensionAttr($attr) and
	    !defined $STYLE_ATTRS{$attr} and
	    !defined $attrs_implemented{$attr}) {
	    if (defined $STYLE_ATTRS_NYI{$attr}) {
		&not_implemented_attr($attr);
	    } else {
		&myWarn ("!!! Unimplemented attribute '$attr' (='$attributes{$attr}') in '$type' $name\n");
	    }
	}
    }
} # end of attrs_implemented

# These hashes contain the number of usage of not implemented attributes and
# the lines on svg source files where a not implemented attributes is used
# so that they can be displayed by the sub &print_warning_for_not_implemented_attr
my %not_implemented_attr;
my %not_implemented_attr_lines;
sub not_implemented_attr {
    my ($attr) = @_;
    $not_implemented_attr{$attr}++;
    if (defined $not_implemented_attr_lines{$attr}) {
	push @{$not_implemented_attr_lines{$attr}},&current_line;
    } else {
	$not_implemented_attr_lines{$attr} = [&current_line];
    }
}

sub print_warning_for_not_implemented_attr {
    foreach my $k (sort keys %not_implemented_attr) {
	print "not implemented/implementable attribute '$k' was used $not_implemented_attr{$k} times in lines ";
	my @lines;
	if ($not_implemented_attr{$k} > 20) {
	    @lines = @{$not_implemented_attr_lines{$k}}[0..19];
	    print join (", ",@lines) ,"...\n";
	} else {
	    @lines = @{$not_implemented_attr_lines{$k}};
	    print join (", ",@lines) ,"...\n";
	}
    }
}


# print a warning for the first use of an attribute of a non-implemented extension to SVG
# return  :
# - true if the attribute belong to an extension of SVG
# - false if its supposed to be a standard SVG attribute (or a non-existing attribute)
sub isAnExtensionAttr {
    my ($attr) = @_;
    if ( $attr =~ /^(.+):.+/ ) {
	my $prefix = $1;
	if (defined $notImplementedExtensionPrefix{$prefix} and
	    $notImplementedExtensionPrefix{$prefix} == 0) {
	    &myWarn ("!! XML EXTENSION '$prefix' IS NOT IMPLEMENTED\n");
	    # we set the value to 1 so that the next time we will not prnt another message
	    $notImplementedExtensionPrefix{$prefix} = 1;
	}
	return 1;
    } else {
	return 0;
    }
} # end of isAnExtensionAttr

{
    my $inMetadata=0;
    sub metadata {
	$inMetadata++;
    }
sub _metadata {
    $inMetadata--;
}

sub inMetadata {
    return $inMetadata;
}
}

sub notYetImplemented {
    my ($elementname) = @_;
    &myWarn ("####### $elementname: Not Yet Implemented\n");
}

{
    my $expat;
sub Init {
    $expat = shift;
}
sub Final {
    undef $expat;
}

## takes 1 arg : 'message'
sub myWarn {
    my ($mess) = @_;
    if (defined $expat) {
	print STDOUT ("at ", $expat->current_line, ": $mess");
    } else {
	print STDOUT $mess;
    }
}

sub current_line {
    if (defined $expat) {
	return $expat->current_line; 
    } else {
	return "_undef_";
    }
}
}

sub display {
    my (@res) = @_;
    $backend->treatLines(@res);
}

sub findINC
{
 my $file = join('/',@_);
 my $dir;
 $file  =~ s,::,/,g;
 foreach $dir (@INC)
  {
   my $path;
   return $path if (-e ($path = "$dir/$file"));
  }
 return undef;
}


###################################################################
### this a slightly different implementation of the subs style as defined in XML::Parser
### Differences are :
# - when an error occure in a callback, the error is handled and a warning is
#     printed with the line number of the SVG source file
# - namespace can be used (this is usefull for example to treat the SVG included in dia data files)
#

package XML::Parser::SVG2zinc;
$XML::Parser::Built_In_Styles{'SVG2zinc'} = 1;


sub Start {
  no strict 'refs';
  my $expat = shift;
  my $tag = shift;
  my $ns = $expat->namespace($tag);
#  print "tag=$tag , ns=",$ns||" ", "\n";
  if (!defined $ns || $ns =~ /\/svg$/) {
      ## the tag is a SVG tag
      ## BUG: we should also get some tags of XML standard used by
      ## the SVG standard. Exemple: xlink:href
      my $sub = $expat->{Pkg} . "::$tag";
#     print "Sub=$sub\n";
      if (defined &$sub) {
	  eval { &$sub($expat, $tag, @_) };
	  if ($@) {
	      $expat->xpcarp("An Error occured while evaluationg $tag {...} :\n$@");
	  }
      } elsif (&SVG::SVG2zinc::inMetadata) {
	  # we do othing, unless tags were treated before!
      }
      else {
	  ## skipping the tag if it is part of not implemented extension
	  my ($extension) = $tag =~ /(\w+):.*/;
	  return if defined $extension &&  defined $notImplementedExtensionPrefix{$extension};
	  warn "## Unimplemented SVG tag: $tag\n";
      }
  }
}

sub End {
  no strict 'refs';
  my $expat = shift;
  my $tag = shift;
  my $ns = $expat->namespace($tag);
  if (!defined $ns || $ns =~ /\/svg$/) {
      my $sub = $expat->{Pkg} . "::${tag}_";
      ## the tag is a SVG tag
      if (defined &$sub) {
	  eval { &$sub($expat, $tag) };
	  if ($@) {
	      $expat->xpcarp("An Error occured while evaluationg ${tag}_ {...}) :\n$@");
	  }
      } else {
	  # the following error message is not usefull, as there were already
	  # an error message at the opening tag
	  # warn "## Unimplemented SVG tag: ${tag}_\n";
      }
  }
}



###################################################################


1;

__END__

=head1 NAME

SVG::SVG2zinc - a module to display or convert svg files in scripts, classes, images...

=head1 SYNOPSIS

 use SVG::SVG2zinc;

 &SVG::SVG2zinc::parsefile('file.svg', 'Backend','file.svg',
			   -out => 'outfile',
			   -verbose => $verbose,
			   -namespace => 0|1,
			   -prefix => 'string', 
			   );

 # to generate a Perl script: 
 &SVG::SVG2zinc::parsefile('file.svg','PerlScript', 
			   -out => 'file.pl');

 # to generate a Perl Class:
 &SVG::SVG2zinc::parsefile('file.svg','PerlClass', 
			   -out => 'Class.pm');

 # to display a svgfile:
 &SVG::SVG2zinc::parsefile('file.svg', 'Display'); 

 #To convert a svgfile in png/jpeg file:
 &SVG::SVG2zinc::parsefile('file.svg', 'Image', 
		           -out => 'file.jpg');

 # to generate a Tcl script: 
 &SVG::SVG2zinc::parsefile('file.svg','TclScript', 
			   -out => 'file.tcl');


=head1 DESCRIPTION

Depending on the used Backend, &SVG::SVG2zinc::parsefile either generates a Perl Class,
Perl script, Tcl Script, bitmap images or displays SVG files inside a Tk::Zinc widget.

SVG::SVG2zinc could be extended to generate Python scripts and/or
classes, or other files, just by sub-classing SVG::SVG2zinc::Backend(3pm)

==head1 HOW IT WORKS

This converter creates some TkZinc items associated to most SVG tags.
For example, <SVG> or <G> tags are transformed in TkZinc groups. <PATH>
are converted in TkZinc curves.... many more to come...

==head2 TkZinc items tags

Every TkZinc item created by the parser get one or more tags. If the
corresponding svg tag has an Id, this Id will be used as a tag, after
some cleaning due to TkZinc limitation on tag values (no dot, star, etc...).
If the corresponding svg tag has no Id, the parser add a tag of the
following form : __<itemtype>__<integer>. If the parser is provided
a B<-prefix> option, the prefix is prepended to the tag:
<prefix>__<itemtype>__<integer>

The TkZinc group associated to the top <SVG> tag has the following tag 'svg_top', as well as 'width=integer' 'heigth=integer' tags if width and height are defined in the top <SVG> tag. These tags can be used to find the group and to get its desired width and height.

==head2 RunTime code

There is currently on new Tk::Zinc method needed when executing perl code generated.
This perl Tk::Zinc::adaptViewport function should be translated and included or
imported in any script generated in an other scripting language (eg. Tcl or Python).

=head1 BUGS and LIMITATIONS

Some limitations are due to differences between Tk::Zinc and SVG graphic models :

=over 2

=item B<Drawing width>

Drawing width are zoomed in SVG but are not in Tk::Zinc where it is constant whatever the zoom factor is.

=item B<Gradient Transformation>

Gradient Transformation is not possible in Tk::Zinc. May be it could be implemented by the converter?

=item B<Rounded Rectangles>

Rectangles cannot have rounded corners in Tk::Zinc. Could be implemented, by producing curve item rather than rectangles in Tk::zinc. Should be implemented in a future release of Tk::Zinc  

=item B<Text and tspan tags>

Text and tspan tags are very complex items in SVG, for example placement can be very precise and complex. Many such features are difficult
to implement in Tk::Zinc and are not currently implemented

=item B<Font>

Font management is still limited. It will be rotatable and zoomable in future release of Tk::Zinc. SVG fonts included in a document are not readed, currently.

=item B<SVG image filtering>

No image filtering functions are (and will be) available with Tk::Zinc, except if YOU want to contribute?

=item B<ClipPath tag>

The SVG ClipPath tag is a bit more powerfull than Tk::Zinc clipping (clipper is limited to one item). So currently this is not implemented at all in SVG::SVG2zinc 

=back


There are also some limitations due to the early stage of the converter:

=over 2

=item B<CSS>

CSS in external url is not yet implemented

=item B<SVG animation and scripting>

No animation is currently available, neither scripting in the SVG file. But Perl or Tcl are scripting languages, are not they?

=item B<switch tag>

The SVG switch tag is only partly implemented, but should work in most situations

=item B<href for images>

href for images can only reference a file in the same directory than the SVG source file.

=back

It was said there is still one hidden bug... but please patch and/or report it to the author! Any (simple ?)

SVG file not correctly rendered by this module (except for limitations
listed previously) could be send to the author with little comments
about the expected rendering and observed differences.

=head1 SEE ALSO

svg2zinc.pl(1) a sample script using and demonstrating SVG::SVG2zinc

SVG::SVG2zinc::Backend(3pm) to defined new backends.

Tk::Zinc(3) TkZinc is available at www.openatc.org/zinc/

=head1 AUTHORS

Christophe Mertz <mertz at intuilab dot com>

many patches and extensions from Alexandre Lemort <lemort at intuilab dot com>

helps from Celine Schlienger <celine at intuilab dot com> and Stéphane Chatty <chatty at intuilab dot com>

=head1 COPYRIGHT

CENA (C) 2002-2004, IntuiLab (C) 2004

This program is free software; you can redistribute it and/or modify it under the term of the LGPL licence.

=cut
