#!/usr/bin/env perl

# Flexbar wrapper for Galaxy tool definition, version 2.5
# Author: Johannes Roehr

use warnings;
use strict;

# this parses the last 4 arguments
# what is id -> xml $output.id
# what is folder -> xml $__new_file_path__
# what is format? -> xml $reads.ext
my ($outFile, $id, $folder, $format) = @ARGV[($#ARGV - 3) .. $#ARGV];

# this parses all but the last four arguments
# contains the call to the flexbar actual
my $call = join " ", @ARGV[0..($#ARGV - 4)];

# this calls flexbar and
# prefix for output files will be "FlexbarTargetFile"
system $call .' --target FlexbarTargetFile > '. $outFile and exit 1;

# now we parse all output files
foreach(<FlexbarTargetFile*>){

	# determine filetype
	my $fileType;

	$fileType = $1         if /\.(\w+)$/;
	$fileType = $format    if /\.\w*fast\w$/;
	$fileType = 'fasta'    if /\.fasta$/;
	$fileType = 'csfasta'  if /\.csfasta$/;
	$fileType = 'tabular'  if /\.lengthdist$/;

	# this is just the filename from the for loop
	my $file = $_;

	# replace underscores by minus in $_
	s/_/-/g;

	# set new name for output files
	my $name = "primary_". $id ."_". $_ ."_visible_". $fileType;

	# rename output file to a pattern recognized by flexbar?
	rename $file, $name;
	# best guess: this seems to move the file into a folder
	# rename behavious is not specified by perl... or implementation differs.
	rename $name, $folder;
}
