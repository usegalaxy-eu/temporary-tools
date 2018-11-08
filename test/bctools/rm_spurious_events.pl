#!/usr/local/perl/bin/perl
use feature ':5.10';
use strict 'vars';
use warnings;
use Getopt::Long;
use Pod::Usage;
use List::Util qw/ min max /;
use POSIX qw/ceil floor/;
use File::Temp qw(tempdir);
use File::Basename;

=head1 NAME

filterSquashedReads.pl --frac_max FLOAT

=head1 SYNOPSIS

this script reads sites after pcr duplicate removal.
for all crosslinking sites sharing start and stop coordinates, the maximum number
of squashed reads is determined.
alignments having less than frac_max * max reads are discarded

assumes bed entries to be sorted chr,strand,start,stop,score with score descending

Options:

    --frac_max  filter out alignments supported by less reads than this fraction of the maximum number of reads per position
    -debug      enable debug output
    -help       brief help message
    -man        full documentation

=head1 DESCRIPTION

=cut

###############################################################################
# parse command line options
###############################################################################
my $frac_max = 0.1;
my $help;
my $man;
my $debug;
my $result = GetOptions (   "frac_max=f" => \$frac_max,
                            "help"	=> \$help,
							"man"	=> \$man,
							"debug" => \$debug);
pod2usage(-exitstatus => 1, -verbose => 1) if $help;
pod2usage(-exitstatus => 0, -verbose => 2) if $man;
($result) or pod2usage(2);

###############################################################################
# main
###############################################################################
my $current_chr = '';
my $current_start = -1;
my $current_stop = -1;
my $current_strand = '';
my $current_max = -1;
my $current_threshold = -1;

while (<>) {
    my ($chr, $start, $stop, $id, $count, $strand) = split("\t");
    # my ($count, undef, $start, $chr, $strand, $stop) = split("\t");

    if ($current_start != $start or
        $current_stop != $stop or
        $current_chr ne $chr or
        $current_strand ne $strand) {
        # if this is the first occourence of these coordinates
        # this must be the new maximum
        # save new state
        $current_start = $start;
        $current_stop = $stop;
        $current_chr = $chr;
        $current_strand = $strand;
        # print record and record maximum
        $current_max = $count;
        $current_threshold = $count*$frac_max;
        print $_;
        $debug and say "new threshold ${current_threshold} @ $chr $start $stop $strand $count";
    } elsif ($count >= $current_threshold) {
        # if it is not the first occourence, evaluate threshold and print if valid
        print $_;
    } else {
        $debug and say "below threshold @ $chr $start $stop $strand $count";
    }
}
