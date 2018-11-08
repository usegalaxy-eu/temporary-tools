#!/usr/local/perl/bin/perl

use strict;
use warnings;
use FindBin qw($Bin);

my $blockbusterPeaks = $ARGV[0]; # blockbuster output type 1

# /scratch/2/wrightp/HFQ_Vogel/scripts/make_gff_from_blockbuster_table.pl c-lib_41_alignments_unique_hits_height_10.blockbuster > c-lib_41_clip_peaks_height_10.gff

open(DATA, $blockbusterPeaks) or die("Error: cannot open file $blockbusterPeaks'\n");
    my @blockbusterLines = <DATA>;
close DATA;

push(@blockbusterLines, ">"); # for the last peak in the file

my @cluster = ();

# load the clusters
for(my $i=0;$i<scalar(@blockbusterLines);$i++) {
  
   if ( $blockbusterLines[$i] =~ m/>/) {
       push (@cluster, $blockbusterLines[$i]);
   } else {
       push (@cluster, $blockbusterLines[$i]);
       if ($blockbusterLines[($i+1)] =~ m/>/) {
           callClusterPeaks(@cluster);
           @cluster = ();
       }
   }
}

system "rm cluster.csv peaks.csv";

sub callClusterPeaks{
    
    my (@cl) = @_; # clusterdesctiption + blocks
    
    my @split = split("\t", $cl[0]);
    my $clustersize = $split[5];
    my $strand = $split[4];
    my $replicon = $split[1];

    if(scalar(@cl) == 2) { # simple case, only one block in cluster

        # print gff line with min and max from the one block you have (not the cluster)
        my $block = $cl[1];
        my @split_block = split(/\t/, $block);
        #my $peakStart = $split_block[2]+1;
        my $peakStart = $split_block[2];
        my $peakEnd = $split_block[3];
 
        # print gff line
        print "$replicon\tblockbuster\tclippeak\t$peakStart\t$peakEnd\t.\t$strand\t.\tpeak_Count_to_be_added_as_number\n"; 

    } else { 
              # non trivial case where blocks have to be split into several peaks
              # select maximum block from cluster and extend it to both sides
              # with all blocks that are at least a certain %age (10% -> 0.1) of the
              # max block 
              # also to be considered as an initail peak block a block must represent at
              # least a certain %age of the clustersize (1% -> 0.01)
              # all blocks that overlap with at least half (0.5) of the initial peaks boundaries
              # are considered for extension to left and right
              # if a block overlaps less (but it actually overlaps) it is discarded and
              # also not used on the next iteration
              # every block that has been added to a peak is removed from the cluster
              # when no more blocks are left to add to one peak, a new maximum is selected.
        
        open(WRITETO, '>cluster.csv');
        foreach(my $i=1; $i<scalar(@cl); $i++) { 
             # don't print first line but extract clustersize to pass on to the R script
             # pass the cluster to the R script
             print WRITETO $cl[$i]; # print to file 
        }
        close WRITETO;
        system "Rscript $Bin/split_cluster_peaks.R $clustersize &> /dev/null";
        
        # use file constructed from R script to build the appropriate gff lines and append them 
        open(DATA, "peaks.csv") or die("Error: cannot open file peaks.csv\n");
            my @peakLines = <DATA>;
        close DATA;
        foreach(my $i=1; $i<scalar(@peakLines); $i++) {
            my @split_peaks = split(",", $peakLines[$i]);
            my $peakStart = $split_peaks[0]+1;
            my $peakEnd = $split_peaks[1];
            chomp $peakEnd;
            # print gff line
            print "$replicon\tblockbuster\tclippeak\t$peakStart\t$peakEnd\t.\t$strand\t.\tpeak_Count_to_be_added_as_number\n"; 
        }
    }
}


