clustersize <- as.integer(commandArgs(1))

data <- read.table("cluster.csv", header=FALSE, row.names=1)
#data <- data[order(data$V3),]

gff_frame <- data.frame()

findClusterPeak <- function(cluster, clustersize, gff_frame) {

    if ( length(cluster$V2) == 0) { # done
        return(gff_frame)
    }

    # get max value
    the_max_block <- which(cluster$V6==max(cluster$V6))
    the_max_block_size <- cluster[the_max_block,]$V6

    if (length(the_max_block_size) > 1) { # for same maxblock case
        the_max_block_size <- the_max_block_size[1] 
        the_max_block <- the_max_block[1] 
    }

    #print(the_max_block)
    #print(the_max_block_size)
    #print(length(cluster$V3))

    if ( (the_max_block_size/clustersize >= 0.01) ) { # tunable parameter
 
    minoverlap <- round( (cluster[the_max_block,]$V4-cluster[the_max_block,]$V3) * 0.5) # tunable parameter / how much does a subblock minimally need to overlap with the top block
    min_perc_of_max_block <- 0.1 # tunable parameter / how much percent of the top block in a peak must a subblock be at minimum
 
    # get the blocks overlapping with the maximum block // specifically it identifies the non overlapping and gets the rest as the overlapping
    peak_blocks_indices <- setdiff(as.integer(rownames(cluster)),c(which(cluster$V4<(cluster[the_max_block,]$V3+minoverlap)),which(cluster$V3>(cluster[the_max_block,]$V4-minoverlap))))
    peak_blocks_indices_strings <- unlist(strsplit(paste(peak_blocks_indices, collapse=","), ","))
    peak_blocks <- cluster[peak_blocks_indices_strings,]

    # prepare for next iteration // all blocks that have an overlap with the_max_block are discarded for the next iteration
    next_cluster_indices <- c(which(cluster$V4<cluster[the_max_block,]$V3),which(cluster$V3>cluster[the_max_block,]$V4))
    next_cluster <- cluster[next_cluster_indices,]
    rownames(next_cluster) <- NULL

    # omit lowly expressed blocks from boundary definition
    peak_blocks<-peak_blocks[which( (peak_blocks$V6/the_max_block_size)>=min_perc_of_max_block),]

    # get the peak boundaries    
    the_min_boundary <- min(peak_blocks$V3)
    the_max_boundary <- max(peak_blocks$V4)
 
    # add boundaries to gff frame
    gff_frame <- rbind(gff_frame, c(the_min_boundary,the_max_boundary))

    # call next iteration
    findClusterPeak(next_cluster, clustersize, gff_frame)

    } else {
        # done
        return(gff_frame)
    }
}

peaks<-findClusterPeak(data, clustersize, gff_frame)

write.csv(file="peaks.csv", peaks, row.names=FALSE)

