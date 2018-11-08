## Setup R error handling to go to stderr
options( show.error.messages=F, error = function () { cat( geterrmessage(), file=stderr() ); q( "no", 1, F ) } )
# we need that to not crash galaxy with an UTF8 error on German LC settings.
Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")

library('getopt');
options(stringAsfactors = FALSE, useFancyQuotes = FALSE)
args <- commandArgs(trailingOnly = TRUE)

#get options, using the spec as defined by the enclosed list.
#we read the options from the default: commandArgs(TRUE).
spec = matrix(c(
    'verbose', 'v', 2, "integer",
    'help' , 'h', 0, "logical",
    'outfile' , 'o', 1, "character",
    'outfilefiltered' , 'f', 1, "character",
    'plots' , 'p', 2, "character",
    'input' , 'i', 1, "character",
    'factors', 'm', 2, "character",
    'fittype', 't', 2, "character",
    'filtermode', 'l', 2, "character",
    'threshold', 'c', 2, "double"
#    'organism', 'g', 2, "character"
), byrow=TRUE, ncol=4);
opt = getopt(spec);

# if help was asked for print a friendly message
# and exit with a non-zero error code
if ( !is.null(opt$help) ) {
    cat(getopt(spec, usage=TRUE));
    q(status=1);
}

if( is.null(opt$fittype))
    opt$fittype = "parametric"

if( is.null(opt$filtermode))
    opt$filtermode = "absolute"

trim <- function (x) gsub("^\\s+|\\s+$", "", x)
opt$samples <- trim(opt$samples)
opt$factors <- trim(opt$factors)

htseqCountTable = read.table(opt$input, sep="\t", comment="", as.is=T,row.names=1)
colnames(htseqCountTable) <- htseqCountTable[2,]
# filter out the count statistics from the count table
htseqCountTable<-htseqCountTable[(grep('no_feature|ambiguous|alignment_not_unique',rownames(htseqCountTable),invert=TRUE)),]
names(htseqCountTable)<-colnames(htseqCountTable)
conditions <- htseqCountTable[1,]
conditions <- unlist(conditions[,-1])
htseqCountTable <- htseqCountTable[-(1:2),]
featurenames<-rownames(htseqCountTable)
htseqCountTable[,1] <- gsub(",", ".", htseqCountTable[,1])
head(htseqCountTable)

library('rjson')
library('DESeq2')
#library('biomaRt')

if ( !is.null(opt$plots) ) {
    pdf(opt$plots)
}

## The following function takes deseq data object, computes, writes and plots the results.
computeAndWriteResults <- function(dds, sampleCols, outputcsv, featurenames_filtered) {
    dds <- DESeq(dds, fitType= opt$fittype)
    print(sizeFactors(dds))
    dds <- nbinomWaldTest(dds, cooksCutoff=FALSE)
    res <- results(dds)
    print(mcols(res)$description)
    resCols <- colnames(res)
    cond <- c(rep(paste(unique(conditions[sampleCols]),collapse="_"), nrow(res)))
    if(missing(featurenames_filtered)){
        res[, "geneIds"] <- featurenames
        title_prefix = "Complete: "
    }else{
        res[, "geneIds"] <- featurenames_filtered
        title_prefix = "Filtered: "
    }
#    if(opt$organism != "other"){
#        dataset = ""
#        if(opt$organism == "mouse")
#            dataset = "mmusculus_gene_ensembl"
#        else if(opt$organism == "human")
#            dataset = "hsapiens_gene_ensembl"
#        else if(opt$organism == "fly")
#            dataset = "dmelanogaster_gene_ensembl"
#        ensembldb = useMart("ensembl",dataset=dataset)

#        annot <- getBM(attributes = c("ensembl_gene_id", "external_gene_name","chromosome_name","start_position","end_position","strand","gene_biotype","description"),
#                    filters = "ensembl_gene_id",
#                    values=res[, "geneIds"],
#                    mart=ensembldb)

#        res <- merge(res, annot,
#                    by.x = "geneIds",
#                    by.y = "ensembl_gene_id",
#                    all.x=TRUE)

#        resCols <- colnames(res)
#        resSorted <- res[order(res$padj),]
#        write.table(as.data.frame(resSorted[,c(resCols)]), file = outputcsv, sep="\t", quote = FALSE, append=TRUE, row.names = FALSE, col.names = FALSE)
#    }else{
        resSorted <- res[order(res$padj),]
        write.table(as.data.frame(resSorted[,c("geneIds", resCols)]), file = outputcsv, sep="\t", quote = FALSE, append=TRUE, row.names = FALSE, col.names = FALSE)
#    }

    if(!is.null(opt$plots)){
        plotDispEsts(dds, main= paste(title_prefix, "Dispersion estimate plot"))
        plotMA(dds, main= paste(title_prefix, "MA-plot"))

        library("RColorBrewer")
        library("gplots")
        select <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)[1:30]
        hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)

        rld <- rlogTransformation(dds, blind=TRUE)
        vsd <- varianceStabilizingTransformation(dds, blind=TRUE)
        distsRL <- dist(t(assay(rld)))
        mat <- as.matrix(distsRL)
        heatmap.2(mat, trace="none", col = rev(hmcol), main = paste(title_prefix, "Sample-to-sample distances"), margin=c(13,13))
    }
    return(res)
}

## This functions filters out the genes with low counts and plots p-value distribution
independentFilter <- function(deseqRes, countTable, sampleColumns, designFormula){
    if(opt$filtermode == "absolute"){
        use = (deseqRes$baseMean > opt$threshold & !is.na(deseqRes$pvalue))
    }
    else if(opt$filtermode == "quantile"){
        opt$threshold = opt$threshold/100
        print(paste("Threshold:", opt$threshold))
        print(quantile(unique(deseqRes$baseMean), probs = opt$threshold))
        use = (deseqRes$baseMean > quantile(unique(deseqRes$baseMean), probs = opt$threshold) & (!is.na(deseqRes$pvalue)))
    }
    else{
        print(paste("Unknown filtermode:", opt$filtermode))
        q("no", 1, FALSE)
    }

    ddsFilt <- DESeqDataSetFromMatrix(countData = countTable[use, ], colData = sampleTable, design = designFormula)
    computeAndWriteResults(ddsFilt, sampleColumns, opt$outfilefiltered, featurenames[use])

    ## p-value distribution
    h1 <- hist(deseqRes$pvalue[!use], breaks=50, plot=FALSE)
    h2 <- hist(deseqRes$pvalue[use], breaks=50, plot=FALSE)
    colori <- c(`filtered out`="khaki", `complete`="powderblue")
    barplot(height = rbind(h1$counts, h2$counts), beside = FALSE,
            col = colori, space = 0, main = "Distribution of p-values", ylab="frequency")
    text(x = c(0, length(h1$counts)), y = 0, label = paste(c(0,1)), adj = c(0.5,1.7), xpd=NA)
    legend("topright", fill=rev(colori), legend=rev(names(colori)))
}

parser <- newJSONParser()
parser$addData( opt$factors )
factorsList <- parser$getObject()
sampleColumns <- c()

for (factor in factorsList){
    factorName<-factor[[1]]
    factorValuesMapList<-factor[[2]]
    for (factorValuesMap in factorValuesMapList){
        for(i in factorValuesMap){
            sampleColumns[[length(sampleColumns) + 1]] <- as.integer(i) - 1
#            sampleTable[sampleNames[i-1],factorName]<-paste(c,names(factorValuesMap),sep="_")
        }
    }
}


#for (i in 1:length(factorsList)){
#    dat2<-unlist(factorsList[i], recursive=F)
#    print(factorsList[i])
#    for (j in 1:length(dat2)){
#        dat3<-unlist(dat2[[2]][j], recursive=F)
#        dat4<-unlist(dat3[[1]])
#        for (k in 1:length(dat4)){
#            sampleColumns[[length(sampleColumns) + 1]] <- as.integer(dat4[k]) - 1
#        }
#    }
#}

sampleColumns<-sort(unique(sampleColumns))
sampleColumns
htseqSubCountTable <- htseqCountTable[,sampleColumns]
head(htseqSubCountTable)
ntabcols <- length(htseqSubCountTable)

htseqSubCountTable <- as.integer(unlist(htseqSubCountTable))
htseqSubCountTable <- matrix(unlist(htseqSubCountTable), ncol = ntabcols, byrow = FALSE)
head(htseqSubCountTable)
colnames(htseqSubCountTable) <- names(htseqCountTable)[sampleColumns]
colnames(htseqSubCountTable)
sampleTable = data.frame(row.names=colnames(htseqSubCountTable))
sampleNames<-colnames(htseqSubCountTable)
names(sampleNames)<-c(1:length(sampleNames))

#factor is a list, where first element is the name of the factor. The second element is list of hashes of factor levels and their sample numbers eg: ["Condition", [{"KO": [4, 5]}, {"WT": [2, 3]}]]
for(factor in factorsList){
    factorName<-factor[[1]]
    factorValuesMapList<-factor[[2]]
    n = length(factorValuesMapList)
    c = 1
    for (factorValuesMap in factorValuesMapList){
        for(i in factorValuesMap){
            sampleTable[sampleNames[i-1],factorName]<-paste(n,names(factorValuesMapList)[c],sep="_")
        }
        n = n-1
        c = c+1
    }
}

sampleTable[is.na(sampleTable)] <- "undefined"
sampleTable

factorNames <- c(names(sampleTable)[-1], names(sampleTable)[1])
designFormula <- as.formula(paste("", paste(factorNames, collapse=" + "), sep=" ~ "))
designFormula
dds = DESeqDataSetFromMatrix(countData = htseqSubCountTable,  colData = sampleTable, design = designFormula)
deseqres <- computeAndWriteResults(dds, sampleColumns, opt$outfile) 

## independent filtering
independentFilter(deseqres, htseqSubCountTable, sampleColumns, designFormula)

mcols(deseqres)$description

dev.off()
sessionInfo()
