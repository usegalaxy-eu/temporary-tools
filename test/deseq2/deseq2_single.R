## Setup R error handling to go to stderr
options( show.error.messages=F, error = function () { cat( geterrmessage(), file=stderr() ); q( "no", 1, F ) } )
# we need that to not crash galaxy with an UTF8 error on German LC settings.
Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")

library('getopt')
library('rjson')
library('DESeq2')
library("RColorBrewer")
library("gplots")

options(stringAsfactors = FALSE, useFancyQuotes = FALSE)
args <- commandArgs(trailingOnly = TRUE)

#get options, using the spec as defined by the enclosed list.
#we read the options from the default: commandArgs(TRUE).
spec = matrix(c(
    'verbose', 'v', 2, "integer",
    'help', 'h', 0, "logical",
    'outfile', 'o', 1, "character",
    "countsfile", "n", 1, "character",
    'plots', 'p', 2, "character",
    'factors', 'm', 2, "character"
), byrow=TRUE, ncol=4);
opt = getopt(spec);


# if help was asked for print a friendly message
# and exit with a non-zero error code
if ( !is.null(opt$help) ) {
    cat(getopt(spec, usage=TRUE));
    q(status=1);
}

trim <- function (x) gsub("^\\s+|\\s+$", "", x)
opt$samples <- trim(opt$samples)
opt$factors <- trim(opt$factors)

parser <- newJSONParser()
parser$addData( opt$factors )
factorsList <- parser$getObject()

sampleTable<-data.frame()
factorNames<-c()
primaryfactor = TRUE
for(factor in factorsList){
    factorName<-factor[[1]]
    factorNames<-append(factorNames, factorName)
    factorValuesMapList<-factor[[2]]
    c = length(factorValuesMapList)
    for (factorValuesMap in factorValuesMapList){
        for(files in factorValuesMap){
            fvc = 0
            for(file in files){
                fvc = fvc+1
                if(primaryfactor) {
                    sampleTable[basename(file),"sampleName"]<-paste(fvc,names(factorValuesMap),sep="_")
                }
                sampleTable[basename(file),"fileName"]<-file
                sampleTable[basename(file),factorName]<-paste(c,names(factorValuesMap),sep="_")
            }
        }
        c = c-1
    }
    primaryfactor = FALSE
}

factorNames<-rev(factorNames)
designFormula <- as.formula(paste("", paste(factorNames, collapse=" + "), sep=" ~ "))
designFormula
ddsHTSeq = DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                      directory = "",
                                      design =  designFormula)

ddsHTSeq <- DESeq(ddsHTSeq) #, minReplicatesForReplace=Inf
deseqRes <- results(ddsHTSeq, alpha=0.05) #, cooksCutoff=FALSE, independentFiltering=FALSE
mcols(deseqRes)$description
summary(deseqRes)
attr(deseqRes,"filterThreshold")
resSorted <- deseqRes[order(deseqRes$padj),]
head(resSorted)
write.table(as.data.frame(resSorted), file = opt$outfile, sep="\t", quote = FALSE, col.names = FALSE)

if (!is.null(opt$countsfile)) {
    factors <- sapply(factorsList, function(x) x[[1]])
    print(head(ddsHTSeq))
    print(factorsList)
    labs <- paste0(seq_len(ncol(ddsHTSeq)), ": ", do.call(paste, as.list(colData(ddsHTSeq)[factors])))
    normalizedCounts<-counts(ddsHTSeq,normalized=TRUE)
    colnames(normalizedCounts)<-labs
    write.table(normalizedCounts, file=opt$countsfile, sep="\t", col.names=NA, quote=FALSE)
}


use <- deseqRes$baseMean > attr(deseqRes,"filterThreshold")
table(use)
deseqResFilt <- results(ddsHTSeq)
pvalFilt = deseqResFilt$pvalue
pvalFilt[!use] = NA
deseqResFilt$padj = p.adjust(pvalFilt, method ="BH")

if(!is.null(opt$plots)){
    pdf(opt$plots)

    # heatmap
    plotDispEsts(ddsHTSeq, main= "Dispersion estimate plot")
    plotMA(ddsHTSeq, main= "MA-plot", ylim=c(-2,2), alpha= 0.05)
    select <- order(rowMeans(counts(ddsHTSeq,normalized=TRUE)),decreasing=TRUE)[1:30]
    hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
    rld <- rlogTransformation(ddsHTSeq, blind=TRUE)
    vsd <- varianceStabilizingTransformation(ddsHTSeq, blind=TRUE)
    distsRL <- dist(t(assay(rld)))
    mat <- as.matrix(distsRL)
    hc = hclust(distsRL)
    heatmap.2(mat, Rowv=as.dendrogram(hc), trace="none", col = rev(hmcol), main = "Sample-to-sample distances", margin=c(13,13))

    # PCA
    factorNames<-rev(factorNames)
    print(factorNames)
    print(plotPCA(rld, intgroup=factorNames))

    ## p-value distribution
    h1 <- hist(deseqRes$pvalue[!use], breaks=0:50/50, plot=FALSE)
    h2 <- hist(deseqRes$pvalue[use], breaks=0:50/50, plot=FALSE)
    colori <- c(`filtered out`="khaki", `complete`="powderblue")
    barplot(height = rbind(h1$counts, h2$counts), beside = FALSE,
            col = colori, space = 0, main = "Distribution of p-values", ylab="frequency")
    text(x = c(0, length(h1$counts)), y = 0, label = paste(c(0,1)), adj = c(0.5,1.7), xpd=NA)
    legend("topright", fill=rev(colori), legend=rev(names(colori)))

    dev.off()
}

sessionInfo()


