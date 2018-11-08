#!/usr/bin/env python
#-*- coding: utf-8 -*-

import itertools
import argparse # to parse command line arguments
import pysam 
import copy as cp
import os
import math;
import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib.font_manager import FontProperties
from pylab import *
import time
import multiprocessing
from collections import OrderedDict
#from bisect import *
# to read bed files
from bx.intervals.io import GenomicIntervalReader

debug = 0
ioff()

def  RdBu():

    cdict = {
    'red'  :  ((0., 1., 1.), (0.5, 1, 1), (1., 0., 0.)),
    'green':  ((0., 0., 0.), (0.5, 1, 1), (1., 0., 0.)),
    'blue' :  ((0., 0., 0.), (0.5, 1, 1), (1., 1., 1.))
    }

    myCmap = matplotlib.colors.LinearSegmentedColormap('my_colormap', cdict, 1024)
    return(myCmap)

def  RdYlBu():

    cdict = {
    'red'  :  ((0., 1., 1.), (0.5, 0.98, 0.98), (1., 0, 0)),
    'green':  ((0., 0.4, 0.4), (0.5, 0.98, 0.98), (1., 0.2, 0.2)),
    'blue' :  ((0., 0.4, 0.4), (0.5, 0.82, 0.82), (1., 0.6, 0.6))
    }

    myCmap = matplotlib.colors.LinearSegmentedColormap('my_colormap', cdict, 1024)
    return(myCmap)

def flattenMatrix(matrixDict):
    """
    concat, flatten and remove nans from matrix
    """
    matrixFlatten = np.concatenate( [ x for x  in matrixDict.values()] ).flatten()
    return matrixFlatten[ logical_not(np.isnan(matrixFlatten)) ]

def plotMatrix( heatmapMatrixDict, regionsDict, windowWidth, 
                beforeRegionStartLength,
                afterRegionStartLength,
                regionBodyLength,
                outFileName, outFileFormat='png', 
                colorMap='binary', missingDataColor='black',
                plotTitle='', 
                xAxisLabel='', yAxisLabel='',  regionsLabel='',  
                zMin=None, zMax=None,
                yMin=None, yMax=None,
                averageType='median',
                referencePointLabel='TSS',
                startLabel='TSS', endLabel="TES", 
                heatmapHeight=25, 
                heatmapWidth=7.5,
                onePlotPerGroup=False, whatToShow='plot, heatmap and scale',
                sortRegions='no'):
#    colorMap = RdYlBu()
    tickPlotAdj = 0.5
    matrixFlatten = None
    if zMin == None:
        matrixFlatten = flattenMatrix( heatmapMatrixDict )    
        zMin = prctile(matrixFlatten, 2.0)
        if np.isnan(zMin):
            zMin = None 

    if zMax == None:
        if matrixFlatten == None:
            matrixFlatten = flattenMatrix( heatmapMatrixDict )    
#        zMax = max(np.median(heatmapMatrix, axis=0)) + np.std(heatmapMatrix) / 2
        zMax = prctile(matrixFlatten, 98.0);
        if np.isnan(zMax):
            zMax = None 

    rcParams['font.size'] = 10.0
#    rcParams['font.size'] = 9.0

    showSummaryPlot = False
    showHeatmap = False
    showColorbar = False

    if whatToShow == 'plot and heatmap':
        showSummaryPlot = True
        showHeatmap = True
    elif whatToShow == 'plot only':
        showSummaryPlot = True
    elif whatToShow == 'heatmap only':
        showHeatmap = True
    elif whatToShow == 'colorbar only':
        showColorbar = True
    elif whatToShow == 'heatmap and colorbar':
        showHeatmap = True
        showColorbar = True
    else:
        showSummaryPlot = True
        showHeatmap = True
        showColorbar = True

    sumFigSpacer = 0
    sumFigHeightInches = 0
    spaceBetweenClusterNames = 0

    if showSummaryPlot:
        spaceBetweenClusterNames = 0.55
        sumFigHeightInches = 1.5

    if showSummaryPlot and showHeatmap:
        sumFigSpacer = 0.3

    # the heatmapHeight value is given in cm
    heatmapHeightInches = float(heatmapHeight) / 2.54 if showHeatmap or showColorbar else 0
    figWidth =  float(heatmapWidth) / 2.54

    # measures are in inches
    topBorder = 0.36
    bottomBorder = 0.1
    smallSpacer = 0.05

#    import pdb;pdb.set_trace()

#    regionGroups.append( (len(regionList),regionsLabel) )
    numGroups = len(heatmapMatrixDict.keys())

    heatmapMatrixDict = mergeSmallGroups(heatmapMatrixDict)

    if numGroups == 1:
        # when the number of regions is just one, is better
        # to use the defaults for one plot per group, because
        # it does not tries to plot a legend below the plot,
        # which for one region is unecesary
        onePlotPerGroup = True

    # heatmapHeightInches  height of heatmap plus bottom spacing
    if onePlotPerGroup:
        # decide the figure size based on the number of groups defined
        # Each summary plot occupies 1 x 1.3 inches plus 0.4 inches for 
        # spacing. The heatmap occupies heatmapHeightInches inches
        figHeight = topBorder + bottomBorder +  numGroups * (sumFigHeightInches + sumFigSpacer) + heatmapHeightInches
    else:
        figHeight = topBorder + bottomBorder + ( numGroups / 2 ) * spaceBetweenClusterNames + sumFigHeightInches + sumFigSpacer + heatmapHeightInches

    sumFracHeight = float(sumFigHeightInches) / figHeight
    # fraction of the height for summary plots

    # fraction of the height for the whole heatmap
    heatFracHeight = heatmapHeightInches / figHeight

    topBorderFrac = topBorder/figHeight
    bottomBorderFrac = bottomBorder/figHeight
    spacerFrac = sumFigSpacer/figHeight
    smallSpacerFrac = smallSpacer/figHeight
    
    # figsize: w,h tuple in inches
    fig = plt.figure(figsize=(figWidth,figHeight), dpi=360)

    b = beforeRegionStartLength
    a = afterRegionStartLength
    m = regionBodyLength
    w = windowWidth

    if b < 1e5:
        quotient = 1000
        symbol = 'Kb'
    if b >= 1e5:
        quotient = 1e6
        symbol = 'Mb'

    if m == 0:
       xTicks = [int(k /float(w))- tickPlotAdj for k in [0, b, b  + a]]
       xTicksLabel = ['{0:.1f}'.format(-(float(b)/quotient)), referencePointLabel, '{0:.1f}{1}'.format(float(a)/quotient, symbol)];
    else:
       xTicks = [int(k /float(w))- tickPlotAdj for k in [0, b, b + m, b + m + a]]
       xTicksLabel = ['{0:.1f}'.format(-(float(b)/quotient)), startLabel, endLabel, '{0:.1f}kb'.format(float(a)/quotient, symbol)];

    fig.suptitle(plotTitle, y=1 - (0.06/figHeight))

    cmap = get_cmap(colorMap) 
    cmap.set_bad(missingDataColor) # nans are printed using this color

    # add_axes( rect ) where rect is [left, bottom, width, height] and
    # all quantities are in fractions of figure width and height.
    # summary plot
    # each group needs its axe
    left = 0.45 / figWidth
    width = 0.6

    if not showSummaryPlot:
        # when only the heatmap is to be printed use more space
        left = 0.28 / figWidth
        width = 0.7 - left
    if not showColorbar:
        width = 0.7


    ###### plot summary plot ###########
    if showSummaryPlot:

        # i.e if multipleLines in one plot:
        if not onePlotPerGroup:
            bottom = 1 - ( topBorderFrac + sumFracHeight)
            if debug:
                print ([left, bottom, width, sumFracHeight, figHeight])
            ax = fig.add_axes([left, bottom, width, sumFracHeight])

        index = -1
        for label, ma in heatmapMatrixDict.iteritems():
            index += 1
            if onePlotPerGroup:
                # create an axe for each sub plot
                bottom = 1 - topBorderFrac - (index + 1)*sumFracHeight - index*spacerFrac
                ax = fig.add_axes([left, bottom, width, sumFracHeight])
                ax.plot(matrixAvg(ma, averageType))
                ax.set_ylim(yMin,yMax)
                ax.axes.set_xticks(xTicks);
                ax.axes.set_xticklabels(xTicksLabel)
                # reduce the number of yticks by half
                if index==0:
#                    ax.axes.set_title(plotTitle, weight='bold');
                    numTicks =  len(ax.get_yticks())
                    yTicks = [ax.get_yticks()[i] for i in range(1, numTicks, 2)]
                ax.set_yticks(yTicks)
                ax.axes.set_ylabel(yAxisLabel);

            else:
                # add new lines to existing plot
                ax.plot(matrixAvg(ma, averageType),
                        label=label)

        # i.e if multipleLines in one plot:
        if not onePlotPerGroup:
            # in the case of one box with all plots the font of the legend and the positions
            # are changed.
            ax.set_ylim(yMin,yMax)
            ax.axes.set_xticks(xTicks);
            ax.axes.set_xticklabels(xTicksLabel)
#            ax.axes.set_title(plotTitle, weight='bold');
            fontP = FontProperties()
            fontP.set_size('small')
            # the legend shows in a box below the plot
            ax.legend(bbox_to_anchor=(-0.15, -1.2, 1.35, 1), loc='upper center',
                      ncol=2, mode="expand", borderaxespad=0., prop=fontP, 
                      frameon=False, markerscale=0.5)


            # reduce the number of yticks by half
            numTicks =  len(ax.get_yticks())
            yTicks = [ax.get_yticks()[i] for i in range(1, numTicks, 2)]
            ax.set_yticks(yTicks)

            ax.axes.set_ylabel(yAxisLabel);

    ###### plot heatmap plot ###########
    if showHeatmap:
#    heatmapLengths = [ float(newRegionGroups[x]-newRegionGroups[x-1])/newRegionGroups[-1] for x in range(1,len(newRegionGroups)) ]
        startHeatmap = heatFracHeight + bottomBorderFrac
        groupLengths = np.array([ len(x) for x in heatmapMatrixDict.values() ], dtype='float64')
        groupLengthsFrac = groupLengths / sum(groupLengths)
        index = -1
        for label, ma in heatmapMatrixDict.iteritems():
#            import pdb;pdb.set_trace()
#        for index in range(0, numGroups):
            index += 1
            # maks nans
            ma = np.ma.array( ma, mask=np.isnan(ma) )
            # the size of the heatmap is proportional to its length
            bottom = startHeatmap - sum(groupLengthsFrac[:index+1])*heatFracHeight 
            height = groupLengthsFrac[index]*heatFracHeight-smallSpacerFrac
            axHeat = fig.add_axes([left, bottom, width, height])
            # PRINT THE HEATMAP
#            SS = np.array([ myAverage( region, averageType ) for region in ma ]).argsort()
#            if sortRegions == 'ascend':
#                ma = ma[SS,:]
#            elif sortRegions == 'descend':
#                ma = ma[SS[::-1],:]
            img = axHeat.imshow(ma, 
                                aspect='auto',
                                interpolation='nearest', 
                                origin='upper',
                                vmin = zMin, 
                                vmax = zMax,
                                cmap = colorMap,
            ) 
            axHeat.axes.get_xaxis().set_visible(False)
            axHeat.axes.set_xlabel(xAxisLabel);
            axHeat.axes.set_ylabel(label);
            axHeat.axes.set_yticks([]);

    ###### add heatmap colorbar ###########
    if showColorbar:
        # this is the case when only the colorbar wants to be printed and nothing else
        if not showHeatmap:
            left = 0.2
            width = 0.6 - ( 0.20/figWidth)
            legend = fig.add_axes([left,bottomBorderFrac,width,heatFracHeight-smallSpacerFrac])
            norm = matplotlib.colors.Normalize(vmin=zMin, vmax=zMax)
            matplotlib.colorbar.ColorbarBase(legend, cmap=colorMap, norm=norm )
        else:
#            left = 1.3 / figWidth
            left = 0.79
            width = 0.05
            legend = fig.add_axes([left,bottomBorderFrac,width,heatFracHeight-smallSpacerFrac])
            if debug:
                print([left,bottomBorderFrac,width,heatFracHeight-smallSpacerFrac])
            fig.colorbar(img, cax=legend)

    if showSummaryPlot and not showHeatmap:
        # for some reason if the bbox_inches is used
        # then the legend for the summary plot is omited
        savefig(outFileName, format=outFileFormat);
    
    else:
        # To remove the white space from the surroundings of the image:
        # bbox_inches = 'tight' 
       # savefig(outFileName, format=outFileFormat, bbox_inches='tight');
        savefig(outFileName, format=outFileFormat);
     
def mergeSmallGroups(heatmapMatrixDict):

    groupLengths = [ len(x) for x in heatmapMatrixDict.values() ] 
    minGroupLength = sum(groupLengths) * 0.01
    toMerge = []
    i = 0
    _mergedHeatMapDict = OrderedDict()
    numGroups = len(groupLengths)
 
    for label, ma in heatmapMatrixDict.iteritems():
        # merge small groups together
        # otherwise visualization is impaired
        if groupLengths[i] > minGroupLength:
            if len(toMerge):
                toMerge.append(label)
                newLabel = " ".join(toMerge)
                newMa = np.concatenate( [ heatmapMatrixDict[item] for item in toMerge ], axis=0 )
            else:
                newLabel = label
                newMa = heatmapMatrixDict[label]

            _mergedHeatMapDict[newLabel] = newMa
            toMerge = []
        else:
            toMerge.append(label)
        i += 1
    if len(toMerge):
        newLabel = " ".join(toMerge)
        newMa = np.array()
        for item in toMerge:
            newMa = np.concatenate( [ newMa, heatmapMatrixDict[item] ] )
        _mergedHeatMapDict[newLabel] = newMa
        
    return _mergedHeatMapDict

def parseArguments(args=None):
    parser = argparse.ArgumentParser(description = 'Creates a heatmap for a score associated to genomic regions. Normally this regions are genes, but any other region defined in a GFF format will work.')

    # define the arguments
    parser.add_argument('--regionsFileName', '-R',
                        metavar = 'File',
                        help = 'File name, in GFF format, containing the regions to plot',
                        required = True)

    parser.add_argument('--scoreFileName', '-S',
                        help = 'File name of either a BigWig file containing a score, usually covering the whole genome or of a BAM file. For this last case, coverage counts will be used for the histogram',
                        metavar = 'File',
                        required = True)

    parser.add_argument('--scoreFileFormat', '-F',
                        help = 'File type, either BAM or BigWig',
                        choices = ["bam",  "bigwig"],
                        required = True)

    parser.add_argument('--regionBodyLength', '-m',
                        default = 0,
                        type = int,
                        help = 'Distance downstream of the start site of the regions defined in the GFF file. Only this section of the region would be plotted.')

    parser.add_argument('--windowWidth', '-w',
                        help = 'Length, in base pairs, of the non-overlapping window for averaging the score over the regions length',
                        type = int,
                        default = 10)

    parser.add_argument('--beforeRegionStartLength', '-b',
                        default = 500,
                        type = int,
                        help = 'Distance upstream of the start site of the regions defined in the GFF file. If the regions are genes, this would be the distance before the transcription start site.')

    parser.add_argument('--afterRegionStartLength', '-a',
                        default = 1500,
                        type = int,
                        help = 'Distance downstream of the start site of the regions defined in the GFF file. Only this section of the region would be plotted.')

    parser.add_argument('--referencePoint',
                        default = 'TSS',
                        choices = ['TSS', 'TES', 'center'],
                        help = '[ONLY] in the case of --regionBodyLength = 0, the reference point for the ploting could be eigther the region start (TSS) the end region end (TES) or the center of the region. The label for this point can be changed using the --startLabel and --endLabel options.')

    parser.add_argument('--startLabel', 
                        default = 'TSS',
                        help = 'Label for the region start. By default is TSS, but could be "start", "gene start", "peak start" etc. Same for the --endLabel option. See below')

    parser.add_argument('--endLabel', 
                        default = 'TES',
                        help = 'Label for the region end. By default is TES.')

    parser.add_argument('--regionsLabel', '-z',
                        default = 'genes',
                        help = 'Description for the regions plotted in the heatmap. This label will be overridden if labels are found in the regions BED file. Default is "genes"')

    parser.add_argument('--plotTitle', '-T',
                        help = 'Title of the plot, to be printed on top of the generated image. Leave blank for no title',
                        default = '')

    parser.add_argument('--sortRegions', 
                        help = 'Whether the heatmap should present the regions sorted. The default is to sort in descending order.',
                        choices = ["descend", "ascend", "no"],
                        default =  'descend')

    parser.add_argument('--missingDataAsZero', 
                        help = '[only for bigwig input] Set to "yes", if missing data should be treated as zeros. Default is to ignore such cases which will be depicted as black areas in the heatmap.',
                        action= 'store_true')


    parser.add_argument('--missingDataColor',
                        default = 'black',
                        help = 'If --missingDataAsZero is not set, such cases will be colored in black by default. Using this parameter a different color can be set. A value between 0 and 1 will be used for a gray scale (black is 0). Also color names can be used, see a list here: http://packages.python.org/ete2/reference/reference_svgcolors.html. Alternatibely colors can be specified using the #rrggbb notation.')


    parser.add_argument('--outFileName', '-out',
                        help = 'File name to save the image. ', 
                        type=writableFile,
                        required=True)

    parser.add_argument('--outFileFormat', '-O',
                        help = 'Output format for the plot. Options are "png", "emf", "eps", "pdf", "svg".',
                        default = "png",
                        choices = ["png", "emf", "eps", "pdf", "svg"])

    parser.add_argument('--outFileNameData',
                        help = 'File name to save the data underlaying data for the summary plot',
                        type=argparse.FileType('w'))

    parser.add_argument('--outFileNameMatrix',
                        help = 'If this option is given, then the matrix of values undelyting the heatmap will be saved using this name.', 
                        metavar = 'FILE',
                        type=writableFile)

    parser.add_argument('--outFileSortedRegions', 
                        help = 'File name in which the regions are saved after skiping zeros or min max threshold values. The order of the regions in the file follows the sorting order selected. This is useful for example to generate other heatmaps keeping the sorting of the first heatmap. ',
                        metavar = 'FILE',
                        type=argparse.FileType('w'))

    parser.add_argument('--averageType', '-at' ,
                        default = 'median',
                        choices = ["mean", "median", "min", "max", "std"],
                        help = 'Define the type of statistic that should be plotted in the summary image above the heatmap. The options are: "mean", "median", "min", "max" and "std". The default is "median". ')

    parser.add_argument('--xAxisLabel', '-x',
                        default = 'gene distance (bp)',
                        help = 'Description for the x-axis label')

    parser.add_argument('--yAxisLabel', '-y',
                        default = '',
                        help = 'Description for the y-axis label for the top panel')


    parser.add_argument('--zMin', '-min',
                        default = None,
                        help = 'Minimun value for the heatmap intensities')

    parser.add_argument('--zMax', '-max',
                        default = None,
                        help = 'Maximun value for the heatmap intensities')

    parser.add_argument('--yMin',
                        default = None,
                        help = 'Minimun value for the Y-axis')

    parser.add_argument('--yMax',
                        default = None,
                        help = 'Maximun value for the Y-axis')

    parser.add_argument('--colorMap',
                        default = 'RdYlBu',
                        help = 'Color map to use for the heatmap. Available values can be seen here: http://www.astro.lsa.umich.edu/~msshin/science/code/matplotlib_cm/',
                        choices = [m for m in matplotlib.cm.datad] # 
                        )

    parser.add_argument('--skipZeros',
                        help = 'Whether regions with only scores of zero should be included or not. Default is to include them. ',
                        action='store_true')

    parser.add_argument('--minThreshold',
                        default = None,
                        type = float, 
                        help = 'Numeric value. Any regions that have a value that is equal or less that this numeric value will be skipped. This is useful to skip, for example, genes where the read count is zero for any of the bins. This could be the result of unmappable areas and can bias the overall results.')

    parser.add_argument('--maxThreshold',
                        default = None,
                        type = float,
                        help = 'Numeric value. Any regions that have a value that is equal or higher that this numeric value will be skipped. The maxThreshold is useful to skip those few regions with very high read counts (e.g. major satellites) that may bias the average values.')


    parser.add_argument('--onePlotPerGroup', 
                        help = 'When the BED file contains groups separated by "#" , the default is to plot the averages for the distinct plots in one plot. If this option is set, each group will get its own plot, stacking on top of each other.',
                        action= 'store_true')

    parser.add_argument('--verbose', 
                        help = 'set to yes, to see Warning messages and other information',
                        action= 'store_true')

    parser.add_argument('--heatmapHeight', 
                        help = 'height in cm. The default for the heatmap height is 25  centimeters. The minimun value is 3cm  and the maximun is 100. The summary plot on top the heatmap has a fixed height of 3.8 cm',
                        type = float,
                        default =  25)

    parser.add_argument('--heatmapWidth',
                        help = 'Width in cm. The default value is 7.5 centimeters. The minimun value is 1cm  and the maximun is 100.',
                        type = float,
                        default =  7.5)

    parser.add_argument('--whatToShow', 
                        help = 'The default is to include a summary plot on top of the heatmap and a heatmap colorbar. Other options are: "plot only", "plot and heatmap", "heatmap only", "colorbar only", "heatmap and colorbar", and the default "plot, heatmap and colorbar" ',
                        choices= ["plot only", "plot, heatmap and colorbar",  "plot and heatmap", "heatmap only", "colorbar only", "heatmap and colorbar"],
                        default =  'plot, heatmap and colorbar')


    parser.add_argument('--scale', 
                        help = 'If set, all values are multiplied by this number',
                        type = float,
                        default =  1)

    parser.add_argument('--numberOfProcessors', '-p',
                        help = 'Number of processors to use. The default is to use the maximun number of procesors',
                        type = int,
                        default = -1,
                        required = False)

    args = parser.parse_args(args)

    # Because of galaxy, the value of this variables is normally
    # set to ''. Therefore this check is needed
    for attr in ['zMin', 'zMax', 'yMax', 'yMin']:
       try:
          args.__setattr__(attr , float(args.__getattribute__(attr)))
#       except ValueError, TypeError:
       except:
          args.__setattr__(attr, None)

    args.heatmapHeight =  args.heatmapHeight if args.heatmapHeight > 3 and args.heatmapHeight <= 100 else 10

    availProc = multiprocessing.cpu_count()
    if args.numberOfProcessors < 1 or args.numberOfProcessors > availProc:
        args.numberOfProcessors = multiprocessing.cpu_count()

    if not matplotlib.colors.is_color_like(args.missingDataColor):
        print "The value {0}  for --missingDataColor is not valid".format(args.missingDataColor)
        exit(1)
    return(args)

def writableFile(string):
    try:
        open(string, 'w').close()
    except:
        msg = "{} file can be opened for writting".format(string)
        raise argparse.ArgumentTypeError(msg)
    return string

def myAverage( valuesArray, avgType='mean' ) :

#    valuesArray[np.isnan(valuesArray)] =  0
    # computes the mean but only for those values that are not nan
    valuesArray = valuesArray[ logical_not(np.isnan(valuesArray)) ]
    if len(valuesArray) > 0:
        if avgType=='mean':
            mean = np.mean(valuesArray)
        elif avgType=='max':
            mean = valuesArray.max()
        elif avgType=='std':
            mean = np.std(valuesArray)
        elif avgType=='min':
            mean = np.min(valuesArray)
        elif avgType=='max':
            mean = np.max(valuesArray)
        else:
            mean = np.median(valuesArray)
    else:
        # if it is not possible to compute the mean, then return zero
        mean = np.nan
    return mean

def matrixAvg(matrix, avgType='mean'):
    avgList = []
    for region in matrix.transpose(): # the default iteration is for the matrix rows
                                  # but I need to iterate over the columns, hence the transpose
        avg =  myAverage(region, avgType)
        avgList.append(avg)
    return avgList

def coverageFromArray(valuesArray, zones, windowWidth):
#    import pdb;pdb.set_trace()
    try:
        valuesArray[0]
    except IndexError:
        print "values array %s, zones %s" % (valuesArray, zones)

    cvgList = []
    start = zones[0][0]
    end = zones[-1][1]
    for zone in zones:
        # the linspace is to get equally spaced positions along the range
        # If the gene is short the sampling regions could overlap,
        # if it is long, the sampling regions would be spaced
        countsList = []
        (posArray, stepSize) = np.linspace(zone[0], zone[1], zone[2], 
                                           endpoint=False,
                                           retstep =True)
        for pos in np.ceil(posArray):
            indexStart = int(pos - start)
            indexEnd   = int(indexStart + windowWidth)
#            indexEnd   = int(indexStart + stepSize)
            countsList.append( myAverage( valuesArray[indexStart:indexEnd] ) )
        cvgList.append( np.array(countsList) )

    return np.concatenate( cvgList )


def coverageFromBam(bamfile, chrom, zones, windowWidth):
    cvgList = []
    start = zones[0][0]
    end = zones[-1][1]
    try:
        valuesArray = np.zeros( end-start )
        for read in bamfile.fetch(chrom, min(0,start), end):
            indexStart = max(read.pos - start, 0)
            indexEnd =   min(read.pos - start + read.qlen, end - start)
            valuesArray[indexStart:indexEnd] += 1
    except ValueError:
        sys.stderr.write("Value out of range for region %s %s %s\n" % (chrom, start, end ) )
        return np.array([0]) # return something inocous

    return coverageFromArray(valuesArray, zones, windowWidth)

def coverageFromBigWig(bigwig, chrom, zones, windowWidth, nansAsZeros = False):
   cvgList = []
   if zones[0][0] < 0:
       valuesArray = np.zeros(zones[-1][1] - zones[0][0])
       valuesArray[:] = np.nan
       valuesArray[abs(zones[0][0]):] = bigwig.get_as_array(chrom, 0, zones[-1][1])
   else:
       valuesArray = bigwig.get_as_array(chrom, zones[0][0], zones[-1][1])
   try:
       valuesArray[0]
   except TypeError :
#       sys.stderr.write("Chromosome '%s' referenced in bed file but absent from bigwig file\n" % (chrom))
      # this error happens when bigwig returns nothing, For example when a chromosome
      # is not nown.
       return None

   except OverflowError as detail:
       print detail
       return None

   # replaces nans for zeros
#   print nansAsZeros
   if nansAsZeros == True:
       valuesArray[np.isnan(valuesArray)] =  0

   return coverageFromArray(valuesArray, zones, windowWidth)

def computeSubMatrix( (regions, matrixCols, args) ):
    # read BAM or scores file
    if args.scoreFileFormat == "bam":
       bamfile = pysam.Samfile(args.scoreFileName)

    if args.scoreFileFormat == "bigwig":
       from bx.bbi.bigwig_file import BigWigFile
       bigwig = BigWigFile( file = open(args.scoreFileName, 'r') )
    # create an empty to store the matrix values
    subMatrix = np.zeros((len(regions) ,matrixCols))
    subMatrix[:] = np.NAN

    j = 0
    subRegions = []
    startTime = time.time()

    for feature in regions:
       # print some information 
       if args.verbose and j % 300 == 0 and j > 0:
          endTime = time.time()
          estimated = ( float( len(regions) - j ) * (endTime - startTime) ) /  j 
          if args.verbose:
              print "%s processing %d (%.1f per sec). Estimated remaining time: %f2" % (multiprocessing.current_process().name, j, j / (endTime - startTime), estimated)
       if args.regionBodyLength > 0 and feature.end - feature.start < args.windowWidth:
           if args.verbose:
               print "Region shorter than window width (%d) %s %s:%s:%s. Skipping..." % ((feature.end - feature.start), feature.name, feature.chrom, feature.start, feature.end)
           continue

       if feature.strand == '-':
          a = args.beforeRegionStartLength / args.windowWidth
          b = args.afterRegionStartLength  / args.windowWidth
          start = feature.end
          end = feature.start
       else:
           b = args.beforeRegionStartLength / args.windowWidth
           a = args.afterRegionStartLength / args.windowWidth
           start = feature.start
           end = feature.end
       # build zones
       if args.regionBodyLength > 0:
          zones = [(feature.start - b * args.windowWidth, feature.start, b ),
                   (feature.start, feature.end -args.regionBodyLength / args.windowWidth, args.regionBodyLength / args.windowWidth),
                   (feature.end, feature.end + a * args.windowWidth, a)]
       elif args.referencePoint == 'TES': # around TES
          zones = [(end - b * args.windowWidth, end, b ),
                   (end, end + a * args.windowWidth, a )]
       elif args.referencePoint == 'center': # at the region center
          middlePoint = feature.start + (feature.end - feature.start)/2
          zones = [(middlePoint - b * args.windowWidth, middlePoint, b ),
                   (middlePoint, middlePoint + a * args.windowWidth, a )]
       else: # around TSS
          zones = [(start - b * args.windowWidth, start, b ),
                   (start, start + a * args.windowWidth, a )]

       if feature.start - b * args.windowWidth < 0:
           if args.verbose:
               print "region too close to chromosome start for %s %s:%s:%s. " % (feature.name, feature.chrom, feature.start, feature.end)

       coverage = None
       if args.scoreFileFormat == "bam":
          if feature.chrom not in bamfile.references:
              if args.verbose:
                  sys.stderr.write("Skipping region located at unknown chromosome for %s %s:%s-%s.\n" % (feature.name, feature.chrom, feature.start, feature.end))
              continue
          coverage = coverageFromBam(bamfile,feature.chrom, zones, args.windowWidth)

       if args.scoreFileFormat == "bigwig":
          coverage = coverageFromBigWig(bigwig, feature.chrom, zones, args.windowWidth, args.missingDataAsZero)
          if coverage == None:
               if args.verbose:
                   sys.stderr.write("No scores defined for region %s %s:%s-%s.\n" % (feature.name, feature.chrom, feature.start, feature.end))
               coverage = np.zeros(matrixCols)
               if not args.missingDataAsZero:
                   coverage[:] = np.nan
               continue
       try:
          temp = coverage.copy()
          temp[np.isnan(temp)] =  0
          totalScore = np.sum( temp )
          
       except:
           if args.verbose:
               sys.stderr.write("No scores defined for region %s %s:%s-%s. Skipping...\n" % (feature.name, feature.chrom, feature.start, feature.end))
           coverage = np.zeros(matrixCols)
           if not args.missingDataAsZero:
               coverage[:] = np.nan
           # to induce skipping if zero regions are omited this is set to zero
           totalScore = 0

       if totalScore == 0:
          if args.skipZeros == True :
              if args.verbose:
                  sys.stderr.write("Skipping region with all scores equal to zero for %s %s:%s-%s.\n" % (feature.name, feature.chrom, feature.start, feature.end))
              continue
          elif args.verbose:
             sys.stderr.write( "Warning: All values are zero for %s %s:%s-%s.\n" % (feature.name, feature.chrom, feature.start, feature.end) )
             sys.stderr.write( "add --skipZeros to exclude such regions\n" )

       if args.minThreshold and coverage.min() <= args.minThreshold:
           continue
       if args.maxThreshold and coverage.max() >= args.maxThreshold:
           continue
       """
       if np.isnan(coverage).any():
           print "All values are zero or some are nan for %s %s:%s-%s. Skipping..." % (feature.name, feature.chrom, feature.start, feature.end)
           continue
       """
       if args.scale !=1:
           coverage = args.scale * coverage

       if feature.strand == "-":
          subMatrix[j,:] = coverage[::-1]
       else:
          subMatrix[j,:] = coverage
       # testing feature, unfinished
       # TODO: sort regions according to feature length
       # add option to argparse
#       if args.skipEndRegion == True:
       if 2 == 1 and args.regionBodyLength == 0:
           regionLength = feature.end - feature.start
           if regionLength < args.afterRegionStartLength:
               toNans = a - int(regionLength / args.windowWidth)
               subMatrix[j,-toNans:] = np.NAN
       subRegions.append(feature)
       j += 1

    # remove empty rows
    subMatrix = subMatrix[0:j,:]
    if len(subRegions) != len(subMatrix[:,0]):
        print "regions lengths do not match"
    return ( subMatrix, subRegions )


class intervalWrapper():
   # class to create a simple object that can be pickled and send to workers
   def __init__(self, genomicInterval):
      self.chrom = genomicInterval.chrom
      self.start = genomicInterval.start
      self.end  = genomicInterval.end
      self.strand = genomicInterval.strand
      try:
         self.name = genomicInterval.fields[3]
      except IndexError:
         self.name = "No name"


def getRegionsAndGroups(regionsFileName, onlyMultiplesOf=1):
    # reads a bed file containing the position
    # of genomic intervals
    # In case is hash sign '#' is found in the
    # file, this is considered as a delimiter
    # to split the heatmap into groups
    
    regions = []
    regionsDict = OrderedDict()
    regionGroups = [(0,'')]

    prevInterval = None
    duplicates = 0
    totalIntervals = 0
    includedIntervals = 0
    # drop some lines
    for ginterval in GenomicIntervalReader( open(regionsFileName, 'r').readlines() ):
       totalIntervals +=1
       if ginterval.__str__()[0] == '#':
           if includedIntervals>1 and  includedIntervals - regionGroups[-1][0] > 1:
               label = ginterval.__str__()[1:]
               newLabel = label
               if label in regionsDict.keys():
                   # loop to find a unique label name
                   i = 0
                   while True:
                       i += 1
                       newLabel = label +"_r" + str(i)
                       if newLabel not in regionsDict.keys():
                           break

               regionsDict[newLabel] = regions[:]
               regions = []
           continue
       # if the list of regions is to big, only consider a fraction of the data
       if totalIntervals % onlyMultiplesOf !=0:
           continue
       # skip regions that have the same position as the previous.
       # This assumes that the regions file given is sorted 
       if prevInterval and prevInterval.chrom == ginterval.chrom and \
               prevInterval.start == ginterval.start and \
               prevInterval.end == ginterval.end:
           if args.verbose:
               print "Gene in same region already included:  %s %s:%s-%s. Skipping" % (ginterval.fields[3], ginterval.chrom, ginterval.start, ginterval.end)

           duplicates += 1
           continue
       else:
           prevInterval = ginterval

       regions.append( intervalWrapper(ginterval) )
       includedIntervals += 1

    if len(regions):
        regionsDict[args.regionsLabel] = regions

    if args.verbose:
        print "%d (%.2f) regions covering the exact same interval were found" % \
            (duplicates, 
             float(duplicates) *100 / totalIntervals)

    return regionsDict


def saveTabulatedValues(heatmapMatrixDict, regionsDict, args):
    regionLength = 0 if args.regionBodyLength == 0 else args.regionBodyLength
    bin = range(-args.beforeRegionStartLength, regionLength + args.afterRegionStartLength, args.windowWidth)

    avgDict = OrderedDict()
    stdDict = OrderedDict()

    for label, heatmapMatrix in heatmapMatrixDict.iteritems():
        avgDict[label] = matrixAvg(heatmapMatrix, args.averageType)
        stdDict[label] = matrixAvg(heatmapMatrix, 'std')
        
    
    args.outFileNameData.write('#bin No.\t{}\n'.format(" mean\t std\t".join(avgDict.keys())))
    for j in range(0, len(avgDict[avgDict.keys()[0]])):
        args.outFileNameData.write('{}\t'.format(bin[j]))
        for label in heatmapMatrixDict.keys():
            args.outFileNameData.write('{}\t{}\t'.format(avgDict[label][j],stdDict[label][j]))
        args.outFileNameData.write('\n')

    args.outFileNameData.close()


################ main ###########


def main(args):
    r"""
    >>> import filecmp
    >>> import os
    >>> args = parseArguments("-R /data/projects/ramirez/tools/heatmapper/test/test2.bed -S /data/projects/ramirez/tools/heatmapper/test/test.bw -F bigwig -b 100 -a 100 -m 0 --outFileName /tmp/_test.png --outFileNameData /tmp/_test.tab -w 1 -at mean ".split())
    >>> main(args)
    >>> filecmp.cmp('/data/projects/ramirez/tools/heatmapper/test/master.tab', '/tmp/_test.tab')
    True
    >>> os.remove('/tmp/_test.png')
    >>> os.remove('/tmp/_test.tab')
    """
    
    # get the number of regions in the BED file. This number is needed to initialize the empty
    # matrix holding the scores.
#    numRegions = sum(1 for line in open(args.regionsFileName, 'r'))

    # fail early
    try:
        for file in [args.regionsFileName, args.scoreFileName ]:
            if os.stat( file ).st_size == 0:
                print "file %s is empty" % ( file );
                exit()
            
    except (IOError, OSError) as detail:
        print detail
        exit()

    if args.regionBodyLength > 0 and args.regionBodyLength % args.windowWidth > 0:
       print "Length of body to print has to be a multiple of --windowWidth"
       exit()

    # the beforeRegionStartLength is extended such that length is a multiple of windowWidth 
    if args.afterRegionStartLength % args.windowWidth > 0:
       print "Length of region after the body has to be a multiple of --windowWidth"
       exit()

    if args.beforeRegionStartLength % args.windowWidth > 0:
       print "Length of region before the body has to be a multiple of --windowWidth"
       exit()

    if args.afterRegionStartLength == 0:
        # this case is when the whole gene wants to be plotted
        
        # we need to find the longest gene to compute the matrix size
        # and we keep the lengths to later order accordingly
        intervalLength = []
        #for ginterval in GenomicIntervalReader( open(args.regionsFileName, 'r').readlines() ):
        

    # determine the number of matrix columns based on the lengths given by the user
    matrixCols = ( (args.afterRegionStartLength + 
                    args.beforeRegionStartLength + args.regionBodyLength) /  
                    args.windowWidth )

    regionsDict = getRegionsAndGroups(args.regionsFileName)
    heatmapMatrixDict = OrderedDict()
    matrixAvgsDict = OrderedDict()

    for label, regions in regionsDict.iteritems():
        # args to pass to the multiprocessing workers 
        mp_args = []

        # prepare chunks of  400 regions to send to workers.
        for index in range(0,len(regions), 400):
            index_end = min(len(regions), index + 400 )
            mp_args.append((regions[index:index_end], matrixCols, args))

        if len(mp_args) > 1 and args.numberOfProcessors > 1:
            if args.verbose:
                print "'{}' total workers: {}, using {} processors ".format(label,len(mp_args), args.numberOfProcessors)
            pool = multiprocessing.Pool(args.numberOfProcessors)
            res = pool.map( computeSubMatrix, mp_args )
    #        res = pool.asyn( computeSubMatrix, mp_args ).get(9999999)
        else:
            res = map( computeSubMatrix, mp_args )
    
        # each worker in the pools returns a tuple containing 
        # the submatrix data and the regions that correspond to the
        # submatrix

        # merge all the submatrices into heatmapMatrix
        heatmapMatrix = np.concatenate( [ r[0] for r in res], axis=0)
        # merge all valid regions
        regionList = np.concatenate( [ r[1] for r in res], axis=0)
        if len(regionList) == 0:
            print "Error: could not compute values for any of the regions"
            exit()

        # sort the matrix using the average of values
        matrixAvgs = np.array([ myAverage( region, args.averageType ) for region in heatmapMatrix ])
        SS = matrixAvgs.argsort()
        if args.sortRegions == 'ascend':
            heatmapMatrix = heatmapMatrix[SS,:]
            regionList = regionList[SS]
            matrixAvgs = matrixAvgs[SS]
        elif args.sortRegions == 'descend':
            heatmapMatrix = heatmapMatrix[SS[::-1],:]
            regionList = regionList[SS[::-1]]
            matrixAvgs = matrixAvgs[SS[::-1]]

        heatmapMatrixDict[label] = heatmapMatrix
        regionsDict[label] = regionList
        matrixAvgsDict[label] = matrixAvgs

    if args.outFileNameMatrix:
        # merge all group matrices
        heatMapMatrix = np.concatenate( [ x for x  in heatmapMatrixDict.values()] )
        np.savetxt(args.outFileNameMatrix, heatmapMatrix, fmt="%.3g")
        if args.verbose:
            print "output matrix file saved as {}".format(args.outFileNameMatrix)

    if args.outFileNameData:
        saveTabulatedValues(heatmapMatrixDict, regionsDict, args)

    if args.outFileSortedRegions:
       for label, regions in regionsDict.iteritems():
           j = 0
           for region in regions:
              args.outFileSortedRegions.write('%s\t%s\t%s\t%s\t%s\t%s\n' % ( region.chrom, 
                                                                             region.start, 
                                                                             region.end, 
                                                                             region.name, 
                                                                             matrixAvgsDict[label][j], 
                                                                             region.strand) )
              j += 1
           args.outFileSortedRegions.write('#{}\n'.format(label))
       args.outFileSortedRegions.close()

    args.referencePoint = args.startLabel if args.referencePoint == "TSS" else args.endLabel

    plotMatrix(heatmapMatrixDict, regionsDict, args.windowWidth,
               args.beforeRegionStartLength,
               args.afterRegionStartLength,
               args.regionBodyLength,
               args.outFileName, args.outFileFormat, 
               args.colorMap, args.missingDataColor, args.plotTitle, 
               args.xAxisLabel, args.yAxisLabel,  args.regionsLabel,  
               args.zMin, args.zMax,
               args.yMin, args.yMax, 
               args.averageType,
               args.referencePoint,
               args.startLabel,
               args.endLabel,
               args.heatmapHeight,
               args.heatmapWidth,
               args.onePlotPerGroup,
               args.whatToShow,
               args.sortRegions)



if __name__ == "__main__":
#   args = parseArguments("-R /data/projects/ramirez/NSL/data/Dm530.genes.gff -S /data/projects/muehlpfordt/2011-Autumn_NSL/results/2012-01-18/log2FC_DESeq_NSL3kd-GFPkd.bg -b 500 -a 1500 --plotTitle deltaNSL3kd --outFileName deltaNSL3kd_autumn.png --outFileNameData deltaNSL3kd --outFileFormat png -w 50 --colorMap hot  -min -2 -max 1.1".split())
#   args = parseArguments("-S ../src/corrected_plust.bw -R ../results/test_input_bias_03_02_2012/Dm530_genes.bed -F bigwig --outFileName corrected_plust.png -a 1500 -b 1500 -w 50 --outFileNameData corrected_plust.tab".split())
   args = parseArguments()

   main(args)
