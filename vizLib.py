import math
import matplotlib.pyplot as pylab
import scipy.cluster.hierarchy as sch
import qMS
import numpy
import matplotlib
from mpl_toolkits.mplot3d import Axes3D

def setRcs(scale=None, legendScale=None, tickScale=1):
    """setRcs sets a series of rc params for matplotlib tomake decent looking plots. Also some useful def
    in there to cutomize things.

    :param scale: the default font scale
    :type scale: float
    
    :returns:  nothing, edits the rcparam file.
            
    """
    if scale is None:
        scale = 12
    if legendScale is None:
        legendScale = scale*0.85
    
    defaultFont = {'family' : 'sans-serif',
                   'variant' : 'normal',
                   'weight' : 400,
                   'size' : scale*1}

    axisFont = {'titlesize' : scale*1.5,
                'labelsize' : scale*1.25}

    xAxisTicks = {'major.size' : 8.0*tickScale,
                  'minor.size' : 4.0*tickScale,
                  'major.width' : 1.0*tickScale,
                  'minor.width' : 1.0*tickScale,
                  'labelsize' : scale*1,
                  'minor.pad' : 3,
                  'major.pad' : 3}

    yAxisTicks = {'major.size' : 8.0*tickScale,
                  'minor.size' : 4.0*tickScale,
                  'major.width' : 1.0*tickScale,
                  'minor.width' : 1.0*tickScale,
                  'labelsize' : scale*1,
                  'minor.pad' : 3,
                  'major.pad' : 3}

    legend = {'fancybox' : True,
              'numpoints' : 1,
              'fontsize' : legendScale,
              'borderaxespad' : 1}

    matplotlib.rc('font', **defaultFont)
    matplotlib.rc('axes', **axisFont)
    matplotlib.rc('xtick', **xAxisTicks)
    matplotlib.rc('ytick', **yAxisTicks)
    matplotlib.rc('legend', **legend)

    matplotlib.rc('lines', linewidth=2)
    matplotlib.rc('axes', linewidth=3)
    ##AXISNAME.xaxis.set_ticks_position('bottom')
    ##AXISNAME.yaxis.set_ticks_position('left')
    ##AXISNAME.xaxis.labelpad = 2

def plotStatsDict(statsDict, name='', proteins=None, offset=0.0, markerSize=12, color='#e31a1c', yMax=1.5, 
                  median=False, figSize = (22,5), noFill=False, mew=1):
    """plotStatsDataStruct plots the contents of a stats dictionary. proteins to be plotted are 
        listed in the non-redundent list, proteins. The data is in statsDict, the name is in in name.
        Decent colors are red (['#ae2221', '#d72c2b', '#e78180']) and blue (['#25557d', '#3170a4', '#5696cc'])

    :param statsDict: a dictionary (easily created by qMS.calcStatsDict)
    :type statsDict: dictionary with keys of proteins names and values of numpy arrays of calcValues
    :param name: the name of the dataset
    :type name: string
    :param proteins: a non-redudent list of the protein names to use (should be the IDs in dataByProtein).
        If none is given, all of the keys will be plotted.
    :type proteins: list
    :param offset: where to center the point (x axis), scales 0 (left edge; default) to 1.0 (right edge)
    :type offset: float
    :param markerSize: size of the dots (default = 12)
    :type markerSize: int
    :param color: color for the dataset (default #e31a1c)
    :type color: color
    :param yMax: the max value for the y axis
    :type yMax: float
    :param median: bool to plot only the median values
    :type median: bool

    :returns:  a pyplot axis with the data plotted
    
        
    """
    
    if proteins is None:
        #proteins = sorted(statsDict.keys())
        proteins = qMS.sort_nicely(statsDict.keys())

    xAxis = range(1,len(proteins)+1)
    fig = pylab.figure(figsize=figSize)
    ax = fig.add_subplot(111)
    xs = []
    ys = []
    for x in xAxis:
        p = proteins[x-1]
        if p in statsDict.keys():
            if median:
                xs.append(x+offset)
                ys.append(numpy.median(statsDict[p]))
            else:
                for v in statsDict[p]:
                    xs.append(x+offset)
                    ys.append(v)

    pylab.grid(b=True, which='major', color='grey', linestyle='--', axis='y', linewidth=1.5, alpha=0.5)
    pylab.grid(b=True, which='major', color='grey', linestyle='-', axis='x', linewidth=1.5, alpha=0.75)
    if noFill:
        ax.plot(xs, ys, 'o', color='none', markeredgecolor=color, mew=mew, markersize=markerSize, label=name)
    else:
        ax.plot(xs, ys, 'o', color=color, markersize=markerSize, label=name)

    pylab.xticks(xAxis, [item for item in proteins], rotation=45, size=15)
    pylab.xlim(1, len(proteins)+1)
    ####################################
    pylab.yticks([0,yMax/5.0, 2*yMax/5.0, 3*yMax/5.0, 4*yMax/5.0, yMax], size=15)
    ####################################
    pylab.ylim(-0.05, yMax)
    ax.set_axisbelow(True)
    return ax

def addStatsDictToPlot(statsDict, ax, name='', offset=0.0, markerSize=12, color='#377db8', median=False, noFill=False, mew=1):
    """addStatsDataStructToPlot adds the contents of a stats dictionary to an existing plot. ONLY PROTEINS
        PRESENT IN THE ORIGINAL PLOT WILL BE PLOTTED IN THE NEW PLOT
        The data is in statsDict, the axis to add to is in ax, the name is in in name.
        Decent colors are red (['#ae2221', '#d72c2b', '#e78180']) and blue (['#25557d', '#3170a4', '#5696cc'])

    :param statsDict: a dictionary (easily created by calcStatsDict)
    :type statsDict: dictionary
    :param ax: the pyplot axis to modify
    :type ax: pyplot axis
    :param name: the name of the dataset
    :type name: string
    :param offset: where to center the point (x axis), scales 0 (left edge; default) to 1.0 (right edge)
    :type offset: float
    :param markerSize: size of the dots (default = 12)
    :type markerSize: int
    :param color: color for the dataset (default #e31a1c)
    :type color: color
    :param median: bool to plot only the median values
    :type median: bool
    :param noFill: a bool if you want open circles
    :type noFill: bool
    :param mew: a float of the marker edge width
    :type mew: float
    
    :returns:  a pyplot axis with the data plotted
        
    """

    a = ax.get_xmajorticklabels()
    proteins = [t.get_text() for t in a]

    xAxis = range(1,len(proteins)+1)
    xs = []
    ys = []
    for x in xAxis:
        p = proteins[x-1]
        if p in statsDict.keys():
            if median:
                xs.append(x+offset)
                ys.append(numpy.median(statsDict[p]))
            else:
                for v in statsDict[p]:
                    xs.append(x+offset)
                    ys.append(v)

    if noFill:
        ax.plot(xs, ys, 'o', color='none', markeredgecolor=color, mew=mew, markersize=markerSize, label=name)
    else:
        ax.plot(xs, ys, 'o', color=color, markersize=markerSize, label=name)

    ax.set_axisbelow(True)
    return ax

def makePlotWithFileList(isoFileList, numerator, denominator, subunits=None, normProtein=None, yMax=1.5, 
                         median=False, names=None, colors=None, figSize=(22,5), markerSize=None, noFill=False, mew=1):
    """makePlotWithFileList is a  helper function that plots massage-style data from a list of files

    :param isoFileList: a list of the files to be ploted (shoudl be full path to _iso.csv files)
    :type isoFileList: list of strings
    :param numerator: strings of the keys in the numerator (ampu, ampl, amps)
    :type numerator: list of strings
    :param denominator: strings of the keys in the denominator (ampu, ampl, amps)
    :type denominator: list of strings
    :param subunits: a list of the proteins to be plotted - strings must match the keys of the 'proteins' field in _iso.csv
    :type subunits: list
    :param normProtein: string of the protein to normalize to for all datasets (defaults to None)
    :type normProtein: string
    :param yMax: float of maximum y value
    :type yMax: float
    :param median: bool to plot only the median values
    :type median: bool
    :param names: a list of the names to be listed in the legend; must be same length as isoFileList
    :type names: list of strings
    :param colors: a list of the colors to be used in plotting, again must be same length as isoFileList
    :type colors: list of strings
    :param figSize: the figure size (list - (10,10))
    :type figSize: list of floats
    :param noFill: a bool if you want open circles
    :type noFill: bool
    :param mew: a float of the marker edge width
    :type mew: float

    :returns: the plotted axis
    
    """
    if names is None:
        names = isoFileList
        
    namesList = [(isoFileList[i], names[i]) for i in range(len(isoFileList))]
    allStats = qMS.multiStatsDict(isoFileList, numerator, denominator, normalization=1.0, offset=0.0, normProtein=normProtein)
    
    return makePlotWithStatsDictDict(allStats, subunits=subunits, yMax=yMax,
                                     median=median, namesList=namesList, colors=colors, figSize=figSize, markerSize=markerSize,
                                     noFill=noFill, mew=mew)

def makePlotWithStatsDictDict(allStats, subunits=None, yMax=1.5, 
                         median=False, namesList=None, colors=None, figSize=(22,5), markerSize=None, noFill=False, mew=1):
    """makePlotWithStatsDictDict is a  helper function that plots massage-style data from a dict of statsDicts

    :param allStats: dictionary of stats dictionaries (returned by multiStatsDict)
    :type allStats: dict of statsDicts
    :param subunits: a list of the proteins to be plotted - strings must match the keys of the 'proteins' field in _iso.csv
    :type subunits: list
    :param yMax: float of maximum y value
    :type yMax: float
    :param median: bool to plot only the median values
    :type median: bool
    :param namesList: a list of lists where each sublist is a pair with the key for the 
        allStats in [0] and the name for the legend in [1]. typically, the [0] would be a full
        path.
    :type namesList: list of pairs of strings
    :param colors: a list of the colors to be used in plotting, again must be same length as allStats
    :type colors: list of strings
    :param figSize: the figure size (list - (10,10))
    :type figSize: list of floats
    :param markerSize: the size of the dots
    :type markerSize: float
    :param noFill: a bool if you want open circles
    :type noFill: bool
    :param mew: a float of the marker edge width
    :type mew: float

    :returns: the plotted axis
    
    """

    if namesList is None:
        namesList = [(k, k) for k in allStats.keys()]
    
    if colors is None:
        colors = [pylab.cm.jet(float(i)/float(len(allStats))) for i in range(len(allStats))]
    
    offsets = float(len(allStats)+1)
    if markerSize is None:
        markerSize = (20.0/offsets)+4

    ax = plotStatsDict( allStats[namesList[0][0]], name=namesList[0][1], proteins=subunits, \
                        offset=1.0/offsets, markerSize=markerSize, yMax=yMax, median=median, 
                        color=colors[0], figSize=figSize, noFill=noFill, mew=mew)
        
    for i in range(1,len(namesList)):
        ax = addStatsDictToPlot(allStats[namesList[i][0]], ax, name=namesList[i][1], \
                                offset=(1.0/offsets)*(i+1.25), markerSize=markerSize, median=median, 
                                color=colors[i], noFill=noFill, mew=mew)
        ax.set_axisbelow(True)
    return ax

def plotCompData(xdat, ydat, proteins, title=None, xlabel='dat1', ylabel='dat2', xMax=1.5, yMax=1.5, figSize=(10,10), saveFile=None):
    """plotCompData is a  makes a scatter plot out of data to be compared. Useful for the pool size stuff. Only the median values are plotted.

    :param xdat: a dictionary of values to be plotted on the x axis (keys are protein names, values are numpy arrays of values)
    :type xdat: a dictionary of numpy arrays (keys=protein names)
    :param ydat: a dictionary of values to be plotted on the y axis (keys are protein names, values are numpy arrays of values)
    :type ydat: a dictionary of numpy arrays (keys=protein names)
    :param proteins: a list of the proteins to be plotted - strings must match the keys of the 'proteins' field in xdat, ydat
    :type proteins: list
    :param title: the plot title
    :type title: string
    :param xlabel: the plot xlabel
    :type xlabel: string
    :param ylabel: the plot ylabel
    :type ylabel: string
    :param xMax: float of maximum x value
    :type xMax: float
    :param yMax: float of maximum y value
    :type yMax: float
    :param figSize: the figure size (list - (10,10))
    :type figSize: list of floats
    :param saveFile: the full path of the file to be saved to
    :type saveFile: string

    :returns: the plotted axis
    
    """
    x = [numpy.median(xdat[i]) for i in proteins]
    y = [numpy.median(ydat[i]) for i in proteins]    
    scat = pylab.figure(figsize=figSize)
    scatAx = scat.add_subplot(111)    
    scatAx.scatter(x,y, c='b', s=150)
    scatAx.set_title(title)
    scatAx.set_xlabel(xlabel)
    scatAx.set_ylabel(ylabel)
    scatAx.set_xlim([-0.1,xMax])
    scatAx.set_ylim([-0.1,yMax])
    scatAx.set_xticks([0,xMax/5,xMax/5*2,xMax/5*3,xMax/5*4,xMax])
    scatAx.set_yticks([0,yMax/5,yMax/5*2,yMax/5*3,yMax/5*4,yMax])
    scatAx.yaxis.tick_left()
    scatAx.xaxis.tick_bottom()
    for prot, xl, yl in zip(proteins, x, y):
        scatAx.annotate(str(prot[4:]), xy = (float(xl), float(yl)), xytext = (15,15), textcoords = 'offset points', arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))
    pylab.tight_layout()
    scatAx.plot(numpy.linspace(0, 10), numpy.linspace(0,10))
    return scatAx

def proteinScatterPlot(yDataDict, xData, xMin=0, xMax=None, yMin=0, yMax=10,
                       title=None, xLabel=None, yLabel=None, colors=None, 
                       figSize=(10,10), markerSize=10, legend=False,
                        linestyle=None, xTicks=None, legendLoc='upper left', legendCols=2):
    """proteinScatterPlot is a  makes a scatter plot out of ydata with a fixed x. Useful for the standard curve stuff or gradients.

    :param yDataDict: a dictionary of values to be plotted on the y axis (keys are protein names, values are numpy arrays of values)
    :type yDataDict: a dictionary of numpy arrays (keys=protein names)
    :param xData: a list of the xvalues (same for all proteins)
    :type xData: a list of floats
    :param xMin: float of minimum x value
    :type xMin: float
    :param yMin: float of minimum y value
    :type yMin: float
    :param xMax: float of maximum x value
    :type xMax: float
    :param yMax: float of maximum y value
    :type yMax: float
    :param title: the plot title
    :type title: string
    :param xlabel: the plot xlabel
    :type xlabel: string
    :param ylabel: the plot ylabel
    :type ylabel: string
    :param colors: a list of colors to be used
    :type colors: list of colors
    :param figSize: the figure size (list - (10,10))
    :type figSize: list of floats
    :param markerSize: the size of the dots
    :type markerSize: float
    :param legend: a bool of whether to have a legend
    :type legend: bool
    :param linestyle: the connector line style
    :type linestyle: the connector line style
    :param xTicks: a list of the xTicks
    :type xTicks: a list of floats
    :param legendloc: a string of where to locate the legend
    :type legendloc: string
    :param legendCols: the number of columns in the legend
    :type legendCols: int

    :returns: the plotted axis
    
    """    
    if xMax is None:
        xMax = max(xData)
    if colors is None:
        colors = [pylab.cm.jet(float(i)/float(len(yDataDict))) for i in range(len(yDataDict))]
    scat = pylab.figure(figsize=figSize)
    scatAx = scat.add_subplot(111)
    for i,p in enumerate(qMS.sort_nicely(yDataDict.keys())):
        scatAx.scatter(xData, yDataDict[p], c=colors[i], s=markerSize, label=p)
        if not (linestyle is None):
            scatAx.plot(xData, yDataDict[p], c=colors[i], markersize=markerSize, linestyle=linestyle)
    scatAx.set_title(title, multialignment='center')
    scatAx.set_xlabel(xLabel)
    scatAx.set_ylabel(yLabel)
    scatAx.set_xlim([xMin,xMax])
    scatAx.set_ylim([-0.1,yMax])
    if xTicks is None:
        scatAx.set_xticks([0,xMax/4,xMax/4*2,xMax/4*3,xMax])
    else:
        scatAx.set_xticks(xTicks)
    scatAx.set_yticks([0,yMax/4,yMax/4*2,yMax/4*3,yMax])
    if legend:
        pylab.legend(loc=legendLoc, ncol=legendCols, scatterpoints=1)
    scatAx.yaxis.tick_left()
    scatAx.xaxis.tick_bottom()
    pylab.tight_layout()
    
    return scatAx

    

def plotMSSpectra3D(listOfFilesToPlot, listOfNames=None, listOfColors=None, gridLines=False, yMin=0.5, yMax=2.5, 
                    legend=True, normalizeToN15=False, subtractRef=None):
    """plotMSSpectra3D is a  makes a 3d plot of MS spectra

    :param listOfFilesToPlot: a list of the spectra to be plotted (full paths)
    :type listOfFilesToPlot: list of strings
    :param listOfNames: a list of the names for each dataset
    :type listOfNames: list of strings
    :param listOfColors: a list of colors to be used
    :type listOfColors: list of colors
    :param gridLines: a bool of whether to draw gridlines
    :type gridLines: bool
    :param yMin: float of minimum y value
    :type yMin: float
    :param yMax: float of maximum y value
    :type yMax: float
    :param legend: a bool of whether to have a legend
    :type legend: bool
    :param normalizeToN15: a bool of whether to nomrmalize each plot the N15 maximum
    :type normalizeToN15: bool
    :param subtractRef: a int pointing to which is the refernece spectra that should be subtracted from each 
    :type subtractRef: int

    :returns: the plotted axis
    
    """  
    if listOfNames==None:
        listOfNames = listOfFilesToPlot
    if listOfColors==None:
        listOfColors = [pylab.cm.jet(float(i)/float(len(listOfFilesToPlot))) for i in range(len(listOfFilesToPlot))]
    
    fig = pylab.figure()
    ax = fig.add_subplot(111, projection='3d')

    yTotal = len(listOfFilesToPlot)
    top = 0.0

    if not (subtractRef is None):
        [bhah, zsRef, blah] = qMS.readMSSpectraFile(listOfFilesToPlot[subtractRef])
        zNorm = max(zsRef[len(zsRef)/2:])
        zsRef = numpy.array(zsRef)/zNorm
    
    for i,f in enumerate(listOfFilesToPlot):
        [xs, zs, name] = qMS.readMSSpectraFile(f)
        ys = [yTotal-i]*len(xs)
        if normalizeToN15:
            zNorm = max(zs[len(zs)/2:])
            zs = numpy.array(zs)/zNorm
        if not (subtractRef is None):
            zNorm = max(zs[len(zs)/2:])
            zs = numpy.array(zs)/zNorm
            zs = zs-zsRef
            zs = zs*zNorm
            xs = xs[:len(xs)/2]
            ys = ys[:len(ys)/2]
            zs = zs[:len(zs)/2]
        ax.plot(numpy.array(xs),numpy.array(ys),numpy.array(zs), color=listOfColors[i], lw=1.5, label=listOfNames[i])
        top = max([top, float(max(zs))])


    ax.w_xaxis.pane.set_visible(False)
    ax.w_yaxis.pane.set_visible(False)
    ax.w_zaxis.pane.set_visible(False)

    if gridLines:        
        ax.w_xaxis.gridlines.set_linewidth(2)
        ax.w_yaxis.gridlines.set_linewidth(2)
        ax.w_zaxis.gridlines.set_linewidth(2)
    
    else:
        ax.w_xaxis.gridlines.set_visible(False)
        ax.w_yaxis.gridlines.set_visible(False)
        ax.w_zaxis.gridlines.set_visible(False)

    [i.set_linewidth(2) for i in ax.w_xaxis.get_ticklines()]
    [i.set_linewidth(2) for i in ax.w_yaxis.get_ticklines()]
    [i.set_linewidth(2) for i in ax.w_zaxis.get_ticklines()]

    ax.w_xaxis.line.set_linewidth(2)
    ax.w_yaxis.line.set_linewidth(2)
    ax.w_zaxis.line.set_linewidth(2)
           
    ax.set_zticks([0, top/3, 2*top/3, top])
    ax.set_zlim3d([0, top])
    ax.set_ylim3d(yMin, yMax)
    ax.set_yticks(range(1,yTotal+1))

    ax.set_xlabel("mass")
    ax.set_zlabel("intensity")

    ax.view_init(10, -67)
    if legend:
        pylab.legend()
    
    return ax

'''
def drawHeatMap(xdat, name="unnamed", colors=pylab.cm.RdBu, dendro=False, protColors=None, 
                cIndex=None, km=None, nameList=None, scale=None, saveName=None):
    """drawHeatMap produces a colored heatmap in a new figure window

    :param xdat: a data object (must contain fields 'data', 'fractions', 'proteins')
    :type xdat: dict
    :param colors: a color scale (a cmap)
    :type colors: cmap
    :param name: figure name and title
    :type name: str.
    :param dendro: a boolean to draw the dendrogram on the left
    :type dendro: bool.
    :param protColors: a color map used to label the protein names with group colors
    :type protColors: cmap
    :param cIndex: a list of groupIds for the proteins
    :type cIndex: list
    :param km: if present, will draw the kmeans cluster profiles at the top of the figure- input is a 2d-matrix - rowVectors for each centroid, each column is a fraction
    :type km: matrix
    :returns:  int -- the return code.
    :raises: AttributeError, KeyError
    :returns: a figure object

    """

    data = xdat['data']
    if nameList is None:
        nameList = xdat['fractions']
    ls = xdat['proteins']
    fig = pylab.figure()

    offset = 0.1
    yStart = 0.1
    yLength = 0.8
    if dendro:
        xStart = 0.35
        xLength = 0.5
    else:
        xStart = 0.1
        xLength = 0.8
    
    figAxes = heatMapAxes(data, dims = [xStart, yStart, xLength, yLength], colors=colors, columns=nameList, rows=ls, protColors=protColors, cIndex=cIndex, fig=fig, colorBar=True)
    figAxes.set_title(name)
    if dendro:
        ax2Data = fig.add_axes([offset, offset, xLength-0.3, yLength])
        sch.dendrogram(xdat['rightDendro'], orientation='right', color_threshold=1)
    if not (scale is None):
        figAxes.text(1.15, 0.5, scale, rotation=270, verticalalignment="center", 
                     horizontalalignment="center", fontsize=18, transform=figAxes.transAxes)
   if not (saveName is None):
        pylab.savefig(saveName)
    return fig
'''

def drawHeatMap(xdat, name="unnamed", colors=pylab.cm.Reds, dendro=False, protColors=None, cIndex=None, km=None,
                nameList=None, scale=None, saveName=None, colorBar=False, figSize=(6,6)):
    """drawHeatMap produces a colored heatmap in a new figure window

    :param xdat: a data object (must contain fields 'data', 'fractions', 'proteins')
    :type xdat: dict
    :param colors: a color scale (a cmap)
    :type colors: cmap
    :param name: figure name and title
    :type name: str.
    :param dendro: a boolean to draw the dendrogram on the left
    :type dendro: bool.
    :param protColors: a color map used to label the protein names with group colors
    :type protColors: cmap
    :param cIndex: a list of groupIds for the proteins
    :type cIndex: list
    :param km: if present, will draw the kmeans cluster profiles at the top of the figure- input is a 2d-matrix - rowVectors for each centroid, each column is a fraction
    :type km: matrix
    :returns:  int -- the return code.
    :raises: AttributeError, KeyError
    :returns: a figure object

    """

    data = xdat['data']
    if nameList is None:
        nameList = xdat['fractions']

    #fractions = xdat['fractions']
    proteins = xdat['proteins']
    fig = pylab.figure(figsize=figSize)
    if not (name is None):
        fig.suptitle(name)
    ##Draw heatmap
    offset = 0.1
    if dendro:
        xStart = 0.35
        xLength = 0.5
    else:
        xStart = 0.1
        xLength = 0.8
    if km is None:
        yStart = 0.1
        yLength = 0.75
    else:
        yStart = 0.05
        yLength = 0.85
    figAxes = heatMapAxes(data, dims = [xStart, yStart, xLength, yLength], columns=nameList, rows=proteins, protColors=protColors, cIndex=cIndex, fig=fig, colors=colors)
    ##Draw colorbar
    if colorBar:
        fig.colorbar(figAxes)

    if dendro:
        ax2Data = fig.add_axes([offset, offset, xLength-0.3, yLength])
        sch.dendrogram(xdat['rightDendro'], orientation='right', color_threshold=0.0)
        ax2Data.set_xticks([])
        ax2Data.set_yticks([])

#    if not (scale is None):
 #       figAxes.text(1.15, 0.5, scale, rotation=270, verticalalignment="center", 
  #                   horizontalalignment="center", fontsize=18, transform=figAxes.transAxes)    
    
    if not km is None:
        small = data.min()
        big = data.max()
        if math.fabs(small) > math.fabs(big):
            big = 0-small
        else:
            small = 0-big
        ax3Data = fig.add_axes([xStart, yLength+offset, xLength-0.1, 0.1])
        ax3Data.matshow(km, aspect='auto', origin='lower', cmap=colors, vmin=small, vmax=big)
        for i in range(len(km)):
            ax3Data.text(-0.75, i, 'clus'+str(i), verticalalignment="center", horizontalalignment="right", fontsize=10, color=cIndex(float(i)/(protColors.max()+1)))
        ax3Data.set_xticks([])
        ax3Data.set_yticks([])
    #fig.tight_layout()
    if not (saveName is None):
        pylab.savefig(saveName)
        
    return fig

'''
def heatMapAxes(data, dims=[0.1, 0.1, 0.7, 0.7], colors=pylab.cm.RdBu, columns=None, rows=None, protColors=None, cIndex=None, fig=None, colorBar=False):
    """heatMapAxes draws a heatmap

    :param data: a datamatrix to draw
    :type xdat: a 2D Matrix
    :param dims: the size of the plot to draw - defaults to full window
    :type dims: list (4 eements)
    :param colors: color index to use - defaults to redblue
    :type colors: cmap
    :param fractions: fraction names
    :type fractions: list
    :param proteins: protein names
    :type proteins: list
    :param protColors: a color map used to label the protein names with group colors
    :type protColors: cmap
    :param cIndex: a list of groupIds for the proteins
    :type cIndex: list
    :param fig: where to plot the axes (which figure); defaults to new figure
    :type fig: matplotlib figure
    :returns:  an axes

    """
    
    if fig is None:
        fig = pylab.figure()
    axData = fig.add_axes(dims)
    for i in range(len(columns)):
        axData.text(i, -0.4 , ' '+str(columns[i]), rotation=270, verticalalignment="top", horizontalalignment="center", fontsize=14)
    if protColors == None:
        for i in range(len(rows)):
            axData.text(-0.6, i, '  '+str(rows[i]), verticalalignment="center", horizontalalignment="right", fontsize=14)
    else:
        for i in range(len(rows)):
            axData.text(-0.6, i, '  '+str(rows[i]), verticalalignment="center", horizontalalignment="right", fontsize=14, color=cIndex(float(protColors[i])/(protColors.max()+1)))
    
    small = data.min()
    big = data.max()
    if math.fabs(small) > math.fabs(big):
        big = 0-small
    else:
        small = 0-big
    
    masked_array = numpy.ma.array (data, mask=numpy.isnan(data))
    colors.set_bad('grey',1.)
    figData = axData.imshow(masked_array, interpolation='nearest', cmap=colors, aspect='auto', origin='lower')
    ##Draw colorbar
    if colorBar:
        fig.colorbar(figData, ax=axData, ticks=[0, 0.25, 0.50, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0], pad=0.01, extend='neither')
    
    axData.set_xticks([])
    axData.set_yticks([])

    return axData
'''

def heatMapAxes(data, dims=[0.1, 0.1, 0.7, 0.7], colors=pylab.cm.autumn, columns=None, rows=None, protColors=None, cIndex=None, fig=None, colorBar=False):
    """heatMapAxes draws a heatmap

    :param data: a datamatrix to draw
    :type xdat: a 2D Matrix
    :param dims: the size of the plot to draw - defaults to full window
    :type dims: list (4 elements)
    :param colors: color index to use - defaults to redblue
    :type colors: cmap
    :param fractions: fraction names
    :type fractions: list
    :param proteins: protein names
    :type proteins: list
    :param protColors: a color map used to label the protein names with group colors
    :type protColors: cmap
    :param cIndex: a list of groupIds for the proteins
    :type cIndex: list
    :param fig: where to plot the axes (which figure); defaults to new figure
    :type fig: matplotlib figure
    :returns:  an axes

    """
    if fig is None:
        fig = pylab.figure()
    axData = fig.add_axes(dims)
    for i in range(len(columns)):
        axData.text(i, -0.5 , ' '+str(columns[i]), rotation=270, verticalalignment="top", horizontalalignment="center", fontsize=12)
    if protColors == None:
        for i in range(len(rows)):
            axData.text(-0.75, i, '  '+str(rows[i]), verticalalignment="center", horizontalalignment="right", fontsize=12)
    else:
        for i in range(len(rows)):
            axData.text(-0.75, i, '  '+str(rows[i]), verticalalignment="center", horizontalalignment="right", fontsize=12, color=cIndex(float(protColors[i])/(protColors.max()+1)))
    small = data.min()
    big = data.max()
    if math.fabs(small) > math.fabs(big):
        big = 0-small
    else:
        small = 0-big
    masked_array = numpy.ma.array (data, mask=numpy.isnan(data))
    colors.set_bad('grey',1.)
    figData = axData.imshow(masked_array, interpolation='nearest', cmap=colors, aspect='auto', origin='lower')
    if colorBar:
        fig.colorbar(figData, ax=axData, ticks=[0, 0.25, 0.50, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0], pad=0.01, extend='neither')
    axData.set_xticks([])
    axData.set_yticks([])
    return figData

def findIndices(g):
    """findIndices is a  helper function that likely should be deleted
    
    """
    change = [0]
    seen = [g[0]]
    for i in range(1, len(g)):
        if not g[i] in seen:
            change.append(i)
            seen.append(g[i])
    return change
    
def mapGroups(groupList, letters):
    """mapGroups is a  helper function that maps groups numbers to letters - likely should be deleted

    :param groupList: a list to be mapped
    :type d: list
    :param letters: a list to be mapped onto
    :type letters: list
    :returns: a list with elements of groupList mapped onto the letters
    
    """
    changeList = findIndices(groupList)
    i = 0
    for index in changeList:
        toReplace = groupList[index]
        groupList = listReplace(groupList, toReplace, letters[i])
        i = i+1
    return list(groupList)