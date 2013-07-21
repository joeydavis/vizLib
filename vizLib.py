import math
import matplotlib.pyplot as pylab
import scipy.cluster.hierarchy as sch
import qMS
import numpy
import matplotlib

def setRcs(scale=None):
    """setRcs sets a series of rc params for matplotlib tomake decent looking plots. Also some useful def
    in there to cutomize things.

    :param scale: the default font scale
    :type scale: float
    
    :returns:  nothing, edits the rcparam file.
            
    """
    if scale is None:
        scale = 12
    
    defaultFont = {'family' : 'serif',
                   'variant' : 'normal',
                   'weight' : 400,
                   'size' : scale*1}

    axisFont = {'titlesize' : scale*1.5,
                'labelsize' : scale*1.25}

    xAxisTicks = {'major.size' : 8,
                  'minor.size' : 4,
                  'major.width' : 1,
                  'minor.width' : 1,
                  'labelsize' : scale*1,
                  'minor.pad' : 3,
                  'major.pad' : 3}

    yAxisTicks = {'major.size' : 8,
                  'minor.size' : 4,
                  'major.width' : 1,
                  'minor.width' : 1,
                  'labelsize' : scale*1,
                  'minor.pad' : 3,
                  'major.pad' : 3}

    legend = {'fancybox' : True,
              'numpoints' : 1,
              'fontsize' : scale*0.85,
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
    
def plotStatsDict(statsDict, name='', proteins=None, offset=0.0, markerSize=12, color='#e31a1c', yMax = 1.5, median=False):
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
        proteins = sorted(statsDict.keys())

    xAxis = range(1,len(proteins)+1)
    fig = pylab.figure(figsize=(22,5))
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
    ax.plot(xs, ys, 'o', color=color, markersize=markerSize, label=name)
    pylab.xticks(xAxis, [item for item in proteins], rotation=45, size=15)
    pylab.xlim(1, len(proteins)+1)
    ####################################
    pylab.yticks([0,yMax/5, 2*yMax/5, 3*yMax/5, 4*yMax/5, yMax], size=15)
    ####################################
    pylab.ylim(0, yMax)
    return ax

def addStatsDictToPlot(statsDict, ax, name='', offset=0.0, markerSize=12, color='#377db8', median=False):
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
    ax.plot(xs, ys, 'o', color=color, markersize=markerSize, label=name)
    return ax

def makePlotWithFileList(isoFileList, numerator, denominator, AllProteins=None, normProtein=None, yMax = 1.5, median=False, names=None, colors=None):
    """makePlotWithFileList is a  helper function that plots massage-style data from a list of files

    :param isoFileList: a list of the files to be ploted (shoudl be full path to _iso.csv files)
    :type isoFileList: list of strings
    :param numerator: strings of the keys in the numerator (ampu, ampl, amps)
    :type numerator: list of strings
    :param denominator: strings of the keys in the denominator (ampu, ampl, amps)
    :type denominator: list of strings
    :param AllProteins: a list of the proteins to be plotted - strings must match the keys of the 'proteins' field in _iso.csv
    :type AllProteins: list
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

    :returns: the plotted axis
    
    """

    allStats = qMS.multiStatsDict(isoFileList, numerator, denominator, normalization=1.0, offset=0.0, normProtein=normProtein)
    
    if names is None:
        names = isoFileList
    
    if colors is None:
        colors= ['#377db8' for i in isoFileList]
    
    offsets = float(len(isoFileList)+1)

    ax = plotStatsDict( allStats[isoFileList[0]], name=names[0], proteins=AllProteins, \
                        offset=1.0/offsets, markerSize=(10.0/offsets)+4, yMax=yMax, median=median, color=colors[0])
    i = 1
    for k in isoFileList[1:]:
        ax = addStatsDictToPlot(allStats[k], ax, name=names[i], \
                                offset=(1.0/offsets)*(i+1.25), markerSize=(10.0/offsets)+4, median=median, color=colors[i])
        i = i + 1
    return ax

def drawHeatMap(xdat, name="unnamed", colors=pylab.cm.RdBu, dendro=False, protColors=None, cIndex=None, km=None):
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
    vecs = xdat['fractions']
    ls = xdat['proteins']
    fig = pylab.figure()
    fig.suptitle(name)
    ##Draw heatmap
    offset = 0.1
    yStart = 0.1
    yLength = 0.8
    if dendro:
        xStart = 0.35
        xLength = 0.5
    else:
        xStart = 0.1
        xLength = 0.8
    
    figData = heatMapAxes(data, dims = [xStart, yStart, xLength, yLength], columns=vecs, rows=ls, protColors=protColors, cIndex=cIndex, fig=fig)
    ##Draw colorbar
    fig.colorbar(figData)

    if dendro:
        ax2Data = fig.add_axes([offset, offset, xLength-0.3, yLength])
        sch.dendrogram(xdat['rightDendro'], orientation='right', color_threshold=8)
        ax2Data.set_xticks([])
        ax2Data.set_yticks([])
        
    return fig

def heatMapAxes(data, dims=[0.1, 0.1, 0.7, 0.7], colors=pylab.cm.RdBu, columns=None, rows=None, protColors=None, cIndex=None, fig=None):
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
        axData.text(i, -0.5 , ' '+str(columns[i]), rotation=270, verticalalignment="top", horizontalalignment="center", fontsize=12)
    if protColors == None:
        for i in range(len(rows)):
            axData.text(-0.525, i, '  '+str(rows[i]), verticalalignment="center", horizontalalignment="right", fontsize=8)
    else:
        for i in range(len(rows)):
            axData.text(-0.525, i, '  '+str(rows[i]), verticalalignment="center", horizontalalignment="right", fontsize=8, color=cIndex(float(protColors[i])/(protColors.max()+1)))
    small = data.min()
    big = data.max()
    if math.fabs(small) > math.fabs(big):
        big = 0-small
    else:
        small = 0-big
    figData = axData.matshow(data, aspect='auto', origin='lower', cmap=colors, vmin=small, vmax=big)
    #fig.colorbar(figData)
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