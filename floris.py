import sympy
import numpy as np
import scipy.optimize
import scipy.stats.distributions as distributions
from matplotlib import patches
import matplotlib.pyplot as plt
import colorline as colorline
from matplotlib.collections import LineCollection
from scipy import signal

class Dummy_Class():
    def __init__(self):
        pass

def atan2(y, x):
    """take arctangent of y/x preserving angle"""
    # fromhttp://en.wikipedia.org/wiki/Atan2
    return 2*sympy.atan( y/(sympy.sqrt(x**2+y**2) + x) )
    
    
def sp2np (A):

    An = np.zeros(A.shape)
    
    for i in range(An.shape[0]):
        for j in range(An.shape[1]):
            An[i,j] = A[i,j]
            
    return An
        
        
def geomean(x):
    return np.exp(np.sum(np.log(x)) / len(x))
def geostd(x):
    G = geomean(x)
    return np.exp(np.sqrt(np.sum(np.log(x/G)**2)/len(x)))
    
        
def dist_point_to_line(pt, linept1, linept2, sign=False):
    # from wolfram mathworld
    x1 = linept1[0]
    x2 = linept2[0]
    y1 = linept1[1]
    y2 = linept2[1]
    x0 = pt[0]
    y0 = pt[1]
    
    if sign:
        d = -1*((x2-x1)*(y1-y0)-(x1-x0)*(y2-y1) )  / np.sqrt( (x2-x1)**2+(y2-y1)**2)
    else:
        d = np.abs( (x2-x1)*(y1-y0)-(x1-x0)*(y2-y1) ) / np.sqrt( (x2-x1)**2+(y2-y1)**2 )
    
    return d
    
        
def binarize(x, threshold=0.00001, val=1):

    def bin(x, threshold, val):
        if x > threshold:
            return val
        if x <= threshold:
            return 0
    
    if type(x) in [int, float]:
        return bin(x,threshold,val)
    
    if type(x) is list:
        return [bin(i,threshold,val) for i in x] 
    
    if type(x) is np.ndarray:
        x[x>threshold] = val
        return x
            
def bootstrap_linear_fit(xdata, ydata, n=None, confidence=0.95):
    if n is None:  
        n = len(xdata)
    print n
    fits = np.zeros([n, 2])
    print fits.shape
    
    for i in range(n):
        # Choose #sample_size members of d at random, with replacement
        choices = np.random.random_integers(0, len(xdata)-1, n)
        xsample = xdata[choices]
        ysample = ydata[choices]
        #fit = linear_fit_type2(xsample, ysample)
        fit = np.polyfit(xsample, ysample, 1)
        fits[i,:] = fit
        
    def get_mean_and_confidence_interval(data, confidence):
        confidence_indices = [int(n*(1-confidence)/2.), int(n-n*(1-confidence)/2.)]
        mean = np.mean(data)
        argsort = np.argsort(data)
        confidence_indices = argsort[confidence_indices]
        confidence_range = data[confidence_indices]
        return mean, confidence_range
        
    slope_mean, slope_confidence_range = get_mean_and_confidence_interval(fits[:,0], 0.95)
    intercept_mean, intercept_confidence_range = get_mean_and_confidence_interval(fits[:,1], 0.95)
        
    #print 'slope: mean: ', slope_mean, ' | confidence: ', slope_confidence_range[0], ',', slope_confidence_range[1]
    #print 'intercept: mean: ', intercept_mean, ' | confidence: ', intercept_confidence_range[0], ',', intercept_confidence_range[1]
        
    return slope_mean, slope_confidence_range, intercept_mean, intercept_confidence_range
    
def bootstrap_histogram(xdata, bins, normed=False, n=None, return_raw=False):
    if type(xdata) is not np.ndarray:
        xdata = np.array(xdata)

    if n is None:  
        n = len(xdata)
    hist_list = np.zeros([n, len(bins)-1])
    
    for i in range(n):
        # Choose #sample_size members of d at random, with replacement
        choices = np.random.random_integers(0, len(xdata)-1, n)
        xsample = xdata[choices]
        hist = np.histogram(xsample, bins, normed=normed)[0].astype(float)
        hist_list[i,:] = hist
        
    hist_mean = np.mean(hist_list, axis=0)
    hist_std = np.std(hist_list, axis=0)
    
    if return_raw:
        return hist_list
    else:
        return hist_mean, hist_std
    
    
    
def linear_fit_confidence_interval(xdata, ydata, alpha=0.05):
    fit = linear_fit_type2(xdata, ydata)    
        
    
    
    
    
def normalize(array):
    normed_array = norm_array(array)
    return array / normed_array
def norm_array(array):
    normed_array = np.zeros_like(array)
    for i in range(len(array)):
        normed_array[i,:] = np.linalg.norm(array[i])
    return normed_array[:,0]
    
def diffa(array):
    if len(array.shape) == 2:
        array = array.reshape(len(array))
    if len(array) < 2:
        return np.array([0])

    d = np.diff(array)
    d = np.hstack( (d[0], d) )
    return d
    
def linear_fit_type2(xdata, ydata, alpha=0.05, full_output=False, weights=None):
    return linear_fit_type(xdata, ydata, alpha=0.05, full_output=full_output, weights=None, fit_type=2)
    
def linear_fit_type1(xdata, ydata, alpha=0.05, full_output=False, weights=None):
    return linear_fit_type(xdata, ydata, alpha=0.05, full_output=full_output, weights=None, fit_type=1)
    
def linear_fit_type(xdata, ydata, alpha=0.05, full_output=False, weights=None, fit_type=1):
    if weights is None:
        weights = np.ones_like(xdata)
    fity = np.polyfit(xdata, ydata, 1)
    #fitx = np.polyfit(ydata, xdata, 1)
    
    def error(params):
        linept1 = [0, params[1]]
        linept2 = [1, params[0]*1+params[1]]
        
        err = 0
        #errarr = []
        for i, x in enumerate(xdata):
            y = ydata[i]
            if fit_type == 2:
                err += dist_point_to_line([x,y], linept1, linept2)**2 *weights[i]
            elif fit_type == 1:
                err += np.abs(y - (x*params[0]+params[1])) *weights[i]
            #errarr.append( dist_point_to_line([x,y], linept1, linept2, sign=True))
        return err
        
    fitfmin = scipy.optimize.fmin(error, fity)
    err = error(fitfmin)
    variance = err / (len(xdata)-2)
    
    # intercept
    SXX = np.sum((xdata - np.mean(xdata))**2)
    SYY = np.sum((ydata - np.mean(ydata))**2)
    SXY = np.sum( (ydata - np.mean(ydata))*(xdata - np.mean(xdata)) )
    Rsq = (SXY)**2 / (SXX*SYY)
    se_intercept = np.sqrt(variance) * (1/len(xdata)+np.mean(xdata)**2/SXX)**(0.5)
    intercept_confidence_interval = [fitfmin[1]-distributions.t.pdf(alpha/2., len(xdata)-2)*se_intercept, fitfmin[1]+distributions.t.pdf(alpha/2., len(xdata)-2)*se_intercept]
    
    # slope
    se_slope = np.sqrt(variance) / (np.sqrt(SXX))
    slope_confidence_interval = [fitfmin[0]-distributions.t.pdf(alpha/2., len(xdata)-2)*se_slope, fitfmin[0]+distributions.t.pdf(alpha/2., len(xdata)-2)*se_slope]
    
    #print fity
    #print fitx
    #print (fity[0]+fitx[0])/2., (fity[1]+fitx[1])/2.
    #print fitfmin
    
    if not full_output:
        return fitfmin
    else:
        return fitfmin, variance, intercept_confidence_interval, slope_confidence_interval, Rsq
    

def ANCOVA(groups):
    # see Sokal and Rolhf Biometry pgs ~518
    a = len(groups)
    GroupResults = [None for i in range(a)]
    for g, group in enumerate(groups):
        G = Dummy_Class()
        G.X = group[0]
        G.Y = group[1]
        G.n = len(G.X)
        G.df = G.n-1
        G.SSY = np.sum(G.Y**2) - (np.sum(G.Y))**2 / G.n
        G.SPXY = np.sum(G.X*G.Y) - (np.sum(G.X)*np.sum(G.Y)) / G.n
        G.SSX = np.sum(G.X**2) - (np.sum(G.X))**2 / G.n
        G.bYX = G.SPXY / G.SSX
        G.SSYhat = G.SPXY**2 / G.SSX
        G.SSYX = G.SSY - G.SSYhat
        G.MSYX = G.SSYX / (G.n-2)
        G.Fs = G.SSYhat / G.MSYX
        GroupResults[g] = G
        
        '''
        print
        print G.SSY
        print G.SPXY
        print G.SSX
        print G.bYX
        print G.SSYhat
        print G.SSYX
        print G.MSYX
        print G.Fs
        print 
        '''
        
    SumofGroups = Dummy_Class()
    SumofGroups.n = np.sum( [G.n for G in GroupResults] )
    SumofGroups.df = SumofGroups.n - 2*a
    SumofGroups.SSYX = np.sum( [G.SSYX for G in GroupResults] )
    SumofGroups.MSYX = SumofGroups.SSYX / SumofGroups.df
    
    PooledWithin = Dummy_Class()
    PooledWithin.df = np.sum( [G.df for G in GroupResults] )
    PooledWithin.SSY = np.sum( [G.SSY for G in GroupResults] )
    PooledWithin.SPXY = np.sum( [G.SPXY for G in GroupResults] )
    PooledWithin.SSX = np.sum( [G.SSX for G in GroupResults] )
    PooledWithin.b = PooledWithin.SPXY/PooledWithin.SSX
    PooledWithin.SSYhat = PooledWithin.SPXY**2 / PooledWithin.SSX
    PooledWithin.SSYX = PooledWithin.SSY - PooledWithin.SSYhat
    PooledWithin.MSYX = PooledWithin.SSYX / (np.sum( [G.n for G in GroupResults] ) - a - 1)
    PooledWithin.FS = PooledWithin.SSY / PooledWithin.MSYX
    
    Amongb = Dummy_Class()
    Amongb.df = a-1
    Amongb.SSYX = PooledWithin.SSYX - SumofGroups.SSYX
    Amongb.MSYX = Amongb.SSYX / (a-1)
    Amongb.FS = Amongb.MSYX / SumofGroups.MSYX
    
    pc1_Y = np.sum( [np.sum(G.Y) for G in GroupResults] )
    pc1_X = np.sum( [np.sum(G.X) for G in GroupResults] )
    pc2_Y = np.sum( [np.sum(G.Y**2) for G in GroupResults] )
    pc2_X = np.sum( [np.sum(G.X**2) for G in GroupResults] )
    pc2_XY = np.sum( [np.sum(G.X*G.Y) for G in GroupResults] )
    pc3_Y = np.sum( np.array([(np.sum(G.Y))**2 / G.n for G in GroupResults]) )
    pc3_X = np.sum( np.array([(np.sum(G.X))**2 / G.n for G in GroupResults]) )
    pc4_Y = pc1_Y**2 / SumofGroups.n
    pc4_X = pc1_X**2 / SumofGroups.n
    
    Total = Dummy_Class()
    Total.SSY = pc2_Y - pc4_Y
    Total.SSX = pc2_X - pc4_X
    Total.SPXY = pc2_XY - pc1_Y*pc1_X / SumofGroups.n
    Total.SSYX = Total.SSY - Total.SPXY**2 / Total.SSX
    
    Groups = Dummy_Class()
    Groups.SSY = pc3_Y - pc4_Y
    Groups.SSX = pc3_X - pc4_X
    Groups.SPXY = Total.SPXY - PooledWithin.SPXY
    Groups.SSYX = Groups.SSY - Groups.SPXY**2 / Groups.SSX
    Groups.MSYX = Groups.SSYX / (a-2)
    
    AdjustedMeans = Dummy_Class()
    AdjustedMeans.df = a-1
    AdjustedMeans.SSYX = Total.SSYX - PooledWithin.SSYX
    AdjustedMeans.MSYX = AdjustedMeans.SSYX / (a-1)
    AdjustedMeans.FS = AdjustedMeans.MSYX / PooledWithin.MSYX
    
    ## calculate significance ##
    alpha_slope = distributions.f.sf(Amongb.FS, Amongb.df, PooledWithin.df)
    alpha_intercept = distributions.f.sf(AdjustedMeans.FS, AdjustedMeans.df, PooledWithin.df-1)
    
    print 'Fs Among b (significance of variation in slopes among data sets)'
    print 'FS: ', Amongb.FS, 'df1: ', Amongb.df, 'df2: ', PooledWithin.df
    print 'alpha (slope): ', alpha_slope
    print
    print 'Fs Among a (significance of variation in intercepts among data sets - only noteworth if alpha for slope not significant)'
    print 'FS: ', AdjustedMeans.FS, 'df1: ', AdjustedMeans.df, 'df2: ', PooledWithin.df
    print 'alpha (intercept): ', alpha_intercept
    
    
def test_ANCOVA():
    # from Sokal and Rolhf Biometry 2nd ed. pgs ~518
    X1 = np.array([-0.31, 0.17, 0.58, 0.81])
    Y1 = np.array([-2.4, 6.3, 15.8, 20.5])
    print 
    print np.sum(X1), np.sum(Y1)
    G1 = [X1, Y1]
        
    X2 = np.array([-1.18, -0.65, .1, .5, .67])
    Y2 = np.array([-7., 2.1, 17.8, 27.3, 32.])
    print 
    print np.sum(X2), np.sum(Y2)
    G2 = [X2, Y2]
    
    X3 = np.array([-1.79, -1.21, -.35, .08, .49, .65])
    Y3 = np.array([-10.8, -2.8, 14.2, 25.5, 35.7, 41.2])
    print 
    print np.sum(X3), np.sum(Y3)
    G3 = [X3, Y3]
    
    X4 = np.array([-1.83, -1.25, -.41, .05, .43, 0.59])
    Y4 = np.array([-5.4, 3., 20.7, 30.5, 39.9, 45])
    print 
    print np.sum(X4), np.sum(Y4)
    G4 = [X4, Y4]
    
    G = [G1, G2, G3, G4]
    
    ANCOVA(G)
        
    
def custom_hist_rectangles(hist, leftedges, width, facecolor='green', edgecolor='none', alpha=1):

    if type(width) is not list:
        width = [width for i in range(len(hist))]
            
    rects = [None for i in range(len(hist))]
    for i in range(len(hist)):
        rects[i] = patches.Rectangle( [leftedges[i], 0], width[i], hist[i], facecolor=facecolor, edgecolor=edgecolor, alpha=alpha)

    return rects




def histogram(ax, data_list, bins=10, bin_width_ratio=0.6, colors='green', edgecolor='none', bar_alpha=0.7, curve_fill_alpha=0.4, curve_line_alpha=0.8, curve_butter_filter=[3,0.3], return_vals=False, show_smoothed=True, normed=False, normed_occurences=False, bootstrap_std=False, exponential_histogram=False):
    
    n_bars = float(len(data_list))
    if type(bins) is int:
    
        mia = np.array([np.min(d) for d in data_list])
        maa = np.array([np.max(d) for d in data_list])
        
        bins = np.linspace(np.min(mia), np.max(maa), bins, endpoint=True)
        
    if type(colors) is not list:
        colors = [colors]
    if len(colors) != n_bars:
        colors = [colors[0] for i in range(n_bars)]
        
    bin_centers = np.diff(bins)/2. + bins[0:-1]
    bin_width = np.mean(np.diff(bins))
    bin_width_buff = (1-bin_width_ratio)*bin_width/2.
    bar_width = (bin_width-2*bin_width_buff)/n_bars
    
    butter_b, butter_a = signal.butter(curve_butter_filter[0], curve_butter_filter[1])
    
    if return_vals:
        data_hist_list = []
        data_curve_list = []
        data_hist_std_list = []
        
    # first get max number of occurences
    max_occur = []
    for i, data in enumerate(data_list):
        data_hist = np.histogram(data, bins=bins, normed=normed)[0].astype(float)
        max_occur.append(np.max(data_hist))
    max_occur = np.max(np.array(max_occur))
        
    for i, data in enumerate(data_list):
        
        if bootstrap_std:
            data_hist, data_hist_std = bootstrap_histogram(data, bins=bins, normed=normed)
        else:
            data_hist = np.histogram(data, bins=bins, normed=normed)[0].astype(float)
            
        if exponential_histogram:
            data_hist = np.log(data_hist)
        
        if normed_occurences is not False:
            if normed_occurences == 'total':
                data_hist /= max_occur 
                if bootstrap_std:
                    data_hist_std /= max_occur
            else:
                div = float(np.max(data_hist))
                print div
                data_hist /= div 
                if bootstrap_std:
                    data_hist_std /= div
                    
        
        rects = custom_hist_rectangles(data_hist, bins[0:-1]+bar_width*i+bin_width_buff, width=bar_width, facecolor=colors[i], edgecolor=edgecolor, alpha=bar_alpha)
        if bootstrap_std:
            for j, s in enumerate(data_hist_std):
                x = bins[j]+bar_width*i+bin_width_buff + bar_width/2.
                #ax.plot([x,x], [data_hist[j], data_hist[j]+data_hist_std[j]], alpha=1, color='w')
                ax.plot([x,x], [data_hist[j], data_hist[j]+data_hist_std[j]], alpha=bar_alpha, color=colors[i])
                
                #ax.plot([x-bar_width/3., x+bar_width/3.], [data_hist[j]+data_hist_std[j],data_hist[j]+data_hist_std[j]], alpha=1, color='w')
                #ax.plot([x-bar_width/3., x+bar_width/3.], [data_hist[j]+data_hist_std[j],data_hist[j]+data_hist_std[j]], alpha=bar_alpha, color=colors[i])
        for rect in rects:
            rect.set_zorder(1)
            ax.add_artist(rect)
        
                
        if show_smoothed:
            data_hist_filtered = signal.filtfilt(butter_b, butter_a, data_hist)
            interped_bin_centers = np.linspace(bin_centers[0]-bin_width/2., bin_centers[-1]+bin_width/2., 100, endpoint=True)
            v = 100 / float(len(bin_centers))
            interped_data_hist_filtered = np.interp(interped_bin_centers, bin_centers, data_hist_filtered)
            interped_data_hist_filtered2 = signal.filtfilt(butter_b/v, butter_a/v, interped_data_hist_filtered)
            #ax.plot(bin_centers, data_hist_filtered, color=facecolor[i])
            if curve_fill_alpha > 0:
                ax.fill_between(interped_bin_centers, interped_data_hist_filtered2, np.zeros_like(interped_data_hist_filtered2), color=colors[i], alpha=curve_fill_alpha, zorder=-100, edgecolor='none')
            if curve_line_alpha:
                ax.plot(interped_bin_centers, interped_data_hist_filtered2, color=colors[i], alpha=curve_line_alpha)
        
        if return_vals:
            data_hist_list.append(data_hist)
            if bootstrap_std:
                data_hist_std_list.append(data_hist_std)
            
            if show_smoothed:
                data_curve_list.append(data_hist_filtered)
                
    if return_vals and bootstrap_std is False:
        return bins, data_hist_list, data_curve_list
    elif return_vals and bootstrap_std is True:
        return bins, data_hist_list, data_hist_std_list, data_curve_list
        

 
 
def boxplot(ax, x_data, y_data_list, nbins=50, colormap='YlOrRd', linewidth=2, boxwidth=1, usebins=None, boxlinewidth=0.5, outlier_limit=0.02):    

    if usebins is None: 
        usebins = nbins

    cl = colorline.Colorline(ax0=ax, hide_colorbar=True, colormap=colormap)
    for i, y_data in enumerate(y_data_list):
        #print len(y_data)
    
        # calc boxplot statistics
        median = np.median(y_data)
        ind = np.where(y_data<=median)[0].tolist()
        first_quartile = np.median(y_data[ind])
        ind = np.where(y_data>=median)[0].tolist()
        last_quartile = np.median(y_data[ind])
        #print first_quartile, median, last_quartile
        
        # find outliers
        ind_sorted = np.argsort(y_data)
        bottom_limit = int(len(ind_sorted)*(outlier_limit/2.))
        top_limit = int(len(ind_sorted)*(1-outlier_limit/2.))
        indices_inrange = ind_sorted[bottom_limit:top_limit]
        outliers = ind_sorted[0:bottom_limit].tolist() + ind_sorted[top_limit:len(ind_sorted)-1].tolist()
        y_data_inrange = y_data[indices_inrange]
        y_data_outliers = y_data[outliers]
    
    
        # plot colorline
        x = x_data[i]
        hist, bins = np.histogram(y_data_inrange, usebins)
        hist = hist.astype(float)
        hist /= np.max(hist)
        x_arr = np.ones_like(bins)*x
        cl.colorline(x_arr, bins, hist, norm=(0,1), linewidth=linewidth)
        
        
        # plot boxplot
        ax.hlines(median, x-boxwidth/2., x+boxwidth/2., color='black', linewidth=boxlinewidth)
        ax.hlines([first_quartile, last_quartile], x-boxwidth/2., x+boxwidth/2., color='black', linewidth=boxlinewidth/2.)
        ax.vlines([x-boxwidth/2., x+boxwidth/2.], first_quartile, last_quartile, color='black', linewidth=boxlinewidth/2.)
        
        # plot outliers
        x_arr_outliers = x*np.ones_like(y_data_outliers)
        ax.plot(x_arr_outliers, y_data_outliers, '.', markerfacecolor='gray', markeredgecolor='none', markersize=1)
        
    #ax.set_xlim(np.min(x_data), np.max(x_data)*1.1)
    #ax.set_ylim(np.min(y_data_list), np.max(y_data_list)*1.1) 
    
    
def test_boxplot():

    from scipy.stats import norm as normal_dist 
    
    x_data = np.arange(5,20,5)
    y_data_list = []
    for i in range(len(x_data)):
        y_data_list.append( normal_dist.rvs(loc=10*np.random.random(), scale=10*np.random.random(), size=1000) )
        
        fig = plt.figure()
        ax = fig.add_subplot(111)
    
    boxplot(ax, x_data, y_data_list)
        
        
    
    fig.savefig('boxplot_test.pdf', format='pdf')        


# plot a line in x and y with changing colors defined by z, and optionally changing linewidths defined by linewidth
def colorline(ax, x,y,z,linewidth=1, colormap='jet', norm=None, zorder=1, alpha=1, linestyle='solid'):
        cmap = plt.get_cmap(colormap)
        
        if type(linewidth) is list or np.array:
            linewidths = linewidth
        else:
            linewidths = np.ones_like(z)*linewidth
        
        if norm is None:
            norm = plt.Normalize(np.min(z), np.max(z))
        else:
            norm = plt.Normalize(norm[0], norm[1])
        
        '''
        if self.hide_colorbar is False:
            if self.cb is None:
                self.cb = matplotlib.colorbar.ColorbarBase(self.ax1, cmap=cmap, norm=norm, orientation='vertical', boundaries=None)
        '''
            
        # Create a set of line segments so that we can color them individually
        # This creates the points as a N x 1 x 2 array so that we can stack points
        # together easily to get the segments. The segments array for line collection
        # needs to be numlines x points per line x 2 (x and y)
        points = np.array([x, y]).T.reshape(-1, 1, 2)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)
        
        # Create the line collection object, setting the colormapping parameters.
        # Have to set the actual values used for colormapping separately.
        lc = LineCollection(segments, linewidths=linewidths, cmap=cmap, norm=norm, zorder=zorder, alpha=alpha, linestyles=linestyle )
        lc.set_array(z)
        lc.set_linewidth(linewidth)
        
        ax.add_collection(lc)

def colorline_example():
    
    def tent(x):
        """
        A simple tent map
        """
        if x < 0.5:
            return x
        else:
            return -1.0*x + 1
    
    pi = np.pi
    t = np.linspace(0, 1, 200)
    y = np.sin(2*pi*t)
    z = np.array([tent(x) for x in t]) 
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    # standard colorline
    colorline(ax,t,y,z)
    
    # colorline with changing widths, shifted in x
    colorline(ax,t+0.5,y,z,linewidth=z*5)
    
    # colorline with points, shifted in x
    colorline(ax,t+1,y,z, linestyle='dotted')
    
    # set the axis to appropriate limits
    ax.set_xlim(0,2)
    ax.set_ylim(0,1.5)
       
    fig.savefig('colorline_example_plot.pdf', format='pdf')

    
