import sympy
import numpy as np
import scipy.optimize
import scipy.stats.distributions as distributions
from matplotlib import patches
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
    # from Sokal and Rolhf Biometry pgs ~518
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




def histogram(ax, data_list, bins=10, bin_width_ratio=0.6, colors='green', edgecolor='none', bar_alpha=0.7, curve_alpha=0.4, curve_butter_filter=[3,0.3], return_vals=False, show_smoothed=True):
    
    n_bars = float(len(data_list))
    if type(bins) is int:
        bins = np.linspace(np.min(data_list), np.max(data_list), bins, endpoint=True)
        
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
    for i, data in enumerate(data_list):
        print data
        data_hist = np.histogram(data, bins=bins)[0]
        print len(data), len(data_hist)
        rects = custom_hist_rectangles(data_hist, bins[0:-1]+bar_width*i+bin_width_buff, width=bar_width, facecolor=colors[i], edgecolor=edgecolor, alpha=bar_alpha)
        for rect in rects:
            rect.set_zorder(1)
            ax.add_artist(rect)
        
        if show_smoothed:
            data_hist_filtered = signal.filtfilt(butter_b, butter_a, data_hist)
            interped_bin_centers = np.linspace(bin_centers[0], bin_centers[-1], 100, endpoint=True)
            v = 100 / float(len(bin_centers))
            interped_data_hist_filtered = np.interp(interped_bin_centers, bin_centers, data_hist_filtered)
            interped_data_hist_filtered2 = signal.filtfilt(butter_b/v, butter_a/v, interped_data_hist_filtered)
            #ax.plot(bin_centers, data_hist_filtered, color=facecolor[i])
            ax.fill_between(interped_bin_centers, interped_data_hist_filtered2, np.zeros_like(interped_data_hist_filtered2), color=colors[i], alpha=curve_alpha, zorder=-100, edgecolor='none')
        
        if return_vals:
            data_hist_list.append(data_hist)
            
            if show_smoothed:
                data_curve_list.append(data_hist_filtered)
                
    if return_vals:
        return bins, data_hist_list, data_curve_list
        

 
    
