import sympy
import numpy as np
import scipy.optimize
import scipy.stats.distributions as distributions

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
    
def linear_fit_type2(xdata, ydata, alpha=0.05, full_output=False):
    
    fity = np.polyfit(xdata, ydata, 1)
    #fitx = np.polyfit(ydata, xdata, 1)
    
    def error(params):
        linept1 = [0, params[1]]
        linept2 = [1, params[0]*1+params[1]]
        
        err = 0
        #errarr = []
        for i, x in enumerate(xdata):
            y = ydata[i]
            err += dist_point_to_line([x,y], linept1, linept2)**2
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
    
    
    
    
    
    
