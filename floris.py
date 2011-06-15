import sympy
import numpy as np



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
    
