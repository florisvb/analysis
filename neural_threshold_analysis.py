import sys
#sys.path.insert(0, '/usr/local/lib/python2.6/dist-packages')
sys.path.append('/home/floris/src/pymovie2')

from matplotlib import rcParams
fig_width = 3.25 # width in inches
fig_height = 3.25  # height in inches
fig_size =  (fig_width, fig_height)

fontsize = 8
params = {'backend': 'Agg',
          'ps.usedistiller': 'xpdf',
          'ps.fonttype' : 3,
          'pdf.fonttype' : 3,
          'font.family' : 'sans-serif',
          'font.serif' : 'Times, Palatino, New Century Schoolbook, Bookman, Computer Modern Roman',
          'font.sans-serif' : 'Helvetica, Avant Garde, Computer Modern Sans serif',
          'font.cursive' : 'Zapf Chancery',
          'font.monospace' : 'Courier, Computer Modern Typewriter',
          'font.size' : fontsize,
          'text.fontsize': fontsize,
          'axes.labelsize': fontsize,
          'axes.linewidth': 1.0,
          'xtick.major.linewidth': 1,
          'xtick.minor.linewidth': 1,
          #'xtick.major.size': 6,
          #'xtick.minor.size' : 3,
          'xtick.labelsize': fontsize,
          #'ytick.major.size': 6,
          #'ytick.minor.size' : 3,
          'ytick.labelsize': fontsize,
          'figure.figsize': fig_size,
          'figure.dpi' : 72,
          'figure.facecolor' : 'white',
          'figure.edgecolor' : 'white',
          'savefig.dpi' : 300,
          'savefig.facecolor' : 'white',
          'savefig.edgecolor' : 'white',
          'figure.subplot.left': 0.2,
          'figure.subplot.right': 0.8,
          'figure.subplot.bottom': 0.25,
          'figure.subplot.top': 0.9,
          'figure.subplot.wspace': 0.0,
          'figure.subplot.hspace': 0.0,
          'lines.linewidth': 1.0,
          'text.usetex': True, 
          }
rcParams.update(params) 

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib

import colorline
import flydra_analysis as fa
import sa1_analysis as sa1
import numpy as np
import flydra_floris_analysis as ffa
import floris

import scipy.optimize

def expansion(vel, a=None, d=None, r=0.009565):
    if d is None:
        d = r / np.sin(a/2.)
    return 2*(r/d**2)*vel / np.sqrt(1-(r/d)**2)
    
def timetocontact(vel, a=None, d=None, r=0.009565):
    exp = expansion(vel=vel, a=a, d=d, r=r)
    return (2*np.sin(a/2.)) / (exp*np.sqrt(1-(np.sin(a/2.))**2))
    
def expansion_from_timetocontact(ttc, a):
    e = (2*np.sin(a/2.)) / (ttc*np.sqrt(1-(np.sin(a/2.))**2))
    e[e<0] = 0
    return e
    

# for match fmin:
# True time-to-contact model
def calc_expthreshold_ttc(param, a):
    ttc, ttc_threshold = param
    e = expansion_from_timetocontact(ttc, a) + 0
    e[e<0] = 0
    return e
def calc_error(param, a):
    expthreshold_ttc = calc_expthreshold_ttc(param, a)
    
    # RSDET model
    m = -0.21
    b = 0.159
    vel = (m*np.log(a)+b)
    expthreshold = expansion(vel, a=a)
    
    a_less = np.where( a < 80*np.pi/180. )[0].tolist()
    err = np.sum( np.abs(expthreshold_ttc[a_less] - expthreshold[a_less]) )
    return err
        
    
def match_ttc_to_rsdet(ttc0=0.13, ttc_threshold0=-0.2):
    
    radius = 0.009565
    a = np.linspace(0.001,2.5,100)
    
    val = scipy.optimize.fmin(calc_error, [ttc0, ttc_threshold0], args=[a], disp=0)
    
    return val
    

def neural_threshold_tti_vs_rsdet_models(dataset_landing, save_plot=True, movie_dataset=None, ttc=None):
    
    distfig = plt.figure()
    distax = distfig.add_subplot(111)
    
    radius = 0.009565
    a = np.linspace(0,2.5,100)
    
    fit, Rsq, x, y, yminus, yplus = fa.get_angle_vs_speed_curve(dataset_landing, plot=False, plot_sample_trajecs=False, post_type=['checkered', 'checkered_angled', 'black', 'black_angled'], filename=None, keys=None, tti=None, color_code_posts=False)
    std = np.mean(yplus - y)
    
    # RSDET model
    m = fit[0]
    b = fit[1]
    vel = (m*np.log(a)+b)
    
    
    print std, fit
    
    expthreshold = expansion(vel, a=a)
    
    expthreshold_plus = expansion(vel+std, a=a)
    expthreshold_minus = expansion(vel-std, a=a)
    distax.plot( np.log(a), expthreshold, color='purple')
    distax.fill_between(np.log(a), expthreshold_plus, expthreshold_minus, color='purple', linewidth=0, alpha=0.3)
    
    # True time-to-contact model
    if ttc is None:
        ttc, ttc_threshold = match_ttc_to_rsdet(ttc0=0.13, ttc_threshold0=0)
    ttc_threshold = 0
    expthreshold_ttc = expansion_from_timetocontact(ttc, a) + ttc_threshold
    distax.plot( np.log(a), expthreshold_ttc, ':', color='purple')
    
    
    
    
    # plot a sample constant velocity trajectory
    vels = [0.2, 0.4, 0.8]
    
    for vel in vels:
        fps = 5000.0
        x = np.arange(.2, 0.0, -vel/fps)
        d = x+radius
        a = 2*np.arcsin(radius / (d))
        #exp = 2/np.sqrt(1-(radius/(d))**2) * (radius/(d)**2) * vel
        exp = expansion(vel, a=a)
        indices = np.where(exp<12)[0].tolist()
        distax.plot( np.log(a[indices]), exp[indices], color='gray', linewidth=0.5)
    
    
    # plot parameters    
    fa.fix_angle_log_spine(distax, histograms=False, set_y=False)
    ylim_max = 1000
    distax.set_ylim(0,ylim_max/180.*np.pi)
    rad_ticks_y = np.linspace(0,ylim_max*np.pi/180.,5,endpoint=True)
    deg_tick_strings_y = [str(s) for s in np.linspace(0,ylim_max,5,endpoint=True)]
    for i, s in enumerate(deg_tick_strings_y):
        deg_tick_strings_y[i] = s.split('.')[0]
    distax.set_yticks(rad_ticks_y)
    distax.set_yticklabels(deg_tick_strings_y)
    distax.set_ylabel('Expansion, deg/s')
    
    if save_plot:
        distfig.savefig('neural_threshold_distance.pdf', format='pdf')
        
    return
