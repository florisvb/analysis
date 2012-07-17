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




def run(neural_delay=0.05):
    
    initial_retinal_size = 15*np.pi/180.
    expansion = np.array([143, 250, 330, 500, 1000, 1430, 2000, 5000])*np.pi/180.
    delay = np.array([0.3, .22, .21, .16, .11, .105, .1, .1])
    angle_subtended = initial_retinal_size + (delay-neural_delay)*expansion
    
    tti = np.abs(np.sin(angle_subtended/2.)-1) / (0.5*1/np.tan(angle_subtended/2.)*expansion)
    print tti
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    ax.plot( angle_subtended[0:5], expansion[0:5], '.', color='blue' ) 
    ax.plot( angle_subtended[5:], expansion[5:], '.', color='red' ) 
    
    deg_ticks = np.array([0, 45, 90, 135])
    deg_tick_strings = [str(d) for d in deg_ticks]
    rad_ticks = deg_ticks*np.pi/180.
    
    
    ## plot paramters    
    ax.set_ylim([0,50])
    ax.set_xlim(rad_ticks[0], rad_ticks[-1])
    ax.set_autoscale_on(False)
    
    fa.adjust_spines(ax, ['left', 'bottom'])
    ax.set_xlabel('position on retina, deg')
    ax.set_ylabel('expansion threshold, rad/s')
    
    ax.set_xticks(rad_ticks)
    ax.set_xticklabels(deg_tick_strings)
    
    filename = 'tammero_data_landing.pdf'
    fig.savefig(filename, format='pdf')

    print np.mean(angle_subtended[0:5])*180/np.pi



