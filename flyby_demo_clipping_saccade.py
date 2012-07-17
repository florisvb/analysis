import sys
#sys.path.insert(0, '/usr/local/lib/python2.6/dist-packages')
sys.path.append('/home/floris/src/pymovie2')

from matplotlib import rcParams
fig_width = 3.6 # width in inches
fig_height = 3.6  # height in inches
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





def plot_demo_flyby(dataset, key):  


    trajec = dataset.trajecs[key]
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_ylim(-.15,.15)
    ax.set_xlim(-.25, .25)
    ax.set_autoscale_on(False)
    ax.set_aspect('equal')
    
    frames = np.arange(0,trajec.frames_of_flyby[-1]).tolist()
    ax.plot(trajec.positions[:,0], trajec.positions[:,1], color='gray')
    ax.plot(trajec.positions[frames,0], trajec.positions[frames,1], color='black')
    
    if len(trajec.all_saccades) > 0:
        for s in trajec.all_saccades:
            if s < trajec.frame_nearest_to_post:
                sac_range = fa.get_saccade_range(trajec, s) 
                ax.plot(trajec.positions[sac_range,0], trajec.positions[sac_range,1], '-', color='green', alpha=1, linewidth=1)
                sac = patches.Circle( (trajec.positions[sac_range[0],0], trajec.positions[sac_range[0],1]), radius=0.002, facecolor='green', edgecolor='none', alpha=1, zorder=100)
                ax.add_artist(sac)
    
    post = patches.Circle( (0, 0), radius=0.009565, facecolor='black', edgecolor='none', alpha=1)
    ax.add_artist(post)
    
    fig.savefig('demo_flyby_flydra.pdf', format='pdf')
