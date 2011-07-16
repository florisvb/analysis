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






def get_slipangles(movie, behavior='straight'):
    fa.prep_movie_trajec(movie)
    trajec = movie.trajec
    
    # look at all frames not in saccade range or within 1cm of post
    closest_dist = 0.01
    
    first_frame = sa1.get_frame_from_timestamp_flydra(movie, movie.timestamps[0])
    last_frame = sa1.get_frame_from_timestamp_flydra(movie, movie.timestamps[-1])
    last_frame = np.min([last_frame, len(trajec.speed)-1]) 
    
    saccade_frames = []
    if len(trajec.saccades) > 0:
        for s in trajec.saccades:
            sac_range = fa.get_saccade_range(trajec, s)
            saccade_frames.extend(sac_range)
    straight_frames = []
    saccade_frames = []
    landing_frames = []
    print first_frame, last_frame
    for f in range(first_frame, last_frame):
        print f
        if f in saccade_frames:
            print 'x'
            saccade_frames.append(f)
        if trajec.dist_to_stim_r_normed[f] < closest_dist and trajec.behavior == 'landing':
            print 'y'
            landing_frames.append(f)
        else:
            straight_frames.append(f)
            
    if behavior == 'straight':
        frames = straight_frames
    elif behavior == 'saccade':
        frames = saccade_frames
    elif behavior == 'landing':
        if 
    
    slipangles = []
    for f in frames:
        sa1f = sa1.get_frame_from_timestamp(movie, trajec.epoch_time[f])
        slipangles.append(movie.scaled.slipangle[sa1f])
        
    for i, s in enumerate(slipangles):
        while slipangles[i] > np.pi:
            slipangles[i] -= np.pi
        while slipangles[i] < -np.pi:
            slipangles[i] += np.pi
        
    return slipangles
        
    
def get_slipangles_for_movie_dataset(movie_dataset):

    keys = movie_dataset.get_movie_keys()

    slipangles = []
    for key in keys:
        movie = movie_dataset.movies[key]
        try:
            tmp = get_slipangles(movie)
            slipangles.extend(tmp)
        except:
            print key
            
    slipangles_in_degrees = np.array(slipangles)*180/np.pi
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    bins = np.linspace(-90, 90, 20, endpoint=True)
    ax.hist(slipangles_in_degrees, bins=bins, edgecolor='none', facecolor='purple')

    ax.set_xlim(-90,90)
    fa.adjust_spines(ax, ['left', 'bottom'])
    xticks = [-90, -45, 0, 45, 90]
    ax.set_xticks(xticks)

    ax.set_xlabel('slipangle, deg')
    ax.set_ylabel('occurences')
    N = len(slipangles)
    n = len(keys)
    string = 'N='+str(N)+'\nn='+str(n)
    ax.text(45,60,string)
        
    fig.savefig('slipangles_nosaccade.pdf', format='pdf')
                
    return slipangles



