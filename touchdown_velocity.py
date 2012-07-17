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
import floris



def plot_touchdown_velocity(movie_dataset):
    
    landing_keys = movie_dataset.get_movie_keys(behavior='landing')
    
    legextension_time = []
    touchdown_velocities = []
    touchdown_position = []
    sitting_position = []
    touchdown_time = []
    angle = []
    for key in landing_keys:
        movie = movie_dataset.movies[key]
        try:
            s = movie.scaled
        except:
            s = None
            print key, ': excepted'
        if movie.landingframe is not None and 'touchandgo' not in movie.subbehavior and movie.trajec is not None and s is not None:
            touchdown_velocities.append(movie.scaled.speed[movie.landingframe_relative])
            touchdown_position.append(movie.scaled.dist_to_post[movie.landingframe_relative][0])
            sitting_position.append( np.min(movie.scaled.dist_to_post) )
            legextensionframe = movie.legextensionrange[0] - movie.firstframe_ofinterest
            touchdown_time.append(movie.timestamps[movie.landingframe_relative] - movie.timestamps[legextensionframe])
            
            legextensionframe = movie.legextensionrange[0] - movie.firstframe_ofinterest
            legextension_time.append(movie.timestamps[legextensionframe])
            
            a = movie.scaled.signed_angletopost[movie.landingframe_relative]
            while a < -np.pi:
                a += np.pi
            while a > np.pi:
                a -= np.pi
            angle.append(a)
        else:
            print key
            
    print 'n: ', len(touchdown_velocities)
    
    time_after_leg_extension = np.array(touchdown_time)*1000
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ma = np.min(time_after_leg_extension)
    mi = np.max(time_after_leg_extension)
    #bins = np.linspace(ma, mi, 10, endpoint=True)
    bins, hists, curves = floris.histogram(ax, [time_after_leg_extension], bins=10, colors='black', bin_width_ratio=0.9, edgecolor='none', bar_alpha=0.8, curve_fill_alpha=0, curve_line_alpha=0, return_vals=True, show_smoothed=False)
    ax.set_ylim(0,np.max(hists))
    ax.set_xlim(0,400)
    ax.set_autoscale_on(False)
    fa.adjust_spines(ax, ['left', 'bottom'])
    ax.set_xlabel('Time to touchdown, ms')
    ax.set_ylabel('Occurences')
    filename = 'time_to_touchdown_at_leg_ext.pdf'
    fig.savefig(filename, format='pdf')
    print 'mean time to touchdown after leg ext: ', np.mean(touchdown_velocities)
            
    touchdown_velocities = np.array(touchdown_velocities)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.hist(touchdown_velocities)
    filename = 'touchdown_velocities.pdf'
    fig.savefig(filename, format='pdf')
    print 'mean touchdown velocity: ', np.mean(touchdown_velocities), '+/-', np.std(touchdown_velocities)

    dist_travelled_in_touchdown = np.array(touchdown_position) - np.array(sitting_position)
    indices = np.where(dist_travelled_in_touchdown > 0)[0].tolist()
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.hist(dist_travelled_in_touchdown[indices])
    filename = 'touchdown_dist_travelled.pdf'
    fig.savefig(filename, format='pdf')
    print 'mean dist travelled: ', np.mean(dist_travelled_in_touchdown)
    
    time_travelled_in_touchdown = np.array(touchdown_time)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.hist(time_travelled_in_touchdown)
    filename = 'time_travelled_in_touchdown.pdf'
    fig.savefig(filename, format='pdf')
    
    angle = np.array(angle)*180/np.pi
    fig = plt.figure()
    ax = fig.add_subplot(111)
    bins = np.linspace(-90,90,20)
    ax.hist(angle, bins=bins)
    filename = 'angle_to_post_at_touchdown.pdf'
    fig.savefig(filename, format='pdf')

    accel = touchdown_velocities[indices]**2 / (2*dist_travelled_in_touchdown[indices])
    mass = 0.001
    force = accel*mass
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.hist(force*10**4)
    filename = 'touchdown_deceleration.pdf'
    fig.savefig(filename, format='pdf')
    print 'mean touchdown deceleration: ', np.mean(accel)
    print 'mean touchdown force: ', np.mean(force)




