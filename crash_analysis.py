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


def make_crash_dataset(movie_dataset, example_dataset):
    dataset_crash = ffa.Dataset(like=example_dataset)
    keys = movie_dataset.get_movie_keys(behavior='landing', crash=True)
    
    for key in keys:
        movie = movie_dataset.movies[key]
        trajec = movie.trajec
        trajec.behavior = movie.behavior
        trajec.key = key
        fa.calc_frame_of_landing(trajec)
        fa.normalize_dist_to_stim_r(trajec)    
        d = np.max(np.max(trajec.dist_to_stim_r_normed, 0.08))
        trajec.frames = np.arange(fa.get_frame_at_distance(trajec, d), trajec.frame_of_landing).tolist()    
        fa.prep_trajectory(trajec)
    
        dataset_crash.trajecs.setdefault(key, trajec)
        
    fa.prep_dataset(dataset_crash)
    
    return dataset_crash


def crash_spagetti_plots(dataset):
    keys = dataset.trajecs.keys()
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    for key in keys:
        trajec = dataset.trajecs[key]
        alpha = 0.5
        linewidth = 0.5
        if key == '20101111_C001H001S0045':
            alpha = 1
            linewidth = 1
            
        ax.plot(trajec.positions[:,0], trajec.positions[:,1], '-', color='black', linewidth=linewidth, alpha=alpha)
        if trajec.frame_at_deceleration is not None:
            d = trajec.frame_at_deceleration   
            dec = patches.Circle( (trajec.positions[d,0], trajec.positions[d,1]), radius=0.002, facecolor='blue', edgecolor='none', alpha=alpha, zorder=100)
            ax.add_artist(dec)
            
    post = patches.Circle( (0, 0), radius=0.009565, facecolor='black', edgecolor='none', alpha=1)
    ax.add_artist(post)
    ax.hlines(-0.1, 0.1, 0.2, color='black', linewidth=1)
    ax.text(.15, -.102, '10cm', horizontalalignment='center', verticalalignment='top')
            
    ax.set_ylim(-.15,.15)
    ax.set_xlim(-.25, .25)
    ax.set_autoscale_on(False)
    fig.set_size_inches(7.2,7.2*3/5.)
    fig.subplots_adjust(bottom=0., top=1, right=1, left=0)
    ax.set_aspect('equal')
    ax.set_axis_off()
    fig.savefig('crashes.pdf', format='pdf')




def crash_analysis(dataset, dataset_landing, keys=None):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    keys = dataset.trajecs.keys()
        
    for key in keys:
        trajec = dataset.trajecs[key]
        ftp = np.arange(trajec.frames[0], trajec.frames[-1]).tolist()
        
        alpha = 0.1
        color = 'blue'
        if key == '20101111_C001H001S0045':
            alpha = 1
            color = 'black'
            
        if trajec.angle_at_deceleration*180/np.pi > 90:
            print key
        
        ax.plot( np.log(trajec.angle_subtended_by_post[ftp]), trajec.speed[ftp], color='black', linewidth=0.5, alpha=alpha)
        ax.plot( np.log(trajec.angle_at_deceleration), trajec.speed_at_deceleration, '.', color=color, alpha=0.8)
        
    fit, Rsq, x, y, yminus, yplus = fa.get_angle_vs_speed_curve(dataset_landing, plot=False)
    ax.plot( x, y, color='blue')
    ax.fill_between(x, yplus, yminus, color='blue', linewidth=0, alpha=0.2)
    
    fa.fix_angle_log_spine(ax, histograms=False)
    fig.savefig('crash_spagetti.pdf', format='pdf')
    
    
def landing_analysis_for_crash_comparison(dataset, keys=None):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    if keys is None:
        classified_keys = fa.get_classified_keys(dataset)
        keys = classified_keys['straight']
        #keys = dataset.trajecs.keys()
        
    for key in keys:
        trajec = dataset.trajecs[key]
        ftp = np.arange(trajec.frames[0]-25, trajec.frames[-1]).tolist()
        ax.plot( np.log(trajec.angle_subtended_by_post[ftp]), trajec.speed[ftp], color='black', linewidth=0.5, alpha=0.05)

    keys_to_plot = ['2_29065', '2_31060', '8_10323', '6_715']
    for key in keys_to_plot:
        trajec = dataset.trajecs[key]
        ftp = np.arange(trajec.frames[0]-25, trajec.frames[-1]).tolist()
        ax.plot( np.log(trajec.angle_subtended_by_post[ftp]), trajec.speed[ftp], color='black', linewidth=0.5)
        ax.plot( np.log(trajec.angle_at_deceleration), trajec.speed_at_deceleration, '.', color='black', alpha=1)
        
    fit, Rsq, x, y, yminus, yplus = fa.get_angle_vs_speed_curve(dataset, plot=False)
    ax.plot( x, y, color='blue')
    ax.fill_between(x, yplus, yminus, color='blue', linewidth=0, alpha=0.2)
        
    fa.fix_angle_log_spine(ax, histograms=False)
    fig.savefig('landing_for_crash_comparison.pdf', format='pdf')


