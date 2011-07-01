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



#20101030_C001H001S0019


def flyby_xy_spagetti(dataset, keys=None, show_all_poi=False, filename='flyby_xy_spagetti.pdf'):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_ylim(-.15,.15)
    ax.set_xlim(-.25, .25)
    ax.set_autoscale_on(False)
    if keys is None:
        classified_keys = fa.get_classified_keys(dataset)
        keys = classified_keys['straight']
        
    keys_with_saccades = []
    keys_without_saccades = []
    for key in keys:
        trajec = dataset.trajecs[key]    
        if len(trajec.saccades) > 0:
            keys_with_saccades.append(key)
        else:
            keys_without_saccades.append(key)
    '''
    key_to_highlight = None
    k = -1
    while key_to_highlight is None:
        k += 1
        trajec = dataset.trajecs[keys[k]]    
        if len(trajec.saccades) > 0:
            key_to_highlight = keys[k]
    '''
    
    key_to_highlight = keys_with_saccades[2]
    print key_to_highlight
    for key in keys_with_saccades:
        trajec = dataset.trajecs[key]    
        
        if key == key_to_highlight:
            alpha = 1
            linewidth = 1
        else:
            alpha = 0.2
            linewidth = 0.5
            
        if show_all_poi:    
            s = trajec.saccades[0]
            sac_range = fa.get_saccade_range(trajec, s) 
            ax.plot(trajec.positions[sac_range,0], trajec.positions[sac_range,1], '-', color='green', alpha=alpha, linewidth=1)
            sac = patches.Circle( (trajec.positions[sac_range[0],0], trajec.positions[sac_range[0],1]), radius=0.002, facecolor='green', edgecolor='none', alpha=alpha, zorder=100)
            ax.add_artist(sac)
            
            if trajec.frame_at_deceleration is not None:
                d = trajec.frame_at_deceleration   
                dec = patches.Circle( (trajec.positions[d,0], trajec.positions[d,1]), radius=0.002, facecolor='blue', edgecolor='none', alpha=alpha, zorder=100)
                ax.add_artist(dec)
            
        
        ax.plot(trajec.positions[trajec.frames_of_flyby,0], trajec.positions[trajec.frames_of_flyby,1], '-', color='black', alpha=alpha, linewidth=linewidth)
        
    
    trajec_to_highlight = dataset.trajecs[key_to_highlight]
    s = trajec_to_highlight.saccades[-1]
    d = trajec_to_highlight.frame_at_deceleration   
    sac_range = fa.get_saccade_range(trajec_to_highlight, s) 
    print 'saccade: ', sac_range, 'decel: ', d
    ax.plot(trajec_to_highlight.positions[sac_range,0], trajec_to_highlight.positions[sac_range,1], '-', color='green', alpha=1, linewidth=1)
    sac = patches.Circle( (trajec_to_highlight.positions[sac_range[0],0], trajec_to_highlight.positions[sac_range[0],1]), radius=0.002, facecolor='green', edgecolor='none', alpha=1, zorder=100)
    ax.add_artist(sac)
    
    dec = patches.Circle( (trajec_to_highlight.positions[d,0], trajec_to_highlight.positions[d,1]), radius=0.002, facecolor='blue', edgecolor='none', alpha=1, zorder=100)
    ax.add_artist(dec)
    
    print 'saccade: ', sac_range, 'decel: ', d
    
    post = patches.Circle( (0, 0), radius=0.009565, facecolor='black', edgecolor='none', alpha=1)
    ax.add_artist(post)
    #ax.text(0,0,'post\ntop view', horizontalalignment='center', verticalalignment='center')
    
    ax.hlines(-0.1, 0.1, 0.2, color='black', linewidth=1)
    ax.text(.15, -.102, '10cm', horizontalalignment='center', verticalalignment='top')
    
    #post = patches.Circle( (0, 0), radius=0.028, facecolor='none', edgecolor='green', alpha=1, linewidth=0.15, linestyle='dashed')
    #ax.add_artist(post)
    fig.set_size_inches(7.2,7.2*3/5.)
    fig.subplots_adjust(bottom=0., top=1, right=1, left=0)
    ax.set_aspect('equal')
    ax.set_axis_off()
    fig.savefig(filename, format='pdf')
    
def landing_xy_spagetti(dataset, keys=None, show_all_poi=False, key_class='straight'):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_ylim(-.15,.15)
    ax.set_xlim(-.25, .25)
    ax.set_autoscale_on(False)
    if keys is None:
        classified_keys = fa.get_classified_keys(dataset)
        keys = classified_keys[key_class]
        #keys = dataset.trajecs.keys()
    
    if key_class == 'straight':
        key_to_highlight = keys[4]
    elif key_class == 'other':
        key_to_highlight = keys[10]
    else:
        key_to_highlight = keys[0]
        
    for key in keys:
        trajec = dataset.trajecs[key]    
        
        if key == key_to_highlight:
            alpha = 1
            linewidth = 1
        else:
            alpha = 0.2
            linewidth = 0.5
            
        if show_all_poi:    
            
            if trajec.frame_at_deceleration is not None:
                d = trajec.frame_at_deceleration   
                dec = patches.Circle( (trajec.positions[d,0], trajec.positions[d,1]), radius=0.002, facecolor='blue', edgecolor='none', alpha=alpha, zorder=100)
                ax.add_artist(dec)
                
            if len(trajec.saccades) > 0:
                for s in trajec.saccades:
                    sac = patches.Circle( (trajec.positions[s,0], trajec.positions[s,1]), radius=0.002, facecolor='green', edgecolor='none', alpha=alpha, zorder=100)
                    ax.add_artist(sac)
            
        
        ax.plot(trajec.positions[:,0], trajec.positions[:,1], '-', color='black', alpha=alpha, linewidth=linewidth)
        
    
    trajec_to_highlight = dataset.trajecs[key_to_highlight]
    d = trajec_to_highlight.frame_at_deceleration   
    dec = patches.Circle( (trajec_to_highlight.positions[d,0], trajec_to_highlight.positions[d,1]), radius=0.002, facecolor='blue', edgecolor='none', alpha=1, zorder=100)
    ax.add_artist(dec)
    if len(trajec_to_highlight.saccades) > 0:
        for s in trajec_to_highlight.saccades:
            #s = trajec_to_highlight.saccades[-1]
            sac = patches.Circle( (trajec_to_highlight.positions[s,0], trajec_to_highlight.positions[s,1]), radius=0.002, facecolor='green', edgecolor='none', alpha=1, zorder=100)
            ax.add_artist(sac)
    
    post = patches.Circle( (0, 0), radius=0.009565, facecolor='black', edgecolor='none', alpha=1)
    ax.add_artist(post)
    #ax.text(0,0,'post\ntop view', horizontalalignment='center', verticalalignment='center')
    
    ax.hlines(-0.1, 0.1, 0.2, color='black', linewidth=1)
    ax.text(.15, -.102, '10cm', horizontalalignment='center', verticalalignment='top')
    
    #post = patches.Circle( (0, 0), radius=0.028, facecolor='none', edgecolor='green', alpha=1, linewidth=0.15, linestyle='dashed')
    #ax.add_artist(post)
    fig.set_size_inches(7.2,7.2*3/5.)
    fig.subplots_adjust(bottom=0., top=1, right=1, left=0)
    ax.set_aspect('equal')
    ax.set_axis_off()
    filename = 'landing_xy_spagetti_' + key_class + '.pdf'
    fig.savefig(filename, format='pdf')
    
def landing_xy_spagetti_types(dataset):
    
    classified_keys = fa.get_classified_keys(dataset)
    
    for key_class in classified_keys.keys():
        landing_xy_spagetti(dataset, keys=None, show_all_poi=False, key_class=key_class)
        
        
        
def flyby_xy_spagetti_types(dataset_flyby):
    classified_keys_flyby = fa.get_classified_keys(dataset_flyby)
    flyby_xy_spagetti(dataset_flyby, keys=classified_keys_flyby['flyby'][100:300], filename='flyby_xy_spagetti_random.pdf')
    flyby_xy_spagetti(dataset_flyby, keys=classified_keys_flyby['straight'], filename='flyby_xy_spagetti_straight.pdf')
