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


def flyby_xy_spagetti(dataset, keys=None, show_all_poi=False, filename='flyby_xy_spagetti.pdf', nkeys=300, alpha2=1, color2='gray', show_saccades=True, keys_to_highlight=None, show_deceleration=False):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_ylim(-.15,.15)
    ax.set_xlim(-.25, .25)
    ax.set_autoscale_on(False)
    if keys is None:
        classified_keys = fa.get_classified_keys(dataset)
        keys = classified_keys['straight']
        keys = dataset.trajecs.keys()
        
    keys_with_saccades = []
    keys_without_saccades = []
    for key in keys:
        trajec = dataset.trajecs[key]    
        if len(trajec.sac_ranges) > 0:
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
    
    if keys_to_highlight is None:
        keys_to_highlight = ['10_80113']
        
    print keys_to_highlight

    if len(keys) > nkeys:
        keys = keys[0:nkeys]
        keys.extend(keys_to_highlight)

    for key in keys:
        trajec = dataset.trajecs[key]    
        
        if key in keys_to_highlight:
            alpha = 1
            linewidth = 1
            color = 'black'
            zorder = 500
            ax.plot(trajec.positions[:,0], trajec.positions[:,1], '-', color='black', alpha=1, linewidth=linewidth, zorder=zorder)
            frames = np.arange(0,trajec.frames_of_flyby[-1]).tolist()
            ax.plot(trajec.positions[frames,0], trajec.positions[frames,1], '-', color='green', alpha=1, linewidth=linewidth, zorder=zorder)
        else:
            linewidth = 0.5
            zorder = 100
            
            frames = np.arange(0,trajec.frames_of_flyby[-1]).tolist()
            ax.plot(trajec.positions[frames,0], trajec.positions[frames,1], '-', color=color2, alpha=alpha2, linewidth=linewidth, zorder=zorder)
            
        if show_saccades:
            if len(trajec.sac_ranges) > 0:
                for sac_range in trajec.sac_ranges:
                
                    if sac_range == trajec.last_saccade_range:
                        color = 'red'
                    else:
                        color = 'blue'
                
                    if sac_range[0] < trajec.frame_nearest_to_post:
                    #if 1:
                        ax.plot(trajec.positions[sac_range,0], trajec.positions[sac_range,1], '-', color=color, alpha=alpha2, linewidth=linewidth, zorder=zorder+1)
                        #sac = patches.Circle( (trajec.positions[sac_range[0],0], trajec.positions[sac_range[0],1]), radius=0.001, facecolor='red', edgecolor='none', alpha=alpha2, zorder=zorder+1)
                        #ax.add_artist(sac)
                        
        if show_deceleration:
            dec = patches.Circle( (trajec.positions[trajec.frame_at_deceleration,0], trajec.positions[trajec.frame_at_deceleration,1]), radius=0.001, facecolor='purple', edgecolor='none', color='purple', linewidth=2, alpha=alpha2, zorder=zorder+1)
            ax.add_artist(dec)
            
    
    
    
    
    post = patches.Circle( (0, 0), radius=0.009565, facecolor='black', edgecolor='none', alpha=1)
    ax.add_artist(post)
    #ax.text(0,0,'post\ntop view', horizontalalignment='center', verticalalignment='center')
    
    
    #post = patches.Circle( (0, 0), radius=0.028, facecolor='none', edgecolor='red', alpha=1, linewidth=0.15, linestyle='dashed')
    #ax.add_artist(post)
    
    prep_xy_spagetti_for_saving(ax)
    fig.savefig(filename, format='pdf')
    
def landing_xy_spagetti(dataset, keys=None, show_all_poi=False, alpha2=1, color2='gray', show_saccades=True):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_ylim(-.15,.15)
    ax.set_xlim(-.25, .25)
    ax.set_autoscale_on(False)
    if keys is None:
        #classified_keys = fa.get_classified_keys(dataset)
        keys = dataset.trajecs.keys()
    
    #keys_to_highlight = [keys[33]]
    # straight: 2_29065, 2_3954, 5_15269
    # curve in: 24_68713
    # multiple saccades: 1_27887
    # hover loop: 10_78813
    
    #keys_to_highlight = ['5_15269', '24_68713', '10_78813', '1_27887']
    #keys_to_highlight = ['5_15269']
    #keys = keys_to_highlight
    keys_to_highlight = []
    print keys_to_highlight
    
    for key in keys:
        trajec = dataset.trajecs[key]    
        
        if key == '24_68713':
            flip = -1
        else:
            flip = 1
            
        initial_frame = fa.get_frame_at_distance(trajec, 0.25)
        frames = np.arange(initial_frame, len(trajec.speed)-1).tolist()
        #frames = np.arange(initial_frame, trajec.frame_of_landing).tolist()
            
        if key in keys_to_highlight:
            alpha = 1
            linewidth = 1
            color = 'black'
            zorder = 500
            ax.plot(flip*trajec.positions[frames,0], flip*trajec.positions[frames,1], '-', color=color, alpha=alpha, linewidth=linewidth, zorder=zorder)
        else:
            linewidth = 0.5
            zorder = 100
            
            ax.plot(flip*trajec.positions[frames,0], flip*trajec.positions[frames,1], '-', color=color2, alpha=alpha2, linewidth=linewidth, zorder=zorder)
            
        if show_saccades:
            if len(trajec.sac_ranges) > 0:
                for sac_range in trajec.sac_ranges:
                
                    if sac_range == trajec.last_saccade_range:
                        color = 'red'
                    else:
                        color = 'blue'    
                    
                    if sac_range[-1] in frames and sac_range[-1] < trajec.frame_of_landing:
                        ax.plot(flip*trajec.positions[sac_range,0], flip*trajec.positions[sac_range,1], '-', color=color, alpha=alpha2, linewidth=linewidth, zorder=zorder+1)
                        #sac = patches.Circle( (flip*trajec.positions[sac_range[0],0], flip*trajec.positions[sac_range[0],1]), radius=0.001, facecolor='red', edgecolor='none', alpha=alpha2, zorder=zorder+1)
                        #ax.add_artist(sac)
        
    '''
    trajec_to_highlight = dataset.trajecs[key_to_highlight]
    d = trajec_to_highlight.frame_at_deceleration   
    dec = patches.Circle( (trajec_to_highlight.positions[d,0], trajec_to_highlight.positions[d,1]), radius=0.002, facecolor='blue', edgecolor='none', alpha=1, zorder=100)
    ax.add_artist(dec)
    if len(trajec_to_highlight.saccades) > 0:
        for s in trajec_to_highlight.saccades:
            #s = trajec_to_highlight.saccades[-1]
            sac = patches.Circle( (trajec_to_highlight.positions[s,0], trajec_to_highlight.positions[s,1]), radius=0.002, facecolor='red', edgecolor='none', alpha=1, zorder=100)
            ax.add_artist(sac)
    '''
    
    post = patches.Circle( (0, 0), radius=0.009565, facecolor='black', edgecolor='none', alpha=1, zorder=1000)
    ax.add_artist(post)
    #ax.text(0,0,'post\ntop view', horizontalalignment='center', verticalalignment='center')
    
    prep_xy_spagetti_for_saving(ax)
    
    filename = 'landing_xy_spagetti' + '.pdf'
    fig.savefig(filename, format='pdf')
    
    
def prep_xy_spagetti_for_saving(ax):
    
    rect = patches.Rectangle( [-.25, -.15], .5, .3, facecolor='none', edgecolor='gray', clip_on=False, linewidth=0.2)
    ax.add_artist(rect)
    
    '''
    offset = 0.00
    dxy = 0.05
    #xarrow = patches.FancyArrowPatch(posA=(-.25+offset, -.15+offset), posB=(-.25+offset+dxy, -.15+offset), arrowstyle='simple') 
    #patches.Arrow( -.25+offset, -.15+offset, dxy, 0, color='black', width=0.002)
    xarrow = patches.FancyArrowPatch((-.25+offset, -.15+offset), (-.25+offset+dxy, -.15+offset), arrowstyle="-|>", mutation_scale=10, color='gray', shrinkA=0, clip_on=False)
    ax.add_patch(xarrow)
    yarrow = patches.FancyArrowPatch((-.25+offset, -.15+offset), (-.25+offset, -.15+offset+dxy), arrowstyle="-|>", mutation_scale=10, color='gray', shrinkA=0, clip_on=False)
    ax.add_artist(yarrow)
    text_offset = -.011
    ax.text(-.25+offset+dxy+text_offset, -.15+offset+.005, 'x', verticalalignment='bottom', horizontalalignment='left', color='gray', weight='bold')
    ax.text(-.25+offset+.005, -.15+offset+dxy+text_offset, 'y', verticalalignment='bottom', horizontalalignment='left', color='gray', weight='bold')
    '''
    
    scale_bar_offset = 0.01
    ax.hlines(-0.15+scale_bar_offset, 0.25-scale_bar_offset-.1, 0.25-scale_bar_offset, linewidth=1, color='gray')
    ax.text(0.25-scale_bar_offset-.1/2., -0.15+scale_bar_offset+.002, '10cm', horizontalalignment='center', verticalalignment='bottom', color='gray')
    
    ax.set_aspect('equal')
    
    scaling = .5/.75
    margin = 0.04
    aspect_ratio = 3/5. # height/width
    
    fig_width = 7.204*scaling
    plt_width = fig_width - 2*margin*(1-aspect_ratio)
    fig_height = plt_width*aspect_ratio + 2*margin
    
    fig = ax.figure
    
    fig.set_size_inches(fig_width,fig_height)
    fig.subplots_adjust(bottom=margin, top=1-margin, right=1, left=0)
    ax.set_axis_off()
    
def prep_radial_spagetti_for_saving(ax):
    
    rect = patches.Rectangle( [0, -.15], .25, .3, facecolor='none', edgecolor='gray', clip_on=False, linewidth=0.2)
    ax.add_artist(rect)
    
    '''
    offset = 0.00
    dxy = 0.05
    #xarrow = patches.FancyArrowPatch(posA=(-.25+offset, -.15+offset), posB=(-.25+offset+dxy, -.15+offset), arrowstyle='simple') 
    #patches.Arrow( -.25+offset, -.15+offset, dxy, 0, color='black', width=0.002)
    xarrow = patches.FancyArrowPatch((0+offset, -.15+offset), (0+offset+dxy, -.15+offset), arrowstyle="-|>", mutation_scale=10, color='gray', shrinkA=0, clip_on=False)
    ax.add_patch(xarrow)
    yarrow = patches.FancyArrowPatch((0+offset, -.15+offset), (0+offset, -.15+offset+dxy), arrowstyle="-|>", mutation_scale=10, color='gray', shrinkA=0, clip_on=False)
    ax.add_artist(yarrow)
    text_offset = -.011
    ax.text(0+offset+dxy+text_offset, -.15+offset+.005, 'r', verticalalignment='bottom', horizontalalignment='left', color='gray', weight='bold')
    ax.text(0+offset+.002, -.15+offset+dxy+text_offset+.004, 'z', verticalalignment='bottom', horizontalalignment='left', color='gray', weight='bold')
    '''
    
    scale_bar_offset = 0.01
    ax.hlines(-0.15+scale_bar_offset, 0.25-scale_bar_offset-.1, 0.25-scale_bar_offset, linewidth=1, color='gray')
    ax.text(0.25-scale_bar_offset-.1/2., -0.15+scale_bar_offset+.002, '10cm', horizontalalignment='center', verticalalignment='bottom', color='gray')
    
    ax.set_aspect('equal')
    
    scaling = .25/.75
    margin = 0.04
    aspect_ratio = 3/2.5 # height/width
    
    fig_width = 7.204*scaling
    plt_width = fig_width - 2*margin*(aspect_ratio)
    fig_height = plt_width*aspect_ratio + 2*margin
    
    fig = ax.figure
    
    fig.set_size_inches(fig_width,fig_height)
    fig.subplots_adjust(bottom=margin, top=1-margin, right=1, left=0)
    ax.set_axis_off()
    
    
    
def landing_xy_spagetti_types(dataset):
    
    classified_keys = fa.get_classified_keys(dataset)
    
    for key_class in classified_keys.keys():
        landing_xy_spagetti(dataset, keys=None, show_all_poi=False, key_class=key_class)
        
        
        
def flyby_xy_spagetti_types(dataset_flyby):
    classified_keys_flyby = fa.get_classified_keys(dataset_flyby)
    flyby_xy_spagetti(dataset_flyby, keys=classified_keys_flyby['flyby'][100:300], filename='flyby_xy_spagetti_random.pdf')
    flyby_xy_spagetti(dataset_flyby, keys=classified_keys_flyby['straight'], filename='flyby_xy_spagetti_straight.pdf')
    
    
    
    
    
    
    
def flyby_xy_spagetti_nopost(dataset, keys=None, filename='flyby_xy_spagetti_nopost.pdf', alpha2=1, color2='gray', show_saccades=True):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_ylim(-.15,.15)
    ax.set_xlim(-.25, .25)
    ax.set_autoscale_on(False)
    if keys is None:
        keys = dataset.trajecs.keys()
        
    
    keys_to_highlight = ['32_8003']
    print keys_to_highlight
    for key in keys:
        trajec = dataset.trajecs[key]    
        frames = np.arange(0,trajec.frames_of_flyby[-1]).tolist()
        
        if key in keys_to_highlight:
            alpha = 1
            linewidth = 1
            color = 'black'
            zorder = 500
            ax.plot(trajec.positions[:,0], trajec.positions[:,1], '-', color='black', alpha=1, linewidth=linewidth, zorder=zorder)
            ax.plot(trajec.positions[frames,0], trajec.positions[frames,1], '-', color='blue', alpha=1, linewidth=linewidth, zorder=zorder)
        else:
            linewidth = 0.5
            color = 'black'
            zorder = 100
            ax.plot(trajec.positions[frames,0], trajec.positions[frames,1], '-', color=color2, alpha=alpha2, linewidth=linewidth, zorder=zorder)
            
        if show_saccades:
            if len(trajec.sac_ranges) > 0:
                for sac in trajec.sac_ranges:
                
                    if sac == trajec.last_saccade_range:
                        color = 'red'
                    else:
                        color = 'blue'    
                
                    if sac[-1] in frames:
                        ax.plot(trajec.positions[sac,0], trajec.positions[sac,1], '-', color=color, alpha=alpha2, linewidth=linewidth, zorder=zorder+1)
                        #sac = patches.Circle( (trajec.positions[sac_range[0],0], trajec.positions[sac_range[0],1]), radius=0.001, facecolor='red', edgecolor='none', alpha=alpha2, zorder=zorder+1)
                        #ax.add_artist(sac)
    
        
        
    
    
    post = patches.Circle( (0, 0), radius=0.009565, facecolor='black', edgecolor='none', alpha=1)
    #ax.add_artist(post)
    #ax.text(0,0,'post\ntop view', horizontalalignment='center', verticalalignment='center')
    
    prep_xy_spagetti_for_saving(ax)
    filename = 'nopost_xy_spagetti' + '.pdf'
    fig.savefig(filename, format='pdf')
    
    
    
    
    
    
    
def radial_trajectory_plot(dataset, filename='radial_trajectory_plot.pdf', keys=None, dataset_type='landing', nkeys=None, highlight_trajecs=False, behavior=['landing', 'flyby'], show_saccades=True, alpha2=1, color2='gray'):
    radius = 0.009565
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_ylim(-.15,.15)
    ax.set_xlim(0, .25)
    ax.set_autoscale_on(False)
    if keys is None:
        #keys = dataset.trajecs.keys()
        keys = fa.get_keys_for_behavior(dataset, behavior=behavior)
        if len(keys) > 300:
            keys = keys[0:300]
        
    if highlight_trajecs:
        if dataset_type == 'landing':
            keys_to_highlight = ['5_15269', '24_68713', '10_78813', '1_27887']
        elif dataset_type == 'flyby':
            keys_to_highlight = ['10_80113']
        else:
            keys_to_highlight = ['32_8003']
    else:
        keys_to_highlight = []
    
    if nkeys is not None:
        if nkeys < len(keys):
            keys = keys[0:nkeys]
            for key in keys_to_highlight:
                if key not in keys:
                    keys.append(key)
    
    print keys_to_highlight
    for key in keys:
        trajec = dataset.trajecs[key]    
        r = trajec.dist_to_stim_r + radius
        
        if trajec.behavior == 'landing':
            initial_frame = fa.get_frame_at_distance(trajec, 0.25)
            frames = np.arange(initial_frame, trajec.frame_of_landing).tolist()
        else:
            #try:
            #    #frames = trajec.frames_below_post
            #    frames = np.arange(0,trajec.frame_of_landing).tolist()
            #except:
            frames = np.arange(0,trajec.frames_of_flyby[-1]).tolist()
        
        if key in keys_to_highlight:
            alpha = 1
            linewidth = 1
            color = 'black'
            zorder = 500
            #ax.plot(r[:], trajec.positions[:,2], '-', color='black', alpha=1, linewidth=linewidth, zorder=zorder)
            ax.plot(r[:], trajec.positions[:,2], '-', color=color, alpha=alpha, linewidth=linewidth, zorder=zorder)
            ax.plot(r[frames], trajec.positions[frames,2], '-', color='blue', alpha=alpha, linewidth=linewidth, zorder=zorder)
            print 'plotted key to highlight'
        else:
            linewidth = 0.5
            zorder = 100
            ax.plot(r[frames], trajec.positions[frames,2], '-', color=color2, alpha=alpha2, linewidth=linewidth, zorder=zorder)
        
        if show_saccades:
            if len(trajec.sac_ranges) > 0:
                for sac_range in trajec.sac_ranges:
                    if sac_range == trajec.last_saccade_range:
                        ax.plot(r[sac_range], trajec.positions[sac_range,2], '-', color='red', alpha=alpha2, linewidth=linewidth, zorder=zorder+1)
                    elif len(trajec.last_saccade_range) > 0:
                        if sac_range[0] in frames:
                            if sac_range[0] < trajec.last_saccade_range[0]:
                                ax.plot(r[sac_range], trajec.positions[sac_range,2], '-', color='blue', alpha=alpha2, linewidth=linewidth, zorder=zorder+1)
                        
    if dataset_type != 'nopost':
        post = patches.Rectangle( (0, -0.15), radius, 0.15, facecolor='black', edgecolor='none', alpha=1)
        ax.add_artist(post)
    
    
    prep_radial_spagetti_for_saving(ax)
    filename = 'radial_spagetti' + dataset_type + '.pdf'
    fig.savefig(filename, format='pdf')
    
    
    
    
