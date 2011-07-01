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

def prep_movie_trajec(movieinfo):
    trajec = movieinfo.trajec
    trajec.key = movieinfo.id
    trajec.frames = np.arange(fa.get_frame_at_distance(trajec, 0.08), trajec.frame_of_landing).tolist()
    fa.prep_trajectory(trajec)

def plot_landing_demo_trajectory(movie_dataset, filename=None):
    movieinfo = movie_dataset.movies['20101110_C001H001S0038']
    nudge = [0, -.003] # for SA1 data
    trajec = movieinfo.trajec
    trajec.key = movieinfo.id
    trajec.frames = np.arange(fa.get_frame_at_distance(trajec, 0.08), trajec.frame_of_landing).tolist()
    fa.prep_trajectory(trajec)
    plot_demo_trajectory(movieinfo, nudge, filename=filename)
    
    legim = movieinfo.frames[movieinfo.legextensionrange[0] - movieinfo.firstframe_ofinterest + 50].uimg
    plt.imsave('legextensionimg', legim, cmap=plt.get_cmap('gray')) 
    
def plot_flyby_demo_trajectory(movie_dataset, filename=None):
    movieinfo = movie_dataset.movies['20101101_C001H001S0024']
    nudge = [0, -.001] # for SA1 data
    trajec = movieinfo.trajec
    trajec.key = movieinfo.id
    frame_nearest_to_post = np.argmin(trajec.dist_to_stim_r)
    fa.calc_radius_at_nearest(trajec)    
    frames = np.arange(0, frame_nearest_to_post).tolist()
    frame_at_distance = fa.get_frame_at_distance(trajec, 0.08, singleframe=True, frames=frames)
    last_frame = np.min( [frame_nearest_to_post+20, len(trajec.speed)-1]) 
    trajec.frames_of_flyby = np.arange(frame_at_distance, last_frame).tolist()
    trajec.frames = trajec.frames_of_flyby
    fa.prep_trajectory(trajec)
    plot_demo_trajectory(movieinfo, nudge, filename=filename)
    
def plot_demo_trajectory(movieinfo, nudge, filename=None):
    trajec = movieinfo.trajec

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_autoscale_on(False)
    ax.set_xlim(-0.02, 0.05)
    ax.set_ylim(-0.02, 0.05)
    ax.set_axis_off()
    cl = colorline.Colorline( xlim=ax.get_xlim(), ylim=ax.get_ylim(), norm=(0,0.8), colormap = 'jet', hide_colorbar=True, ax0=ax)
    # plot trajectory
    frames = np.arange(50, trajec.frame_of_landing).tolist()
    cl.colorline( trajec.positions[frames,0], trajec.positions[frames,1], trajec.speed[frames],linewidth=1, alpha=0.5)

    # parameters
    radius = trajec.radius_at_nearest
    
    if trajec.behavior == 'landing':
        # show pt of deceleration
        f = trajec.frame_at_deceleration
        x = trajec.positions[f,0]
        y = trajec.positions[f,1]
        pt_of_deceleration = patches.Circle( (x,y), radius=0.0005, facecolor='blue', edgecolor='none')
        ax.add_artist(pt_of_deceleration)
        
        # anotate deceleration point
        string = 'deceleration initiation'
        string_position = (x-0.005, y+.015)
        arrow = {'facecolor':'blue', 'arrowstyle':"->", 'edgecolor':'blue', 'linewidth':1}
        ax.annotate(string, (x-.0001,y+.0001),
            xytext=string_position,
            arrowprops=arrow,
            horizontalalignment='right', verticalalignment='top', color='blue')
        
    # show angle subtended
    if trajec.behavior == 'landing':
        f = trajec.frames[20]
        x = np.mean(trajec.positions[f-1:f+1,0])
        y = np.mean(trajec.positions[f-1:f+1,1])
    else:   
        f = trajec.frames[17]
        x = trajec.positions[f,0]
        y = trajec.positions[f,1]
    d = np.linalg.norm([x,y])
    half_angle_to_post = np.arcsin( radius / d )
    world_angle = np.arctan2(y,x)
    
    a = half_angle_to_post - world_angle
    visual_intercept_1 = [0+np.cos(np.pi/2.-a)*radius, 0+np.sin(np.pi/2.-a)*radius]
    
    a = half_angle_to_post + world_angle
    visual_intercept_2 = [0+np.cos(np.pi/2.-a)*radius, 0-np.sin(np.pi/2.-a)*radius]
    
    xy = np.vstack( (visual_intercept_1, visual_intercept_2, [x,y]) )
    triangle = patches.Polygon( xy, facecolor='purple', edgecolor='none', zorder=-10, alpha=0.2 )
    ax.add_artist(triangle)
    
    # dist to post
    ox = np.cos(world_angle)*radius
    oy = np.sin(world_angle)*radius
    ax.plot([ox,x], [oy,y], color='brown', linewidth=1, alpha=0.5)
    
        
    # angle arc
    if trajec.behavior == 'landing':
        arc = patches.Arc( (x,y), .03, .03, 180, np.abs((half_angle_to_post - world_angle))*180/np.pi, (half_angle_to_post + world_angle)*180/np.pi, edgecolor='purple', linewidth=0.5)
    else:
        arc = patches.Arc( (x,y), .05, .05, 180, (world_angle - half_angle_to_post)*180/np.pi, (half_angle_to_post + world_angle)*180/np.pi, edgecolor='purple', linewidth=1)
    ax.add_artist(arc)
    if trajec.behavior == 'landing':
        ax.text(.02, .006, 'retinal size', color='purple', horizontalalignment='center', verticalalignment='center')
        ax.text(.02, 0, 'distance\nto post', color='brown', horizontalalignment='center', verticalalignment='center')
    else:
        ax.text(.02, -0.01, 'retinal size', color='purple', horizontalalignment='center', verticalalignment='center')
        ax.text(.02, -0.015, 'distance\nto post', color='brown', horizontalalignment='center', verticalalignment='center')
    
    if trajec.behavior == 'landing': 
        # show pt of leg extension
        legextensionframe = movieinfo.legextensionrange[0] - movieinfo.firstframe_ofinterest
        time_of_leg_ext = movieinfo.timestamps[legextensionframe]
        flydra_frame_of_leg_ext = fa.get_flydra_frame_at_timestamp(trajec, time_of_leg_ext)+1
        x = trajec.positions[flydra_frame_of_leg_ext, 0]
        y = trajec.positions[flydra_frame_of_leg_ext, 1]
        pt_of_leg_ext = patches.Circle( (x,y), radius=0.0005, facecolor='red', edgecolor='none')
        ax.add_artist(pt_of_leg_ext)

        # anotate leg extension point
        string = 'leg extension'
        string_position = (x, y+.01)
        arrow = {'facecolor':'red', 'arrowstyle':"->", 'edgecolor':'red', 'linewidth':1}
        ax.annotate(string, (x-.0001,y+.0001),
            xytext=string_position,
            arrowprops=arrow,
            horizontalalignment='center', verticalalignment='top', color='red')
            
                    
    if trajec.behavior == 'flyby':
        '''
        s = trajec.saccades[0]
        x = trajec.positions[s, 0]
        y = trajec.positions[s, 1]
        pt_of_saccade = patches.Circle( (x,y), radius=0.0005, facecolor='green', edgecolor='none')
        ax.add_artist(pt_of_saccade)
        '''
        
        # show saccade duration
        saccade_frames = fa.get_saccade_range(trajec, trajec.saccades[-1])
        x = trajec.positions[saccade_frames, 0]
        y = trajec.positions[saccade_frames, 1]
        ax.plot(x,y,'-', color='green')
        
        s = saccade_frames[0]
        x = trajec.positions[s, 0]
        y = trajec.positions[s, 1]
        pt_of_saccade = patches.Circle( (x,y), radius=0.0005, facecolor='green', edgecolor='none')
        ax.add_artist(pt_of_saccade)
        
        # anotate saccade point
        string = 'saccade initiation'
        string_position = (x+0.01, y+.01)
        arrow = {'facecolor':'green', 'arrowstyle':"->", 'edgecolor':'green', 'linewidth':1}
        ax.annotate(string, (x-.0001,y+.0001),
            xytext=string_position,
            arrowprops=arrow,
            horizontalalignment='right', verticalalignment='top', color='green')
            
            
            
        if trajec.frame_at_deceleration is not None:
            x = trajec.positions[trajec.frame_at_deceleration, 0]
            y = trajec.positions[trajec.frame_at_deceleration, 1]
            pt_of_deceleration = patches.Circle( (x,y), radius=0.0005, facecolor='blue', edgecolor='none')
            ax.add_artist(pt_of_deceleration)
            
            # anotate deceleration point
            string = 'deceleration initiation'
            string_position = (x-0.005, y+.015)
            arrow = {'facecolor':'blue', 'arrowstyle':"->", 'edgecolor':'blue', 'linewidth':1}
            ax.annotate(string, (x-.0001,y+.0001),
                xytext=string_position,
                arrowprops=arrow,
                horizontalalignment='right', verticalalignment='top', color='blue')
    
    post = patches.Circle( (0,0), radius=radius, facecolor='black', edgecolor='none')
    ax.add_artist(post)
    
    # annotate post
    ax.text( 0,0, 'post\n(top view)', color='gray', verticalalignment='center', horizontalalignment='center', fontsize=9)
    
    if movieinfo is not None:
        strobe = sa1.strobe_from_movieinfo(movieinfo, interval=200)
        ax.imshow(strobe.T, plt.get_cmap('gray'), origin='lower', extent = [-1*movieinfo.post_pos[0]*movieinfo.scale+nudge[0], (1024-movieinfo.post_pos[0])*movieinfo.scale+nudge[0], -1*movieinfo.post_pos[1]*movieinfo.scale+nudge[1], (1024-movieinfo.post_pos[1]+nudge[1])*movieinfo.scale], zorder=-20)
    
    
    
    
    # scale bar
    ax.hlines(-.012, -.005+.04, .005+.04, color='black')
    ax.vlines(-.005+.04, -.012+.0005, -.012-.0005, color='black')
    ax.vlines(.005+.04, -.012+.0005, -.012-.0005, color='black')
    ax.text( .04, -.0125, '1 cm', horizontalalignment='center', verticalalignment='top', fontsize=9)
    
    
    ax.set_aspect('equal')
    fig.subplots_adjust(bottom=0., top=1, right=1, left=0)
    #fig.set_size_inches(6.5,3.25)
    fig.savefig(filename, format='pdf')
        
    return 
    
    
    
    
def plot_demo_trajectory_speed_vs_angle(movie_dataset=None, movie_landing=None, movie_flyby=None, filename=None):
    
    if movie_dataset is not None:
        movie_landing = movie_dataset.movies['20101110_C001H001S0038']
        movie_flyby = movie_dataset.movies['20101101_C001H001S0024']
    
    fig = plt.figure()
    angleax = fig.add_subplot(111)

    
    ## landing    
    trajec = movie_landing.trajec
    ftp = np.arange(trajec.frames[0]-25, trajec.frames[-1]).tolist()
    angleax.plot( np.log(trajec.angle_subtended_by_post[ftp]), trajec.speed[ftp], color='black')
    angleax.plot( np.log(trajec.angle_at_deceleration), trajec.speed_at_deceleration, '.', markerfacecolor='blue', markeredgecolor='blue')
    
    
    movieinfo = movie_landing
    legextensionframe = movieinfo.legextensionrange[0] - movieinfo.firstframe_ofinterest
    time_of_leg_ext = movieinfo.timestamps[legextensionframe]
    flydra_frame_of_leg_ext = fa.get_flydra_frame_at_timestamp(trajec, time_of_leg_ext)+1
    angleax.plot( np.log(trajec.angle_subtended_by_post[flydra_frame_of_leg_ext]), trajec.speed[flydra_frame_of_leg_ext], '.', markerfacecolor='red', markeredgecolor='red')
    
    ## flyby
    trajec = movie_flyby.trajec
    ftp = np.arange(trajec.frames[0]-25, trajec.frames[-1]).tolist()
    angleax.plot( np.log(trajec.angle_subtended_by_post[ftp]), trajec.speed[ftp], '-', color='black')

    sf = fa.get_saccade_range(trajec, trajec.saccades[-1]) 
    angleax.plot( np.log(trajec.angle_subtended_by_post[sf[0]]), trajec.speed[sf[0]], '.', markerfacecolor='green', markeredgecolor='green')
    angleax.plot( np.log(trajec.angle_subtended_by_post[sf]), trajec.speed[sf], '-', color='green')
    angleax.plot( np.log(trajec.angle_at_deceleration), trajec.speed_at_deceleration, '.', markerfacecolor='blue', markeredgecolor='blue')
    
    fa.fix_angle_log_spine(angleax, histograms=False) 
    
    
    if filename is not None:
        fig.subplots_adjust(right=0.9, bottom=0.3)
        fig.savefig(filename, format='pdf')
    return angleax
    
    
    
    
    
    
    
def plot_demo_trajectory_speed_vs_dist(movie_dataset=None, movie_landing=None, movie_flyby=None):
    
    if movie_dataset is not None:
        movie_landing = movie_dataset.movies['20101110_C001H001S0038']
        movie_flyby = movie_dataset.movies['20101101_C001H001S0024']
    
    fig = plt.figure()
    fig.set_facecolor('white')
    distax = fig.add_subplot(111)

    ## landing    
    trajec = movie_landing.trajec
    ftp = np.arange(trajec.frames[0]-25, trajec.frames[-1]).tolist()
    distax.plot( -1*(trajec.dist_to_stim_r_normed[ftp]), trajec.speed[ftp], '-', color='black')
    distax.plot( -1*(trajec.dist_at_deceleration), trajec.speed_at_deceleration, 'o', markerfacecolor='blue', markeredgecolor='none')
    
    movieinfo = movie_landing
    legextensionframe = movieinfo.legextensionrange[0] - movieinfo.firstframe_ofinterest
    time_of_leg_ext = movieinfo.timestamps[legextensionframe]
    flydra_frame_of_leg_ext = fa.get_flydra_frame_at_timestamp(trajec, time_of_leg_ext)+1
    distax.plot( -1*(trajec.dist_to_stim_r_normed[flydra_frame_of_leg_ext]), trajec.speed[flydra_frame_of_leg_ext], 'o', markerfacecolor='red', markeredgecolor='none')
    
    ## flyby
    trajec = movie_flyby.trajec
    ftp = np.arange(trajec.frames[0]-25, trajec.frames[-1]).tolist()
    distax.plot( -1*(trajec.dist_to_stim_r_normed[ftp]), trajec.speed[ftp], '--', color='black')

    s = trajec.saccades[0]    
    distax.plot( -1*(trajec.dist_to_stim_r_normed[s]), trajec.speed[s], 'o', markerfacecolor='green', markeredgecolor='none')
    
    ## plot paramters    
    adjust_spines(distax, ['left', 'bottom'], color={'bottom':'brown'})
    distax.set_xlabel('distance to post, cm', color='brown')
    
    distax.set_ylim([0,0.8])
    distax.set_xlim([-0.1,0])
    
    dist_ticks = np.array([-0.1, -0.08, -0.06, -0.04, -0.02, 0])
    distax.set_xticks(dist_ticks.tolist())
    ticklabels = [str(d) for d in (dist_ticks*100)]
    distax.set_xticklabels(ticklabels, color='brown') 
    
    distax.set_ylabel('speed, m/s')
    
    return distax
    
    
    
    
    
def save(fig, filename):
    fig.savefig(filename, format='pdf', bbox='tight')
    
    
    
    
    
    
    
    
    
    
################
def adjust_spines(ax,spines, color={}, spine_locations={}, smart_bounds=False):
    if type(spines) is not list:
        spines = [spines]
    spine_locations_dict = {'top': 10, 'right': 10, 'left': 10, 'bottom': 10}
    for key in spine_locations.keys():
        spine_locations_dict[key] = spine_locations[key]
        
    if 'none' in spines:
        for loc, spine in ax.spines.iteritems():
            spine.set_color('none') # don't draw spine
        ax.yaxis.set_ticks([])
        ax.xaxis.set_ticks([])
        return
    
    for loc, spine in ax.spines.iteritems():
        if loc in spines:
            spine.set_position(('outward',spine_locations_dict[loc])) # outward by x points
            
            if loc in color.keys():
                c = color[loc]
            else:
                c = 'black'
            
            spine.set_color(c)
            #spine.set_smart_bounds(True)
            if loc == 'bottom' and smart_bounds:
                spine.set_smart_bounds(True)
        else:
            spine.set_color('none') # don't draw spine

    # turn off ticks where there is no spine
    if 'left' in spines:
        ax.yaxis.set_ticks_position('left')
    elif 'right' in spines:
        ax.yaxis.set_ticks_position('right')
    else:
        # no yaxis ticks
        ax.yaxis.set_ticks([])

    if 'bottom' in spines:
        ax.xaxis.set_ticks_position('bottom')
    else:
        # no xaxis ticks
        ax.xaxis.set_ticks([])
        
        
        
def example():
    
    # example sin function for plotting
    x = np.linspace(0, np.pi, 100)
    y = 0.9*np.sin(x)
    
    # set up a figure and axes instance to use for plotting
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    # set limits, and aspect of axis
    ax.set_autoscale_on(False)
    ax.set_xlim(0,3)
    ax.set_ylim(0,1)
    ax.set_aspect('equal') # make the x and y axis with the same scale
    
    ######################
    #    basic plot      #
    ax.plot(x,y, color='black')
    ######################
    
    ######################
    # example annotation #
    string = 'look, the top\nof the graph'
    string_position = (np.pi/2., 0.5)
    arrow_position = (np.pi/2,0.9)
    arrow = {'facecolor':'red', 'arrowstyle':"->", 'edgecolor':'red'}
    ax.annotate(string, arrow_position,
        xytext=string_position,
        arrowprops=arrow,
        horizontalalignment='center', verticalalignment='center', color='red')
    ######################
    
    ######################
    # example shape      #
    pts = np.vstack( ([.2,.5], [.5,.8], [2,.1]) )
    triangle = patches.Polygon( pts, facecolor='purple', edgecolor='none', zorder=-10, alpha=0.2 )
    circle = patches.Circle( (2.5,.2), 0.1, facecolor='green', edgecolor='green', alpha=1)
    ax.add_artist(triangle)
    ax.add_artist(circle)
    ######################
    
    #######################
    # example text        #
    string = 'hi, i am a\npurple triangle'
    ax.text(0.9, 0.2, string, horizontalalignment='center', verticalalignment='center', color='purple')
    #######################
    
    
    #######################
    # fix and label stuff #
    
    # note, this call to adjust_spines needs to happen before you label the axes or edit the ticklabels
    # the color and spine_locations dicts allow you to adjust the color and location of individual spines. if you leave them blank they default to black and 10. 
    adjust_spines(ax, ['left', 'bottom'], color={'left':'purple'}, spine_locations={'bottom':20}) 
    
    # x labels and ticks
    ax.set_xlabel('time, s')

    # y labels and ticks
    yticks = ax.get_yticks()
    # make custom ticks for y axis
    nticks = 3
    newyticks = np.linspace(yticks[0], yticks[-1], nticks, endpoint=True)
    ax.set_yticks(newyticks)
    ax.set_yticklabels([str(s) for s in newyticks], color='purple')
    ax.set_ylabel('sin(t)', color='purple')
    #######################
    
    
    #######################
    # save                #
    fig.savefig('example_figure.pdf', format='pdf') 
    # use pdf, as eps has issues with text output, and cannot handle opacity (alpha) settings. these pdfs are layered and perfectly editable in adobe illustrator.. from there you can output an eps file.
    #######################
    
    return fig
