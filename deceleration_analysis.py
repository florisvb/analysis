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

def landing_spagetti_plots(dataset, gain=[170000], keys=None):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    fig2 = plt.figure()
    ax2 = fig2.add_subplot(111)
    
    fig3 = plt.figure()
    ax3 = fig3.add_subplot(111)
    
    if keys is None:
        classified_keys = fa.get_classified_keys(dataset)
        keys = classified_keys['straight']
        #keys = dataset.trajecs.keys()
        
    #keys_to_highlight = [keys[2], keys[0], keys[4], keys[6]]
    keys_to_highlight = ['2_29065', '2_31060', '8_10323', '6_715']
    points_at_deceleration = True
        
    for key in keys:
        trajec = dataset.trajecs[key]
        ftp = np.arange(trajec.frames[0]-25, trajec.frame_of_landing-1).tolist()
        #ftp = np.arange(trajec.frame_at_deceleration-10, trajec.frame_of_landing-1).tolist()
        if key in keys_to_highlight:
            color = 'black'
            zorder = 10
        else:
            color = 'blue'
            zorder = 1
        ax.plot( np.log(trajec.angle_subtended_by_post[ftp]), trajec.speed[ftp], color=color, linewidth=0.5, alpha=1, zorder=zorder)
        if points_at_deceleration and key in keys_to_highlight:
            ax.plot( np.log(trajec.angle_at_deceleration), trajec.speed_at_deceleration, '.', color='black', alpha=1, zorder=zorder)
        
        ax3.plot( np.log(trajec.angle_subtended_by_post[ftp]), trajec.speed[ftp], color='black', linewidth=0.5, alpha=1)
        
        if 0:
            # simulate
            angle, speed, expansion = sim_deceleration(trajec, gain)
            if speed[0] > 0.1:        
                ax2.plot( np.log(angle), speed, color='black', linewidth=0.5, alpha=1)
                ax3.plot( np.log(angle), speed, color='black', linewidth=0.5, alpha=1)
                f = np.where( np.abs(speed-speed[0])>0 )[0][0]
                
                #ax2.plot( np.log(angle[f]), speed[f], '.', color='blue')
            
    fit, Rsq, x, y, yminus, yplus = fa.get_angle_vs_speed_curve(dataset, plot=False)
    ax2.plot( x, y, color='purple', alpha=1)
    ax2.fill_between(x, yplus, yminus, color='purple', linewidth=0, alpha=0.2)

    ax.plot( x, y, color='purple', alpha=1)
    ax.fill_between(x, yplus, yminus, color='purple', linewidth=0, alpha=0.2)
        
    fa.fix_angle_log_spine(ax, histograms=False)
    fig.savefig('deceleration_real.pdf', format='pdf')
    
    fa.fix_angle_log_spine(ax2, histograms=False)
    fig2.savefig('deceleration_sim.pdf', format='pdf')
    
    fa.fix_angle_log_spine(ax3, histograms=False)
    fig3.savefig('deceleration_comparison.pdf', format='pdf')
    
    
#### linear scale to show not the same as srini paper

def landing_srini_plots(dataset, keys=None):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    
    if keys is None:
        classified_keys = fa.get_classified_keys(dataset)
        keys = classified_keys['straight']
        #keys = dataset.trajecs.keys()
        
    #keys = [keys[2], keys[0], keys[4], keys[6]]
        
    for key in keys:
        trajec = dataset.trajecs[key]
        ftp = np.arange(trajec.frame_at_deceleration-10, trajec.frame_of_landing-1).tolist()
        ax.plot( trajec.dist_to_stim_r_normed[ftp], trajec.speed[ftp], color='black', linewidth=0.5, alpha=0.2)
        # simulate
        
    ax.set_xlim([0,.1])
            
    fig.savefig('deceleration_linear_scale.pdf', format='pdf')
    
    
    
def sim_deceleration(trajec, gain, constant_vel=False):  
    # plot constant velocity comparison
    fps = 2000.
    dt = 1/fps
    r = 0.009565
    r = 0.009565
    radius = r
    
    
    nf = 5000
    positions = np.zeros([nf])
    speed = np.zeros([nf])
    distance = np.zeros([nf])
    angle = np.zeros([nf])
    
    for i in range(2):
        positions[i] = -0.22#-1*(trajec.dist_to_stim_r[trajec.frame_at_deceleration]+radius)
        speed[i] = trajec.speed[trajec.frame_at_deceleration]
        #angle[i] = trajec.angle_subtended_by_post[trajec.frame_at_deceleration]
        distance[i] = np.linalg.norm(positions[i]) - radius
        angle[i] = 2*np.arcsin( radius / (distance[i]+radius) )
    frames = [0,2]
    
    
    for f in range(frames[-1],nf-1): 
        if np.linalg.norm(positions[f-1])-radius <= 0.001:
            landed = True
            break
        else:
            landed = False
            
            
        if not landed:
            frames.append(f)
            positions[f] = positions[f-1] + speed[f-1]*dt
            distance[f] = np.linalg.norm(positions[f]) - radius
            angle[f] = 2*np.arcsin( radius / (distance[f]+radius) )
            
            a = angle
            exp0 = (a[f]-a[f-1])/dt #/ (-2.*np.tan(a[f]/2.))
            exp1 = (a[f-1]-a[f-2])/dt #/ (-2.*np.tan(a[f-1]/2.))
            
            
            m = -0.21/radius
            b = 0.159/radius
            expthreshold = (m*np.log(a[f])+b)*(2*np.tan(a[f]/2.)*np.sin(a[f]/2.))
            #tti_threshold = 0.15
            #expthreshold = 2*np.tan(af/2.) / tti_threshold
    
            exp0 -= expthreshold
            exp1 -= expthreshold
            
            #tti = 2*np.tan(af/2.) / exp0
            #tti = np.max([tti, 0])
            
            
            exp0 = np.max([exp0, 0])
            exp1 = np.max([exp1, 0])
            
            c = -1*exp0 / 3500.
            
            dda = -1*np.max([(exp0-exp1)/dt, 0])
            #dda = -1*(exp0-exp1)/dt
            c = dda / gain[0]
            
            if c > 0:
                print c, dda, np.max((exp0-exp1)/dt, 0)
            
            
            
            
            #c = np.min([c,0])
            #max_accel = -.007
            #c = np.max([c,max_accel])
            if constant_vel is False:
                v = np.max([speed[f-1] + c, 0.07])
            else:
                v = speed[f-1]    
            speed[f] = v
        
            
    expansion = floris.diffa(angle)/dt
    return angle[frames], speed[frames], expansion[frames]
    
    
    
    
def get_error_between_sim_and_trajec(gain, trajec):
    
    angle_sim, speed_sim = sim_deceleration(trajec, gain)
    
    angles_to_check_over = np.linspace(trajec.angle_at_deceleration, trajec.angle_subtended_by_post[trajec.frame_of_landing],20)
    errors = np.zeros(len(angles_to_check_over))
    for i, a in enumerate(angles_to_check_over):
        s_real = np.interp(a, trajec.angle_subtended_by_post[trajec.frames], trajec.speed[trajec.frames])
        s_sim = np.interp(a, angle_sim, speed_sim)
        errors[i] = np.abs(s_real-s_sim)
        
    return np.sum(errors) / float(len(angles_to_check_over))        
            
            
def get_error_between_sim_and_dataset(gain, *args):
    dataset = args[0]
    sum_error = 0
    
    keys = None
    if keys is None:
        classified_keys = fa.get_classified_keys(dataset)
        keys = classified_keys['straight']
    #keys = [keys[2], keys[0], keys[4], keys[6]]
    #nkeys = 5
    #keys = [keys[i] for i in range(20, nkeys+20)]
    n = 0
    for key in keys:
        trajec = dataset.trajecs[key]
        print trajec.obj_id
        e = get_error_between_sim_and_trajec(gain, trajec)
        if not np.isnan(e):
            sum_error += e
            n += 1             
    
    mean_error = sum_error / float(n)
    
    return mean_error
            
            
def optimize_gain_for_dataset(dataset):
    
    gain0 = 150000.
    results = scipy.optimize.fmin(get_error_between_sim_and_dataset, gain0, args=(dataset,))
    
    return results
   
            
            
def test_neural_threshold(save_plot=True, tti=0.12, tti_thresh=0.7, movie_dataset=None):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    distfig = plt.figure()
    distax = distfig.add_subplot(111)
    
    radius = 0.009565
    a = np.linspace(0,2.5,100)
    m = -0.21/radius
    b = 0.159/radius
    expthreshold = (m*np.log(a)+b)*(2*np.tan(a/2.)*np.sin(a/2.))
    
    vel = (expthreshold*radius/(2*np.tan(a/2.)*np.sin(a/2.)))[1:]
    print vel
    
    f = np.where(vel<0.07)[0][0]+1
    print f
    expthreshold_clipped = expthreshold[0:f]
    a_clipped = a[0:f]
    a_clipped_flipped = (a_clipped[::-1]-a_clipped[-1])[::-1]
    
    a_clipped_mirrored = np.hstack( (a_clipped_flipped, a_clipped) )
    expthreshold_clipped_mirrored = np.hstack( (expthreshold_clipped[::-1], expthreshold_clipped) )
    
    ax.plot( a_clipped_mirrored, expthreshold_clipped_mirrored, color='blue' )
    ax.fill_between(a_clipped_mirrored, expthreshold_clipped_mirrored, np.ones_like(a_clipped_mirrored)*30, facecolor='blue', edgecolor=None, alpha=0.1 )
    #ax.plot( a_clipped, expthreshold_clipped, color='blue' )
    #ax.plot( a_clipped[::-1]-a_clipped[-1], expthreshold_clipped, color='blue')
    #ax.fill_between( a_clipped, expthreshold_clipped, np.ones_like(a_clipped)*30, facecolor='blue', edgecolor=None, alpha=0.2 )
    #ax.fill_between( a_clipped[::-1]-a_clipped[-1], expthreshold_clipped, np.ones_like(a_clipped)*30, facecolor='blue', edgecolor=None, alpha=0.2 )
    
    # distance plot
    distax.plot( np.log(a), expthreshold, color='blue')
    
    
    
    
    
    
    
    
    
    
    
    
    # for true time to impact model
    #tti = 0.05
    expthreshold_tti = 2*np.tan(a/2.) / tti - tti_thresh
    vel = (expthreshold_tti*radius/(2*np.tan(a/2.)*np.sin(a/2.)))[1:]
    #f = np.where(vel<0.07)[0][0]+1
    expthreshold_tti_clipped = expthreshold_tti#[0:f]
    a_clipped = a#[0:f]
    
    ax.plot( a_clipped, expthreshold_tti_clipped, color='red' )
    ax.plot( a_clipped[::-1]-a_clipped[-1], expthreshold_tti_clipped, color='red')
    
    deg_ticks = np.array([-90, -45, 0, 45, 90])
    deg_tick_strings = [str(d) for d in deg_ticks]
    rad_ticks = [-np.pi, -np.pi/2., 0, np.pi/2., np.pi]
    
    distax.plot( np.log(a), expthreshold_tti, color='red')
    
    ## plot paramters    
    ax.set_ylim([0,1000*np.pi/180.])
    rad_ticks_y = np.linspace(0,1000*np.pi/180.,5,endpoint=True)
    deg_tick_strings_y = [str(s) for s in np.linspace(0,1000,5,endpoint=True)]
    
    for i, s in enumerate(deg_tick_strings_y):
        deg_tick_strings_y[i] = s.split('.')[0]
    
    ax.set_xlim(rad_ticks[0], rad_ticks[-1])
    ax.set_autoscale_on(False)
    
    fa.adjust_spines(ax, ['left', 'bottom'])
    ax.set_xlabel('position on retina, deg')
    ax.set_ylabel('expansion threshold, deg/s')
    
    ax.set_xticks(rad_ticks)
    ax.set_xticklabels(deg_tick_strings)
    ax.set_yticks(rad_ticks_y)
    ax.set_yticklabels(deg_tick_strings_y)
    
    
    if save_plot:
        fig.savefig('neural_threshold.pdf', format='pdf')
        
    if movie_dataset is not None:
        angle_at_leg_extension, bins, data_filtered, xvals = fa.leg_extension_angle_histogram(movie_dataset, plot=False)
        #ax2.plot(xvals, data_filtered, color='green', alpha=0.3)
        data_filtered /= np.max(data_filtered)
        data_filtered *= 7
        distax.fill_between(xvals, data_filtered, np.zeros_like(xvals), color='green', linewidth=0, alpha=0.2)
    
    fa.fix_angle_log_spine(distax, histograms=False, set_y=False)
    ylim_max = 1000
    distax.set_ylim(0,ylim_max/180.*np.pi)
    rad_ticks_y = np.linspace(0,ylim_max*np.pi/180.,5,endpoint=True)
    deg_tick_strings_y = [str(s) for s in np.linspace(0,ylim_max,5,endpoint=True)]
    for i, s in enumerate(deg_tick_strings_y):
        deg_tick_strings_y[i] = s.split('.')[0]
    distax.set_yticks(rad_ticks_y)
    distax.set_yticklabels(deg_tick_strings_y)
    distax.set_ylabel('expansion threshold, deg/s')
    
    if save_plot:
        distfig.savefig('neural_threshold_distance.pdf', format='pdf')
        
    return a, expthreshold, expthreshold_tti
    
    
def neural_threshold_tti_vs_rsdet_models(save_plot=True, tti=0.12, tti_thresh=0.7, movie_dataset=None):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    distfig = plt.figure()
    distax = distfig.add_subplot(111)
    
    radius = 0.009565
    a = np.linspace(0,2.5,100)
    m = -0.21/radius
    b = 0.159/radius
    expthreshold = (m*np.log(a)+b)*(2*np.tan(a/2.)*np.sin(a/2.))
    
    vel = (expthreshold*radius/(2*np.tan(a/2.)*np.sin(a/2.)))[1:]
    print vel
    
    f = np.where(vel<0.07)[0][0]+1
    print f
    expthreshold_clipped = expthreshold[0:f]
    a_clipped = a[0:f]
    #a_clipped_flipped = (a_clipped[::-1]-a_clipped[-1])[::-1]
    
    #a_clipped_mirrored = np.hstack( (a_clipped_flipped, a_clipped) )
    #expthreshold_clipped_mirrored = np.hstack( (expthreshold_clipped[::-1], expthreshold_clipped) )
    
    ax.plot( a_clipped, expthreshold_clipped, color='purple' )
    ax.fill_between(a_clipped, expthreshold_clipped, np.ones_like(a_clipped)*30, facecolor='purple', edgecolor=None, alpha=0.1 )
    #ax.plot( a_clipped, expthreshold_clipped, color='blue' )
    #ax.plot( a_clipped[::-1]-a_clipped[-1], expthreshold_clipped, color='blue')
    #ax.fill_between( a_clipped, expthreshold_clipped, np.ones_like(a_clipped)*30, facecolor='blue', edgecolor=None, alpha=0.2 )
    #ax.fill_between( a_clipped[::-1]-a_clipped[-1], expthreshold_clipped, np.ones_like(a_clipped)*30, facecolor='blue', edgecolor=None, alpha=0.2 )
    
    # distance plot
    distax.plot( np.log(a), expthreshold, color='purple')
    
    # for true time to impact model
    #tti = 0.05
    expthreshold_tti = 2*np.tan(a/2.) / tti - tti_thresh
    vel = (expthreshold_tti*radius/(2*np.tan(a/2.)*np.sin(a/2.)))[1:]
    #f = np.where(vel<0.07)[0][0]+1
    expthreshold_tti_clipped = expthreshold_tti#[0:f]
    a_clipped = a#[0:f]
    
    ax.plot( a_clipped, expthreshold_tti_clipped, color='black' )
    #ax.plot( a_clipped[::-1]-a_clipped[-1], expthreshold_tti_clipped, color='red')
    
    #deg_ticks = np.array([-90, -45, 0, 45, 90])
    #deg_tick_strings = [str(d) for d in deg_ticks]
    #rad_ticks = [-np.pi, -np.pi/2., 0, np.pi/2., np.pi]
    
    distax.plot( np.log(a), expthreshold_tti, color='black')
    
    
    
    
    ######## plot a sample constant velocity trajectory ##############
    
    vel = 0.4
    fps = 200.0
    x = np.arange(.2, 0, -vel/fps)
    d = x+radius
    a = 2*np.arcsin(radius / (d))
    #exp = 2*vel*np.tan(a/2.)*np.sin(a/2.)
    exp = floris.diffa(a)*fps
    exp = 2/np.sqrt(1-(radius/(d))**2) * (radius/(d)**2) * vel
    distax.plot( np.log(a), exp, color='gray')
    
    
    #return x, a
    
    
    
    ## plot paramters    
    ax.set_ylim([0,1000*np.pi/180.])
    rad_ticks_y = np.linspace(0,1000*np.pi/180.,5,endpoint=True)
    deg_tick_strings_y = [str(s) for s in np.linspace(0,1000,5,endpoint=True)]
    
    for i, s in enumerate(deg_tick_strings_y):
        deg_tick_strings_y[i] = s.split('.')[0]
    
    ax.set_autoscale_on(False)
    
    fa.fix_angle_log_spine(ax, set_y=False, histograms=False)
    
    #fa.adjust_spines(ax, ['left', 'bottom'])
    #ax.set_xlabel('Retinal size')
    ax.set_ylabel('expansion threshold, deg/s')
    
    ax.set_yticks(rad_ticks_y)
    ax.set_yticklabels(deg_tick_strings_y)
    
    
    if save_plot:
        fig.savefig('neural_threshold.pdf', format='pdf')
        
    if movie_dataset is not None:
        angle_at_leg_extension, bins, data_filtered, xvals = fa.leg_extension_angle_histogram(movie_dataset, plot=False)
        #ax2.plot(xvals, data_filtered, color='green', alpha=0.3)
        data_filtered /= np.max(data_filtered)
        data_filtered *= 7
        distax.fill_between(xvals, data_filtered, np.zeros_like(xvals), color='green', linewidth=0, alpha=0.2)
    
    fa.fix_angle_log_spine(distax, histograms=False, set_y=False)
    ylim_max = 1000
    distax.set_ylim(0,ylim_max/180.*np.pi)
    rad_ticks_y = np.linspace(0,ylim_max*np.pi/180.,5,endpoint=True)
    deg_tick_strings_y = [str(s) for s in np.linspace(0,ylim_max,5,endpoint=True)]
    for i, s in enumerate(deg_tick_strings_y):
        deg_tick_strings_y[i] = s.split('.')[0]
    distax.set_yticks(rad_ticks_y)
    distax.set_yticklabels(deg_tick_strings_y)
    distax.set_ylabel('expansion threshold, deg/s')
    
    if save_plot:
        distfig.savefig('neural_threshold_distance.pdf', format='pdf')
        
    return a, expthreshold, expthreshold_tti
    
def plot_expansion_for_sim_trajec(dataset, gain=[170000], keys=None):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    
    if keys is None:
        classified_keys = fa.get_classified_keys(dataset)
        keys = classified_keys['straight']
        #keys = dataset.trajecs.keys()
        
    keys = keys[0:10]
        
    for key in keys:
        trajec = dataset.trajecs[key]
        ftp = np.arange(trajec.frame_at_deceleration-10, trajec.frame_of_landing-1).tolist()
        # simulate
        angle, speed, expansion = sim_deceleration(trajec, gain)
        #ax.plot( np.log(angle), expansion, color='red', linewidth=0.5, alpha=0.5)
        ax.plot( np.log(trajec.angle_subtended_by_post[ftp]), trajec.expansion[ftp], color='black', linewidth=0.5, alpha=0.1)
        ax.plot( np.log(trajec.angle_at_deceleration), trajec.expansion_at_deceleration, '.', color='blue') 
        angle, speed, expansion = sim_deceleration(trajec, gain, constant_vel=True)
        #ax.plot( np.log(angle), expansion, color='green', linewidth=1, alpha=0.5)
        
    angle, expthreshold, expthreshold_tti = test_neural_threshold(save_plot=False, tti=0.05, tti_thresh=6)
    
    ax.plot(np.log(angle), expthreshold, color='blue')
    ax.plot(np.log(angle), expthreshold_tti, color='green')
    
    angle, expthreshold, expthreshold_tti = test_neural_threshold(save_plot=False, tti=0.12, tti_thresh=0.7)
    ax.plot(np.log(angle), expthreshold_tti, color='red')
        
    fa.fix_angle_log_spine(ax, histograms=False, set_y=False)
    exp_limit = 1500
    ax.set_ylim(0,exp_limit/180.*np.pi)
    rad_ticks_y = np.linspace(0,exp_limit*np.pi/180.,5,endpoint=True)
    deg_tick_strings_y = [str(s) for s in np.linspace(0,exp_limit,5,endpoint=True)]
    for i, s in enumerate(deg_tick_strings_y):
        deg_tick_strings_y[i] = s.split('.')[0]
    ax.set_yticks(rad_ticks_y)
    ax.set_yticklabels(deg_tick_strings_y)
    ax.set_ylabel('expansion threshold, deg/s')
    fig.savefig('deceleration_real.pdf', format='pdf')




def plot_tti_wrt_retinal_size(dataset, keys=None, tti_thresh=0.05, plot=False):
    
    if plot:
        fig = plt.figure()
        ax = fig.add_subplot(111)
    
    
    if keys is None:
        classified_keys = fa.get_classified_keys(dataset)
        keys = classified_keys['straight']
        keys = dataset.trajecs.keys()
        
    #keys = keys[0:10]
    
    tti_below_thresh = []
        
    for key in keys:
        trajec = dataset.trajecs[key]
        ftp = np.arange(trajec.frame_at_deceleration-10, trajec.frame_of_landing-1).tolist()
        # simulate
        if plot:
            ax.plot( np.log(trajec.angle_subtended_by_post[ftp]), trajec.time_to_impact[ftp], '.', color='black', alpha=0.3)
            
        #tti_check = (np.sin(angle_subtended/2.)-1) / (np.sin(angle_subtended/2.) + 0.5*1/np.tan(angle_subtended/2.)*expansion)
            
        try:
            tti_below_thresh.append( trajec.angle_subtended_by_post[ np.where( (trajec.time_to_impact < tti_thresh)*(trajec.angle_to_post < 45*np.pi/180.) )[0][0] ])
        except:
            pass
            
    print 
    print 'time to impact below threshold of ', str(tti_thresh), ' sec: '
    print 'mean: ', np.mean(tti_below_thresh)*180/np.pi
    print 'std dev: ', np.std(tti_below_thresh)*180/np.pi
    print 'n: ', len(tti_below_thresh)
    print 'percent: ', len(tti_below_thresh) / float(len(keys))
        
    if plot:
        fa.fix_angle_log_spine(ax, histograms=False, set_y=False)
        ax.set_ylim(0,0.4)
        ticks_y = np.linspace(0,0.4,5,endpoint=True)
        tick_strings_y = [str(s) for s in np.linspace(0,0.4,5,endpoint=True)]
        for i, s in enumerate(tick_strings_y):
            tick_strings_y[i] = s
        ax.set_yticks(ticks_y)
        ax.set_yticklabels(tick_strings_y)
        ax.set_ylabel('time to impact, sec')
        fig.savefig('time_to_impact_vs_angle.pdf', format='pdf')
        
    
    
    
def plot_leg_ext_exp(movie_dataset, delay=0.05):
    
    keys = movie_dataset.movies.keys()
    landing_keys = movie_dataset.get_movie_keys(behavior='landing')
    
    angle = []
    expansion = []
    dist = []
    speed = []
    n = 0
    for key in keys:
        movie = movie_dataset.movies[key]
        if movie.legextensionrange is not None:
            try:
                s = movie.scaled
            except:
                s = None
                print key, ': excepted'
            
            legextensionframe = movie.legextensionrange[0] - movie.firstframe_ofinterest
            if movie.landingframe is not None and movie.trajec is not None and s is not None and 'crash' not in movie.subbehavior and legextensionframe < movie.landingframe_relative:
                
                frame = legextensionframe - int(delay*movie.framerate)
                
                a = movie.scaled.angle_subtended_by_post[frame]
                '''
                while a < -np.pi:
                    a += np.pi
                while a > np.pi:
                    a -= np.pi
                '''
                angle.append(a)
                
                
                e = floris.diffa(movie.scaled.angle_subtended_by_post)
                e *= movie.framerate
                expansion.append(e[frame])
                
                dist.append(movie.scaled.dist_head_to_post[frame])
                speed.append(movie.scaled.speed[frame])
                n += 1
            else:
                print key
            
            
    angle = np.array(angle)
    expansion = np.array(expansion)
    dist = np.array(dist)
    speed = np.array(speed)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    ax.plot(angle*180/np.pi, expansion*180/np.pi, '.')
    
    
    filename = 'leg_extension_angle_expansion.pdf'
    fig.savefig(filename, format='pdf')
    
    print 'mean angle subtended: ', np.mean(angle)*180/np.pi, '+/-', np.std(angle)*180/np.pi, ' deg'
    print 'distance head to post: ', np.mean(dist), '+/-', np.std(dist)
    print 'time to impact: ', np.mean(dist/speed), '+/-', np.std(dist/speed)
    print 'n = ', n
    
    
    
    
def decleration_color_coded_for_post(dataset, plot=True):

    fa.calc_stats_at_deceleration(dataset, keys=None, time_offset=0, return_vals=False)
    angleok = np.where( (dataset.angle_at_deceleration < 2)*(dataset.angle_at_deceleration > .01) )[0].tolist()
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    black_post = []
    checkered_post = []
    for i, p in enumerate(dataset.post_type_at_deceleration):
        if i in angleok:
            if 'black' in p:
                black_post.append(i)
            if 'checkered' in p:
                checkered_post.append(i)
            
    print '*'*40
    print len(black_post), len(checkered_post)
    
    bins, hists, hist_std, curves = floris.histogram(ax, [np.log(dataset.angle_at_deceleration[black_post]), np.log(dataset.angle_at_deceleration[checkered_post])], bins=16, colors=['black', 'teal'], bin_width_ratio=0.9, edgecolor='none', bar_alpha=0.2, curve_fill_alpha=0, curve_line_alpha=0.8, return_vals=True, normed_occurences='total', bootstrap_std=True)
    
    ax.plot(np.log(dataset.angle_at_deceleration[black_post]), dataset.speed_at_deceleration[black_post], '.', color='black', alpha=1)
    ax.plot(np.log(dataset.angle_at_deceleration[checkered_post]), dataset.speed_at_deceleration[checkered_post], '.', color='teal', alpha=1)
    
    #ax.set_ylim(0,1.2)
    fa.fix_angle_log_spine(ax, histograms=True)
    
    fig.savefig('deceleration_color_coded.pdf', format='pdf')









    
