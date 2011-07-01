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
        
    #keys = [keys[2], keys[0], keys[4], keys[6]]
        
    for key in keys:
        trajec = dataset.trajecs[key]
        ftp = np.arange(trajec.frames[0]-25, trajec.frames[-1]).tolist()
        ax.plot( np.log(trajec.angle_subtended_by_post[ftp]), trajec.speed[ftp], color='black', linewidth=0.5, alpha=0.1)
        ax3.plot( np.log(trajec.angle_subtended_by_post[ftp]), trajec.speed[ftp], color='black', linewidth=0.5, alpha=0.1)
        # simulate
        angle, speed = sim_deceleration(trajec, gain)
        if speed[0] > 0.1:        
            ax2.plot( np.log(angle), speed, color='red', linewidth=0.5, alpha=0.3)
            ax3.plot( np.log(angle), speed, color='red', linewidth=0.5, alpha=0.3)
            f = np.where( np.abs(speed-speed[0])>0 )[0][0]
            
            #ax2.plot( np.log(angle[f]), speed[f], '.', color='blue')
            
    x, y, yminus, yplus = fa.get_angle_vs_speed_curve(dataset, plot=False)
    ax2.plot( x, y, color='blue', alpha=0.3)
    ax2.fill_between(x, yplus, yminus, color='blue', linewidth=0, alpha=0.2)

        
    fa.fix_angle_log_spine(ax, histograms=False)
    fig.savefig('deceleration_real.pdf', format='pdf')
    
    fa.fix_angle_log_spine(ax2, histograms=False)
    fig2.savefig('deceleration_sim.pdf', format='pdf')
    
    fa.fix_angle_log_spine(ax3, histograms=False)
    fig3.savefig('deceleration_comparison.pdf', format='pdf')
    
    
def sim_deceleration(trajec, gain):  
    # plot constant velocity comparison
    fps = 2000.
    dt = 1/fps
    r = 0.009565
    radius = r
    
    
    nf = 5000
    positions = np.zeros([nf])
    speed = np.zeros([nf])
    distance = np.zeros([nf])
    angle = np.zeros([nf])
    
    for i in range(2):
        positions[i] = -0.2#-1*(trajec.dist_to_stim_r[trajec.frame_at_deceleration]+radius)
        speed[i] = trajec.speed[trajec.frame_at_deceleration]
        #angle[i] = trajec.angle_subtended_by_post[trajec.frame_at_deceleration]
        distance[i] = np.linalg.norm(positions[i]) - radius
        angle[i] = 2*np.arcsin( radius / (distance[i]+radius) )
    frames = [0,1]
    
    
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
            af = np.min([angle[f], 3])
            af1 = np.min([angle[f-1], 3])
            af2 = np.min([angle[f-2], 3])
            exp0 = (a[f]-a[f-1])/dt #/ (-2.*np.tan(a[f]/2.))
            exp1 = (a[f-1]-a[f-2])/dt #/ (-2.*np.tan(a[f-1]/2.))
            
            
            m = -0.21/radius
            b = 0.159/radius
            expthreshold = (m*np.log(af)+b)*(2*np.tan(af/2.)*np.sin(af/2.))
            
            exp0 -= expthreshold
            exp1 -= expthreshold
            
            exp0 = np.max([exp0, 0])
            exp1 = np.max([exp1, 0])
            
            
            dda = -1*(exp0-exp1)/dt
            
            dda = np.min([dda, 0])
            
            c = dda / gain[0]
            
            
            c = np.min([c,0])
            v = np.max([speed[f-1] + c, 0])
            speed[f] = v
            
    return angle[frames], speed[frames]
    
    
    
    
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
            
            
def test_neural_threshold():
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
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
    
    ax.plot( a_clipped, expthreshold_clipped, color='blue' )
    ax.plot( a_clipped[::-1]-a_clipped[-1], expthreshold_clipped, color='blue')
    
    # for true time to impact model
    tti = 0.15
    expthreshold = 2*np.tan(a/2.) / tti
    vel = (expthreshold*radius/(2*np.tan(a/2.)*np.sin(a/2.)))[1:]
    f = np.where(vel<0.07)[0][0]+1
    expthreshold_clipped = expthreshold[0:f]
    a_clipped = a[0:f]
    
    ax.plot( a_clipped, expthreshold_clipped, color='red' )
    ax.plot( a_clipped[::-1]-a_clipped[-1], expthreshold_clipped, color='red')
    
    deg_ticks = np.array([-90, -45, 0, 45, 90])
    deg_tick_strings = [str(d) for d in deg_ticks]
    rad_ticks = [-np.pi, -np.pi/2., 0, np.pi/2., np.pi]
    
    
    ## plot paramters    
    ax.set_ylim([0,20])
    ax.set_xlim(rad_ticks[0], rad_ticks[-1])
    ax.set_autoscale_on(False)
    
    fa.adjust_spines(ax, ['left', 'bottom'])
    ax.set_xlabel('position on retina, deg')
    ax.set_ylabel('expansion threshold, rad/s')
    
    ax.set_xticks(rad_ticks)
    ax.set_xticklabels(deg_tick_strings)
    
    
    
    
    fig.savefig('neural_threshold.pdf', format='pdf')
            



            
