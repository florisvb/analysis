import sys
#sys.path.insert(0, '/usr/local/lib/python2.6/dist-packages')
sys.path.append('/home/floris/src/pymovie2')
sys.path.append('/home/floris/src/floris_functions')
import floris_plot_lib as fpl
import matplotlib.cm as cm
import colorgrid

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

import copy

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib

import scipy.signal as signal
import scipy.optimize

import colorline
import flydra_analysis as fa
import sa1_analysis as sa1
import numpy as np
import floris

import trajectory_plots as tp
import numpyimgproc as nim


def saccade_figure(dataset_all, dataset_flyby):
    trajec = dataset_flyby.trajecs['10_80113']
    sac.calc_saccades2(trajec, plot=True)
    plot_single_trajec_stuff(trajec)
    
    histogram_angular_velocities(dataset_all)
    plot_turn_angle_vs_threshold(dataset_all, keys=None, plot=True)
    plot_hists_change_in_heading(dataset_all)
    plot_histogram_ratio_saccade_to_straight(dataset)

    
    plot_hists_change_in_heading_straight_vs_saccade(dataset_all)
    
    plot_saccade_speed_vs_amplitude(dataset_all, thresh=300)
    plot_turn_angle_vs_threshold(dataset_all, plot=True)
    plot_saccade_snips(dataset)

def plot_single_trajec_stuff(trajec):
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(trajec.positions[:,0], trajec.positions[:,1], '-', color='black', alpha=1, linewidth=1, zorder=1)
    post = patches.Circle( (0, 0), radius=0.009565, facecolor='black', edgecolor='none', alpha=1)
    ax.add_artist(post)
    tp.prep_xy_spagetti_for_saving(ax)   
    fig.savefig('single_trajectory.pdf', format='pdf')
    
    fig = plt.figure()
    ax = fig.add_axes([0.25,0.2,0.7,0.2])
    ax.plot(trajec.fly_time, floris.diffa(trajec.heading)*100.*180/np.pi, '-', color='blue', alpha=1, linewidth=1, zorder=1)
    ax.plot(trajec.fly_time, trajec.heading_smooth_diff*100.*180/np.pi, '-', color='black', alpha=1, linewidth=1, zorder=1)
    
    # fill in saccades:
    for sac in trajec.sac_ranges:
        angle_vel = trajec.heading_smooth_diff[sac]*100.*180/np.pi
        if np.mean(angle_vel) > 0:
            ax.fill_between(trajec.fly_time[sac], 300*np.ones_like(angle_vel), angle_vel, edgecolor='none', facecolor='red', alpha=0.3)
        else:
            ax.fill_between(trajec.fly_time[sac], angle_vel, -300*np.ones_like(angle_vel), edgecolor='none', facecolor='red', alpha=0.3)
        
        ax.plot(trajec.fly_time[sac], trajec.heading_smooth_diff[sac]*100.*180/np.pi, '-', color='red', alpha=1, linewidth=1, zorder=2)
        
        ax.hlines([-300,300], 0,3, color='red', linewidth=0.5, alpha=1)
        
    yticks = [-1000, 0, 1000]
    xticks = [0, 1, 2, 3]
    ax.set_ylim(-1200, 1200)
    fpl.adjust_spines(ax, ['left', 'bottom'], yticks=yticks, xticks=xticks, smart_bounds=True)
    ax.set_xlabel('Time, s')
    ax.set_ylabel('Angular velocity, deg/s')
    fig.savefig('single_trajectory_angular_velocity.pdf', format='pdf')
    

def plot_histogram_ratio_saccade_to_straight(dataset):
    
    if type(dataset) is not list:
        
        ratios_list = []
        threshes = [200,300,400]
        
        for thresh in threshes:
            
            ratios = []
            for k, trajec in dataset.trajecs.items():
                calc_saccades2(trajec, threshold_lo=thresh)
                frames_sac = np.sum([len(sac) for sac in trajec.sac_ranges])
                all_frames = len(trajec.speed)
                ratios.append(float(frames_sac) / float(all_frames))
            ratios = np.array(ratios)
            ratios_list.append(ratios)
        

        fig = plt.figure()
        ax = fig.add_axes([.2,.2,.7,.4])
        fpl.histogram(ax, ratios_list, bins=20, bin_width_ratio=0.8, colors=['blue', 'red', 'green'], edgecolor='none', bar_alpha=1, curve_fill_alpha=0.3, curve_line_alpha=0, curve_butter_filter=[3,0.3], return_vals=False, show_smoothed=True, normed=False, normed_occurences=False, bootstrap_std=False, exponential_histogram=False, smoothing_range=(0.02, 0.6))
        ax.set_autoscale_on(False)
        ax.set_xlim([0,1])
        ax.set_ylim([0,300])
        
        fpl.adjust_spines(ax, ['left', 'bottom'], yticks = [0, 100, 200, 300])
        ax.set_ylabel('Trajectories')
        ax.set_xlabel('Ratio saccades to straight flight')
        fig.savefig('ratio_saccade_to_straight_hist.pdf', format='pdf')
    
    if type(dataset) is list:
    
        ratios_list = []
        
        gaussian_stats = []
        
        for data in dataset:
            
            ratios = []
            for k, trajec in data.trajecs.items():
                frames_sac = np.sum([len(sac) for sac in trajec.sac_ranges])
                all_frames = len(trajec.speed)
                ratios.append(float(frames_sac) / float(all_frames))
            ratios = np.array(ratios)
            ratios_list.append(ratios)
            
            # some stats:
            indices = np.where(ratios>0)[0].tolist()
            gstats = [np.mean(ratios[indices]), np.std(ratios[indices])]
            gaussian_stats.append(gstats)

        fig = plt.figure()
        ax = fig.add_subplot(111)
        bins, data_hist_list, data_curve_list = floris.histogram(ax, ratios_list, bins=40, bin_width_ratio=0.8, colors=['green', 'gray'], edgecolor='none', bar_alpha=1, curve_fill_alpha=0.3, curve_line_alpha=0, curve_butter_filter=[3,0.3], return_vals=True, show_smoothed=False, normed=False, normed_occurences=False, bootstrap_std=False, exponential_histogram=False)
        
        ax.set_autoscale_on(False)
        ax.set_xlim([0,1])
        ax.set_ylim([0,100])
        fa.adjust_spines(ax, ['left', 'bottom'])
        ax.set_ylabel('Trajectories')
        ax.set_xlabel('Ratio saccades to straight flight')
        
        #plot_stats
        if 0:
            xdata = np.linspace(0, 1, 100)
            for i, gstats in enumerate(gaussian_stats):
                print gstats
                occurences = np.max(data_hist_list[i])
                gaus = np.exp(-1*(xdata-gstats[0])**2 / (2*gstats[1]**2))*occurences
                ax.plot(xdata, gaus)
            
        fig.savefig('ratio_saccade_to_straight_hist_multi.pdf', format='pdf')
            
        
        
        
        
        
        
        
        
        
        
        
        

def saccade_vs_speed_hist(dataset):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    speeds = []
    for k, trajec in dataset.trajecs.items():
        foi = trajec.frame_of_landing
        for sac_range in trajec.sac_ranges:
            if sac_range[-1] < foi:
                speeds.append(np.mean(trajec.speed[sac_range]))
    speeds = np.array(speeds)
            
    floris.histogram(ax, [speeds], bins=20)
    
    fig.savefig('saccade_vs_speed_hist.pdf', format='pdf')
    
def saccade_amp_vs_speed(dataset):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    speeds = []
    amps = []
    for k, trajec in dataset.trajecs.items():
        foi = trajec.frame_of_landing
        for sac_range in trajec.sac_ranges:
            if sac_range[-1] < foi:
                speeds.append(np.mean(trajec.speed[sac_range]))
                amps.append(get_angle_of_saccade(trajec, sac_range)*180./np.pi)
    speeds = np.array(speeds)
    amps = np.abs(np.array(amps))
    
    ax.plot(speeds, amps, '.', color='black')
    
    fig.savefig('saccade_amp_vs_speed.pdf', format='pdf')
    

def plot_heading_to_post(trajec):
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    frames = np.arange(0, trajec.frame_of_landing).tolist()
    
    ax.plot(trajec.fly_time[frames], trajec.heading_smooth[frames], color='black')
    
    for sac_range in trajec.sac_ranges:
        if sac_range[-1] < trajec.frame_of_landing:
            ax.plot(trajec.fly_time[sac_range], trajec.heading_smooth[sac_range], color='red')
        
    fig.savefig('heading_over_time.pdf', format='pdf')
    

def get_trajec_simplified(threshold, trajec, rand_snip=False, plot=False):


    calc_saccades2(trajec, threshold_lo=threshold)

    if rand_snip is False:
        straight_snips, saccade_snips = snip_saccades(trajec, plot=False, return_frames=True, return_pos=False)
    else:
        straight_snips, saccade_snips = snip_saccades(trajec, plot=False, return_frames=True, return_pos=False)
        fake_snips_to_gen = [len(snip) for snip in saccade_snips]
        #for lsnip in fake_snips_to_gen:
            
            
        
    if plot:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(trajec.positions[:,0], trajec.positions[:,1], '-', color='black', alpha=1, linewidth=1, zorder=1)
    
    error_in_straight_snips = 0
    for snip in straight_snips:
        if plot:
            ax.plot(trajec.positions[[snip[0],snip[-1]],0], trajec.positions[[snip[0],snip[-1]],1], '-', color='blue', alpha=1, linewidth=1, zorder=1)
        linept1 = trajec.positions[snip[0],0:2]
        linept2 = trajec.positions[snip[-1],0:2]
        for s in snip:
            error_in_straight_snips += floris.dist_point_to_line(trajec.positions[s,0:2], linept1, linept2, sign=False)
            
    dist_travelled = np.linalg.norm(linept2 - linept1)
    error_in_straight_snips_normed = error_in_straight_snips / dist_travelled
    
    if plot:
        for snip in saccade_snips:
            ax.plot(trajec.positions[snip,0], trajec.positions[snip,1], '-', color='red', alpha=1, linewidth=1, zorder=1)
        post = patches.Circle( (0, 0), radius=0.009565, facecolor='black', edgecolor='none', alpha=1)
        ax.add_artist(post)
        tp.prep_xy_spagetti_for_saving(ax)  
        fig.savefig('simplified_trajec.pdf', format='pdf')
        
    return error_in_straight_snips_normed
    
    
def calc_windowed_heading_changes(trajec, window=10, plot=False):

    fps = 100.
    w = 10
    
    f1 = 0
    f2 = f1+w
    
    trajec.heading_window_diff = np.zeros_like(trajec.speed)
    trajec.heading_window_change = np.zeros_like(trajec.speed)
    
    while f2<len(trajec.speed):
        heading_change_in_window = np.abs(np.sum(trajec.heading_smooth_diff[f1:f2]))*180./np.pi
        #heading_change_per_second_in_window = heading_change_in_window / float(w) * fps
        heading_change_per_second_in_window = np.max(np.abs(trajec.heading_smooth_diff[f1:f2]))*180./np.pi
        f = int((f2-f1)/2.)+f1
        trajec.heading_window_diff[f] = heading_change_per_second_in_window
        trajec.heading_window_change[f] = heading_change_in_window
        
        f1+=1
        f2 = f1+w
        
    if plot:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        
        ax.plot(trajec.heading_window_diff, trajec.heading_window_change, '.', color='black', alpha=0.5)
        
        fig.savefig('heading_change_vs_diff.pdf', format='pdf')
        
def plot_windowed_heading_changes(dataset):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for k, trajec in dataset.trajecs.items():
        calc_windowed_heading_changes(trajec, window=10, plot=False)
        ax.plot(trajec.heading_window_diff, trajec.heading_window_change, '.', color='black', alpha=0.1)
    
    fig.savefig('heading_change_vs_diff.pdf', format='pdf')
        

def histogram_turnsize(dataset):
    
    turns = []
    for k, trajec in dataset.trajecs.items():   
        switch_direction = floris.diffa(np.sign(trajec.heading_smooth_diff[0:trajec.frame_of_landing]))
        switch_indices = np.where( np.abs(switch_direction) > 1)[0].tolist()
        
        
        for i, s in enumerate(switch_indices[0:-2]):
            turn_angle = np.sum(trajec.heading_smooth_diff[switch_indices[i]: switch_indices[i+1]])
            turns.append(turn_angle)
    
    turns = np.abs(np.array(turns))*180./np.pi
        
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    bins = np.linspace(0, 180, 400)
    floris.histogram(ax, [turns], bins=bins, bin_width_ratio=0.8, colors=['gray'], edgecolor='none', bar_alpha=1, curve_fill_alpha=0.3, curve_line_alpha=0, curve_butter_filter=[3,0.3], return_vals=False, show_smoothed=False, normed=False, normed_occurences=False, bootstrap_std=False, exponential_histogram=True)
    
    # show a gaussian:
    mean = np.mean(turns)
    std = np.std(turns)
    std = 100
    xdata = np.linspace(-4000, 4000, 1000)
    gaus = np.exp(-1*(xdata-mean)**2 / (2*std**2)) * np.exp(9) #(1/(std*np.sqrt(2*np.pi)))
    lam = 0.0023
    exp = lam*np.exp(-1*np.abs(xdata)*lam) * np.exp(12.5)
    
    gaus_exp = gaus*.25+exp*.75
    
    dx = np.mean(np.diff(xdata))
    
    gaus_exp_integral = np.sum(gaus_exp)*dx
    print 'gaus exp ingegral: ', gaus_exp_integral
    
    gaus_integral = np.sum(gaus)*dx
    print 'gauss integral: ', gaus_integral
    
    exp_integral = np.sum(exp)*dx
    print 'exp integral: ', exp_integral
    
    ax.plot(xdata, np.log(gaus), linewidth=0.5, color='red')
    ax.plot(xdata, np.log(exp), linewidth=0.5, color='red')
    ax.plot(xdata, np.log(gaus)+np.log(exp), linewidth=0.5, color='blue')
    
    fa.adjust_spines(ax, ['left','bottom'])
    ax.set_ylim([1,8])
    ax.set_xlim([0,180])
    ax.set_xlabel('Turn size, deg')
    ax.set_ylabel('Occurences')
    
    fig.savefig('turn_size_histogram.pdf', format='pdf')
    
    
def plot_continuous_snip_correlation(dataset, behavior='straight'):    

    angle_vel = []
    post_pos = []
    for k, trajec in dataset.trajecs.items():
        straight_snips, saccade_snips = snip_saccades(trajec, plot=False, return_frames=True, return_pos=False)
        if behavior == 'straight':
            snips = straight_snips
        elif behavior == 'saccade':
            snips = saccade_snips
        
        for snip in snips:
            post_pos.extend((trajec.angle_to_post[snip]*180/np.pi).tolist())
            angle_vel.extend((trajec.heading_smooth_diff[snip]*100*180/np.pi).tolist())
    angle_vel = np.array(angle_vel)
    post_pos = np.array(post_pos)
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    hist,x,y = np.histogram2d(post_pos, angle_vel, 300, normed=False)
    colorgrid.colorgrid(ax, np.log(hist.T+1), xlim=(x[0],x[-1]), ylim=(y[0],y[-1]))
    
    fitfmin, variance, intercept_confidence_interval, slope_confidence_interval, Rsq = floris.linear_fit_type2(post_pos, angle_vel, full_output=True, weights=None)
    x = np.linspace(-180,180,100)
    y = np.polyval(fitfmin, x)
    ax.plot(x,y,color='white', linewidth=0.5)
    
    ax.set_aspect('auto')
    ax.set_xlim(-180,180)
    
    if behavior == 'straight':
        ax.set_ylim(-500,500)
    elif behavior == 'saccade':
        ax.set_ylim(-3000, 3000)
        
    fa.adjust_spines(ax, ['left', 'bottom'])
    ax.set_xlabel('Angle to post, deg')
    ax.set_ylabel('Angular velocity, deg/s')
    
    print 'R sqaured: ', Rsq
    
    filename = 'continuous_snip_correlation_' + behavior + '.pdf'
    fig.savefig(filename, format='pdf')
    
    
    
    
    

def histogram_angular_velocities(dataset):
    
    angular_vels = []
    for k, trajec in dataset.trajecs.items():   
        angular_vels.extend(trajec.heading_smooth_diff[0:trajec.frame_of_landing])
    angular_vels = np.array(angular_vels)*180./np.pi*100.
        
    fig = plt.figure()
    ax = fig.add_axes([0.25, 0.2, 0.7, 0.7])
    
    bins = np.linspace(-4000, 4000, 1000)
    fpl.histogram(ax, [angular_vels], bins=bins, bin_width_ratio=0.8, colors=['gray'], edgecolor='none', bar_alpha=1, curve_fill_alpha=0.3, curve_line_alpha=0, curve_butter_filter=[3,0.3], return_vals=False, show_smoothed=False, normed=False, normed_occurences=False, bootstrap_std=False, exponential_histogram=True)
    
    xdata = bins
    dx = np.mean(np.diff(xdata))
    
    # make a gaussian:
    mean = 0
    std = 85
    gaus = np.exp(-1*(xdata-mean)**2 / (2*std**2)) * (1/(std*np.sqrt(2*np.pi)))
    
    # make a lognorm
    mean_lognorm = np.log(300)
    std_lognorm = np.log(2.3)
    lognorm = np.exp(-1*( np.log( np.abs(xdata) )-mean_lognorm)**2 / (2*std_lognorm**2)) * (1/(np.abs(xdata)*np.sqrt(2*np.pi*std_lognorm**2)))
    
    # gaussian + lognorm    
    alpha = 0.2
    gaus_lognorm = gaus*(1-alpha) + lognorm*alpha
    
    # scaling to match occurences histogram
    mult = 15
    
    ax.plot(xdata, np.log(gaus*(1-alpha)*np.exp(mult)), linewidth=1, color='black')
    ax.plot(xdata, np.log(lognorm*(alpha)*np.exp(mult)), linewidth=1, color='red')
    ax.plot(xdata, np.log(gaus_lognorm*np.exp(mult)), linewidth=1, color='magenta')
    
    ax.set_ylim([0,10])
    ax.set_xlim([-3000,3000])
    
    actual_vals = 10**np.arange(0, 5)
    yticks = np.log(actual_vals).tolist()
    xticks = [-3000, -1500, 0, 1500, 3000]
    ytick_labels = [str(a) for a in actual_vals]
    
    fpl.adjust_spines(ax, ['left','bottom'], yticks=yticks, xticks=xticks)
    
    ax.set_yticklabels(ytick_labels)
    ax.set_xlabel('Angular velocity, deg/s')
    ax.set_ylabel('Occurences (log scale)')
    
    fig.savefig('angular_velocity_histogram.pdf', format='pdf')
    
    
    
    
    
    
    # probability of angular velocity falling in one distribution over the other:
    fig2 = plt.figure()
    ax2 = fig2.add_subplot(111)
    ax2.plot(xdata, lognorm*alpha / gaus_lognorm, color='black')
    fa.adjust_spines(ax2, ['left','bottom'])
    ax2.vlines(300, 0, 1)
    print 'probability of being saccade > 95%: ', np.interp(300, xdata, lognorm*alpha / gaus_lognorm)
    
    ax2.set_ylim([0,1.2])
    ax2.set_xlim([-1000,1000])
    fig2.savefig('prob_of_being_saccade.pdf', format='pdf')
    
    ax.plot(xdata, np.log(gaus_lognorm*np.exp(mult)), linewidth=0.5, color='blue')
    
    #ax.plot(xdata, np.log(exp* np.exp(9)), linewidth=0.5, color='red')
    alpha = 0.01
    
    
    

def calc_saccades2(trajec, threshold_lo=300, threshold_hi=100000000, plot=False):

        
    #mean = np.mean( np.abs(trajec.heading_smooth_diff) )*180/np.pi
    #std = np.std( np.abs(trajec.heading_smooth_diff) )*180/np.pi
    
    
    fps = 100.
        
    possible_saccade_array = (np.abs(trajec.heading_smooth_diff)*fps*180/np.pi > threshold_lo)*(np.abs(trajec.heading_smooth_diff)*180/np.pi < threshold_hi)
    possible_saccades = nim.find_blobs(possible_saccade_array, [3,100])
    
    if len(possible_saccades) == 1:
        if np.sum(possible_saccades[0]) == 0:
            possible_saccades = []
    
    trajec.all_saccades = None
    trajec.saccades = None
    trajec.sac_ranges = []
    
    
    if len(possible_saccades) > 0:
        for sac in possible_saccades:
            indices = np.where(sac==1)[0].tolist()
            if len(indices) > 0:
                # expand saccade range to make sure we get full turn
                #lo = np.max([indices[0]-5, 0])
                #hi = np.min([indices[-1]+5, len(trajec.speed)-2])
                #new_indices = np.arange(lo, hi).tolist()
                #tmp = np.where( np.abs(trajec.heading_diff_window[new_indices])*180/np.pi > 350 )[0].tolist()
                #indices = np.array(new_indices)[ tmp ].tolist()
                angle_of_saccade = np.abs(get_angle_of_saccade(trajec, indices)*180./np.pi)
                mean_speed = np.mean(trajec.speed[indices])
                if len(indices) > 3: #angle_of_saccade > 1 and mean_speed > 0.005:
                    trajec.sac_ranges.append(indices)
                    s_rel = np.argmax(trajec.heading_smooth_diff[indices])
                    s = indices[s_rel]
                
                #trajec.all_saccades.append(s)
            
    if plot:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        
        ax.plot(trajec.positions[:,0], trajec.positions[:,1], '-', color='black', alpha=1, linewidth=1, zorder=1)
            
        for sac in trajec.sac_ranges:        
            ax.plot(trajec.positions[sac,0], trajec.positions[sac,1], '-', color='red', alpha=1, linewidth=1, zorder=1+1)
            
                        
        post = patches.Circle( (0, 0), radius=0.009565, facecolor='black', edgecolor='none', alpha=1)
        ax.add_artist(post)
        tp.prep_xy_spagetti_for_saving(ax)   
                        
        fig.savefig('saccade_trajectory.pdf', format='pdf')
        
        
    
    
        
        
def plot_saccade_speed_vs_amplitude(dataset, keys=None, thresh=300):
    if keys is not None:
        if type(keys) is not list:
            keys = [keys]
            
    speeds = []
    angles = []
    duration = []
    forward_speeds = []
    
    for k, trajec in dataset.trajecs.items():
        calc_saccades2(trajec, threshold_lo=thresh)
        if keys is not None:
            if k not in keys:
                continue
                
        foi = trajec.frame_of_landing
    
        for sac_range in trajec.sac_ranges:
            if sac_range[-1] < foi:
                angle_prior = trajec.angle_to_post[sac_range[0]]*-1
                if np.abs(angle_prior*180./np.pi) < 180:
                    if trajec.angle_subtended_by_post[sac_range[0]] < 180:
                        angles.append( get_angle_of_saccade(trajec, sac_range)*180./np.pi )
                        speeds.append( np.max( np.abs(trajec.heading_smooth_diff[sac_range]) )*100.*180./np.pi )
                        duration.append( float(len(sac_range))/trajec.fps )
                        forward_speeds.append( trajec.speed[sac_range[0]] )
    
    
    
    # try with histograms
    fig = plt.figure()
    ax = fig.add_subplot(111)
    hist,x,y = np.histogram2d(np.abs(angles), np.abs(speeds), 300, normed=False)
    colorgrid.colorgrid(ax, np.log(hist.T+1), xlim=(x[0],x[-1]), ylim=(y[0],y[-1]))
    fa.adjust_spines(ax, ['left', 'bottom'])
    ax.set_ylim([0,3000])
    ax.set_xlim([0,180])
    ax.set_aspect('auto')
    ax.set_xlabel('Amplitude, deg')
    ax.set_ylabel('Peak angular velocity, deg/s')
    filename = 'saccade_angle_vs_speed_heatmap' + str(thresh) + '.pdf'
    fig.savefig(filename, format='pdf')
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    hist,x,y = np.histogram2d(np.abs(angles), np.abs(duration)*1000., (300,50), normed=False)
    colorgrid.colorgrid(ax, np.log(hist.T+1), xlim=(x[0],x[-1]), ylim=(y[0],y[-1]))
    fa.adjust_spines(ax, ['left', 'bottom'])
    ax.set_ylim([0,300])
    ax.set_xlim([0,180])
    ax.set_aspect('auto')
    ax.set_xlabel('Amplitude, deg')
    ax.set_ylabel('Duration, ms')
    filename = 'saccade_angle_vs_duration_' + str(thresh) + '.pdf'
    fig.savefig(filename, format='pdf')
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    hist,x,y = np.histogram2d(np.abs(angles), np.abs(forward_speeds), (300,100), normed=False)
    colorgrid.colorgrid(ax, np.log(hist.T+1), xlim=(x[0],x[-1]), ylim=(y[0],y[-1]))
    fa.adjust_spines(ax, ['left', 'bottom'])
    ax.set_ylim([0,1])
    ax.set_xlim([0,180])
    ax.set_aspect('auto')
    ax.set_xlabel('Amplitude, deg')
    ax.set_ylabel('Forward speed, ms')
    filename = 'saccade_angle_vs_forwardspeed_' + str(thresh) + '.pdf'
    fig.savefig(filename, format='pdf')

    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    floris.histogram(ax, [angles], bins=100, bin_width_ratio=0.8, colors=['gray'], edgecolor='none', bar_alpha=1, curve_fill_alpha=0.3, curve_line_alpha=0, curve_butter_filter=[3,0.3], return_vals=False, show_smoothed=False, normed=False, normed_occurences=False, bootstrap_std=False)
    fa.adjust_spines(ax, ['left','bottom'])
    ax.set_ylim([0,700])
    ax.set_xlim([0,180])
    ax.set_xlabel('Amplitude, deg')
    ax.set_ylabel('Duration, ms')
    filename = 'saccade_angle_histogram' + str(thresh) + '.pdf'
    fig.savefig(filename, format='pdf')


def snip_saccades(trajec, plot=False, return_pos=True, return_frames=False):
    
    straight_snips = []
    saccade_snips = []    
    
    frame_counter = 0
    
    for sac_range in trajec.sac_ranges:
        if sac_range[-1] < trajec.frame_of_landing:
            frames_up_to_saccade = np.arange(frame_counter, sac_range[0]).tolist()
            if len(frames_up_to_saccade) > 2:
                straight_snips.append( frames_up_to_saccade )
            saccade_snips.append( sac_range )
            
            frame_counter = sac_range[-1]
        
    # catch remaining frames
    if frame_counter < len(trajec.speed):
        frames_after_last_saccade = np.arange(frame_counter, trajec.frame_of_landing).tolist()
        if len(frames_after_last_saccade) > 2:
            straight_snips.append( frames_after_last_saccade ) 
            
    def rotate_snip(snip):
        pos = trajec.positions[snip,0:2] - trajec.positions[snip[0],0:2]
        
        vec = pos[2]-pos[0]
        
        
        #vec = pos[-1]-pos[0]
        
        angle_to_rotate = np.arctan2(vec[1], vec[0])
        rotation_matrix = np.array([[np.cos(angle_to_rotate), -np.sin(angle_to_rotate)],[np.sin(angle_to_rotate), np.cos(angle_to_rotate)]])
        pos_rot = np.dot(pos, rotation_matrix)
        
        if 0:
            # now rotate to match heading:
            heading_change = -1*np.sum(trajec.heading_smooth_diff[snip])
            heading_rotation_matrix = np.array([[np.cos(heading_change), -np.sin(heading_change)],[np.sin(heading_change), np.cos(heading_change)]])
            print heading_change
            pos_rot = np.dot(pos_rot, heading_rotation_matrix)
        
        
        # flip
        #pos_rot[:,1] = np.abs(pos_rot[:,1])

        return pos_rot
        
        
    straight_snips_pos_rotated = [rotate_snip(snip) for snip in straight_snips]
    saccade_snips_pos_rotated = [rotate_snip(snip) for snip in saccade_snips]
        
    if plot:
                        
        ###### deconstructed
        
        fig = plt.figure()
        ax = fig.add_subplot(111)
        
        for snip in straight_snips_pos_rotated:
            ax.plot(snip[:,0], snip[:,1], '-', color='black', alpha=1, linewidth=0.5, zorder=1)
        ax.hlines(0, 0, 0.1, color='blue', linewidth=0.5)
        
        ax.set_aspect('equal')

        fig.savefig('straight_snips_single_trajec.pdf', format='pdf')       
        
    elif return_pos:
        return straight_snips_pos_rotated, saccade_snips_pos_rotated
    elif return_frames:
        return straight_snips, saccade_snips
                        
def plot_saccade_snips(dataset):

    fig_straight = plt.figure()
    ax_straight = fig_straight.add_subplot(111)
    
    fig_sac = plt.figure()
    ax_sac = fig_sac.add_subplot(111)
    
    fig_tracking_correlation = plt.figure()
    

    for k, trajec in dataset.trajecs.items():
        #print k
        straight_snips, saccade_snips = snip_saccades(trajec, plot=False)
        for snip in straight_snips:
            ax_straight.plot(snip[:,0], snip[:,1], '-', color='black', alpha=0.5, linewidth=0.5, zorder=1)
        for snip in saccade_snips:
            ax_sac.plot(snip[:,0], snip[:,1], '-', color='red', alpha=0.5, linewidth=0.5, zorder=1)
            

    ax_straight.hlines(0, 0, 0.4, color='blue', linewidth=0.5, zorder=100)
    ax_sac.hlines(0, 0, 0.06, color='blue', linewidth=0.5, zorder=100)  
    
    ax_straight.set_autoscale_on(False)
    ax_straight.set_xlim([0,0.6])
    ax_straight.set_ylim([-.2,0.2])
    ax_straight.set_aspect('equal')
    
    ax_sac.set_autoscale_on(False)
    ax_sac.set_xlim([-.01,0.08])
    ax_sac.set_ylim([-0.06, .06])
    ax_sac.set_aspect('equal')

    fig_straight.savefig('straight_snips_trajectory.pdf', format='pdf')       
    fig_sac.savefig('saccade_snips_trajectory.pdf', format='pdf')       
    
    

def get_change_in_heading(trajec, return_lists=False):
    #for k, trajec in dataset.trajecs.items():
    #print k
    straight_snips, saccade_snips = snip_saccades(trajec, plot=False, return_frames=True, return_pos=False)

    if return_lists is False:
        change_in_course_straight = np.sum([ np.abs(np.sum( trajec.heading_smooth_diff[snip] )) for snip in straight_snips])
        change_in_course_saccade = np.sum([ np.abs(np.sum( trajec.heading_smooth_diff[snip] )) for snip in saccade_snips])
        change_in_course_all = np.sum( np.abs(trajec.heading_smooth_diff[0:trajec.frame_of_landing] ) )
        return change_in_course_straight, change_in_course_saccade, change_in_course_all
    else:
        change_in_course_straight = [ np.abs(np.sum( trajec.heading_smooth_diff[snip] )) for snip in straight_snips]
        change_in_course_saccade = [ np.abs(np.sum( trajec.heading_smooth_diff[snip] )) for snip in saccade_snips]
        change_in_course_all = np.sum(np.abs(trajec.heading_smooth_diff[0:trajec.frame_of_landing]))
        return change_in_course_straight, change_in_course_saccade, change_in_course_all
    
def get_proportion_saccade_of_trajec(trajec):
    
    straight_snips, saccade_snips = snip_saccades(trajec, plot=False, return_frames=True, return_pos=False)
    
    total_saccade_frames = np.sum( [np.sum(sac_snip) for sac_snip in saccade_snips] )
    total_straight_frames = np.sum( [np.sum(straight_snip) for straight_snip in straight_snips] )
    
    #ratio_saccade_frames = float(total_saccade_frames) / float(total_saccade_frames + total_straight_frames)
    
    ratio_saccade_frames = len(saccade_snips)
    
    return ratio_saccade_frames
    

    
def plot_hists_num_saccades_per_meter(dataset_list, thresh=300):
    
    sac_per_meter_list = []
    for dataset in dataset_list:
        fa.calc_func(dataset, calc_saccades2, thresh)

        sac_per_meter = []
        n_trajecs = 0
        for k, trajec in dataset.trajecs.items():
            dist_travelled = np.sum(trajec.speed*trajec.fps)
            n_saccades = len(trajec.sac_ranges)
            sac_per_meter.append( float(n_saccades) / dist_travelled )
            n_trajecs += 1
    
        print 'n, trajecs: ', n_trajecs
        sac_per_meter = np.array(sac_per_meter)
        sac_per_meter_list.append(sac_per_meter)
    
        
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    nbins = 50
    bins = np.linspace(0, 0.006, nbins, endpoint=True)
    bins, data_hist_list, data_curve_list = floris.histogram(ax, sac_per_meter_list, bins=bins, bin_width_ratio=0.8, colors=['blue', 'red', 'green'], edgecolor='none', bar_alpha=1, curve_fill_alpha=0.3, curve_line_alpha=0, curve_butter_filter=[3,0.3], return_vals=True, show_smoothed=False, normed=False, normed_occurences=True, bootstrap_std=False)
    
    ax.set_autoscale_on(False)
    ax.set_xlim([0,0.01])
    ax.set_ylim([0,1])
    
    fa.adjust_spines(ax, ['left', 'bottom'])
    
    ax.set_xlabel('Number of saccades per meter')
    ax.set_ylabel('Occurences (normalized)')
    
    fig.savefig('num_sac_per_meter_histograms.pdf', format='pdf')
    
    

def plot_hists_change_in_heading_straight_vs_saccade(dataset, thresh=300):

    fig_frames = plt.figure()
    ax_frames = fig_frames.add_subplot(111)

    
    fa.calc_func(dataset, calc_saccades2, thresh)

    straight = []
    saccade = []
    full = []
    n_saccades = []
    
    n_trajecs = 0
    for k, trajec in dataset.trajecs.items():
        change_in_course_straight, change_in_course_saccade, change_in_course_all = get_change_in_heading(trajec, return_lists=True)
        straight.extend(change_in_course_straight)
        saccade.extend(change_in_course_saccade)
        full.append(change_in_course_all)
        n_trajecs += 1
        
    print 'n, straight: ', len(straight)
    print 'n, saccade: ', len(saccade)
    print 'n, trajecs: ', n_trajecs
    
    straight = np.array(straight)*180./np.pi
    saccade = np.array(saccade)*180./np.pi
    full = np.asarray(full)*180./np.pi
    
        
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    nbins = 50
    bins = np.linspace(0, 180.+180./float(nbins), nbins, endpoint=True)
    bins, data_hist_list, data_curve_list = floris.histogram(ax, [straight, saccade], bins=bins, bin_width_ratio=0.8, colors=['black', 'red'], edgecolor='none', bar_alpha=1, curve_fill_alpha=0.3, curve_line_alpha=0, curve_butter_filter=[3,0.3], return_vals=True, show_smoothed=True, normed=False, normed_occurences=True, bootstrap_std=False)
    
    ax.set_autoscale_on(False)
    ax.set_xlim([0,180])
    ax.set_ylim([0,1])
    
    fa.adjust_spines(ax, ['left', 'bottom'])
    
    ax.set_xlabel('Net heading change in flight segment')
    ax.set_ylabel('Occurences (normalized)')
    
    fig.savefig('heading_histograms.pdf', format='pdf')
    

def plot_hists_change_in_heading(dataset):


    saccade_ratio_list = []
    threshes = [200,300,400]
    
    for thresh in threshes:
    
        fa.calc_func(dataset, calc_saccades2, thresh)
    
        straight = []
        saccade = []
        full = []
        n_saccades = []
        
        for k, trajec in dataset.trajecs.items():
            change_in_course_straight, change_in_course_saccade, change_in_course_all = get_change_in_heading(trajec)
            straight.append(change_in_course_straight)
            saccade.append(change_in_course_saccade)
            full.append(change_in_course_all)
            
            straight_snips, saccade_snips = snip_saccades(trajec, plot=False, return_frames=True, return_pos=False)
            n_saccades.append( len(saccade_snips)  )
            
        straight = np.array(straight)
        saccade = np.array(saccade)
        full = np.array(full)
        n_saccades = np.array(n_saccades, dtype=float)
        
        where_zeros = np.where( n_saccades==0 )[0].tolist()
        n_saccades[where_zeros] = 1
        
        heading_ratio = saccade/(saccade+straight)
        saccade_ratio_list.append( heading_ratio )
        
    fig = plt.figure()
    ax = fig.add_axes([0.2, 0.2, 0.7, 0.4])
    
    #bins, data_hist_list, data_curve_list = floris.histogram(ax, [straight, saccade], bins=20, bin_width_ratio=0.6, colors=['black', 'red'], edgecolor='none', bar_alpha=0.7, curve_fill_alpha=0, curve_line_alpha=0, curve_butter_filter=[3,0.3], return_vals=True, show_smoothed=False, normed=False, normed_occurences=False, bootstrap_std=False)
    
    bins, data_hist_list, data_curve_list = fpl.histogram(ax, saccade_ratio_list, bins=25, bin_width_ratio=0.8, colors=['blue', 'red', 'green'], edgecolor='none', bar_alpha=1, curve_fill_alpha=0.3, curve_line_alpha=0, curve_butter_filter=[3,0.3], return_vals=True, show_smoothed=True, normed=False, normed_occurences=False, bootstrap_std=False, smoothing_range=(0.1, 1.0))
    
    print bins
    
    ax.set_autoscale_on(False)
    ax.set_xlim([0,1])
    ax.set_ylim([0,300])
    
    fpl.adjust_spines(ax, ['left', 'bottom'], yticks = [0, 100, 200, 300])
    
    ax.set_xlabel('Perc Heading Changes During Saccade')
    ax.set_ylabel('Occurences')
    
    fig.savefig('heading_histograms.pdf', format='pdf')

def get_angle_of_saccade(trajec, sac_range, method='integral', smoothed=True):
    
    if method != 'integral':
        f0 = sac_range[0]
        f1 = sac_range[-1]

        obj_ori_0 = trajec.velocities[f0] / np.linalg.norm(trajec.velocities[f0])   
        obj_ori_1 = trajec.velocities[f1] / np.linalg.norm(trajec.velocities[f1])  

        obj_ori_0_3vec = np.hstack( ( obj_ori_0, 0) ) 
        obj_ori_1_3vec = np.hstack( (obj_ori_1, 0 ) ) 

        sign_of_angle_of_saccade = np.sign( np.sum(np.cross( obj_ori_0, obj_ori_1 ) ) )

        cosangleofsaccade = np.dot(obj_ori_0, obj_ori_1)
        angleofsaccade = np.arccos(cosangleofsaccade)
         
        signed_angleofsaccade = -1*angleofsaccade*sign_of_angle_of_saccade
        
        return signed_angleofsaccade
    
    else:
        if smoothed is False:
            return np.sum(trajec.heading_smooth_diff[sac_range])*-1
        else:
            return np.sum( floris.diffa(trajec.heading)[sac_range])*-1
            

def plot_turn_angle_vs_threshold(dataset, keys=None, plot=False):

    threshes = np.hstack( ([300], np.linspace(150, 450, 15, endpoint=True)) )
    sac_angles = [[] for i in threshes]  
    
    for k, trajec in dataset.trajecs.items():
        if keys is not None:
            if trajec.key not in keys:
                continue
        calc_saccades2(trajec)
        
        for sac in trajec.sac_ranges:
            if sac[-1] > trajec.frame_of_landing:
                continue
            peak = np.argmax(trajec.heading_smooth_diff[sac])+sac[0]
            
            # for each threshold walk forwards and backwards to find where saccade would start/end:
            for i, thresh in enumerate(threshes):
                # find upper limit
                topindex = peak
                angvel = np.abs(trajec.heading_smooth_diff[topindex]*100.*180/np.pi)
                while angvel > thresh:    
                    topindex += 1
                    if topindex >= len(trajec.speed)-1:
                        topindex = len(trajec.speed)-1
                        angvel = np.abs(trajec.heading_smooth_diff[topindex]*100.*180/np.pi)
                        break
                    if topindex-peak > 20:
                        angvel = np.abs(trajec.heading_smooth_diff[topindex]*100.*180/np.pi)
                        break
                    angvel = np.abs(trajec.heading_smooth_diff[topindex]*100.*180/np.pi)
                    
                # find lower limit
                botindex = peak
                angvel = np.abs(trajec.heading_smooth_diff[botindex]*100.*180/np.pi)
                while angvel > thresh:  
                    botindex -= 1
                    if botindex <= 0:
                        botindex = 0
                        angvel = np.abs(trajec.heading_smooth_diff[botindex]*100.*180/np.pi)
                        break
                    if peak-botindex > 20:
                        angvel = np.abs(trajec.heading_smooth_diff[botindex]*100.*180/np.pi)
                        break
                    angvel = np.abs(trajec.heading_smooth_diff[botindex]*100.*180/np.pi)
                sac_range = np.arange(botindex, topindex).tolist()
                
                if len(sac_range) > 0:
                    sac_angle = np.abs(get_angle_of_saccade(trajec, sac_range))
                else:
                    sac_angle = np.nan
                    
                sac_angles[i].append(sac_angle*180./np.pi)
                
                #print thresh, sac_range, sac_angle
                
    sac_angles = [np.array(sac_angle) for sac_angle in sac_angles]

    master_sac = 0
    diffs = []
    for i, sac_angle in enumerate(sac_angles):
        if i == 0:
            continue  
        tmp = sac_angle-sac_angles[master_sac]
        ind_not_nan = np.where(np.isnan(sac_angle)==0)[0].tolist() # make sure to exclude saccades if they outright dissappeared!
        diffs.append(tmp[ind_not_nan])
        
    means = np.array([np.mean(diff) for diff in diffs])
    stds = np.array([np.std(diff) for diff in diffs])
    
    #means = np.array([floris.geomean(diff) for diff in diffs])
    #stds = np.array([floris.geostd(diff) for diff in diffs])
    
    
    if plot:
        fig = plt.figure()
        ax = fig.add_axes([0.2,0.2,0.7,0.2])
        
        fpl.boxplot(ax, threshes[1:], diffs, colormap=None, boxwidth=10, show_outliers=False, show_whiskers=False, boxlinewidth=2)
        
        # show 10deg
        #ax.hlines([-10,10], 125, 400, linestyle='solid', linewidth=0.5, color='blue')
        
        
        ax.set_xlim(140, 460)
        ax.set_ylim(-15, 15)
        ax.set_autoscale_on(False)
        
        yticks = [-10, 0, 10]
        xticks = [150, 300, 450]
        fpl.adjust_spines(ax, ['left', 'bottom'], smart_bounds=True, yticks=yticks, xticks=xticks)
        
        #ax.set_xticks(xticks)
        #ax.set_yticks(yticks)
        
        
        
        ax.set_xlabel('Threshold, deg/s')
        ax.set_ylabel('Change in turn angle, deg')
        
        #fig_size =  (3.25, 2)
        #fig.set_size_inches(fig_size)
        
        fig.savefig('turnangle_vs_threshold.pdf', format='pdf')
            
    #return sac_angles
            
    
    
def calc_saccade_stats(dataset, keys=None, dist_thresh=0.2):
    if keys is None:
        keys = dataset.trajecs.keys()
        
    speed_at_saccade_array = []
    angle_of_saccade_array = []
    flipped_angle_of_saccade_array = []
    angle_prior_array = []
    flipped_angle_prior_array = []
    angle_subtended_array = []
    attractiveness_array = []
    dist_at_saccade_array = []
    post_types_saccade_array = []
    keys_saccade_array = []
            
    for key in keys:
        trajec = dataset.trajecs[key]
        #sa1.calc_post_dynamics_for_flydra_trajectory(trajec)
        
        if len(trajec.sac_ranges) > 0:
            try:
                frames = trajec.frames_below_post
            except:
                frames = np.arange(0, trajec.frame_nearest_to_post).tolist()
        
            for sac_range in trajec.sac_ranges:
                if sac_range[0] in frames:
                    if sac_range[-1] <= trajec.frame_of_landing and trajec.dist_to_stim_r_normed[sac_range[0]] < 0.2:
                        angle_prior = trajec.angle_to_post[sac_range[0]]*-1
                        signed_angleofsaccade = get_angle_of_saccade(trajec, sac_range)
                        if signed_angleofsaccade is None:
                            continue
                        attractiveness = np.abs(angle_prior - signed_angleofsaccade)
                        
                        if signed_angleofsaccade < 0:
                            flipped_signed_angleofsaccade = -1*signed_angleofsaccade
                            flipped_angle_prior = -1*angle_prior
                        else:
                            flipped_signed_angleofsaccade = signed_angleofsaccade
                            flipped_angle_prior = angle_prior
                        
                        speed_at_saccade_array.append(trajec.speed[sac_range[0]])
                        angle_of_saccade_array.append(signed_angleofsaccade)
                        flipped_angle_of_saccade_array.append(flipped_signed_angleofsaccade)
                        angle_prior_array.append(angle_prior)
                        flipped_angle_prior_array.append(flipped_angle_prior)
                        angle_subtended_array.append( trajec.angle_subtended_by_post[sac_range[0]])
                        attractiveness_array.append(attractiveness)
                        dist_at_saccade_array.append(trajec.dist_to_stim_r_normed[sac_range[0]])
                        post_types_saccade_array.append(trajec.post_type)
                        keys_saccade_array.append(trajec.key)
                            
    dataset.speed_at_saccade_array = np.array(speed_at_saccade_array)
    dataset.angle_of_saccade_array = np.array(angle_of_saccade_array)
    dataset.flipped_angle_of_saccade_array = np.array(flipped_angle_of_saccade_array)
    dataset.angle_prior_array = np.array(angle_prior_array)
    dataset.flipped_angle_prior_array = np.array(flipped_angle_prior_array)
    dataset.angle_subtended_array = np.array(angle_subtended_array)
    dataset.attractiveness_array = np.array(attractiveness_array)
    dataset.dist_at_saccade_array = np.array(dist_at_saccade_array)
    dataset.post_types_saccade_array = post_types_saccade_array
    dataset.keys_saccade_array = keys_saccade_array
    
def get_saccade_angle_mean_and_variance(dataset, bins=None, plot=True, angular_range=[25,120]):
    
    angular_range = np.array(angular_range)*np.pi/180.
    
    colormap_norm = matplotlib.colors.Normalize(0, 0.1, clip=True)
    colormap_norm_dist = matplotlib.colors.Normalize(0, 0.2, clip=True)
    colormap_speed_norm = matplotlib.colors.Normalize(0, 0.4, clip=True)
    colormap_angle_norm = matplotlib.colors.Normalize(0, np.pi/2., clip=True)
    cmap = plt.get_cmap('jet')
    
    in_range = (dataset.angle_subtended_array>angular_range[0])*(dataset.angle_subtended_array<angular_range[1])
    in_range_indices = np.where(in_range)[0].tolist()
    
    
    fit_angle_prior = []
    fit_angle_saccade = []
        
    var = 30*np.pi/180.
        
    if plot:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        colormap_norm_angle_subtended = matplotlib.colors.Normalize(angular_range[0], 90*np.pi/180, clip=True)
        
        for i in range(len(dataset.angle_of_saccade_array)):
            if i in in_range_indices: 
            
                c = cmap(colormap_norm(0.1 - dataset.dist_at_saccade_array[i]))
                alpha = 1 #colormap_norm_dist(0.2 - dataset.dist_at_saccade_array[i])
                r = colormap_speed_norm(dataset.speed_at_saccade_array[i])/100.*4 + 0.03
                
                angle_prior = dataset.angle_prior_array[i] + np.sign(dataset.angle_prior_array[i])*dataset.angle_subtended_array[i]/2.
            
                sac = patches.Circle( (angle_prior, dataset.angle_of_saccade_array[i]), radius=r, facecolor=c, edgecolor='none', alpha=alpha, zorder=100, linewidth=0.3)
                ax.add_artist(sac)
                
                if dataset.angle_of_saccade_array[i] > angle_prior:
                    fit_angle_prior.append(angle_prior)
                    fit_angle_saccade.append(dataset.angle_of_saccade_array[i])
                
                mirror = True
                if mirror:
                    mirrored_signed_angleofsaccade = -1*dataset.angle_of_saccade_array[i]
                    mirrored_angle_prior = -1*angle_prior
                    sac = patches.Circle( (mirrored_angle_prior, mirrored_signed_angleofsaccade), radius=r, facecolor=c, edgecolor='none', alpha=alpha, zorder=100, linewidth=0.3)
                    ax.add_artist(sac)
                    
                    if mirrored_signed_angleofsaccade > mirrored_angle_prior:
                        fit_angle_prior.append(mirrored_angle_prior)
                        fit_angle_saccade.append(mirrored_signed_angleofsaccade)
                    
                
        fit_angle_prior = np.array(fit_angle_prior)
        fit_angle_saccade = np.array(fit_angle_saccade)
                
        fit, variance, intercept_confidence_interval, slope_confidence_interval, Rsq = floris.linear_fit_type( fit_angle_prior, fit_angle_saccade, weights=None, alpha=0.05, full_output=True, fit_type=2)
        
        print 'linear fit: '
        print fit
        print Rsq
        print variance
        
        ytanh = np.tanh(fit_angle_prior-var)*np.pi + np.pi
        err_ytanh = np.mean(np.abs((fit_angle_saccade - ytanh)))
        print 'tanh fit: '
        print err_ytanh
        
        if bins is None:
            bins = np.log(np.linspace(5*np.pi/180., np.pi, 20, endpoint=True))
        bin_centers = np.diff(bins)/2.+bins[0:-1]
                    
        x = np.linspace(-180*np.pi/180., 180*np.pi/180., 100)
        #yp = np.tanh(x-var)/2.*np.pi*m + (1/2.*np.pi*m + var )
        #yp = np.tanh(x-var)/2.*np.pi*2 + (1/2.*np.pi*2)
        #ym = np.tanh(x+var)/2.*np.pi*2 - (1/2.*np.pi*2)
        #ax.fill_between(x,ym,yp, facecolor='red', edgecolor='none', alpha=0.2, zorder=-100)
        
        safety = 90*np.pi/180.
        #yp = x+safety
        #ym = x-safety
        
        yp = np.polyval(fit, x)
        ym = yp - 2*fit[1]
        
        ax.plot(x,yp,color='black')
        ax.plot(x,ym,color='black')
        #ym = np.tanh(x+var)/2.*np.pi*m - (1/2.*np.pi*m + var )
        
        ax.set_ylim(-np.pi, np.pi)
        ax.set_xlim(-np.pi, np.pi)
        ax.set_aspect('equal')
        ax.set_autoscale_on(False)
        
        rad_ticks_x = [-np.pi, -np.pi/2., 0, np.pi/2., np.pi]
        deg_tick_strings_x = ['-180', '-90', '0', '90', '180']
        
        rad_ticks_y = [0, np.pi/2., np.pi]
        deg_tick_strings_y = ['0', '90', '180']
        
        fa.adjust_spines(ax, ['left', 'bottom'])
        ax.set_xticks(rad_ticks_x)
        ax.set_xticklabels(deg_tick_strings_x) 
        ax.set_yticks(rad_ticks_x)
        ax.set_yticklabels(deg_tick_strings_x) 
        
        ax.set_xlabel('retinal position of post, deg')
        ax.set_ylabel('angle of saccade, deg')
        
        fig.subplots_adjust(bottom=0.2, top=0.9, right=0.9, left=0.2)
        fig.savefig('fitted_saccades.pdf', format='pdf')
        
        
    return
    
    
def plot_saccade_angle_for_different_threshes(dataset):

    threshes = [200,300,400]
    
    for thresh in threshes:
        fa.calc_func(dataset, calc_saccades2, thresh)

        filename = 'saccade_angle_' + str(thresh) + '.pdf'
        plot_saccade_angle(dataset, keys=None, angle_subtended_range=[0,180], retinal_position_range=[-160, 160], show_evasive_thresholds=False, show_tanh_evasive_thresholds=False, edge='far', colorscale=True, speedscale=False, regress=None, plot_regress=False, plot=True, plot_color_regress=False, mirror=False, postcolorscale=False, return_left_right=False, recalc_saccade_stats=True, filename=filename)
        
        
    
    
def plot_saccade_angle(dataset, keys=None, angle_subtended_range=[0,180], retinal_position_range=[-160, 160], show_evasive_thresholds=False, show_tanh_evasive_thresholds=False, edge='far', colorscale=False, speedscale=False, regress=None, plot_regress=False, plot=True, plot_color_regress=False, mirror=False, postcolorscale=False, return_left_right=False, recalc_saccade_stats=True, filename=None):
    if keys is None:
        keys = dataset.trajecs.keys()
    
    if recalc_saccade_stats:
        calc_saccade_stats(dataset, keys=keys)
    
    angle_subtended_range = np.array(angle_subtended_range)*np.pi/180.
    retinal_position_range = np.array(retinal_position_range)*np.pi/180.
        
    if plot:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        
        colormap_norm = matplotlib.colors.Normalize(0, 0.1, clip=True)
        colormap_norm_dist = matplotlib.colors.Normalize(0, 0.2, clip=True)
        colormap_speed_norm = matplotlib.colors.Normalize(0, 0.4, clip=True)
        colormap_angle_norm = matplotlib.colors.Normalize(0, np.pi/2., clip=True)
        cmap = plt.get_cmap('jet')

    if edge == 'middle':
        angle_prior = dataset.angle_prior_array
    if edge == 'near':
        angle_prior = dataset.angle_prior_array - np.sign(dataset.angle_prior_array)*dataset.angle_subtended_array/2.
    if edge == 'far':
        angle_prior = dataset.angle_prior_array + np.sign(dataset.angle_prior_array)*dataset.angle_subtended_array/2.

    in_range = (dataset.angle_subtended_array>angle_subtended_range[0])*(dataset.angle_subtended_array<angle_subtended_range[1])*(angle_prior>retinal_position_range[0])*(angle_prior<retinal_position_range[1])
    in_range_indices = np.where(in_range)[0].tolist()
        
        
    if plot:
        for i in range(len(dataset.angle_of_saccade_array)):
            if i in in_range_indices: 
            
                if colorscale is True:
                    c = cmap(colormap_norm(0.1 - dataset.dist_at_saccade_array[i]))
                    alpha = 0.8#colormap_norm_dist(0.2 - dataset.dist_at_saccade_array[i])
                elif postcolorscale is True:
                    if 'checkered' in dataset.post_types_saccade_array[i]:
                        c = 'teal'
                        alpha = 0.8
                    elif 'black' in dataset.post_types_saccade_array[i]:
                        c = 'black'
                        alpha = 0.8
                    else:
                        c = 'black'
                        alpha = 0.8
                else:
                    c = 'black'
                    alpha = 0.8
                
                if speedscale is True:    
                    r = colormap_speed_norm(dataset.speed_at_saccade_array[i])/100.*4 + 0.03
                else:
                    r = .5/100.*4 + 0.03
                
                
                if edge == 'middle':
                    angle_prior = dataset.angle_prior_array[i]# + np.sign(dataset.angle_prior_array[i])*dataset.angle_subtended_array[i]/2.
                if edge == 'near':
                    angle_prior = dataset.angle_prior_array[i] - np.sign(dataset.angle_prior_array[i])*dataset.angle_subtended_array[i]/2.
                if edge == 'far':
                    angle_prior = dataset.angle_prior_array[i] + np.sign(dataset.angle_prior_array[i])*dataset.angle_subtended_array[i]/2.
            
            
                if alpha > 0.1:
                    sac = patches.Circle( (angle_prior, dataset.angle_of_saccade_array[i]), radius=r, facecolor='white', edgecolor='none', alpha=1, zorder=90, linewidth=0.3)
                    ax.add_artist(sac)
                                
                    sac = patches.Circle( (angle_prior, dataset.angle_of_saccade_array[i]), radius=r, facecolor=c, edgecolor='none', alpha=alpha, zorder=100, linewidth=0.3)
                    ax.add_artist(sac)
                    
                    if mirror:
                        sac = patches.Circle( (angle_prior*-1, dataset.angle_of_saccade_array[i]*-1), radius=r, facecolor='white', edgecolor='none', alpha=1, zorder=90, linewidth=0.3)
                        ax.add_artist(sac)
                                    
                        sac = patches.Circle( (angle_prior*-1, dataset.angle_of_saccade_array[i]*-1), radius=r, facecolor=c, edgecolor='none', alpha=alpha, zorder=100, linewidth=0.3)
                        ax.add_artist(sac)
                
    if regress is not None:
        linear_fit_type = 2
        
        if edge == 'middle':
            angle_prior = regress.angle_prior_array
        if edge == 'near':
            angle_prior = regress.angle_prior_array - np.sign(regress.angle_prior_array)*regress.angle_subtended_array/2.
        if edge == 'far':
            angle_prior = regress.angle_prior_array + np.sign(regress.angle_prior_array)*regress.angle_subtended_array/2.

        in_range = (regress.angle_subtended_array>angle_subtended_range[0])*(regress.angle_subtended_array<angle_subtended_range[1])*(angle_prior>retinal_position_range[0])*(angle_prior<retinal_position_range[1])
        in_range_indices = np.where(in_range)[0].tolist()    
    
        if regress.trajecs[regress.trajecs.keys()[0]].behavior == 'landing':
            fit, variance, intercept_confidence_interval, slope_confidence_interval, Rsq = floris.linear_fit_type(angle_prior[in_range_indices], regress.angle_of_saccade_array[in_range_indices], full_output=True, fit_type=linear_fit_type)
            std = np.sqrt(variance)
            
        if regress.trajecs[regress.trajecs.keys()[0]].behavior == 'flyby':
        
            neg_indices = np.where( regress.angle_of_saccade_array < angle_prior )
            flipped_angle_prior = copy.copy(angle_prior)
            flipped_angle_prior[neg_indices] *= -1
            flipped_angle_of_saccade = copy.copy(regress.angle_of_saccade_array)
            flipped_angle_of_saccade[neg_indices] *= -1
        
            fit, variance, intercept_confidence_interval, slope_confidence_interval, Rsq = floris.linear_fit_type(flipped_angle_prior[in_range_indices], flipped_angle_of_saccade[in_range_indices], full_output=True, fit_type=linear_fit_type)
            std = np.sqrt(variance)
            
        # make left / right turn datasets based on saccades within one standard deviation of the fit
        left_turn_angle_prior = []
        right_turn_angle_prior = []
        for i, angle in enumerate(angle_prior):
            if i in in_range_indices:
                if 1:
                    if regress.angle_of_saccade_array[i] > angle_prior[i]:
                        right_turn_angle_prior.append(angle_prior[i])
                    elif regress.angle_of_saccade_array[i] < angle_prior[i]:
                        left_turn_angle_prior.append(angle_prior[i])
              
        x = np.linspace(-np.pi, np.pi, 100)
        y = np.polyval(fit, x)
        
        print 'regression results: '
        print 'variance: ', variance
        print 'std dev: ', std
        print 'slope confidence interval: ', slope_confidence_interval
        print 'intercept confidence interval: ', intercept_confidence_interval
        
        if plot_regress:
            ax.plot(x,y,'-', color='black')
            ax.fill_between(x,y+std, y-std, facecolor='black', edgecolor='none', alpha=0.15, zorder=-100)
            ax.plot(x,x,'-', color='red', linewidth=0.5)
            
            if dataset.trajecs[keys[0]].behavior == 'flyby':
                fit_neg = copy.copy(fit)
                fit_neg[1] *= -1
                y = np.polyval(fit_neg, x)
                ax.plot(x,y,'-', color='black')
                ax.fill_between(x,y+std, y-std, facecolor='black', edgecolor='none', alpha=0.15, zorder=-100)
            
            s1 = 'y = ' + str(fit[0])[0:4] + 'x + ' + str(fit[1])[0:4]
            s2 = 'Rsq = ' + str(Rsq)
            ax.text(0, 0, s1)
            ax.text(0, 1, s2)
            
        if plot_color_regress:
            ax.plot(x,x,'-', color='red', linewidth=0.5)
            
            s1 = 'y = ' + str(fit[0])[0:4] + 'x + ' + str(fit[1])[0:4]
            s2 = 'Rsq = ' + str(Rsq)
            ax.text(0, 0, s1)
            ax.text(0, 1, s2)
            
            ax.fill_between(x,y+std, y-std,facecolor='blue', edgecolor='none', alpha=0.15, zorder=-100)
            ax.fill_between(x,-np.pi*np.ones_like(x), y-std, facecolor='orange', edgecolor='none', alpha=0.2, zorder=-100)
            ax.fill_between(x, y+std, np.ones_like(x)*np.pi, facecolor='orange', edgecolor='none', alpha=0.2, zorder=-100)
            
        
    if show_evasive_thresholds and regress is False:            
        x = np.linspace(-np.pi, np.pi, 100)
        yp = x + 50*np.pi/180.
        ym = x - 50*np.pi/180.
        #ax.plot(x,yp,'--',color='black', zorder=1000, linewidth=0.5)
        #ax.plot(x,ym,'--',color='black', zorder=1000, linewidth=0.5)
        ax.fill_between(x,ym,yp,facecolor='blue', edgecolor='none', alpha=0.15, zorder=-100)
        ax.fill_between(x,-np.pi*np.ones_like(yp), ym, facecolor='orange', edgecolor='none', alpha=0.2, zorder=-100)
        ax.fill_between(x, yp, np.ones_like(yp)*np.pi, facecolor='orange', edgecolor='none', alpha=0.2, zorder=-100)
        
    if show_tanh_evasive_thresholds:
        x = np.linspace(-100*np.pi/180., 100*np.pi/180., 100)
        var = 50*np.pi/180.
        #yp = np.tanh(x-var)/2.*np.pi*m + (1/2.*np.pi*m + var )
        yp = np.tanh(x-var)/2.*np.pi*2 + (1/2.*np.pi*2 + var + var)
        ym = np.tanh(x-var)/2.*np.pi*2 + (1/2.*np.pi*2 + var - var)
        ax.fill_between(x,ym,yp, facecolor='orange', edgecolor='none', alpha=0.2, zorder=-100)
        #ym = np.tanh(x+var)/2.*np.pi*m - (1/2.*np.pi*m + var )
    
    if plot:
        ax.set_ylim(-np.pi, np.pi)
        ax.set_xlim(-np.pi, np.pi)
        ax.set_autoscale_on(False)
        ax.set_aspect('equal')
        rad_ticks = [-np.pi, -np.pi/2., 0, np.pi/2., np.pi]
        deg_tick_strings = ['-180', '-90', '0', '90', '180']
        fa.adjust_spines(ax, ['left', 'bottom'])
        ax.set_xticks(rad_ticks)
        ax.set_xticklabels(deg_tick_strings) 
        ax.set_yticks(rad_ticks)
        ax.set_yticklabels(deg_tick_strings) 
        ax.set_xlabel('retinal position of post, deg')
        ax.set_ylabel('angle of saccade, deg')
        fig.subplots_adjust(bottom=0.2, top=0.9, right=0.9, left=0.2)
        
        dataset_type = dataset.trajecs[keys[0]].behavior
        if dataset.trajecs[keys[0]].post_type == 'none':
            dataset_type = 'nopost'
        if filename is None:
            filename = 'saccade_angle_' + dataset_type + '_' + edge + 'edge' + '.pdf'
        fig.savefig(filename, format='pdf')
        
    if return_left_right:
        return left_turn_angle_prior, right_turn_angle_prior
                        
    if plot is False and regress is not None:
        return angle_prior[in_range_indices], dataset.angle_of_saccade_array[in_range_indices], dataset.angle_subtended_array[in_range_indices], fit, variance
    if plot is False and regress is None:
        return angle_prior[in_range_indices], dataset.angle_of_saccade_array[in_range_indices], dataset.angle_subtended_array[in_range_indices]
    

def saccade_angle(dataset, keys=None, show_evasive_thresholds=False, flip=False, near=0*np.pi/180., far=180*np.pi/180., plot=True):
        
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    fig_attractiveness = plt.figure()
    ax_attractiveness = fig_attractiveness.add_subplot(111)
    
    colormap_norm = matplotlib.colors.Normalize(0, 0.1, clip=True)
    colormap_norm_dist = matplotlib.colors.Normalize(0, 0.2, clip=True)
    colormap_speed_norm = matplotlib.colors.Normalize(0, 0.4, clip=True)
    colormap_angle_norm = matplotlib.colors.Normalize(0, np.pi/2., clip=True)
    #colormap_norm = matplotlib.colors.Normalize(0, 0.7, clip=True)
    cmap = plt.get_cmap('jet')
    
    if keys is None:
        classified_keys = fa.get_classified_keys(dataset)
        keys = classified_keys['straight']
        keys = dataset.trajecs.keys()
        
    keys_with_saccades = []
    keys_without_saccades = []
    for key in keys:
        trajec = dataset.trajecs[key]    
        if len(trajec.saccades) > 0:
            keys_with_saccades.append(key)
        else:
            keys_without_saccades.append(key)
            
    speed_at_saccade_array = []
    angle_of_saccade_array = []
    angle_prior_array = []
    angle_subtended_array = []
    attractiveness_array = []
            
    for key in keys_with_saccades:
        trajec = dataset.trajecs[key]
        #sa1.calc_post_dynamics_for_flydra_trajectory(trajec)
        
        try:
            frames = trajec.frames_below_post
        except:
            frames = np.arange(0, len(trajec.speed)-1).tolist()
        
        for s in trajec.all_saccades:
            if s in frames:
            #if 1:
            #    s = trajec.saccades[-1]
                sac_range = fa.get_saccade_range(trajec, s)
                if sac_range is None:
                    continue
                if trajec.angle_subtended_by_post[sac_range[0]] < far and trajec.angle_subtended_by_post[sac_range[0]] > near:
                    if sac_range[0] <= trajec.frame_nearest_to_post:
                        
                        angle_prior = trajec.angle_to_post[sac_range[0]]*-1
                        
                        #angle_of_saccade = (trajec.worldangle[sac_range[-1]] - trajec.worldangle[sac_range[0]])
                        
                        f0 = sac_range[0]
                        f1 = sac_range[-1]

                        obj_ori_0 = trajec.velocities[f0] / np.linalg.norm(trajec.velocities[f0])   
                        obj_ori_1 = trajec.velocities[f1] / np.linalg.norm(trajec.velocities[f1])  

                        obj_ori_0_3vec = np.hstack( ( obj_ori_0, 0) ) 
                        obj_ori_1_3vec = np.hstack( (obj_ori_1, 0 ) ) 

                        sign_of_angle_of_saccade = np.sign( np.sum(np.cross( obj_ori_0, obj_ori_1 ) ) )

                        cosangleofsaccade = np.dot(obj_ori_0, obj_ori_1)
                        angleofsaccade = np.arccos(cosangleofsaccade)
                         
                        signed_angleofsaccade = -1*angleofsaccade*sign_of_angle_of_saccade
                        if flip:
                            if signed_angleofsaccade < 0:
                                signed_angleofsaccade *= -1
                                angle_prior *= -1
                        
                        '''
                        while angle_of_saccade > np.pi:
                            angle_of_saccade -= np.pi
                        while angle_of_saccade < -1*np.pi:
                            angle_of_saccade += np.pi
                        '''
                        #angle_after = trajec.angle_to_post[sac_range[0]] + angle_of_saccade
                        
                        dist_at_saccade = trajec.dist_to_stim_r_normed[s]
                        c = cmap(colormap_norm(0.1 - dist_at_saccade))
                        alpha = colormap_norm_dist(0.2 - dist_at_saccade)
                        #c = cmap(colormap_speed_norm(trajec.speed[s]))
                        #c = cmap(colormap_norm( np.log(trajec.angle_subtended_by_post[s])))
                        
                        speed_at_saccade_array.append(trajec.speed[sac_range[0]])
                        
                        if 'checkered' in trajec.post_type:
                            edgecolor = 'none'
                        elif 'black' in trajec.post_type:
                            edgecolor = 'none'
                        else:
                            edgecolor = 'none'
                            
                        r = colormap_speed_norm(trajec.speed[s])/100.*4 + 0.03
                        
                        angle_subtended_array.append( trajec.angle_subtended_by_post[sac_range[0]])
                        angle_prior_array.append(angle_prior)
                        angle_of_saccade_array.append(signed_angleofsaccade)
                        
                        if alpha > 0.1:
                            sac = patches.Circle( (angle_prior, signed_angleofsaccade), radius=r, facecolor='white', edgecolor=edgecolor, alpha=1, zorder=99, linewidth=0.3)
                            ax.add_artist(sac)
                        
                        sac = patches.Circle( (angle_prior, signed_angleofsaccade), radius=r, facecolor=c, edgecolor=edgecolor, alpha=alpha, zorder=100, linewidth=0.3)
                        ax.add_artist(sac)
                        
                        #attractiveness = floris.dist_point_to_line((angle_prior, signed_angleofsaccade), [-np.pi,-np.pi], [np.pi,np.pi], sign=False)
                        attractiveness = np.abs(angle_prior - signed_angleofsaccade)
                        attractiveness_array.append(attractiveness)
                        
                        c = cmap(colormap_angle_norm( np.pi-np.abs(angle_prior) ))
                        #pt = patches.Circle( ( np.log(trajec.angle_subtended_by_post[f0]), attractiveness), radius=0.05, facecolor='green', edgecolor='none', alpha=0.2, zorder=100)
                        #ax_attractiveness.add_artist(pt)
                        
                            
    angle_prior_array = np.array(angle_prior_array)
    angle_of_saccade_array = np.array(angle_of_saccade_array)
    angle_subtended_array = np.array(angle_subtended_array)
    speed_at_saccade_array = np.array(speed_at_saccade_array)
    attractiveness_array = np.array(attractiveness_array)
    
    
          
    if plot:
                        
        if show_evasive_thresholds:            
            x = np.linspace(-np.pi, np.pi, 100)
            yp = x + 50*np.pi/180.
            ym = x - 50*np.pi/180.
            #ax.plot(x,yp,'--',color='black', zorder=1000, linewidth=0.5)
            #ax.plot(x,ym,'--',color='black', zorder=1000, linewidth=0.5)
            ax.fill_between(x,ym,yp,facecolor='blue', edgecolor='none', alpha=0.15, zorder=-100)
            ax.fill_between(x,-np.pi*np.ones_like(yp), ym, facecolor='orange', edgecolor='none', alpha=0.2, zorder=-100)
            ax.fill_between(x, yp, np.ones_like(yp)*np.pi, facecolor='orange', edgecolor='none', alpha=0.2, zorder=-100)
            
        # regular ax
        #ax_hist = ax.twinx()
        ax.set_ylim(-np.pi, np.pi)
        ax.set_xlim(-np.pi, np.pi)
        ax.set_autoscale_on(False)
        
        rad_ticks = [-np.pi, -np.pi/2., 0, np.pi/2., np.pi]
        deg_tick_strings = ['-180', '-90', '0', '90', '180']
        
        fa.adjust_spines(ax, ['left', 'bottom'])
        ax.set_xticks(rad_ticks)
        ax.set_xticklabels(deg_tick_strings) 
        ax.set_yticks(rad_ticks)
        ax.set_yticklabels(deg_tick_strings) 
        
        ax.set_xlabel('retinal position of post, deg')
        ax.set_ylabel('angle of saccade, deg')
        
        
        #ax_hist.hist(angle_of_saccade_array, bins=20, facecolor='black', edgecolor='none', normed=True, orientation='horizontal', alpha=0.2)
        
        fig.subplots_adjust(bottom=0.2, top=0.9, right=0.9, left=0.2)
        ax.set_aspect('equal')
        fig.savefig('saccades.pdf', format='pdf')
        
        
        if 1:
            
            
            # attractiveness ax
            ax_attractiveness.set_ylim(0, 1)
            ax_attractiveness.set_xlim( np.log(5*np.pi/180.), np.log(np.pi))
            ax_attractiveness.set_autoscale_on(False)
            bins = 25
            ax_attractiveness.hist( angle_of_saccade_array, alpha=0.2, facecolor='green', edgecolor='none', normed=True, bins=bins)
            
            
            fa.fix_angle_log_spine(ax_attractiveness, set_y=False, histograms=False)
            ax_attractiveness.set_yticks([0, 0.5, 1])
            #ax_attractiveness.set_xlabel('retinal size, deg')
            ax_attractiveness.set_ylabel('attractiveness')
            
            fig_attractiveness.savefig('attractiveness.pdf', format='pdf')


            # histogram stuff
            '''
            speed_during_saccade_array = np.array(speed_during_saccade_array)
            angle_of_saccade_array = np.array(angle_of_saccade_array)
            
            fig_hist_speed = plt.figure()
            ax_hist_speed = fig_hist_speed.add_subplot(111)
            ax_hist_speed.hist(speed_during_saccade_array)
            fig_hist_speed.savefig('speed_during_saccade_hist.pdf', format='pdf')
            
            fig_hist_angle = plt.figure()
            ax_hist_angle = fig_hist_angle.add_subplot(111)
            ax_hist_angle.hist(angle_of_saccade_array)
            fig_hist_angle.savefig('angle_of_saccade_hist.pdf', format='pdf')

            print np.mean(angle_of_saccade_array)
            
            '''
            
            ## fit saccade turn magnitude
            fig_fit = plt.figure()
            ax_fit = fig_fit.add_subplot(111)
            
            colormap_norm_angle_subtended = matplotlib.colors.Normalize(near, 90*np.pi/180, clip=True)
            
            for i in range(len(angle_prior_array)):
                alpha = colormap_norm_angle_subtended(angle_subtended_array[i])
                ax_fit.plot( (angle_prior_array[i]), angle_of_saccade_array[i], '.', color='black', alpha=alpha)
                
                mirror = True
                if mirror:
                    mirrored_signed_angleofsaccade = -1*angle_of_saccade_array[i]
                    mirrored_angle_prior = -1*angle_prior_array[i]
                    ax_fit.plot( mirrored_angle_prior, mirrored_signed_angleofsaccade, '.', color='black', alpha=alpha)
                    
                r = np.random.random()
                if r < p:
                    sign = np.sign(angle_prior_array[i])
                else:
                    sign = -1*np.sign(angle_prior_array[i])
                
                
                var = 60*np.pi/180.
                rand = ((np.random.random()*2)-1)*var
                
                if sign > 0:
                    #tmp = np.max([5*np.pi/180., tmp])
                    m = 2
                    tmp = np.tanh(angle_prior_array[i]-var)/2.*np.pi*m + (1/2.*np.pi*m + var + rand)
                    #tmp = -1/np.arctan2(angle_prior_array[i])-0.5
                    #m = 1.5
                    #tmp = m*angle_prior_array[i]+safety*m
                    tmp = np.min([np.pi, tmp])
                else:
                    #tmp = np.min([-5*np.pi/180., tmp])
                    m = 2
                    tmp = np.tanh(angle_prior_array[i]+var)/2.*np.pi*m - (1/2.*np.pi*m + var + rand)
                    tmp = np.max([-1*np.pi, tmp])
                    
                ax_fit.plot( angle_prior_array[i], tmp, '.', color='red', alpha=alpha)
                
                
            x = np.linspace(-np.pi, np.pi, 100)
            yp = np.tanh(x-var)/2.*np.pi*m + (1/2.*np.pi*m + var )
            ym = np.tanh(x+var)/2.*np.pi*m - (1/2.*np.pi*m + var )
            
            ax_fit.plot(x,ym,color='red')
            ax_fit.plot(x,yp,color='red')
            
            #print fit
            ax_fit.set_xlim(-np.pi, np.pi)
            ax_fit.set_ylim(-np.pi,np.pi)
            ax_fit.set_aspect('equal')
            ax_fit.set_autoscale_on(False)
            
            rad_ticks_x = [-np.pi, -np.pi/2., 0, np.pi/2., np.pi]
            deg_tick_strings_x = ['-180', '-90', '0', '90', '180']
            
            rad_ticks_y = [0, np.pi/2., np.pi]
            deg_tick_strings_y = ['0', '90', '180']
            
            fa.adjust_spines(ax_fit, ['left', 'bottom'])
            ax_fit.set_xticks(rad_ticks_x)
            ax_fit.set_xticklabels(deg_tick_strings_x) 
            ax_fit.set_yticks(rad_ticks_x)
            ax_fit.set_yticklabels(deg_tick_strings_x) 
            
            ax_fit.set_xlabel('retinal position of post, deg')
            ax_fit.set_ylabel('angle of saccade, deg')
            
            fig_fit.subplots_adjust(bottom=0.2, top=0.9, right=0.9, left=0.2)
            
            fig_fit.savefig('saccade_turns_fit.pdf', format='pdf')
            
            
            ## speed
            fig_speed = plt.figure()
            ax_speed = fig_speed.add_subplot(111)
            
            for i in range(len(angle_prior_array)):
                
                c = cmap(colormap_speed_norm(speed_at_saccade_array[i]))
                ax_speed.plot( np.log(angle_subtended_array[i]), np.abs(angle_prior_array[i])/np.pi, '.', color=c)
            
            fa.fix_angle_log_spine(ax_speed, set_y=True, histograms=False)    
            
            ax_speed.set_xlabel('retinal position of post, deg')
            ax_speed.set_ylabel('speed of saccade, deg')
            
            fig_speed.subplots_adjust(bottom=0.2, top=0.9, right=0.9, left=0.2)
            
            fig_speed.savefig('saccade_speed.pdf', format='pdf')
        
    return angle_prior_array, angle_of_saccade_array, angle_subtended_array, speed_at_saccade_array, attractiveness_array
    

def saccade_distance_histograms(dataset, dataset_landing, keys=None):
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    bins = np.linspace( np.log(5*np.pi/180.), np.log(180*np.pi/180.), 20)
    
    #angle_prior, angle_of_saccade, angle_subtended, speed_at_saccade, attractiveness_array = saccade_angle(dataset, keys=None, show_evasive_thresholds=False, flip=False, near=0*np.pi/180., far=180*np.pi/180., plot=False)
    
    angle_prior_landing, angle_of_saccade_landing, angle_subtended_landing, fit, variance = plot_saccade_angle(dataset_landing, angle_subtended_range=[0,180], retinal_position_range=[-160, 160], show_evasive_thresholds=False, show_tanh_evasive_thresholds=False, edge='far', colorscale=False, speedscale=False, regress=dataset_landing, plot_regress=False, plot=False)
    std = np.sqrt(variance)
    
    angle_prior, angle_of_saccade, angle_subtended = plot_saccade_angle(dataset, keys=keys, angle_subtended_range=[0,180], retinal_position_range=[-160, 160], show_evasive_thresholds=False, show_tanh_evasive_thresholds=False, edge='far', colorscale=False, speedscale=False, regress=None, plot_regress=False, plot=False)
    
    #attractiveness_array -= angle_subtended/2.
    
    predicted_turn = np.polyval(fit, angle_prior)
    nonevasive = np.where( (angle_of_saccade < predicted_turn+std)*(angle_of_saccade > predicted_turn-std) )[0].tolist()
    evasive = np.where( (angle_of_saccade > predicted_turn+std) + (angle_of_saccade < predicted_turn-std) )[0].tolist()
    print len(evasive), len(nonevasive), len(angle_prior)
    '''
    threshold = 50*np.pi/180.
    evasive = np.where(attractiveness_array>threshold)[0].tolist()
    nonevasive = np.where(attractiveness_array<=threshold)[0].tolist()
    '''
    
    bins, hists, hist_std, curves = floris.histogram(ax, [np.log(angle_subtended[evasive]), np.log(angle_subtended[nonevasive])], bins=bins, colors=['orange', 'blue'], return_vals=True, bootstrap_std=True, normed_occurences='total', curve_line_alpha=0)
    
    #ax.hist( np.log(angle_subtended[evasive]), bins=bins, facecolor='orange', edgecolor='none', alpha=0.5, normed=False)
    #ax.hist( np.log(angle_subtended[nonevasive]), bins=bins, facecolor='blue', edgecolor='none', alpha=0.5, normed=False)
    
    evasive_hist = hists[0]
    nonevasive_hist = hists[1]
    
    ymax = 1
    
    deg_ticks = np.array([5, 10, 30, 60, 90, 180])
    deg_tick_strings = [str(d) for d in deg_ticks]
    rad_ticks = deg_ticks*np.pi/180.
    rad_ticks_log = np.log(rad_ticks)
    
    dist_tick_strings = ['(21)', '(10)', '(2.7)', '(0.9)', '(0.4)', '(0)']
    x_tick_strings = []
    for i, d in enumerate(dist_tick_strings):
        x_tick_strings.append( deg_tick_strings[i] + '\n' + dist_tick_strings[i] )
    
    ax.set_xlim(rad_ticks_log[0], rad_ticks_log[-1])
    ax.set_ylim(0,1.2)
    
    fa.adjust_spines(ax,['left', 'bottom'])
    ax.set_xlabel('retinal size, deg\n(distance, cm)', labelpad=10)
    ax.set_ylabel('occurances')
    
    ax.set_xticks(rad_ticks_log.tolist())
    ax.set_xticklabels(x_tick_strings) 
    
    
    fig.subplots_adjust(bottom=0.25, top=0.9, right=0.9, left=0.2)
    
    dataset_type = dataset.trajecs[dataset.trajecs.keys()[0]].behavior
    if dataset.trajecs[dataset.trajecs.keys()[0]].post_type == 'none':
        dataset_type = 'nopost'
    filename = 'saccade_histogram_' + dataset_type + '.pdf'
        
    fig.savefig(filename, format='pdf')


def saccade_distance_histograms_post_type_coloring(dataset, dataset_landing, keys=None, saccade_type='evasive'):
    calc_saccade_stats(dataset, keys=keys)
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    bins = np.linspace( np.log(5*np.pi/180.), np.log(180*np.pi/180.), 20)
    
    #angle_prior, angle_of_saccade, angle_subtended, speed_at_saccade, attractiveness_array = saccade_angle(dataset, keys=None, show_evasive_thresholds=False, flip=False, near=0*np.pi/180., far=180*np.pi/180., plot=False)
    
    angle_prior_landing, angle_of_saccade_landing, angle_subtended_landing, fit, variance = plot_saccade_angle(dataset_landing, angle_subtended_range=[0,180], retinal_position_range=[-160, 160], show_evasive_thresholds=False, show_tanh_evasive_thresholds=False, edge='far', colorscale=False, speedscale=False, regress=dataset_landing, plot_regress=False, plot=False)
    std = np.sqrt(variance)
    
    
    predicted_turn = np.polyval(fit, dataset.angle_prior_array)
    
    black_post = np.zeros_like(dataset.angle_subtended_array)
    checkered_post = np.zeros_like(dataset.angle_subtended_array)
    for i in range(len(black_post)):
        if 'black' in dataset.post_types_saccade_array[i]:
            black_post[i] = 1
        if 'checkered' in dataset.post_types_saccade_array[i]:
            checkered_post[i] = 1
    
    evasive = np.zeros_like(dataset.angle_subtended_array)
    nonevasive = np.zeros_like(dataset.angle_subtended_array)
    for i in range(len(evasive)):
        if (dataset.angle_of_saccade_array[i] < predicted_turn[i]+std)*(dataset.angle_of_saccade_array[i] > predicted_turn[i]-std):
            nonevasive[i] = 1
        else:
            evasive[i] = 1
            
    if saccade_type == 'evasive':
        black_keys = np.where( black_post*evasive == 1 )[0].tolist()
        checkered_keys = np.where( checkered_post*evasive == 1 )[0].tolist()
    if saccade_type == 'nonevasive':
        black_keys = np.where( black_post*nonevasive == 1 )[0].tolist()
        checkered_keys = np.where( checkered_post*nonevasive == 1 )[0].tolist()

    
    bins, hists, hist_std, curves = floris.histogram(ax, [np.log(dataset.angle_subtended_array[black_keys]), np.log(dataset.angle_subtended_array[checkered_keys])], bins=bins, colors=['black', 'teal'], return_vals=True, normed_occurences='total', curve_line_alpha=0, bootstrap_std=True)
    
    #ax.hist( np.log(angle_subtended[evasive]), bins=bins, facecolor='orange', edgecolor='none', alpha=0.5, normed=False)
    #ax.hist( np.log(angle_subtended[nonevasive]), bins=bins, facecolor='blue', edgecolor='none', alpha=0.5, normed=False)
    
    black_hist = hists[0]
    checkered_hist = hists[1]
    
    ymax = 1
    
    deg_ticks = np.array([5, 10, 30, 60, 90, 180])
    deg_tick_strings = [str(d) for d in deg_ticks]
    rad_ticks = deg_ticks*np.pi/180.
    rad_ticks_log = np.log(rad_ticks)
    
    dist_tick_strings = ['(21)', '(10)', '(2.7)', '(0.9)', '(0.4)', '(0)']
    x_tick_strings = []
    for i, d in enumerate(dist_tick_strings):
        x_tick_strings.append( deg_tick_strings[i] + '\n' + dist_tick_strings[i] )
    
    ax.set_xlim(rad_ticks_log[0], rad_ticks_log[-1])
    ax.set_ylim(0,1.2)
    
    fa.adjust_spines(ax,['left', 'bottom'])
    ax.set_xlabel('Retinal size\n(Distance, cm)', labelpad=10)
    ax.set_ylabel('Occurences (normalized)')
    
    ax.set_xticks(rad_ticks_log.tolist())
    ax.set_xticklabels(x_tick_strings) 
    
    
    fig.subplots_adjust(bottom=0.25, top=0.9, right=0.9, left=0.2)
    
    dataset_type = dataset.trajecs[dataset.trajecs.keys()[0]].behavior
    if dataset.trajecs[dataset.trajecs.keys()[0]].post_type == 'none':
        dataset_type = 'nopost'
    filename = 'saccade_histogram_checker_vs_black' + dataset_type + saccade_type + '.pdf'
        
    fig.savefig(filename, format='pdf')



def colorbar():
    fig = plt.figure()
    cax = fig.add_subplot(111)
    
    cticks = np.linspace(0, 0.1, 5, endpoint=True)
    cb = matplotlib.colorbar.ColorbarBase(cax, cmap=plt.get_cmap('jet'), norm=plt.Normalize(0, 0.1), orientation='horizontal', boundaries=None, ticks=cticks, drawedges=False)
    
    fig.savefig('colorbar_for_saccades.pdf', format='pdf')
    
def radius_scale():
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    def get_r(s):
        colormap_speed_norm = matplotlib.colors.Normalize(0, 0.4, clip=True)
        r = colormap_speed_norm(s)/100.*4 + 0.03
        return r
        
    speeds = [0, 0.2, 0.4, 0.6, 0.8]
    
    for s in speeds:
        r = get_r(s)
        circle = patches.Circle( (0, s), radius=r, facecolor='black', edgecolor='none', alpha=1)
        ax.add_artist(circle)
                            
    ax.set_aspect('equal')
    ax.set_xlim(-np.pi, np.pi)
    ax.set_ylim(-.2, 1.0)
    ax.set_autoscale_on(False)
    
    fa.adjust_spines(ax,['left', 'bottom'])

    ax.set_yticks(speeds)

    ax.set_ylabel('speed scale, m/s')    
    
    fig.savefig('radius_scale_for_saccades.pdf', format='pdf')
    
    
def is_evasive(trajec, s):
    sac_range = fa.get_saccade_range(trajec, s)
    if sac_range is None:
        return None
    angle_prior = trajec.angle_to_post[sac_range[0]]*-1
    
    f0 = sac_range[0]
    f1 = sac_range[-1]

    obj_ori_0 = trajec.velocities[f0] / np.linalg.norm(trajec.velocities[f0])   
    obj_ori_1 = trajec.velocities[f1] / np.linalg.norm(trajec.velocities[f1])  

    obj_ori_0_3vec = np.hstack( ( obj_ori_0, 0) ) 
    obj_ori_1_3vec = np.hstack( (obj_ori_1, 0 ) ) 

    sign_of_angle_of_saccade = np.sign( np.sum(np.cross( obj_ori_0, obj_ori_1 ) ) )

    cosangleofsaccade = np.dot(obj_ori_0, obj_ori_1)
    angleofsaccade = np.arccos(cosangleofsaccade)
     
    signed_angleofsaccade = -1*angleofsaccade*sign_of_angle_of_saccade
    attractiveness = np.abs(angle_prior - signed_angleofsaccade)
    
    if np.abs(attractiveness) < 45*np.pi/180.:
        return False
    else:
        if np.abs(signed_angleofsaccade) < 45*np.pi/180.:
            return True
        else:
            return False    
        
        

                            
                            
        
        
def save_landing_edge_option_plots(dataset_landing):
    
    plot_saccade_angle(dataset_landing, show_evasive_thresholds=False, edge='near', retinal_position_range=[-160,160], regress=dataset_landing, plot_regress=True)
    plot_saccade_angle(dataset_landing, show_evasive_thresholds=False, edge='middle', retinal_position_range=[-160,160], regress=dataset_landing, plot_regress=True)
    plot_saccade_angle(dataset_landing, show_evasive_thresholds=False, edge='far', retinal_position_range=[-160,160], regress=dataset_landing, plot_regress=True)    
    
    
def save_color_far_edge_all_behaviors(dataset_landing, dataset_flyby, dataset_nopost):
    
    plot_saccade_angle(dataset_landing, show_evasive_thresholds=False, edge='far', retinal_position_range=[-160,160], regress=dataset_landing, plot_regress=False, colorscale=True, speedscale=False, plot_color_regress=True)
    plot_saccade_angle(dataset_flyby, show_evasive_thresholds=False, edge='far', retinal_position_range=[-160,160], regress=dataset_landing, plot_regress=False, colorscale=True, speedscale=False, plot_color_regress=True)
    plot_saccade_angle(dataset_nopost, show_evasive_thresholds=False, edge='far', retinal_position_range=[-160,160], regress=dataset_landing, plot_regress=False, colorscale=True, speedscale=False, plot_color_regress=True)    
    
def save_post_type_far_edge_all_behaviors(dataset_landing, dataset_flyby, dataset_nopost):
    
    plot_saccade_angle(dataset_landing, show_evasive_thresholds=False, edge='far', retinal_position_range=[-160,160], regress=dataset_landing, plot_regress=False, colorscale=False, speedscale=False, plot_color_regress=True, postcolorscale=True)
    
    plot_saccade_angle(dataset_flyby, show_evasive_thresholds=False, edge='far', retinal_position_range=[-160,160], regress=dataset_landing, plot_regress=False, colorscale=False, speedscale=False, plot_color_regress=True, postcolorscale=True)
    
    plot_saccade_angle(dataset_nopost, show_evasive_thresholds=False, edge='far', retinal_position_range=[-160,160], regress=dataset_landing, plot_regress=False, colorscale=False, speedscale=False, plot_color_regress=True, postcolorscale=True)    
    
def save_flyby_evasive_fit(dataset_flyby):

    plot_saccade_angle(dataset_flyby, show_evasive_thresholds=False, edge='far', retinal_position_range=[-160,160], angle_subtended_range=[20,180], regress=dataset_flyby, plot_regress=True, colorscale=True, speedscale=False, plot_color_regress=False, mirror=False)


def save_left_right_turn_histogram(dataset_flyby, mirror=False):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    left, right = plot_saccade_angle(dataset_flyby, show_evasive_thresholds=False, edge='far', retinal_position_range=[-160,160], angle_subtended_range=[20,180], regress=dataset_flyby, plot_regress=True, colorscale=True, speedscale=True, plot_color_regress=False, mirror=True, return_left_right=True)
    
    print '*** left ***'
    print left
    
    bins = np.linspace(-160*np.pi/180., 160*np.pi/180., 14, endpoint=True)
    bins, hists, hist_std, curves = floris.histogram(ax, [left, right], bins=bins, colors=['black', 'gray'], return_vals=True, show_smoothed=False, bootstrap_std=True)
    
    ### 
    left_raw_hist = floris.bootstrap_histogram(left, bins=bins, return_raw=True)
    right_raw_hist = floris.bootstrap_histogram(right, bins=bins, return_raw=True)

    perc_left_turn = None

    for i in range(left_raw_hist.shape[0]):
        for j in range(right_raw_hist.shape[0]):
            p = left_raw_hist[i] / (left_raw_hist[i]+right_raw_hist[j])
            p = np.nan_to_num(p)
            
            if perc_left_turn is not None:
                perc_left_turn = np.vstack( (perc_left_turn, p) )
            else:
                perc_left_turn = p
                
    perc_left_turn_mean = np.mean(perc_left_turn, axis=0)
    perc_left_turn_std = np.std(perc_left_turn, axis=0)
    ###
    
    
    left_hist = np.array(hists[0], dtype=float)
    right_hist = np.array(hists[1], dtype=float)
    
    '''
    left_std = np.array(hist_std[0], dtype=float)
    right_std = np.array(hist_std[1], dtype=float)
    
    if mirror:
        percent_turn_left = (left_hist + right_hist[::-1])  / (left_hist+right_hist+right_hist[::-1]+left_hist[::-1])
        percent_turn_right = (right_hist + left_hist[::-1])  / (left_hist+right_hist+right_hist[::-1]+left_hist[::-1])
    else:
        percent_turn_left = left_hist / (left_hist+right_hist)
        percent_turn_left_std = np.abs( (left_hist + left_std) / (left_hist+right_hist+left_std+right_std) - percent_turn_left)
        percent_turn_right = right_hist / (left_hist+right_hist)
    '''
    
    ymax = np.max([np.max(left_hist), np.max(right_hist)])
    ax.set_ylim(0, 70)
    ax.set_xlim(-np.pi, np.pi)
    ax.set_autoscale_on(False)
    rad_ticks = [-np.pi, -np.pi/2., 0, np.pi/2., np.pi]
    deg_tick_strings = ['-180', '-90', '0', '90', '180']
    fa.adjust_spines(ax, ['left', 'bottom'])
    ax.set_xticks(rad_ticks)
    ax.set_xticklabels(deg_tick_strings) 
    ax.set_xlabel('angle to far edge')
    ax.set_ylabel('occurrences')
    fig.subplots_adjust(bottom=0.2, top=0.9, right=0.9, left=0.2)
    fig.savefig('turn_direction_histograms.pdf', format='pdf')
    
    fig = plt.figure()
    percax = fig.add_subplot(111)
    bincenters = np.diff(bins)/2. + bins[0:-1]
    percax.plot(bincenters, perc_left_turn_mean, color='black')
    percax.plot(bincenters, perc_left_turn_mean, '.', color='black')
    #percax.plot(bincenters, percent_turn_right, color='green')
    percax.vlines(bincenters, perc_left_turn_mean-perc_left_turn_std, perc_left_turn_mean+perc_left_turn_std, color='black')
    
    percax.set_ylim(0, 1)
    percax.set_xlim(-np.pi, np.pi)
    percax.set_autoscale_on(False)
    rad_ticks = [-np.pi, -np.pi/2., 0, np.pi/2., np.pi]
    deg_tick_strings = ['-180', '-90', '0', '90', '180']
    fa.adjust_spines(percax, ['left', 'bottom'])
    percax.set_xticks(rad_ticks)
    percax.set_xticklabels(deg_tick_strings) 
    percax.set_xlabel('angle to far edge')
    percax.set_ylabel('percent turns to left')
    fig.subplots_adjust(bottom=0.2, top=0.9, right=0.9, left=0.2)
    fig.savefig('turn_direction_probability.pdf', format='pdf')
    
    
    
def get_percent_same_side_turn(dataset_flyby):

    left, right = plot_saccade_angle(dataset_flyby, show_evasive_thresholds=False, edge='far', retinal_position_range=[-160,160], angle_subtended_range=[20,180], regress=dataset_flyby, plot_regress=True, colorscale=True, speedscale=True, plot_color_regress=False, mirror=True, return_left_right=True)
    
    bins = np.linspace(-np.pi, np.pi, 14, endpoint=True)
    bins, hists, curves = floris.histogram(ax, [left, right], bins=bins, colors=['crimson', 'green'], return_vals=True, show_smoothed=False)
    

def save_histograms(dataset_landing, dataset_flyby, dataset_nopost):

    post_type = ['black', 'black_angled']
    keys = None

    if post_type is not None:
        keys = fa.get_keys_for_post_type(dataset_landing, post_type=post_type)
    saccade_distance_histograms(dataset_landing, dataset_landing, keys=keys)
    if post_type is not None:
        keys = fa.get_keys_for_post_type(dataset_flyby, post_type=post_type)
    saccade_distance_histograms(dataset_flyby, dataset_landing, keys=keys)
    saccade_distance_histograms(dataset_nopost, dataset_landing)
    
    
    


def plot_residuals(dataset_landing, edge='far'):

    angle_prior_landing, angle_of_saccade_landing, angle_subtended_landing, fit, variance = plot_saccade_angle(dataset_landing, angle_subtended_range=[0,180], retinal_position_range=[-160, 160], show_evasive_thresholds=False, show_tanh_evasive_thresholds=False, edge=edge, colorscale=False, speedscale=False, regress=dataset_landing, plot_regress=False, plot=False)
    std = np.sqrt(variance)
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    predicted_turns = np.polyval(fit, angle_prior_landing)
    residuals = (angle_of_saccade_landing - predicted_turns) 
    
    indices = np.where( np.abs(residuals) <= np.abs(np.pi - np.abs(predicted_turns)) )[0].tolist()
        
    ax.plot(angle_prior_landing[indices], residuals[indices], '.', color='black')
    
    fit = floris.linear_fit_type2(angle_prior_landing[indices], residuals[indices], alpha=0.05, full_output=False, weights=None)
    x = np.linspace(-np.pi, np.pi, 100)
    y = np.polyval(fit, x)
    ax.plot(x,y,'-', color='black')
    
    ax.set_ylim(-np.pi, np.pi)
    ax.set_xlim(-np.pi, np.pi)
    ax.set_autoscale_on(False)
    ax.set_aspect('equal')
    rad_ticks = [-np.pi, -np.pi/2., 0, np.pi/2., np.pi]
    deg_tick_strings = ['-180', '-90', '0', '90', '180']
    fa.adjust_spines(ax, ['left', 'bottom'])
    ax.set_xticks(rad_ticks)
    ax.set_xticklabels(deg_tick_strings) 
    ax.set_yticks(rad_ticks)
    ax.set_yticklabels(deg_tick_strings) 
    ax.set_xlabel('retinal position of post, deg')
    ax.set_ylabel('residual, deg')
    fig.subplots_adjust(bottom=0.2, top=0.9, right=0.9, left=0.2)
    
    dataset_type = dataset_landing.trajecs[dataset_landing.trajecs.keys()[0]].behavior
    if dataset_landing.trajecs[dataset_landing.trajecs.keys()[0]].post_type == 'none':
        dataset_type = 'nopost'
    filename = 'saccade_angle_fit_residual' + dataset_type + '_' + edge + 'edge' + '.pdf'
    fig.savefig(filename, format='pdf')
    
    

def save_landing_residuals(dataset_landing):

    plot_residuals(dataset_landing, 'far')
    plot_residuals(dataset_landing, 'middle')
    plot_residuals(dataset_landing, 'near')




def get_retinal_size_at_saccade_with_time_delay(trajec, s, delay):
    
    sac_range = fa.get_saccade_range(trajec, s)
    if sac_range is None:
        return None
    else:

        if sac_range[0] <= trajec.frame_nearest_to_post:
            t_sac = trajec.epoch_time[sac_range[0]]
            t_sac_delayed = t_sac - delay
            angle_subtended_delayed = np.interp(t_sac_delayed, trajec.epoch_time, trajec.angle_subtended_by_post)
            return angle_subtended_delayed
            
        else:
            return None
    
    
    
def plot_delay_effect_on_average_saccade_stats():
    
    
    speed = 0.3
    dist = 0.0369
    
    delays = np.linspace(0,0.15,50)
    delayed_dist = dist + delays*speed
    delayed_retinal_size = 2*np.sin(0.009565 / delayed_dist)
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    ax.plot(delays, delayed_retinal_size*180./np.pi)
    ax.vlines(0.049, (delayed_retinal_size*180./np.pi)[0], (delayed_retinal_size*180./np.pi)[-1], linestyle=':')
    
    
    fa.adjust_spines(ax, ['left', 'bottom'])
    
    filename = 'delayed_saccade_retinal_size' + '.pdf'
    fig.savefig(filename, format='pdf')
    
    
    
    
    
    
        
    
    
    
    
