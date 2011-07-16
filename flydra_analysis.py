import sys
#sys.path.insert(0, '/usr/local/lib/python2.6/dist-packages')
sys.path.append('/home/floris/src/pymovie2')

from matplotlib import rcParams
fig_width = 3.6 # width in inches
fig_height = 3.6  # height in inches
fig_size =  (fig_width, fig_height)

fontsize = 9
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
          'figure.subplot.right': 0.9,
          'figure.subplot.bottom': 0.3,
          'figure.subplot.top': 0.9,
          'figure.subplot.wspace': 0.0,
          'figure.subplot.hspace': 0.0,
          'lines.linewidth': 1.0,
          'text.usetex': True, 
          }
rcParams.update(params) 


import numpy as np
#from scipy.optimize import curve_fit
import scipy.signal as signal

import tables
import time
import datetime
import numpy as np

import pickle
import copy
import floris

import analysis_plot as ap
import sa1_analysis as sa1 
import numpyimgproc as nim
import colorline
import flydra_floris_analysis as ffa

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf as pdf
import matplotlib.colorbar
import matplotlib.patches as patches
#from mpl_toolkits.axes_grid1 import make_axes_locatable

import saccade_analysis
import numpyimgproc as nim

REQUIRED_LENGTH = 30
REQUIRED_DIST = 0.1



def dist_pt_to_line(p, p0, p1):

    dy = p1[1]-p0[1]
    dx = p1[0]-p0[0]
    
    return np.abs( ( dy*(p[0]-p0[0]) - dx*(p[1]-p0[1]) ) / np.sqrt( dx**2 + dy**2 ) )


def diff_windowed(arr, window):
    diff = np.zeros_like(arr)
    
    for i in range(window, len(diff)-window-1):
        v1 = np.mean(arr[i-window:i])
        v2 = np.mean(arr[i:window+i])
        
        #v1 = arr[i-window]
        #v2 = arr[window+i]
        
        diff[i] = (v2-v1) / 2.
        
    diff[0:window] = diff[window]*np.ones_like(diff[0:window])
    diff[len(diff)-window-1:len(diff)-1] = diff[len(diff)-window-2]*np.ones_like(diff[len(diff)-window-1:len(diff)-1])
        
    return diff
    
    
def find_extrema_points(arr, window=5, sign=None):
    arr_diff = diff_windowed(arr, window)
    return find_critical_points(arr_diff, sign=sign)
        
def find_critical_points(arr, sign=None):
    if sign is None:
        return np.where( diff_windowed( np.sign(arr), 1 ) != 0 )[0].tolist()
    elif sign is 'positive':
        return np.where( diff_windowed( np.sign(arr), 1 ) > 0 )[0].tolist()
    elif sign is 'negative':
        return np.where( diff_windowed( np.sign(arr), 1 ) < 0 )[0].tolist()
    
    
def cumsum(array):
    cumarr = copy.copy(array)
    for i in range(1, len(cumarr)):
        cumarr[i] += cumarr[i-1]
    return cumarr

def load(filename, prep=False):
    fname = (filename)
    fd = open( fname, mode='r')
    print 'loading data... '
    dataset = pickle.load(fd)
    
    if prep is True:
        prep_dataset(dataset)
        
    return dataset

def save(dataset, filename):
    print 'saving data to file: ', filename
    fname = (filename)  
    fd = open( fname, mode='w' )
    pickle.dump(dataset, fd)
    return 1
    
def show_dataset_boundary(dataset):
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    for k, trajec in dataset.trajecs.items():
        xmax = np.max(trajec.positions[:,0])
        xmin = np.min(trajec.positions[:,0])
        ymax = np.max(trajec.positions[:,1])
        ymin = np.min(trajec.positions[:,1])
        zmax = np.max(trajec.positions[:,2])
        zmin = np.min(trajec.positions[:,2])
    
        ax.plot(zmax, ymax, '.', color='black')
        ax.plot(zmin, ymin, '.', color='black')
        
    ax.set_xlim(-1.5, 1.5)
    ax.set_ylim(-1.5, 1.5)
    ax.set_autoscale_on(False)
        
    fig.savefig('trajectory_boundaries.pdf', format='pdf')
    
def make_no_post_dataset(dataset):
    #get_post_type_for_dataset(dataset)
    new_dataset = ffa.Dataset(like=dataset)
    for k,trajec in dataset.trajecs.items():
        trajec.key = k
        if trajec.post_type == 'none' and len(trajec.speed) > 30:
            trajec.behavior = 'flyby'
            new_dataset.trajecs.setdefault(k, trajec)
            trajec.calc_dist_to_stim_r()
    return new_dataset
    
def make_post_dataset(dataset):
    #get_post_type_for_dataset(dataset)
    new_dataset = ffa.Dataset(like=dataset)
    for k,trajec in dataset.trajecs.items():
        trajec.key = k
        if trajec.post_type != 'none' and len(trajec.speed) > 30:
            new_dataset.trajecs.setdefault(k, trajec)
    return new_dataset
    
    
def make_behavior_dataset(dataset, filename='landing_dataset_10cm', behavior='landing'):
    new_dataset = ffa.Dataset(like=dataset)
    if type(behavior) is not list:
        behavior = [behavior]
    for k,trajec in dataset.trajecs.items():
        trajec.key = k
        if trajec.behavior in behavior:
            
            calc_frame_of_landing (trajec, threshold = 0.0005)
            normalize_dist_to_stim_r(trajec)
            
            if trajec.behavior == 'landing' and trajec.post_type != 'none':
                if trajec.dist_to_stim_r_normed[0] >= REQUIRED_DIST:
                    #if np.max(trajec.positions[:,2]) < 0 and np.min(trajec.positions[:,2]) > -0.15:
                    if np.min(trajec.positions[:,2]) > -0.15:
                        trajec.frames = np.arange(get_frame_at_distance(trajec, REQUIRED_DIST), trajec.frame_of_landing).tolist()
                        if trajec.frame_of_landing > REQUIRED_LENGTH:
                            if np.max(trajec.positions[trajec.frames,2]) < 0: # altitude check: needs to be below post top
                                #classify(trajec, dfar=REQUIRED_DIST, dnear=0.005)
                                new_dataset.trajecs.setdefault(k, trajec)
                                
                                frames_below_post = np.where(trajec.positions[:,2]<0)[0]
                                
                                continuous_blocks_of_frames = np.abs(floris.diffa(frames_below_post))
                                break_pts = np.where( continuous_blocks_of_frames != 1)[0].tolist()
                                continuous_blocks_of_frames[break_pts] = 0
                                continuous_block_before_landing = nim.find_blob_nearest_to_point(continuous_blocks_of_frames, trajec.frame_of_landing)
                                first_frame = frames_below_post[np.where(continuous_block_before_landing == 1)[0][0]]
                                trajec.frames_below_post = np.arange(first_frame, trajec.frame_of_landing+1).tolist()

                                if np.max(trajec.positions[trajec.frames_below_post,2]) > 0:
                                    print k, np.max(trajec.positions[trajec.frames_below_post,2])
                                
            elif trajec.behavior == 'flyby' and trajec.post_type != 'none':
                frame_nearest_to_post = np.argmin(trajec.dist_to_stim_r)
                print k
                if frame_nearest_to_post > 10 and np.max(trajec.dist_to_stim_r[0:frame_nearest_to_post]) > REQUIRED_DIST:
                    if np.max(trajec.positions[:,2]) < 0 and np.min(trajec.positions[:,2]) > -0.15:
                        if trajec.dist_to_stim_r[frame_nearest_to_post] < 0.1:
                            fs = np.arange(frame_nearest_to_post,len(trajec.speed)).tolist()
                            try:
                                last_frame = get_frame_at_distance(trajec, REQUIRED_DIST, frames=fs)
                            except:
                                last_frame = len(trajec.speed)-1
                            first_frame = get_frame_at_distance(trajec, REQUIRED_DIST, frames=np.arange(0,frame_nearest_to_post).tolist())
                            
                            trajec.frames = np.arange(first_frame, last_frame).tolist()
                            
                            # get frame at 8cm away, prior to nearest approach
                            frame_nearest_to_post = np.argmin(trajec.dist_to_stim_r)
                            trajec.frame_nearest_to_post = frame_nearest_to_post
                            frames = np.arange(0, frame_nearest_to_post).tolist()
                            trajec.frames_of_flyby = frames
                            frame_at_distance = get_frame_at_distance(trajec, 0.08, singleframe=True, frames=frames)
                            
                            last_frame = np.min( [frame_nearest_to_post+20, len(trajec.speed)-1]) 
                            calc_heading(trajec)
                            calc_saccades(trajec)
                            
                            sac_sgns = np.array(trajec.all_saccades) - frame_at_distance
                            sac_negs = np.where(sac_sgns<0)[0]
                            if len(sac_negs) > 0:
                                sac_neg = sac_negs[0]
                            else:
                                sac_neg = 0
                            first_frame = sac_neg + 1
                            
                            try:
                                trajec.frames_of_flyby = np.arange(first_frame, last_frame).tolist()
                                new_dataset.trajecs.setdefault(k, trajec)
                            except:
                                print 'ignored key: ', k, first_frame, last_frame
                            
                            new_dataset.trajecs.setdefault(k, trajec)
                            
            elif trajec.behavior == 'flyby' and trajec.post_type == 'none':
                frame_nearest_to_post = np.argmin(trajec.dist_to_stim_r)
                #print frame_nearest_to_post, np.max(trajec.dist_to_stim_r[0:frame_nearest_to_post]), trajec.positions[frame_nearest_to_post,2], trajec.positions[frame_nearest_to_post,2], trajec.dist_to_stim_r[frame_nearest_to_post]
                if frame_nearest_to_post > 10 and np.max(trajec.dist_to_stim_r[0:frame_nearest_to_post]) > REQUIRED_DIST and np.mean(trajec.speed) > 0.05:
                    if trajec.positions[frame_nearest_to_post,2] < 0.1 and trajec.positions[frame_nearest_to_post,2] > -0.15:
                        if trajec.dist_to_stim_r[frame_nearest_to_post] < 0.1:
                            fs = np.arange(frame_nearest_to_post,len(trajec.speed)).tolist()
                            try:
                                last_frame = get_frame_at_distance(trajec, REQUIRED_DIST, frames=fs)
                            except:
                                last_frame = len(trajec.speed)-1
                            first_frame = get_frame_at_distance(trajec, REQUIRED_DIST, frames=np.arange(0,frame_nearest_to_post).tolist())
                            
                            trajec.frames = np.arange(first_frame, last_frame).tolist()
                            
                            # get frame at 8cm away, prior to nearest approach
                            frame_nearest_to_post = np.argmin(trajec.dist_to_stim_r)
                            trajec.frame_nearest_to_post = frame_nearest_to_post
                            frames = np.arange(0, frame_nearest_to_post).tolist()
                            trajec.frames_of_flyby = frames
                            frame_at_distance = get_frame_at_distance(trajec, 0.08, singleframe=True, frames=frames)
                            
                            last_frame = np.min( [frame_nearest_to_post+20, len(trajec.speed)-1]) 
                            calc_heading(trajec)
                            calc_saccades(trajec)
                            
                            
                            sac_sgns = np.array(trajec.all_saccades) - frame_at_distance
                            sac_negs = np.where(sac_sgns<0)[0]
                            if len(sac_negs) > 0:
                                sac_neg = sac_negs[0]
                            else:
                                sac_neg = 0
                            first_frame = sac_neg + 1
                            
                            try:
                                trajec.frames_of_flyby = np.arange(first_frame, last_frame).tolist()
                                new_dataset.trajecs.setdefault(k, trajec)
                            except:
                                print 'ignored key: ', k, first_frame, last_frame
                            
                            
                
              
            
    save(new_dataset, filename)

    return new_dataset
    
def add_datasets(dataset_list):
    
    new_dataset = dataset_list[0]
    
    for d, dataset in enumerate(dataset_list[1:]):
        for k,trajec in dataset.trajecs.items():
            if k not in new_dataset.trajecs.keys():
                new_dataset.trajecs.setdefault(k, trajec)
            else:
                print 'duplicate id: ', k, new_dataset.trajecs[k].epoch_time[0] - trajec.epoch_time[0]
                new_k = str(time.time()).split('.')[0] + '_' + k
                trajec.key = new_k
                new_dataset.trajecs.setdefault(new_k, trajec)
                
    return new_dataset
    
def get_min_rrev(dataset, keys):
    min_rrev = []
    speed = []
    delay = []
    rrev = []
    for key in keys:
        trajec = dataset.trajecs[key]
        if trajec.min_rrev is not None:
            min_rrev.append(trajec.min_rrev)
            speed.append(trajec.speed[trajec.min_rrev_index])
            delay.append(trajec.min_rrev_delay)
            rrev.append(trajec.rrev[trajec.frame_at_deceleration])          
    
    speed = np.nan_to_num(speed)
    min_rrev = np.array(min_rrev)
    indices = np.where( speed>0.2)[0].tolist()
    dataset.min_rrev_vs_speed_fit = np.polyfit(speed[indices], min_rrev[indices], 1, rcond=None, full=False)
  
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    ax.plot(speed, min_rrev, 'o', color='red')
    x = np.linspace(np.min(speed), np.max(speed), 100)
    ymin = x*dataset.min_rrev_vs_speed_fit[0] + dataset.min_rrev_vs_speed_fit[1]
    ax.plot(x,ymin,color='red')
    
    ydec = x*dataset.rrev_vs_speed_fit[0] + dataset.rrev_vs_speed_fit[1] 
    ax.plot(x,ydec,color='blue')
    
    return min_rrev, speed, delay, rrev
        
        
def prep_dataset(dataset, distance=REQUIRED_DIST, do_classification=True):
    
    for k,trajec in dataset.trajecs.items():
        prep_trajectory(trajec, distance)
    calc_stats_at_deceleration(dataset)
            
            
def prep_trajectory(trajec, distance=REQUIRED_DIST, do_classification=True):
    print trajec.key
    calc_frame_of_landing(trajec)
    calc_radius_at_nearest(trajec)    
    calc_heading(trajec)
    
    #calc_heading_cumsum(trajec)
    trajec.frame_nearest_to_post = np.argmin(trajec.dist_to_stim_r_normed)
    frame_nearest_to_post = trajec.frame_nearest_to_post
    calc_saccades(trajec)
    
    if trajec.behavior == 'landing':
        trajec.frames = np.arange(get_frame_at_distance(trajec, distance), trajec.frame_of_landing+1).tolist()
    if trajec.behavior == 'flyby':
        try:
            frame_at_distance = get_frame_at_distance(trajec, distance, frames=np.arange(0,frame_nearest_to_post).tolist())
            sac_sgns = np.array(trajec.all_saccades) - frame_at_distance
            sac_neg_ind = np.where(sac_sgns<0)[0][-1]
            sac_neg = trajec.all_saccades[sac_neg_ind]
        except:
            sac_neg = 0
        print 'saccade: ', sac_neg
        first_frame = sac_neg + 1
        last_frame = np.min( [frame_nearest_to_post + 20, len(trajec.speed)-1] )
        trajec.frames_of_flyby = np.arange(first_frame, last_frame).tolist()

    calc_time_to_impact(trajec)
    ##calc_wv_ratio_cumsum(trajec)
    calc_dist_travelled(trajec)
    #calc_dist_travelled_in_window(trajec, window=0.1, window_units='sec')
    sa1.calc_post_dynamics_for_flydra_trajectory(trajec)
    trajec.expansion = sa1.diffa(trajec.angle_subtended_by_post)*trajec.fps
    if do_classification:
        classify(trajec, dfar=distance, dnear=0.005)
    
    if trajec.behavior == 'landing':
        calc_deceleration_initiation(trajec, plot=False)

def prep_movie_trajec(movieinfo):
    trajec = movieinfo.trajec
    trajec.key = movieinfo.id
    trajec.behavior = movieinfo.behavior
    calc_frame_of_landing(trajec)
        
    trajec.frames = np.arange(get_frame_at_distance(trajec, 0.08), trajec.frame_of_landing).tolist()
    prep_trajectory(trajec)
    
###
def normalize_dist_to_stim_r(trajec):
    if trajec.behavior == 'landing':
        trajec.dist_to_stim_r_normed = trajec.dist_to_stim_r - trajec.dist_to_stim_r[trajec.frame_of_landing]
    else:
        trajec.dist_to_stim_r_normed = trajec.dist_to_stim_r
def get_frame_at_distance(trajec, distance, singleframe=True, frames=None):
    normalize_dist_to_stim_r(trajec)
    if frames is None:
        frames = np.arange(0, trajec.frame_of_landing).tolist()
    dist_to_post = trajec.dist_to_stim_r_normed[frames]
    #print trajec.key, 'frames: ', frames
    dist_crossovers = np.where( sa1.diffa(np.sign(dist_to_post - distance)) != 0 )[0]
    
    if len(dist_crossovers) > 0:
        if singleframe:
            return dist_crossovers[-1]+frames[0]
        else:
            return dist_crossovers+frames[0]
    else:
        return frames[1]
        
def get_speed_at_distance(trajec, distance, singleframe=False):
    frames = get_frame_at_distance(trajec, distance, singleframe=singleframe)
    speed = np.max(trajec.speed[frames])
    return speed
    
def classify_dataset(dataset):
    for k,trajec in dataset.trajecs.items():
        classify(trajec)
    
def classify(trajec, dfar=REQUIRED_DIST, dnear=0.005):
    speed_hi_threshold = 0.18
    speed_lo_threshold = 0.18
    calc_dist_travelled(trajec)
    calc_deceleration_initiation(trajec, plot=False)
    
    '''
    if trajec.time_at_deceleration is not None:
        if trajec.speed_at_deceleration < 0.19:
            trajec.classification = 'slow'
        elif trajec.speed_at_deceleration < 0.5:
            trajec.classification = 'mid'
        else:
            trajec.classification = 'fast'
    else:
        trajec.classification = 'unkown'
    
    '''
    if trajec.behavior == 'landing' and trajec.dist_to_stim_r_normed[0] >= REQUIRED_DIST:
        trajec.speed_far = get_speed_at_distance(trajec, dfar, singleframe=True)
        trajec.speed_near = get_speed_at_distance(trajec, dnear, singleframe=True)
        
        change_in_heading = np.linalg.norm( trajec.heading_smooth[trajec.frames[0]] - trajec.heading_smooth[trajec.frame_of_landing] )
        projected_min_dist_to_post = dist_pt_to_line( (0,0), trajec.positions[trajec.frames[0],0:2], trajec.positions[trajec.frames[1],0:2] )
        
        #if trajec.speed_far < speed_lo_threshold:
        #    trajec.classification = 'slow'
        
            
        if len(trajec.saccades) == 0 and trajec.tortuosity < 1.1:
            trajec.classification = 'straight'
            
        #elif len(trajec.saccades) == 1:
        #    trajec.classification = 'single_saccade'
            
        else:
            trajec.classification = 'other'
            
    elif trajec.behavior == 'flyby':
        
        # get frame at 8cm away, prior to nearest approach
        frame_nearest_to_post = np.argmin(trajec.dist_to_stim_r)
        frames = np.arange(0, frame_nearest_to_post).tolist()
        frame_at_distance = get_frame_at_distance(trajec, 0.08, singleframe=True, frames=frames)
        #last_frame = np.min( [frame_nearest_to_post+20, len(trajec.speed)-1]) 
        #trajec.frames_of_flyby = np.arange(frame_at_distance, last_frame).tolist()
        
        if np.abs(trajec.angle_to_post[frame_at_distance]) < 5*np.pi/180. and trajec.tortuosity < 1.1: # and len(trajec.saccades)==1:  
            trajec.classification = 'straight'
        else:
            trajec.classification = 'flyby'
  
    else:
        
        trajec.speed_far = trajec.speed[0]
        
        change_in_heading = np.linalg.norm( trajec.heading_smooth[trajec.frames[0]] - trajec.heading_smooth[trajec.frame_of_landing] )
        projected_min_dist_to_post = dist_pt_to_line( (0,0), trajec.positions[trajec.frames[0],0:2], trajec.positions[trajec.frames[1],0:2] )
        
        if trajec.speed_far < speed_lo_threshold:
            trajec.classification = 'slow'
        
            
        elif len(trajec.saccades) == 0 and trajec.tortuosity < 1.1:
            trajec.classification = 'straight'
            
        elif len(trajec.saccades) == 1 and projected_min_dist_to_post > 0.03:
            trajec.classification = 'single_saccade'
        
        else:
            trajec.classification = 'mid'
            
def calc_frame_of_landing (trajec, threshold = 0.0005):
        # search forward in time until we get below threshold, then check to make sure fly stays below threshold for three more frames
        # if it's a flyby, find nearest point to post?
        if trajec.behavior == 'flyby':
            trajec.frame_of_landing = len(trajec.speed)-1
            return trajec.frame_of_landing
        if trajec.behavior == 'landing':
            diffdist = np.abs(sa1.diffa(trajec.dist_to_stim_r))
            #print 'calculating for landing'
            frame_of_landing = 0
            counter = 0
            for i in range(trajec.length):
                if frame_of_landing == 0:
                    if diffdist[i] < threshold and trajec.dist_to_stim_r[i] < 0.005:
                        frame_of_landing = i
                if frame_of_landing > 0:
                    if counter >= 3:
                        trajec.frame_of_landing = frame_of_landing
                        trajec.time_of_landing = trajec.fly_time[frame_of_landing]
                        #print 'frame of landing: ', frame_of_landing
                        return frame_of_landing
                    elif diffdist[i] < threshold:
                        counter = counter + 1
                    else:
                        counter = 0
                        frame_of_landing = 0
                        
            trajec.frame_of_landing = len(trajec.speed)-1
            #print 'frame of landing: ', frame_of_landing
            return len(trajec.speed)-1
        else:
            #print trajec.behavior
            trajec.frame_of_landing = len(trajec.speed)-1
            return len(trajec.speed)-1
            
def get_classified_keys(dataset):
    classification_dict = {}    
    for k,trajec in dataset.trajecs.items():
        if trajec.classification not in classification_dict.keys():
            classification_dict.setdefault( trajec.classification, [] )
        classification_dict[trajec.classification].append(k)
    return classification_dict
    
    
def calc_heading(trajec):

    ## kalman
    data = trajec.speed.reshape([len(trajec.speed),1])
    ss = 2 # state size
    os = 1 # observation size
    F = np.array([   [1,1], # process update
                     [0,1]],
                    dtype=np.float)
    H = np.array([   [0,1]], # observation matrix
                    dtype=np.float)
    Q = 1*np.eye(ss) # process noise
    R = 25*np.eye(os) # observation noise
    init_vel = data[0]
    initx = np.array([0, init_vel], dtype=np.float)
    initv = 0*np.eye(ss)
    xsmooth,Vsmooth = sa1.kalman_smoother(data, F, H, Q, R, initx, initv, plot=False)
    ## 
    
    #plt.plot(xsmooth[:,1])
    
    trajec.smooth_accel = sa1.diffa(xsmooth[:,1])
    
    #plt.plot(trajec.smooth_accel*20.)
    #plt.plot(trajec.dist_to_stim_r)
    
    
    trajec.heading = sa1.remove_angular_rollover(np.arctan2(trajec.velocities[:,1], trajec.velocities[:,0]), 3)
    ## kalman
    data = trajec.heading.reshape([len(trajec.heading),1])
    ss = 3 # state size
    os = 1 # observation size
    F = np.array([   [1,1,0], # process update
                     [0,1,1],
                     [0,0,1]],
                    dtype=np.float)
    H = np.array([   [1,0,0]], # observation matrix
                    dtype=np.float)
    Q = np.eye(ss) # process noise
    Q[0,0] = .01
    Q[1,1] = .01
    Q[2,2] = .01
    R = 1*np.eye(os) # observation noise
    initx = np.array([data[0,0], data[1,0]-data[0,0], 0], dtype=np.float)
    initv = 0*np.eye(ss)
    xsmooth,Vsmooth = sa1.kalman_smoother(data, F, H, Q, R, initx, initv, plot=False)
    ## 
    trajec.heading_smooth = xsmooth[:,0]
    trajec.heading_smooth_diff = xsmooth[:,1]
    trajec.heading_smooth_diff2 = xsmooth[:,2]
    
def calc_heading_cumsum(trajec, initial_dist=REQUIRED_DIST, plot=False):
    initial_frame = get_frame_at_distance(trajec, initial_dist)
    shifted_angles = np.abs(trajec.heading[initial_frame:trajec.frame_of_landing] - trajec.heading[initial_frame])
    cumsum_angles = cumsum(shifted_angles)
    trajec.heading_cumsum = cumsum_angles
    
    trajec.heading_cumsum_diff = sa1.diffa(trajec.heading_cumsum)
    trajec.heading_cumsum_diff2 = sa1.diffa(trajec.heading_cumsum_diff)
    trajec.heading_cumsum_pts_of_inflection = np.where(   np.diff( np.sign(trajec.heading_cumsum_diff2) ) != 0 )[0].tolist()
    
    
    trajec.heading_cumsum_pos_pts_of_inflection = []
    trajec.heading_cumsum_neg_pts_of_inflection = []
    
    window = 10
    for pt in trajec.heading_cumsum_pts_of_inflection:
        if pt-window >= 0 and pt+window < len(trajec.heading_cumsum):
            change_diff = trajec.heading_cumsum_diff[pt+window] - trajec.heading_cumsum_diff[pt-window]
            #print change_diff
            if np.abs(change_diff > .2):
                    trajec.heading_cumsum_pos_pts_of_inflection.append(pt)
                    
    
    if plot:
        plt.plot(trajec.heading_cumsum)
        plt.plot(trajec.heading_cumsum_pos_pts_of_inflection, trajec.heading_cumsum[trajec.heading_cumsum_pos_pts_of_inflection], 'o')
        plt.plot(trajec.heading_cumsum_diff)
        
        
def calc_time_to_impact(trajec):
    trajec.time_to_impact = trajec.dist_to_stim_r_normed / trajec.speed
    
def calc_dist_travelled(trajec):
    if trajec.behavior == 'landing':
        vt = sa1.norm_array(trajec.positions[trajec.frames, 0:2] - trajec.positions[ trajec.frames[0]-1:trajec.frames[-1], 0:2])  
        trajec.dist_travelled = cumsum(vt)
        vt -= vt[0]
        euclidean_dist_travelled = np.linalg.norm( trajec.positions[trajec.frames[0], 0:2] - trajec.positions[trajec.frames[-1], 0:2] )
        trajec.tortuosity = trajec.dist_travelled[-1] / euclidean_dist_travelled
        trajec.mean_speed = cumsum(trajec.speed[trajec.frames]) / cumsum(np.ones_like(trajec.speed[trajec.frames]))
    if trajec.behavior == 'flyby':
        frame_nearest_post = np.argmin(trajec.dist_to_stim_r)
        if len(trajec.saccades) > 0 and trajec.saccades[-1] > trajec.frames_of_flyby[0]+10:
            frames = np.arange(trajec.frames_of_flyby[0],trajec.saccades[-1]+1).tolist()
        else:
            frames = np.arange(trajec.frames_of_flyby[0],frame_nearest_post).tolist()
        print 'herehere'
        print trajec.frames_of_flyby
        print trajec.saccades
        print frames[0], frames[-1]
        vt = sa1.norm_array(trajec.positions[frames, 0:2] - trajec.positions[ frames[0]-1:frames[-1], 0:2])  
        trajec.dist_travelled = cumsum(vt)
        vt -= vt[0]
        euclidean_dist_travelled = np.linalg.norm( trajec.positions[frames[0], 0:2] - trajec.positions[frames[-1], 0:2] )
        trajec.tortuosity = trajec.dist_travelled[-1] / euclidean_dist_travelled
        trajec.mean_speed = cumsum(trajec.speed[frames]) / cumsum(np.ones_like(trajec.speed[frames]))
        
def calc_tortuosity_for_frame_range(trajec, framerange):
    frames = np.arange(framerange[0], framerange[-1]).tolist()
    vt = sa1.norm_array(trajec.positions[frames, 0:2] - trajec.positions[ frames[0]-1:frames[-1], 0:2])  
    dist_travelled = cumsum(vt)
    vt -= vt[0]
    euclidean_dist_travelled = np.linalg.norm( trajec.positions[frames[0], 0:2] - trajec.positions[frames[-1], 0:2] )
    tortuosity = dist_travelled[-1] / euclidean_dist_travelled
    return tortuosity
def calc_tortuosity_btwn_saccades(trajec, s1, s2):
    f1 = get_saccade_range(trajec, s1)[-1]
    f2 = get_saccade_range(trajec, s2)[0]
    return calc_tortuosity_for_frame_range(trajec, [f1,f2])
    
def calc_saccades(trajec, magnitude=351):
    raw_saccades = find_extrema_points(trajec.heading_smooth_diff)
    
        
    
    #windowed_heading_diff = diff_windowed(trajec.heading, 5)
    
    angular_vel = np.abs( trajec.heading_smooth_diff )*100.*180/np.pi
    
    trajec.saccades = []
    trajec.all_saccades = []
    for saccade in raw_saccades:
        if trajec.speed[saccade-10] > 0.01:
            #if np.abs(windowed_heading_diff[saccade]) > magnitude:
            
            if angular_vel[saccade] > magnitude:
                if trajec.dist_to_stim_r_normed[saccade] > 0.005:
                
                    sac_range = get_saccade_range(trajec, saccade)
                    if len(sac_range) < 3:
                        continue
                    for s in trajec.all_saccades:
                        s_range = get_saccade_range(trajec, s)
                        if saccade in s_range:
                            continue
                                    
                    if trajec.behavior == 'flyby':
                        if saccade > trajec.frames_of_flyby[0] and saccade < trajec.frame_nearest_to_post:
                            trajec.saccades.append(saccade)
                    else:
                        if saccade in trajec.frames:
                            trajec.saccades.append(saccade)
                    trajec.all_saccades.append(saccade)
                            
def calc_wv_ratio_cumsum(trajec, plot=False): 

    ratio = trajec.heading_smooth[trajec.frames] / trajec.speed[trajec.frames]
    shifted_ratio = np.abs(ratio - ratio[0])
    cumsum_ratio = cumsum(shifted_ratio)
    trajec.wv_ratio = cumsum_ratio
    
    if plot:
        plt.plot(trajec.time_to_impact[trajec.frames]-trajec.time_to_impact[trajec.frames][-1], trajec.wv_ratio)
        for saccade in trajec.saccades:
            if saccade in trajec.frames:
                #print saccade, trajec.dist_to_stim_r_normed[saccade-trajec.frames[0]], trajec.wv_ratio[saccade-trajec.frames[0]]
                plt.plot(trajec.time_to_impact[saccade]-trajec.time_to_impact[trajec.frames][-1], trajec.wv_ratio[saccade-trajec.frames[0]], 'o')
        plt.show()
        
        
def plot_wv_ratio_trajecs(dataset, keys):

    for key in keys:
        trajec = dataset.trajecs[key]
        plt.plot(trajec.time_to_impact[trajec.frames]-trajec.time_to_impact[trajec.frames][-1], trajec.wv_ratio)
        for saccade in trajec.saccades:
            if saccade in trajec.frames:
                #print saccade, trajec.dist_to_stim_r_normed[saccade-trajec.frames[0]], trajec.wv_ratio[saccade-trajec.frames[0]]
                plt.plot(trajec.time_to_impact[saccade]-trajec.time_to_impact[trajec.frames][-1], trajec.wv_ratio[saccade-trajec.frames[0]], 'o')
    plt.show()
    
def get_first_saccade(trajec):
    for saccade in trajec.saccades:
        if saccade > trajec.frames[0]:
            return saccade
            
            
def calc_dist_travelled_in_window(trajec, window=0.1, window_units='sec'):

    # give window in terms of full window, calculated with half on each side

    if window_units == 'sec':
        window = trajec.fps*window
    else:
        window = window
    window /= 2.
    window = int(window)

    frames = np.arange(window, trajec.frame_of_landing-window)
    trajec.dist_travelled_windowed = np.zeros_like(trajec.speed)
    for i, frame in enumerate(frames):
        frames_in_window = np.arange(frame-window, frame+window).tolist()
        vt = trajec.speed[frames_in_window] / trajec.fps
        straight_line_dist = np.linalg.norm( trajec.positions[frame-window,0:2] - trajec.positions[frame+window,0:2] )
        trajec.dist_travelled_windowed[i] = cumsum(vt)[0] / straight_line_dist
    
    trajec.dist_travelled_windowed[0:frames[0]] = trajec.dist_travelled_windowed[frames[0]]*np.ones_like( trajec.dist_travelled_windowed[0:window] )
    trajec.dist_travelled_windowed[frames[-1]:trajec.frame_of_landing] = trajec.dist_travelled_windowed[frames[-1]]*np.ones_like( trajec.dist_travelled_windowed[frames[-1]:trajec.frame_of_landing] )
    
################################################################################################
## PLOTTING
################################################################################################

###
def pdf_flydra_trajecs_of_classifications(dataset, filename='classified_flydra_trajectories', pages=None, classification=None):
    
    if classification is None:
        classified_keys = get_classified_keys(dataset)
    else:
        classified_keys = get_classified_keys(dataset)[classification]
        classified_keys = {classification: classified_keys}
    
    if pages is None:
        pages = np.inf
    
    page = -1
    plt.close('all')
    pp =  pdf.PdfPages(filename)
    
    for classification, keys in classified_keys.items():
        numfliestoplot = 5
        firstfly = 0
        while firstfly < len(keys):
            page += 1
            if page > pages:
                break
            cl = ap.xy_trajectories(dataset, trajectory=keys, show_saccades=True, trajec_alpha=0.8, firstfly=firstfly, numfliestoplot=numfliestoplot, print_obj_ids=True, frames_to_plot='flyby_slice')
            title = str(classification) + ' course'
            cl.ax0.set_title(title)
            cl.ax0.figure.set_size_inches(2*10,1*10)
            cl.ax0.figure.set_dpi(72)
            pp.savefig()
            plt.close('all')
            
            firstfly += numfliestoplot
    
    pp.close()
    print 'closed'

###
def plot_classification(dataset, behavior='landing'):
    if type(behavior) is not list:
        behavior = [behavior]
    for k,trajec in dataset.trajecs.items():
        if trajec.behavior in behavior:


            if trajec.behavior == 'landing' and trajec.dist_to_stim_r[0] >= REQUIRED_DIST:
                
                if trajec.classification == 'fast':
                    color = 'red'
                elif trajec.classification == 'slow':
                    color = 'blue'
                elif trajec.classification == 'mid':
                    color = 'black'
                else:
                    color = 'green'
                    
                plt.plot(trajec.speed_far, trajec.speed_near, 'o', color=color, )
                
        plt.xlabel('speed far (0.06m)')
        plt.ylabel('speed near (o.01m)')
        
    
    
    
    
    
    
######################## MISC #######################


def pdf_plot_saccade_histogram(dataset, filename='saccade_histogram'):
    
    classified_keys = get_classified_keys(dataset)
    
    page = -1
    plt.close('all')
    pp =  pdf.PdfPages(filename)
    
    for classification, keys in classified_keys.items():
        
        time_to_impact_of_saccades = []        
        
        for key in keys:
            trajec = dataset.trajecs[key]
            
            for saccade in trajec.saccades:
                time_to_impact_of_saccades.append( trajec.time_to_impact[saccade] )
            
            
        plt.hist(time_to_impact_of_saccades, bins=40)            
    
        title = classification
        plt.title(title)
        #plt.fig.set_size_inches(2*10,1*10)
        #ax.figure.set_dpi(72)
        pp.savefig()
        plt.close('all')
    
    pp.close()
    print 'closed'
    
    
    
    
    
    
    
def convert_stringtime_to_epochtime(stringtime, dst=1):
    # stringtime of form '20101103_102345'
    
    year = int(stringtime[0:4])
    month = int(stringtime[4:6])
    day = int(stringtime[6:8])
    hour = int(stringtime[9:11])
    minute = int(stringtime[11:13])
    second = int(stringtime[13:15])
    
    dt = datetime.datetime(year,month,day,hour,minute,second)
    cal = dt.isocalendar()
    dayofyear = cal[1]*7+cal[2]
    dayofweek = cal[2]
    
    struct_time = time.struct_time((year,month,day,hour,minute,second,dayofyear,dayofweek,dst))
    epochtime = time.mktime(struct_time)
    return epochtime

def get_post_type_for_dataset(dataset):
    
    epochtimes = []
    
    post_types_black = ['black' for i in range(8)]
    post_types_black_angled = ['black_angled' for i in range(7)]
    post_types_checkered = ['checkered' for i in range(2)]
    post_types_checkered_angled = ['checkered_angled' for i in range(13)]
    post_types_none = ['none' for i in range(2)]
    post_types = post_types_black + post_types_black_angled + post_types_checkered + post_types_checkered_angled + post_types_none
    
    print len(post_types)
    
    '''
    for f in dataset.filename:
        stringtime = f[4:19]
        epochtime = convert_stringtime_to_epochtime(stringtime)
        epochtimes.append(epochtime)
    epochtimes = np.array(epochtimes)
    '''
    
    for key, trajec in dataset.trajecs.items():
        trajec.key = key
        '''
        time_diff = trajec.epoch_time[0]*np.ones_like(epochtimes) - epochtimes
        time_diff[time_diff<0] = np.inf
        f = np.argmin(time_diff)
        trajec.post_type = post_types[ f ]
        '''
        try:
            s = int(trajec.key.split('_')[0])
            trajec.post_type = post_types[ s-1 ]
            
            calc_radius_at_nearest(trajec)
        except:
            print key
            
def print_key_prefix(dataset, post_type):
    
    for key in dataset.trajecs.keys():
        if dataset.trajecs[key].post_type == post_type:
            print dataset.trajecs[key].key[0:3]

    
def print_posttype_for_classification(dataset, c):
    classified_keys = get_classified_keys(dataset)
    keys = classified_keys[c]
    for key in keys:
        print dataset.trajecs[key].post_type
        

    
def calc_deceleration_initiation(trajec, plot=False):
    trajec.frame_at_deceleration = None

    x = trajec.dist_to_stim_r_normed[0:trajec.frame_of_landing]
    #a = trajec.accel_1d[0:trajec.frame_of_landing]
    t = trajec.time_to_impact[0:trajec.frame_of_landing]
    s = trajec.speed[0:trajec.frame_of_landing]
    ang = trajec.angle_subtended_by_post[0:trajec.frame_of_landing]
    
    a = sa1.diffa(s) / np.abs(sa1.diffa(ang))
    '''
    steady_state_acceleration = a[np.where( (t<0.4)*(t>0.2) )[0].tolist()]
    std = np.std(steady_state_acceleration)
    mean = np.mean(steady_state_acceleration)
    
    # find where acceleration starts to go down outside of normal flight
    a_past_2std_pts = np.where( (a < mean-2*std)*(sa1.diffa(a)<0) )[0].tolist()
    '''
    try:
        if trajec.behavior == 'landing':
            f = np.where( sa1.diffa( np.sign(a) ) < 0 )[0][-1]
            trajec.frame_at_deceleration = f
            frames = np.arange(trajec.frame_at_deceleration-1,trajec.frame_at_deceleration+1).tolist()
            # interpolate to get best estimate of time
            trajec.time_at_deceleration = np.interp(0, a[frames], trajec.epoch_time[frames]) 
            trajec.angle_at_deceleration = np.interp(trajec.time_at_deceleration, trajec.epoch_time[frames], trajec.angle_subtended_by_post[frames])
            trajec.expansion_at_deceleration = np.interp(trajec.time_at_deceleration, trajec.epoch_time[frames], trajec.expansion[frames])
            trajec.speed_at_deceleration = np.interp(trajec.time_at_deceleration, trajec.epoch_time[frames], trajec.speed[frames])
            trajec.dist_at_deceleration = np.interp(trajec.time_at_deceleration, trajec.epoch_time[frames], trajec.dist_to_stim_r[frames])
        elif trajec.behavior == 'flyby':
            
            if len(trajec.saccades) > 0:
                sf = get_saccade_range(trajec, trajec.saccades[-1])
                #frame_of_interest = np.min([trajec.frame_nearest_to_post, sf[0]])
                frame_of_interest = trajec.frame_nearest_to_post
            else:
                frame_of_interest = trajec.frame_nearest_to_post
            indices = np.arange(trajec.frames_of_flyby[0],frame_of_interest).tolist()
            print 'INDICES: ', indices
            f = np.where( sa1.diffa( np.sign(a[indices]) ) < 0 )[0][-1]
            print 'DECEL: ', f
            trajec.frame_at_deceleration = f + indices[0]
            frames = np.arange(trajec.frame_at_deceleration-1,trajec.frame_at_deceleration+1).tolist()
            # interpolate to get best estimate of time
            trajec.time_at_deceleration = np.interp(0, a[frames], trajec.epoch_time[frames]) 
            trajec.angle_at_deceleration = np.interp(trajec.time_at_deceleration, trajec.epoch_time[frames], trajec.angle_subtended_by_post[frames])
            trajec.expansion_at_deceleration = np.interp(trajec.time_at_deceleration, trajec.epoch_time[frames], trajec.expansion[frames])
            trajec.speed_at_deceleration = np.interp(trajec.time_at_deceleration, trajec.epoch_time[frames], trajec.speed[frames])
            trajec.dist_at_deceleration = np.interp(trajec.time_at_deceleration, trajec.epoch_time[frames], trajec.dist_to_stim_r[frames])
            print 'decel calculated for flyby'            
            
        if plot:
            #plt.plot(trajec.epoch_time[0:trajec.frame_of_landing], a)
            #plt.plot(trajec.epoch_time[0:trajec.frame_of_landing], trajec.speed[0:trajec.frame_of_landing])
            #plt.plot(trajec.time_at_deceleration, 0, 'o', color='red')
            #plt.plot(trajec.epoch_time[trajec.frame_at_deceleration], a[trajec.frame_at_deceleration], 'o', color='blue')
            plt.plot(trajec.angle_subtended_by_post[trajec.frames_of_flyby], trajec.speed[trajec.frames_of_flyby])
            plt.plot(trajec.angle_subtended_by_post[trajec.frame_at_deceleration], trajec.speed[trajec.frame_at_deceleration], 'o')
            plt.show()
        #print trajec.key, ' processed', trajec.frame_at_deceleration
    #if 0:
    except:
        print 'exception'
        pt_at_deceleration = None
        trajec.frame_at_deceleration = None
        trajec.time_at_deceleration = None
        trajec.angle_at_deceleration = None
        trajec.expansion_at_deceleration = None
        trajec.speed_at_deceleration = None
        trajec.dist_at_deceleration = None
        
        
        #print trajec.key, ' unprocessed'
        
def plot_deceleration_vs_time_to_impact(dataset, keys, show_legs=True):
    figure = None
    ax0_size=[0.1,0.1,0.7,0.7]
    norm=[0,0.6]
    trajec_alpha = 0.9
    
    cl = colorline.Colorline(xlim=[-.4,0.1], ylim =[-3,3], norm=norm, colormap = 'jet', figure=figure, hide_colorbar=False, ax0_size=ax0_size)
    
    for key in keys:
        trajec = dataset.trajecs[key]
        #calc_deceleration_initiation(trajec)
        if trajec.time_at_deceleration is not None:
             
            a = trajec.accel_1d
            t = trajec.epoch_time - trajec.time_at_deceleration
            s = trajec.speed
            cl.colorline(t[trajec.frames], a[trajec.frames], s[trajec.frames],linewidth=1, norm=norm, alpha=trajec_alpha)
            
            cl.ax0.plot(t[trajec.frame_of_landing], a[trajec.frame_of_landing], '.', color='red')
            
            if show_legs:
                acc_at_legext = np.interp(trajec.legextension_time, trajec.epoch_time, trajec.accel_1d)
                timpact_at_legext = np.interp(trajec.legextension_time, trajec.epoch_time[0:trajec.frame_of_landing], t)
                cl.ax0.plot( timpact_at_legext, acc_at_legext, 'o', color='red')
    cl.ax0.set_xlabel('time, sec')
    cl.ax0.set_ylabel('acceleration, m/s^2')
    cl.ax1.set_ylabel('speed, m/s')
            
def plot_exp_vs_angle(dataset, keys=None, show_legs=False, time_offset=0):
    # point of deceleration initiation figure
    if keys is None:
        keys = dataset.trajecs.keys()
        
    fig = plt.figure()
    ax = fig.add_subplot(111)
    colormap_norm = matplotlib.colors.Normalize(0, .6, clip=True)
    cmap = plt.get_cmap('jet')
    
    for key in keys:
        trajec = dataset.trajecs[key]
        calc_deceleration_initiation(trajec)
        if trajec.frame_at_deceleration is not None:
            frame_offset = int(time_offset / trajec.fps)

            c = cmap(colormap_norm(trajec.speed[trajec.frame_at_deceleration+frame_offset]))
            
            frames = np.arange(trajec.frame_at_deceleration-1+frame_offset,trajec.frame_at_deceleration+1+frame_offset).tolist()
            angle_at_deceleration = np.interp(trajec.time_at_deceleration+time_offset, trajec.epoch_time[frames], trajec.angle_subtended_by_post[frames])
            expansion_at_deceleration = np.interp(trajec.time_at_deceleration+time_offset, trajec.epoch_time[frames], trajec.expansion[frames])
            
            speed_at_deceleration = np.interp(trajec.time_at_deceleration+time_offset, trajec.epoch_time[frames], trajec.speed[frames])
            dist_at_deceleration = np.interp(trajec.time_at_deceleration+time_offset, trajec.epoch_time[frames], trajec.dist_to_stim_r[frames])
            
            ax.plot( angle_at_deceleration, expansion_at_deceleration, 'o', color=c)
            #ax.plot( dist_at_deceleration, speed_at_deceleration, 'o', color=c)
            
            if show_legs:
                if trajec.legextension_time is not None:    
                    angle_at_legext = np.interp(trajec.legextension_time, trajec.epoch_time, trajec.angle_subtended_by_post)
                    expansion_at_legext = np.interp(trajec.legextension_time, trajec.epoch_time, trajec.expansion)
                    
                    speed_at_legext = np.interp(trajec.legextension_time, trajec.epoch_time, trajec.speed)
                    dist_at_legext = np.interp(trajec.legextension_time, trajec.epoch_time, trajec.dist_to_stim_r)
                    
                    ax.plot( angle_at_legext, expansion_at_legext, 'o', color='red')
                    ax.plot( [angle_at_legext, angle_at_deceleration], [expansion_at_legext, expansion_at_deceleration], color='red')
            
    ax.set_ylabel('expansion at deceleration, rad/s')
    ax.set_xlabel('angle at deceleration, rad')
    ax.set_title('retinal expansion vs angular size at initiation of deceleration')
    
    colorbar_pad = 0
    colorbar_size = "3%"
    divider = make_axes_locatable(ax)
    divider.set_anchor('E')
    cax = divider.append_axes("right", size=colorbar_size, pad=colorbar_pad)
    cb = matplotlib.colorbar.ColorbarBase(cax, cmap=cmap, norm=colormap_norm, orientation='vertical', boundaries=None)
    cax.set_ylabel('speed, m/s')
    
def plot_exp_vs_angle_for_saccades(dataset, keys=None, show_legs=False, time_offset=0):
    # point of deceleration initiation figure
    expansion = []
    rrev = []
    speed = []
    angle = []
    
    if keys is None:
        keys = dataset.trajecs.keys()
        
    fig = plt.figure()
    ax = fig.add_subplot(111)
    colormap_norm = matplotlib.colors.Normalize(0, .6, clip=True)
    cmap = plt.get_cmap('jet')
    
    for key in keys:
        trajec = dataset.trajecs[key]
        calc_deceleration_initiation(trajec)
        if trajec.frame_at_deceleration is not None:
            frame_offset = int(time_offset / trajec.fps)
            frame_of_interest = trajec.saccades[-1] 
            time_of_interest = trajec.epoch_time[frame_of_interest] # use interpolated values for better plots
            
            c = cmap(colormap_norm(trajec.speed[frame_of_interest+frame_offset]))
            
            frames = np.arange(frame_of_interest-1+frame_offset,frame_of_interest+1+frame_offset).tolist()
            angle_at_poi = trajec.angle_subtended_by_post[frame_of_interest]
            expansion_at_poi = trajec.expansion[frame_of_interest]
            speed_at_poi = trajec.speed[frame_of_interest]
            dist_at_poi = trajec.dist_to_stim_r[frame_of_interest]
            rrev_at_poi = expansion_at_poi / angle_at_poi
            
            ax.plot( angle_at_poi, expansion_at_poi, 'o', color=c)
            #ax.plot( dist_at_deceleration, speed_at_deceleration, 'o', color=c)

            expansion.append(expansion_at_poi)
            rrev.append(rrev_at_poi)
            speed.append(speed_at_poi)
            angle.append(angle_at_poi)
            
    ax.set_ylabel('expansion at deceleration, rad/s')
    ax.set_xlabel('angle at deceleration, rad')
    ax.set_title('retinal expansion vs angular size at initiation of deceleration')
    
    colorbar_pad = 0
    colorbar_size = "3%"
    divider = make_axes_locatable(ax)
    divider.set_anchor('E')
    cax = divider.append_axes("right", size=colorbar_size, pad=colorbar_pad)
    cb = matplotlib.colorbar.ColorbarBase(cax, cmap=cmap, norm=colormap_norm, orientation='vertical', boundaries=None)
    cax.set_ylabel('speed, m/s')
    
    return np.array(expansion), np.array(rrev), np.array(speed), np.array(angle)
    
    
    
        
def calc_stats_at_deceleration(dataset, keys=None, time_offset=0, return_vals=False):
    if keys is None:
        keys = dataset.trajecs.keys()
    
    expansion = []
    speed = []
    angle = []
    distance = []
    trajectories = []
    radius = []
    time_to_impact = []
    for key in keys:
        trajec = dataset.trajecs[key]
        if trajec.frame_at_deceleration is not None:
            print trajec.frame_of_landing - trajec.frame_at_deceleration 
            frame_offset = int(time_offset / trajec.fps)
            frames = np.arange(trajec.frame_at_deceleration-1+frame_offset,trajec.frame_at_deceleration+1+frame_offset).tolist()
            angle_at_deceleration = np.interp(trajec.time_at_deceleration+time_offset, trajec.epoch_time[frames], trajec.angle_subtended_by_post[frames])
            expansion_at_deceleration = np.interp(trajec.time_at_deceleration+time_offset, trajec.epoch_time[frames], trajec.expansion[frames])
            speed_at_deceleration = np.interp(trajec.time_at_deceleration+time_offset, trajec.epoch_time[frames], trajec.speed[frames])
            
            distance_at_deceleration = np.interp(trajec.time_at_deceleration+time_offset, trajec.epoch_time[frames], trajec.dist_to_stim_r_normed[frames])
            v_over_d_at_deceleration = speed_at_deceleration / distance_at_deceleration
            
            radius_at_deceleration = trajec.radius_at_nearest
            
            if np.sum(np.isnan(expansion_at_deceleration)) > 0:
                continue
            else:
                expansion.append(expansion_at_deceleration)
                speed.append(speed_at_deceleration)
                angle.append(angle_at_deceleration)
                distance.append(distance_at_deceleration)
                trajectories.append(key)
                radius.append(radius_at_deceleration)
                time_to_impact.append( distance_at_deceleration / speed_at_deceleration )
        else:
            print 'no deceleration'
    expansion = np.array(expansion)
    speed = np.array(speed)
    angle = np.array(angle)
    distance = np.array(distance)
    trajectories = np.array(trajectories)
    radius = np.array(radius)
    time_to_impact = np.array(time_to_impact)
        
    dataset.speed_at_deceleration = speed
    dataset.angle_at_deceleration = angle
    dataset.expansion_at_deceleration = expansion
    dataset.trajectories_at_deceleration = trajectories
    dataset.distance_at_deceleration = distance
    dataset.radius_at_deceleration = radius
    dataset.time_to_impact_at_deceleration = time_to_impact
        
    if return_vals:
        return expansion, speed, angle, distance, trajectories 
    else:
        return
        
def calc_stats_at_saccade(dataset, keys=None, time_offset=0, return_vals=False):
    if keys is None:
        classified_keys = get_classified_keys(dataset)
        keys = classified_keys['straight']
    
    expansion = []
    speed = []
    angle = []
    distance = []
    trajectories = []
    radius = []
    time_to_impact = []
    for key in keys:
        trajec = dataset.trajecs[key]
        
        try: 
            s = trajec.saccades[-1]
            docalc = True
        except:
            docalc = False
            
        if docalc:
            angle_at_saccade = trajec.angle_subtended_by_post[s]
            expansion_at_saccade = trajec.expansion[s]
            speed_at_saccade = trajec.speed[s]
            distance_at_saccade = trajec.dist_to_stim_r_normed[s]
            radius_at_saccade = trajec.radius_at_nearest
            
            if np.sum(np.isnan(expansion_at_saccade)) > 0:
                continue
            else:
                expansion.append(expansion_at_saccade)
                speed.append(speed_at_saccade)
                angle.append(angle_at_saccade)
                distance.append(distance_at_saccade)
                trajectories.append(key)
                radius.append(radius_at_saccade)
                time_to_impact.append( distance_at_saccade / speed_at_saccade )
        else:
            print 'no saccade'
            
    expansion = np.array(expansion)
    speed = np.array(speed)
    angle = np.array(angle)
    distance = np.array(distance)
    trajectories = np.array(trajectories)
    radius = np.array(radius)
    time_to_impact = np.array(time_to_impact)
        
    dataset.speed_at_saccade = speed
    dataset.angle_at_saccade = angle
    dataset.expansion_at_saccade = expansion
    dataset.trajectories_at_saccade = trajectories
    dataset.distance_at_saccade = distance
    dataset.radius_at_saccade = radius
    dataset.time_to_impact_at_saccade = time_to_impact
        
    if return_vals:
        return expansion, speed, angle, distance, trajectories 
    else:
        return
        
        
    
def plot_exp_vs_angle_for_speed_range(dataset, keys, velrange=[0.4,0.8], time_offset=0, plot=False):
    expansion, rrev, speed, angle, distance = get_rrev_and_speed(dataset, keys, time_offset=0, plot=False)
    indices = np.where( (speed>velrange[0])*(speed<velrange[1]) )[0].tolist()
    
    result = np.polyfit(angle[indices], expansion[indices], 1, rcond=None, full=False)
    
    xpts = np.linspace(np.min(angle[indices]), np.max(angle[indices]), 100)
    ypts = result[0]*xpts + result[1]
    
    
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(angle[indices], expansion[indices], 'o')
    ax.plot(xpts, ypts)
    ax.set_ylabel('expansion, rad/s')
    ax.set_xlabel('angle, rad')
    plt.show()
    
def get_rrev_vs_speed(dataset, keys, time_offset=0, plot=False):
    min_num_data_pts = 1
    expansion, rrev, speed, angle, distance = get_rrev_and_speed(dataset, keys, time_offset=0)
    speed = np.nan_to_num(speed)
    
    if plot:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        fig2 = plt.figure()
        ax2 = fig2.add_subplot(111)
        fig3 = plt.figure()
        ax3 = fig3.add_subplot(111)
        
        colormap_norm = matplotlib.colors.Normalize(0, .6, clip=True)
        cmap = plt.get_cmap('jet')
    
    # fit expansion vs angle according to speed: (and fix (0,0) as a point)
    sorted_speed = np.argsort(speed)
    window = 10
    fit_array = np.zeros(len(sorted_speed)-window)
    fit_speed = np.zeros(len(sorted_speed)-window)
    for i in range( len(sorted_speed)-window ):
        indices = sorted_speed[i:i+window]
        angle_zeros = np.zeros_like(angle)
        angle_ext = np.hstack( (angle[indices], angle_zeros) )
        expansion_zeros = np.zeros_like(expansion)
        expansion_ext = np.hstack( (expansion[indices], expansion_zeros) )
        
        fit = np.polyfit(angle_ext**2, expansion_ext, 1, rcond=None, full=False) # only fit square: expansion = K*angle^2
        fit = fit[0]
        fit_array[i] = fit
        fit_speed[i] = np.mean(speed[indices])
        
        if plot:        
            x = (np.linspace(0, np.max(angle)))**2
            y = fit*x
            c = cmap(colormap_norm( np.mean(speed[indices]) ))
            ax2.plot( x,y, color=c, alpha=0.2)
            
    dataset.expansion_vs_angle_fit_raw = fit_array
    
    # now fit K vs speed
    dataset.rrev_vs_speed_fit = np.polyfit(fit_speed, fit_array, 1, rcond=None, full=False)[0]
    
    if plot:
        x = np.linspace(0, np.max(speed))
        y = dataset.rrev_vs_speed_fit*x
        ax.plot(x,y,color='blue')
        ax.plot(fit_speed, fit_array, 'o', color='gray')
        ax.set_xlabel('speed, m/s')
        ax.set_ylabel('RREV quadratic fit')
        
        
    if plot:
        for i in range(len(angle)):
            c = cmap(colormap_norm(speed[i]))
            ax2.plot( (angle[i])**2, expansion[i], 'o', color=c)
            ax2.set_xlabel('retinal size squared, deg')
            ax2.set_ylabel('retinal expansion velocity, deg/s')
        
    plt.show()
    return
    

def calc_rrev(dataset, keys=None):
    if keys is None:
        keys = dataset.trajecs.keys()
    for key in keys:
        trajec = dataset.trajecs[key]
        coeff = dataset.rrev_vs_speed_fit*trajec.speed
        trajec.rrev = trajec.expansion / (trajec.angle_subtended_by_post**2)
    

def calc_coeff_of_var_rrev(dataset, keys, variable='deceleration'):
    fig2 = plt.figure()
    ax2 = fig2.add_subplot(111)
    
    time = np.linspace(-.5, .15, 100)
    
    for t in time:
        
        rrev_list = []
        for key in keys:
            trajec = dataset.trajecs[key] 
            if trajec.frame_at_deceleration is not None:
                if variable == 'deceleration':
                    time_at_poi = trajec.time_at_deceleration
                elif variable == 'saccade':
                    try:
                        time_at_poi = trajec.epoch_time[trajec.saccades[-1]]
                    except:
                        time_at_poi = None
                if time_at_poi is not None:
                    rrev_at_poi = np.interp(t, trajec.epoch_time-time_at_poi, trajec.rrev)
                    rrev_list.append(rrev_at_poi)
                
        coeff_var = np.std(rrev_list) / np.abs(np.mean(rrev_list))
            
        ax2.plot(t, np.mean(rrev_list), '.', color='red')
        ax2.plot(t, coeff_var*20, '.', color='blue')
    ax2.set_xlabel('time, s')
    ax2.set_ylabel('cv = std/mean')
        
        
def plot_dist_travelled_histogram(dataset):
    
    classified_keys = get_classified_keys(dataset)
    
    color_dict = {'fast': 'red', 'slow': 'blue', 'mid': 'green'}
    
    dist_travelled = []
    colors = []
    for classification, keys in classified_keys.items():
        dist_travelled_for_class = []
        for key in keys:
            trajec = dataset.trajecs[key] 
            dist_travelled_for_class.append(trajec.dist_travelled[-1])
        dist_travelled.append(dist_travelled_for_class)
        
        try:
            color = color_dict[classification]
        except:
            color = 'black'
            
        colors.append(color)
        
    bins = np.linspace(0.07,0.16,30)
    plt.hist(dist_travelled, bins=bins, color=colors, histtype='barstacked')
    
    return dist_travelled


    
def plot_angle_vs_saccade(dataset):
    classified_keys = get_classified_keys(dataset)
    keys = classified_keys['single_saccade']
    
    norm = matplotlib.colors.Normalize(.02, .5, clip=True)
    cmap = plt.get_cmap('jet')
    
    for key in keys:
        trajec = dataset.trajecs[key]
            
        c = cmap(norm(trajec.speed[trajec.saccades[-1]-5]))
        plt.plot(trajec.angle_to_post[trajec.saccades[-1]-5]*180/np.pi, trajec.dist_to_stim_r_normed[trajec.saccades[-1]-5], 'o', color=c)
        

def simulate_rrev_for_fake_trajectory(dataset, pos0, vel0, time, behavior='landing'):
    trajec_example = dataset.trajecs[dataset.trajecs.keys()[0]]

    trajec = copy.copy(trajec_example)
    dt = 0.01
    radius = trajec.stimulus.radius
    trajec.epoch_time = np.arange(0, time, dt)
    trajec.fly_time = trajec.epoch_time
    nf = len(trajec.epoch_time)
    trajec.length = nf
    trajec.velocities = np.zeros([nf, 3])
    trajec.positions = np.zeros([nf, 3])
    trajec.dist_to_stim_r_normed = np.zeros([nf])
    #print trajec.dist_to_stim_r_normed.shape
    trajec.velocities[0] = vel0
    trajec.positions[0] = pos0
    trajec.dist_to_stim_r_normed[0] = np.linalg.norm(pos0)
    for f in range(1,nf): 
        trajec.positions[f] = trajec.positions[f-1] + trajec.velocities[f-1]*dt
        trajec.dist_to_stim_r_normed[f] = np.linalg.norm(trajec.positions[f]) - radius
        
        if trajec.dist_to_stim_r_normed[f] > 0:
            trajec.velocities[f] = trajec.velocities[f-1]
        else:
            trajec.velocities[f] = 0
            
            
            
    trajec.behavior = behavior
    trajec.speed = sa1.norm_array(trajec.velocities)
    trajec.calc_accel_1d()
    trajec.dist_to_stim_r = sa1.norm_array(trajec.positions[:,0:2]) - radius
    
    prep_trajectory(trajec)
    coeff = dataset.rrev_vs_speed_fit*trajec.speed
    trajec.rrev = trajec.expansion / trajec.angle_subtended_by_post**2

    return trajec
    
    
def plot_simulated_rrevs(dataset, keys):
    trajec_example = dataset.trajecs[dataset.trajecs.keys()[0]]
    
    figure=None
    fig = plt.figure(figure)
    ax = fig.add_subplot(111)
    colormap_norm = matplotlib.colors.Normalize(0, .01, clip=True)
    cmap = plt.get_cmap('jet')
    arrow = dict(arrowstyle="->")
    
    angle_subtended_norm = (5, 40)
    cl = colorline.Colorline(xlim=[0,1], ylim =[0,20], norm=angle_subtended_norm, colormap = 'jet', figure=figure, hide_colorbar=True, ax0=ax)   
    
    speed = np.linspace(0.25, 0.6, 5)
    delays = np.linspace(0.001, 0.01, 10)
    
    for delay in delays:
        rrev_for_speed_arr = []
        t_at_delay_arr = []
        for i, s in enumerate(speed):
            # center, slow
            pos0 = [-0.2, 0, 0]
            vel0 = [s, 0, 0]
            time = 2
            
            trajec = simulate_rrev_for_fake_trajectory(dataset, pos0, vel0, time)
            #post_first_visible = np.where( trajec.angle_subtended_by_post > 5*np.pi/180. )[0][0]
            frames = np.arange(0,trajec.frame_of_landing-1).tolist()
            
            rrev_for_speed = dataset.rrev_vs_speed_fit*s
            
            t = np.interp(rrev_for_speed, trajec.rrev[frames[10:]], trajec.epoch_time[frames[10:]])
            
            '''
            tdelayed = t-delay
            tminangle = np.interp(min_angle_to_see_post, trajec.angle_subtended_by_post[frames], trajec.epoch_time[frames])
            tthreshold = np.interp(rrev_threshold, trajec.rrev[frames], trajec.epoch_time[frames])
            
            ttrigger = np.max([tminangle, tthreshold])+delay
            
            rrev_at_time = np.interp(ttrigger, trajec.epoch_time[frames], trajec.rrev[frames])
            rrev_for_speed_arr.append(rrev_at_time)
            t_at_delay_arr.append(ttrigger)
            '''
            tdelayed = t-delay
            rrev_at_tdelay = np.interp(tdelayed, trajec.epoch_time[frames], trajec.rrev[frames])
            c = cmap(colormap_norm(delay))
            #ax.plot(tdelayed, rrev_at_tdelay, 'o', color=c)
            
            if delay==delays[0]:
                cl.colorline(trajec.epoch_time[frames], trajec.rrev[frames], trajec.angle_subtended_by_post[frames]*180./np.pi,linewidth=1, norm=angle_subtended_norm, alpha=1)
                ax.plot(t, rrev_for_speed, 'o', color='black')
                
            if delay==delays[0] and s==speed[-1]:   
                
                #ax.plot(trajec.epoch_time[frames], trajec.rrev[frames], color='gray')
                string = 'RREV at initiation of deceleration'
                string_position = (t-0.05, rrev_for_speed+1)
                '''
                ax.annotate(string, (t, rrev_for_speed),
                    xytext=string_position,
                    arrowprops=arrow,
                    horizontalalignment='right', verticalalignment='top')
                ''' 
        #c = cmap(colormap_norm(delay))
        c = 'red'
        
        cl.ax0.set_ylim([0,20])
        
        #ax.plot(t_at_delay_arr, rrev_for_speed_arr, '--', color=c, linewidth=2)
        #ax.plot(t_at_delay_arr, rrev_for_speed_arr, 'o', color=c, markeredgecolor=c)
        
    # plot angle grid lines
    angles = [5*np.pi/180.]
    for min_angle_to_see_post in angles:
        rrev_for_angles = []
        t_for_angles = []
        for i, s in enumerate(speed):
            # center, slow
            pos0 = [-.2, 0, 0]
            vel0 = [s, 0, 0]
            time = 2
            trajec = simulate_rrev_for_fake_trajectory(dataset, pos0, vel0, time)
            frames = np.arange(0,trajec.frame_of_landing-1).tolist()
        
            rrev_for_angles.append( np.interp(min_angle_to_see_post, trajec.angle_subtended_by_post[frames], trajec.rrev[frames]) )
            t_for_angles.append( np.interp(min_angle_to_see_post, trajec.angle_subtended_by_post[frames], trajec.epoch_time[frames]) )
    
        c = 'black' #cmap(angles_norm(min_angle_to_see_post*180/np.pi))
        ax.plot(t_for_angles, rrev_for_angles, '-', color=c)
    #ax.hlines(rrev_threshold, 0, 1, color='black')
    #ax.text(.1, rrev_threshold, 'RREV threshold')
    
    '''
    string = '20 deg threshold'
    string_position = (.1, 10)
    ax.annotate(string, (.314, 7.6),
        xytext=string_position,
        arrowprops=arrow,
        horizontalalignment='center', verticalalignment='top')
    '''         
        
    ax.set_title('delay in RREV calculation: 60 ms')
    ax.set_xlabel('time, sec')
    ax.set_ylabel('simulated RREV (1/s) for constant velocity trajectories')
    
   
    colorbar_pad = 0
    colorbar_size = "3%"
    divider = make_axes_locatable(ax)
    divider.set_anchor('E')
    cax = divider.append_axes("right", size=colorbar_size, pad=colorbar_pad)
    angle_norm = matplotlib.colors.Normalize(angle_subtended_norm[0], angle_subtended_norm[1], clip=True)
    cb = matplotlib.colorbar.ColorbarBase(cax, cmap=cmap, norm=colormap_norm, orientation='vertical', boundaries=None)
    cax.set_ylabel('angle subtended, deg')
    
    plt.show()
   
    
def get_leg_extension_from_movie_dataset(dataset, movie_dataset, objid_dict):
    
    for key, trajec in dataset.trajecs.items():
        try:
            movieid = objid_dict[trajec.obj_id][0][1]
            trajec.movieid = movieid
            movie = movie_dataset.movies[movieid]
            trajec.legextensionrange = movie.legextensionrange
            
            if trajec.legextensionrange is not None:
                trajec.legextension_time = sa1.frame_to_timestamp(movie, movie.legextensionrange)[0]
                print movieid, movie.landingtime - trajec.legextension_time
            else:  
                trajec.legextension = None
        except:
            trajec.legextension = None
            
def create_dataset_from_movie_dataset(dataset, movie_dataset):
    hs_dataset = ffa.Dataset(like=dataset)
    landing_keys = movie_dataset.get_movie_keys('landing')
    
    for key in landing_keys:
        movie = movie_dataset.movies[key]
        if movie.trajec is not None:
            hs_dataset.trajecs.setdefault(key, movie.trajec)
            trajec = movie.trajec
            trajec.legextensionrange = movie.legextensionrange
            trajec.behavior = movie.behavior
            trajec.key = key
            
            if trajec.legextensionrange is not None:
                trajec.legextension_time = sa1.frame_to_timestamp(movie, movie.legextensionrange)[0]
                print movie.landingtime - trajec.legextension_time
            else:  
                trajec.legextension = None
            
    return hs_dataset
    
    
def get_flydra_frame_at_timestamp(trajec, t):
    return np.argmin(np.abs(trajec.epoch_time - t))
            
                
def plot_leg_extension_from_movie_dataset(dataset, movie_dataset, objid_dict, keys):
    
    for key in keys:
        trajec = dataset.trajecs[key]
        
        movie = None
        try:
            movieid = objid_dict[trajec.obj_id][0][1]
            movie = movie_dataset.movies[movieid]
        except:
            pass
            
        if movie is not None:
            if movie.legextensionrange is not None:
                
                legextensionframe = movie.legextensionrange[0] - movie.firstframe_ofinterest
                #legextension_time = sa1.frame_to_timestamp(movie, movie.legextensionrange)[0]
                angles_subtended = movie.scaled.angle_subtended_by_post[legextensionframe - 1: legextensionframe+1]
                exp_at_legextension = angles_subtended[1]-angles_subtended[0]
                angle_at_legextension = angles_subtended[1]
                plt.plot( angle_at_legextension, exp_at_legextension / angle_at_legextension, 'o' )
                
                
            else:
                trajec.legextension = None
                
def get_radius_from_altitude(a):
    m = (.0095965 - 0.005) / 0.15
    return m*np.abs(a) + 0.005
def calc_radius_at_nearest(trajec):
    if 'angled' in trajec.post_type:
        frame_nearest_post = np.argmin(trajec.dist_to_stim_r)
        a = trajec.positions[frame_nearest_post, 2]
        trajec.radius_at_nearest = get_radius_from_altitude(a)
    else:
        trajec.radius_at_nearest = 0.009565
    
def test_legs(movie_dataset, dataset_landing=None):
    
    keys = movie_dataset.movies.keys()
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    angle_at_leg_extension = []
    for movieid in keys:
        movie = movie_dataset.movies[movieid]
        if movie.behavior == 'landing' and 'crash' not in movie.subbehavior and 'wingcrash' not in movie.subbehavior:        
            try:
                tmp = movie.scaled
                tmp = True
            except:
                tmp = False
            if tmp:
                if movie.legextensionrange is not None:
                    legextensionframe = movie.legextensionrange[0] - movie.firstframe_ofinterest
                    print movie.scaled.dist_to_post[legextensionframe] - get_radius_from_altitude(movie.trajec.positions[movie.trajec.frame_of_landing,2]), movie.scaled.dist_to_post[movie.landingframe_relative] - get_radius_from_altitude(movie.trajec.positions[movie.trajec.frame_of_landing,2])
                    
                    radius = get_radius_from_altitude(movie.trajec.positions[movie.trajec.frame_of_landing,2])
                    
                    dleg = movie.scaled.dist_head_to_post[legextensionframe] - radius
                    dland = movie.scaled.dist_head_to_post[movie.landingframe_relative] - radius
                    
                    angle_leg_ext = 2*np.arcsin( radius / movie.scaled.dist_head_to_post[legextensionframe] )
                    speed_leg_ext = movie.scaled.speed[legextensionframe]
                    
                    ax.plot( angle_leg_ext, speed_leg_ext, 'o', color='red' ) 
                    
                    t_decel = movie.trajec.time_at_deceleration
                    sa1_frame_of_deceleration = sa1.get_frame_from_timestamp(movie, t_decel)
                    angle_decel = 2*np.arcsin( radius / movie.scaled.dist_head_to_post[sa1_frame_of_deceleration] )
                    
                    
                    
                    d = angle_leg_ext - angle_decel
                    if d > 0:
                        c = 'black'
                    else:
                        c = 'blue'
                            
                    ax.plot( angle_decel, movie.scaled.speed[sa1_frame_of_deceleration], 'o', color=c)
                    #plt.plot( [angle_leg_ext, angle_decel], [speed_leg_ext, movie.scaled.speed[sa1_frame_of_deceleration]], '--', color=c)
                        
                        
                    
                    angle_at_leg_extension.append( angle_leg_ext )
            else:
                #sa1.process_movieinfo(movie)
                print movieid
                
    if dataset_landing is not None:
        curve, mean_error, std_error, errors = plot_angle_vs_speed_landing(dataset_landing, plot=False)
        ax.plot(curve[:,0],curve[:,1],color='black')
        
    fig2 = plt.figure()
    ax2 = fig2.add_subplot(111)
    ax2.hist( np.log(angle_at_leg_extension), bins=10 )
                
    return angle_at_leg_extension
                
def plot_movie_dataset_leg_extension(movie_dataset, dataset_landing=None):
    
    #straight_ids = ['20101110_C001H001S0038', '20101110_C001H001S0032', '20101110_C001H001S0008', '20101113_C001H001S0020', '20101101_C001H001S0009', '20101030_C001H001S0035', '20101111_C001H001S0035', '20101101_C001H001S0001', '20101101_C001H001S0009', '20101110_C001H001S0035', '20101111_C001H001S0034', '20101110_C001H001S0004', '20101111_C001H001S0001', '20101101_C001H001S0020', '20101101_C001H001S0002', '20101110_C001H001S0039', '20101111_C001H001S0054', '20101111_C001H001S0058', '20101110_C001H001S0027', '20101110_C001H001S0014', '20101111_C001H001S0005']
    
    keys = movie_dataset.movies.keys()
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    dist_at_leg_extension = []
    
    for movieid in keys:
        movie = movie_dataset.movies[movieid]
        if movie.behavior == 'landing':        
            if movie.legextensionrange is not None:
                if 1: #try:
                    legextensionframe = movie.legextensionrange[0] - movie.firstframe_ofinterest
                    sa1_time = movie.timestamps[legextensionframe]
                    flydra_frame = get_flydra_frame_at_timestamp(movie.trajec, sa1_time)
                    
                    
                    #movie.trajec.behavior = movie.behavior
                    #movie.trajec.key = movie.id
                    #prep_trajectory(movie.trajec, distance = 0.04)                    
                    
                    try:
                        tmp = movie.scaled
                        tmp = True
                    except:
                        tmp = False
                    
                    if tmp: #movie.trajec.classification == 'straight':
                        
                        ax.plot( movie.scaled.angle_subtended_by_post[legextensionframe], movie.scaled.speed[legextensionframe], 'o', color='red')
                        t_decel = movie.trajec.time_at_deceleration
                        sa1_frame_of_deceleration = sa1.get_frame_from_timestamp(movie, t_decel)
                        
                        ax.plot( movie.scaled.angle_subtended_by_post[sa1_frame_of_deceleration], movie.scaled.speed[sa1_frame_of_deceleration], 'o', color='black')
                        

                        d = movie.scaled.angle_subtended_by_post[legextensionframe] - movie.scaled.angle_subtended_by_post[sa1_frame_of_deceleration]
                        if d > 0:
                            c = 'black'
                        else:
                            c = 'red'
                            

                        ax.plot( [movie.scaled.angle_subtended_by_post[legextensionframe], movie.scaled.angle_subtended_by_post[sa1_frame_of_deceleration]], [movie.scaled.speed[legextensionframe], movie.scaled.speed[sa1_frame_of_deceleration]], '--', color=c)
                        
                        dist_at_leg_extension.append(movie.scaled.dist_head_to_post[legextensionframe])
                        
                        print movie.trajec.frame_at_deceleration
                    
                    #legextension_time = sa1.frame_to_timestamp(movie, movie.legextensionrange)[0]
                    #angles_subtended = movie.scaled.angle_subtended_by_post[legextensionframe - 1: legextensionframe+1]
                    #exp_at_legextension = (angles_subtended[1]-angles_subtended[0])*movie.framerate
                    #angle_at_legextension = angles_subtended[1]
                    #plt.plot( angle_at_legextension*180/np.pi, (exp_at_legextension / angle_at_legextension)*180/np.pi, 'o', color='red' )
                    
                    #plt.text( angle_at_legextension, exp_at_legextension / angle_at_legextension, movieid)
                #except:
                #    pass
                
            
    if dataset_landing is not None:
        curve, mean_error, std_error, errors = plot_angle_vs_speed_landing(dataset_landing, plot=False)
        ax.plot(curve[:,0],curve[:,1],color='black')
                
    return np.array(dist_at_leg_extension)
                
def test(dataset, movie_dataset, objid_dict):
    
    for movieid, movie in movie_dataset.movies.items():
        if movie.trajec is not None:
            print movie.trajec.obj_id
            movie.trajec.key = find_trajec_key_from_obj_id(dataset, movie.trajec.obj_id)
            if movie.trajec.key is None:
                print 'no key: ', len(movie.trajec.speed), movie.trajec.behavior
            
        else:
            print 'no trajec: ', movieid
            
def find_trajec_key_from_obj_id(dataset, obj_id):
    
    obj_id_list = []
    for key, trajec in dataset.trajecs.items():
        if trajec.obj_id == obj_id:
            print key
            obj_id_list.append(key)
    if len(obj_id_list) == 1:
        return obj_id_list[0]
    elif len(obj_id_list) > 1:
        print 'too many matches'
        return None
    else:
        return None
        
        
            
def assign_post_type(dataset, post_type):
    
    for k, trajec in dataset.trajecs.items():
        trajec.post_type = post_type
            
            
def plot_angle_vs_speed_landing(dataset, keys_to_plot_frames=None, keys=None, plot=True):
    
    classified_keys = get_classified_keys(dataset)
    if keys_to_plot_frames is None:
        keys_to_plot_frames = classified_keys['straight']
    if keys is None:
        keys = dataset.trajecs.keys()#classified_keys['straight']
    
    speed_norm = (0,0.6)
    colormap_norm = matplotlib.colors.Normalize(speed_norm[0], speed_norm[1], clip=True)
    cmap = plt.get_cmap('jet')
    
    if plot:
        figure=None
        fig1 = plt.figure(figure)
        ax_accel_vs_speed = fig1.add_subplot(111)

        figure=None
        fig2 = plt.figure(figure)
        ax_speed_vs_angle = fig2.add_subplot(111)
        cl = colorline.Colorline(xlim=[0, np.pi], ylim =[0, 0.8], norm=speed_norm, colormap = 'RdBu', hide_colorbar=True, figure=ax_speed_vs_angle.figure.number, ax0_size=ax_speed_vs_angle)
    
    speed = dataset.speed_at_deceleration
    trajectories = dataset.trajectories_at_deceleration
    angle = dataset.angle_at_deceleration
    expansion = dataset.expansion_at_deceleration
    distance = dataset.distance_at_deceleration
    
    sorted_speed = np.argsort(speed)
    accel_vs_speed_individual_fit = np.zeros([len(sorted_speed),2])
    
    n_trajecs_to_show = 5
    interval = float(int(len(keys_to_plot_frames) / float(n_trajecs_to_show)))
    print interval
    
    for i in range(len(sorted_speed)):
        key = trajectories[i]
        trajec = dataset.trajecs[key]
        
        frames = np.arange(trajec.frame_at_deceleration, trajec.frame_of_landing).tolist()
        accel_vs_speed = np.polyfit(trajec.angle_subtended_by_post[frames], trajec.speed[frames],1)
        accel_vs_speed_individual_fit[i] = accel_vs_speed
        
        '''
        if plot:
            if key in keys_to_plot_frames:
                for f in trajec.frames:
                    c = cmap(colormap_norm(trajec.speed[f]))
                    ax_speed_vs_angle.plot(trajec.angle_subtended_by_post[f], trajec.speed[f], 'o', color=c, markeredgecolor=c, alpha=0.1)
            
        
        '''
        
    
            
            
        #accel_vs_speed_fit = np.polyfit(speed, accel_vs_speed_individual_fit[:,0], 1)
        accel_vs_speed_fit, variance, intercept_confidence_interval, slope_confidence_interval = floris.linear_fit_type2(speed, accel_vs_speed_individual_fit[:,0], full_output=True)
        x = np.linspace(0, np.max(speed), 50)
        y = np.polyval(accel_vs_speed_fit, x)
        if plot:
            ax_accel_vs_speed.plot(speed, accel_vs_speed_individual_fit[:,0], 'o', color='gray')
            ax_accel_vs_speed.set_xlabel('speed at point of deceleration, m/s')
            ax_accel_vs_speed.set_ylabel('"acceleration", m / s*rad')
            ax_accel_vs_speed.plot(x,y,color='black')
    
    
    indices = np.where( (np.isnan(distance) == False)*(np.isnan(expansion) == False)*(distance>0)*(angle<3)*(angle>5*np.pi/180.) )[0].tolist()
    post_radius = 0.009565
    def speed_vs_distance_func(x,a, b):
        return a*x**2 + b*x + post_radius
    #fit_dist_vs_speed = curve_fit(dist_vs_speed_func, distance[indices]+post_radius, speed[indices])[0]
    fit_speed_vs_distance = curve_fit(speed_vs_distance_func, speed[indices], distance[indices]+post_radius)[0]
    
    # get std dev of errors
    #fitted_dist_for_data = speed_vs_distance_func(fitted_speed,fit_speed_vs_distance[0], fit_speed_vs_distance[1])
    
    #
    
    fitted_speed = np.linspace(0.0001, np.max(speed), 1000)
    fitted_dist = speed_vs_distance_func(fitted_speed,fit_speed_vs_distance[0], fit_speed_vs_distance[1])
    fitted_angle = 2*np.arcsin( post_radius / fitted_dist )
    
    errors = np.zeros_like(speed)
    curve = np.zeros([len(fitted_angle), 2])
    curve[:,0] = fitted_angle
    curve[:,1] = fitted_speed
    for i in indices:
        pt = np.array([ angle[i], speed[i] ])
        y_from_curve = np.interp( pt[0], curve[::-1,0],curve[::-1,1])
        sign = np.sign( pt[1]-y_from_curve )
        errors[i] = np.min( sa1.norm_array( (curve-pt) ) )*sign
    mean_error = np.mean( errors[indices] )
    std_error = np.std( errors[indices] )
    print 'mean error of fit: ', mean_error    
    print 'std dev of error of fit: ', std_error  
    
    if plot:
        ax_speed_vs_angle.plot( angle[indices], speed[indices], 'o', color='gray', markeredgecolor='gray', alpha=1)
        ax_speed_vs_angle.plot(fitted_angle,fitted_speed,'--', color='black')
        
        #ax_speed_vs_angle.plot(fitted_angle+std_error/np.sqrt(2),fitted_speed+std_error/np.sqrt(2),'--', color='black')
        #ax_speed_vs_angle.plot(fitted_angle-std_error/np.sqrt(2),fitted_speed-std_error/np.sqrt(2),'--', color='black')
        
        ax_speed_vs_angle.set_xlabel('angle, rad')
        ax_speed_vs_angle.set_ylabel('speed, m/s')
        ax_speed_vs_angle.semilogx()
    
    if plot:
        # get evenly spaced keys in decreasing speed order
        p_speed = []
        for key in keys_to_plot_frames:
            trajec = dataset.trajecs[key]
            p_speed.append( trajec.speed_at_deceleration )
        p_speed_sorted = np.argsort( np.array( p_speed ))
        '''
        n_trajecs_to_plot = 5
        indices_to_plot_tmp = np.linspace(0, len(p_speed_sorted)-1, n_trajecs_to_plot).tolist()
        indices_to_plot = [int(i) for i in indices_to_plot_tmp]
        keys_to_plot = [keys_to_plot_frames[i] for i in indices_to_plot]
        print keys_to_plot
        '2_29065', '2_31060', '1_19435', '8_10323'
        '''
        
        keys_to_plot = ['2_29065', '2_31060', '8_10323', '6_715']
        print keys_to_plot
        
        
        for key in keys_to_plot:
            trajec = dataset.trajecs[key]
            ftp = np.arange(trajec.frames[0]-25, trajec.frames[-1]).tolist()
            cl.colorline(trajec.angle_subtended_by_post[ftp], trajec.speed[ftp], -1*trajec.accel_1d[ftp], norm=(-2,2)) 
            ax_speed_vs_angle.plot( trajec.angle_at_deceleration, trajec.speed_at_deceleration, 'o', color='black', markeredgecolor='black', alpha=1)
        
    fig3 = plt.figure()
    ax = fig3.add_subplot(111)
    ax.plot( angle, expansion, '.')
    adot = 2*(post_radius*fitted_dist**(-2))*fitted_speed / np.sqrt(1-(post_radius/fitted_dist)**2)
    ax.plot( fitted_angle, adot)
    ax.set_xlabel('retina size, rad')
    ax.set_ylabel('retinal expansion velocity, rad/s')
    
    
    plt.show()
    
    return curve, mean_error, std_error, errors
    
        
def plot_angle_vs_speed_flyby(dataset, keys=None, dataset_landing=None):
    if keys is None:
        classified_keys = get_classified_keys(dataset)
        keys = classified_keys['straight']

    fig = plt.figure()
    ax = fig.add_subplot(111)

    colormap_norm = matplotlib.colors.Normalize(0.003, 0.01, clip=True)
    cmap = plt.get_cmap('cool')
        
    if dataset_landing is not None:
        curve, mean_error, std_error, errors = plot_angle_vs_speed_landing(dataset_landing, plot=False)
        ax.plot(curve[:,0],curve[:,1],color='black')
        ax.plot( dataset_landing.angle_at_deceleration, dataset_landing.speed_at_deceleration, '.', color='black', markeredgecolor='black', alpha=1)
        
    saccades = []

    for key in keys:
        trajec = dataset.trajecs[key]
        frame_nearest_post = np.argmin(trajec.dist_to_stim_r)
        
        c = cmap(colormap_norm( np.min(trajec.dist_to_stim_r) ))
        
        if 1:#trajec.post_type == 'checkered':
                    
            if len(trajec.saccades) > 0:
                f = trajec.saccades[-1]
                ax.plot(trajec.angle_subtended_by_post[f], trajec.speed[f], 'o', color=c, markeredgecolor=c)
                
                saccades.append( trajec.angle_subtended_by_post[f] )
                    
            elif trajec.frame_at_deceleration is not None:
                f = trajec.frame_at_deceleration
                ax.plot(trajec.angle_subtended_by_post[f], trajec.speed[f], '*', color=c, markeredgecolor=c)
                
        
    
        
    plt.show()
        
    return saccades






def saccade_angle_histogram(dataset, keys=None, plot=False, plot_sample_trajecs=False, post_type=['checkered', 'checkered_angled', 'black', 'black_angled'], filename=None, all_saccades=False):
    if keys is None:
        classified_keys = get_classified_keys(dataset)
        #keys = classified_keys['straight']
        keys = dataset.trajecs.keys()
        
    angle = []
    speed = []
    for key in keys:
        trajec = dataset.trajecs[key]
        
        if trajec.post_type in post_type:
            if len(trajec.saccades) > 0:
                if all_saccades is False:
                    saccade_frames = get_saccade_range(trajec, trajec.saccades[-1])
                    #f = trajec.saccades[-1]
                    f = saccade_frames[0]
                    if trajec.angle_subtended_by_post[f] > 0:
                        angle.append( trajec.angle_subtended_by_post[f] )
                        speed.append( trajec.speed[f] )
                elif all_saccades is True:
                    for s in trajec.saccades:
                        if trajec.behavior == 'flyby':
                            if saccade_analysis.is_evasive(trajec, s):
                                saccade_frames = get_saccade_range(trajec, s)
                                #f = trajec.saccades[-1]
                                f = saccade_frames[0]
                                if trajec.angle_subtended_by_post[f] > 0:
                                    angle.append( trajec.angle_subtended_by_post[f] )
                                    speed.append( trajec.speed[f] )
                        if trajec.behavior == 'landing':
                            if not saccade_analysis.is_evasive(trajec, s):
                                saccade_frames = get_saccade_range(trajec, s)
                                #f = trajec.saccades[-1]
                                f = saccade_frames[0]
                                if trajec.angle_subtended_by_post[f] > 0:
                                    angle.append( trajec.angle_subtended_by_post[f] )
                                    speed.append( trajec.speed[f] )
    print 'N = ', len(speed)
    speed = np.array(speed)
    angle = np.log(np.array(angle))
    
    data, bins = np.histogram(angle, bins=16, normed=True)
    xvals = np.diff(bins) + bins[0:-1]
    
    butter_b, butter_a = signal.butter(3, 0.3)
    data_filtered = signal.filtfilt(butter_b, butter_a, data)

    if plot:
        fig = plt.figure()
        fig.set_facecolor('white')
        ax = fig.add_subplot(111)
        ax.hist( angle, bins=bins, normed=True, facecolor='green', alpha=0.1, edgecolor='green')
        ax.plot( xvals, data_filtered, color='green' )  
        ax.plot( angle, speed, '.', color='green', alpha=0.8)   
        
        #fit = np.polyfit( angle, speed, 1)
        fit, variance, intercept_confidence_interval, slope_confidence_interval, Rsq = floris.linear_fit_type2(angle, speed, full_output=True)
        x = np.linspace(np.min(angle), np.max(angle), 100)
        y = np.polyval(fit, x)
        ax.plot(x,y,color='green')
        
        yplus = np.polyval(fit, x+np.sqrt(variance))+np.sqrt(variance)
        yminus = np.polyval(fit, x-np.sqrt(variance))-np.sqrt(variance)
        ax.fill_between(x, yplus, yminus, color='green', linewidth=0, alpha=0.2)
    
        slope = str(fit[0])[0:5]
        intercept = str(fit[1])[0:5]
        string = 'y = ' + slope + 'x + ' + intercept
        ax.text(0, 0.5, string)
        string = 'Rsq=' + str(Rsq)[0:4]
        ax.text(0,0.4, string)
        
        if plot_sample_trajecs:
            k = [40]
            ks = [keys[i] for i in k]
            print ks
            for key in ks:
                trajec = dataset.trajecs[key]
                ax.plot( np.log(trajec.angle_subtended_by_post[trajec.frames_of_flyby]), trajec.speed[trajec.frames_of_flyby], '--', color='black', linewidth=2)
                ax.plot( np.log(trajec.angle_subtended_by_post[trajec.saccades[-1]]), trajec.speed[trajec.saccades[-1]], '.', color='black')
        
        
        fix_angle_log_spine(ax)    
        string = 'N=' + str(len(speed))
        ax.text( 10.*np.pi/180., 0.8, string, color='green')
        if filename is not None:
            fig.savefig(filename, format='pdf', bbox='tight')
    return angle, speed, data_filtered, xvals
    
def deceleration_angle_histogram_flyby(dataset, keys=None, plot=False, saccades=None, filename=None, angle_to_post_range=[0,5*np.pi/180.]):
    if keys is None:
        classified_keys = get_classified_keys(dataset)
        keys = classified_keys['straight']
        keys = dataset.trajecs.keys()
        
    angle = []
    speed = []
    for key in keys:
        trajec = dataset.trajecs[key]
        if trajec.frame_at_deceleration is not None:
            if np.abs(trajec.angle_to_post[trajec.frame_at_deceleration]) < angle_to_post_range[-1] and np.abs(trajec.angle_to_post[trajec.frame_at_deceleration]) > angle_to_post_range[0]:
                
                if len(trajec.saccades) > 0:
                    tmp = np.array(trajec.saccades) - trajec.frame_at_deceleration
                    try:
                        s = np.where(tmp > 0)[0][0]
                    except:
                        s = None
                        
                    if s is not None:
                        s = trajec.saccades[s]
                        dt = trajec.epoch_time[s] - trajec.epoch_time[trajec.frame_at_deceleration]
                    else:
                        dt = np.inf
                else:
                    dt = np.inf
                    
                if dt < 0.05:
                    continue
                
                if saccades is False:
                    if len(trajec.saccades) == 0:
                        angle.append( trajec.angle_at_deceleration )
                        speed.append( trajec.speed_at_deceleration )
                elif saccades is True:
                    if len(trajec.saccades) > 0:
                        angle.append( trajec.angle_at_deceleration )
                        speed.append( trajec.speed_at_deceleration )
                elif saccades is None:
                    angle.append( trajec.angle_at_deceleration )
                    speed.append( trajec.speed_at_deceleration )
            
    speed = np.array(speed)
    angle = np.log(np.array(angle))
            
    data, bins = np.histogram(angle, bins=25, normed=True)
    xvals = np.diff(bins) + bins[0:-1]
    
    butter_b, butter_a = signal.butter(3, 0.3)
    data_filtered = signal.filtfilt(butter_b, butter_a, data)

    #fit = np.polyfit( angle, speed, 1 )
    fit, variance, intercept_confidence_interval, slope_confidence_interval, Rsq = floris.linear_fit_type2(angle, speed, full_output=True)
    x = np.linspace( np.min(angle), np.max(angle), 20)
    y = np.polyval(fit, x)
    curve = np.zeros([len(x), 2])
    curve[:,0] = x
    curve[:,1] = y
    yplus = np.polyval(fit, x+np.sqrt(variance))+np.sqrt(variance)
    yminus = np.polyval(fit, x-np.sqrt(variance))-np.sqrt(variance)

    if plot:
        fig = plt.figure()
        fig.set_facecolor('white')
        ax = fig.add_subplot(111)
        ax.hist( angle, bins=bins, normed=True, facecolor='blue', alpha=0.1, edgecolor='blue')
        ax.plot( xvals, data_filtered, '-', color='blue' )     
        ax.plot( angle, speed, '.', color='blue', alpha=0.8)
        
        ax.plot(x,y, '-', color='blue')
        slope = str(fit[0])[0:5]
        intercept = str(fit[1])[0:5]
        string = 'y = ' + slope + 'x + ' + intercept
        ax.text(0, 0.5, string)
        string = 'Rsq=' + str(Rsq)[0:4]
        ax.text(0,0.4, string)
        
        
        ax.fill_between(x, yplus, yminus, color='blue', linewidth=0, alpha=0.2)
        
        fix_angle_log_spine(ax)    
        
        string = 'N=' + str(len(speed))
        ax.text( 10.*np.pi/180., 0.8, string, color='blue')
        if filename is None:
            filename = 'deceleration_histogram_flyby_saccades_' + str(saccades) + '.pdf'
        fig.savefig(filename, format='pdf')
    return angle, bins, data_filtered, xvals, curve, yplus, yminus
    
    
def deceleration_saccade_angle_flyby(dataset, keys=None, plot=True, filename=None, angle_to_post_range=[0,5*np.pi/180.]):
    if keys is None:
        classified_keys = get_classified_keys(dataset)
        keys = classified_keys['straight']
        keys = dataset.trajecs.keys()
        
    fig = plt.figure()
    fig.set_facecolor('white')
    ax = fig.add_subplot(111)
    
    keys = keys[0:50]
        
    for key in keys:
        trajec = dataset.trajecs[key]
        if trajec.frame_at_deceleration is not None:
            if 1: #np.abs(trajec.angle_to_post[trajec.frame_at_deceleration]) < angle_to_post_range[-1] and np.abs(trajec.angle_to_post[trajec.frame_at_deceleration]) > angle_to_post_range[0]:
            
                # find saccade after frame of deceleration:
                if len(trajec.saccades) > 0:
                    tmp = np.array(trajec.saccades) - trajec.frame_at_deceleration
                    try:
                        s = np.where(tmp > 0)[0][0]
                    except:
                        s = None
                        
                    if s is not None:
                        s = trajec.saccades[s]
                        print trajec.epoch_time[s] - trajec.epoch_time[trajec.frame_at_deceleration]
                        ax.plot( np.log(trajec.angle_subtended_by_post[trajec.frame_at_deceleration]), trajec.speed[trajec.frame_at_deceleration], '.', color='blue')
                        ax.plot( np.log(trajec.angle_subtended_by_post[s]), trajec.speed[s], '.', color='green')
                        ax.plot( np.log(trajec.angle_subtended_by_post[trajec.frame_at_deceleration:s+1]), trajec.speed[trajec.frame_at_deceleration:s+1], '-', color='black', linewidth=0.3)
                
        if filename is None:
            filename = 'deceleration_saccade_angle_flyby.pdf'
        fig.savefig(filename, format='pdf')
    return
    
    
def leg_extension_angle_histogram(movie_dataset, plot=False, behavior='landing'):
    keys = movie_dataset.movies.keys()
    
    angle_at_leg_extension = []
    speed_at_leg_extension = []
    n = 0
    for movieid in keys:
        movie = movie_dataset.movies[movieid]
        if movie.behavior == behavior:        
            if movie.legextensionrange is not None:
                legextensionframe = movie.legextensionrange[0] - movie.firstframe_ofinterest
                #sa1_time = movie.timestamps[legextensionframe]
                #flydra_frame = get_flydra_frame_at_timestamp(movie.trajec, sa1_time)
                
                try:
                    tmp = movie.scaled
                    tmp = True
                except:
                    tmp = False
                
                if tmp: #movie.trajec.classification == 'straight':
                    angle_at_leg_extension.append(movie.scaled.angle_subtended_by_post[legextensionframe][0])
                    speed_at_leg_extension.append(movie.scaled.speed[legextensionframe])
                    n += 1
              
    angle_at_leg_extension = np.array(angle_at_leg_extension)
    speed_at_leg_extension = np.array(speed_at_leg_extension)
    
    data, bins = np.histogram( np.log(angle_at_leg_extension), bins=16, normed=True)
    xvals = np.diff(bins) + bins[0:-1]
    
    butter_b, butter_a = signal.butter(3, 0.3)
    data_filtered = signal.filtfilt(butter_b, butter_a, data)
    
    print 'N = ', n
    
    if plot:
        fig = plt.figure()
        fig.set_facecolor('white')
        ax = fig.add_subplot(111) 
        ax.hist( np.log(angle_at_leg_extension), bins=bins, normed=True, facecolor='red', alpha=0.1, edgecolor='red')
        ax.plot( xvals, data_filtered, color='red' )     
        ax.plot( np.log(angle_at_leg_extension), speed_at_leg_extension, '.', color='red')
        
        #fit = np.polyfit( np.log(angle_at_leg_extension), speed_at_leg_extension, 1)
        print np.log(angle_at_leg_extension).shape, speed_at_leg_extension.shape
        fit, variance, intercept_confidence_interval, slope_confidence_interval, Rsq = floris.linear_fit_type2(np.log(angle_at_leg_extension), speed_at_leg_extension, full_output=True)
        x = np.linspace(np.min(np.log(angle_at_leg_extension)), np.max(np.log(angle_at_leg_extension)), 100)
        y = np.polyval(fit, x)
        ax.plot(x,y,color='red')
        
        yplus = np.polyval(fit, x+np.sqrt(variance))+np.sqrt(variance)
        yminus = np.polyval(fit, x-np.sqrt(variance))-np.sqrt(variance)
        ax.fill_between(x, yplus, yminus, color='red', linewidth=0, alpha=0.2)
    
        slope = str(fit[0])[0:5]
        intercept = str(fit[1])[0:5]
        string = 'y = ' + slope + 'x + ' + intercept
        ax.text(0, 0.5, string)
        string = 'Rsq=' + str(Rsq)[0:4]
        ax.text(0,0.4, string)
        
        fix_angle_log_spine(ax)
        plt.show()
        string = 'N=' + str(n)
        ax.text( -2, 0.8, string, color='red')
        fig.savefig('leg_ext_histogram.pdf', format='pdf')
    return angle_at_leg_extension, bins, data_filtered, xvals
    
    
def summary_figure(dataset_landing, dataset_flyby, movie_dataset):
    fig = plt.figure()
    fig.set_facecolor('white')
    ax = fig.add_subplot(111)
    
    # landing stuff
    x, y, yminus, yplus = get_angle_vs_speed_curve(dataset_landing, plot=False)
    ax.plot( x, y, color='blue')
    
    angle_at_leg_extension, bins, data_filtered, xvals = leg_extension_angle_histogram(movie_dataset, plot=False)
    ax.plot(xvals, data_filtered, color='red')
    
    
    # flyby stuff
    angle, bins, data_filtered, xvals = saccade_angle_histogram(dataset_flyby, keys=None, plot=False)
    ax.plot(xvals, data_filtered, ':', color='green')
    
    angle, bins, data_filtered, xvals, curve, yplus, yminus = deceleration_angle_histogram_flyby(dataset_flyby, keys=None, plot=False)
    ax.plot(xvals, data_filtered, ':', color='blue')
    ax.plot(curve[:,0], curve[:,1], ':', color='blue')
    
    fix_angle_log_spine(ax)    
    fig.savefig('summary_fig.pdf', format='pdf')
    
    
def get_angle_vs_speed_curve(dataset, plot=False, plot_sample_trajecs=False, post_type=['checkered', 'checkered_angled', 'black', 'black_angled'], filename=None, keys=None, tti=None):
    angleok = np.where( (dataset.angle_at_deceleration < 2)*(dataset.angle_at_deceleration > .01) )[0].tolist()
    
    if keys is None:
        classified_keys = get_classified_keys(dataset)
        keys = dataset.trajecs.keys()
    
    indices = []
    for i, key in enumerate(dataset.trajectories_at_deceleration):
        if key in keys:
            if i in angleok:
                if dataset.trajecs[key].post_type in post_type:
                    indices.append(i)
                else:
                    pass#print dataset.trajecs[key].post_type
    
    #fit = np.polyfit( np.log(dataset.angle_at_deceleration[indices]), dataset.speed_at_deceleration[indices], 1, )
    fit, variance, intercept_confidence_interval, slope_confidence_interval, Rsq  = floris.linear_fit_type2(np.log(dataset.angle_at_deceleration[indices]), dataset.speed_at_deceleration[indices], full_output=True)
    print fit
    print 'variance: ', variance, np.sqrt(variance)
    x = np.linspace(np.min(np.log(dataset.angle_at_deceleration)), np.max(np.log(dataset.angle_at_deceleration)), 100)
    y = np.polyval(fit, x)
    yplus = np.polyval(fit, x+np.sqrt(variance))+np.sqrt(variance)
    yminus = np.polyval(fit, x-np.sqrt(variance))-np.sqrt(variance)
        
    #slope_mean, slope_confidence_range, intercept_mean, intercept_confidence_range = floris.bootstrap_linear_fit(np.log(dataset.angle_at_deceleration[indices]), dataset.speed_at_deceleration[indices], n=1000)
    
    
    data, bins = np.histogram( np.log(dataset.angle_at_deceleration[indices]), bins=16, normed=True)
    xvals = np.diff(bins) + bins[0:-1]
    
    butter_b, butter_a = signal.butter(3, 0.3)
    data_filtered = signal.filtfilt(butter_b, butter_a, data)
    
    if plot:
        fig = plt.figure()
        fig.set_facecolor('white')
        ax = fig.add_subplot(111)
        ax.hist( np.log(dataset.angle_at_deceleration[indices]), bins=bins, normed=True, facecolor='blue', alpha=0.1, edgecolor='blue')
        ax.plot( xvals, data_filtered, color='blue' )    
        
        ax.plot( np.log(dataset.angle_at_deceleration[indices]), dataset.speed_at_deceleration[indices], '.', color='blue', alpha=0.8)
        ax.plot(x,y, color='blue', linewidth=1)
        
        if tti is not None:
            #tti = 0.15
            a = np.linspace(np.min((dataset.angle_at_deceleration)), np.max((dataset.angle_at_deceleration)), 100)
            exp = 2*np.tan(a/2.) / tti
            s = exp*0.009565/(2.*np.tan(a/2.)*np.sin(a/2.))
            ax.plot(np.log(a),s, color='red', linewidth=1)
        
        ax.fill_between(x, yplus, yminus, color='blue', linewidth=0, alpha=0.2)
        #ax.plot(x+np.sqrt(variance), y+np.sqrt(variance), ':')
        
        '''
        ax.plot(x,slope_confidence_range[0]*x+intercept_confidence_range[0], ':')
        ax.plot(x,slope_confidence_range[0]*x+intercept_confidence_range[1], ':')
        ax.plot(x,slope_confidence_range[1]*x+intercept_confidence_range[0], ':')
        ax.plot(x,slope_confidence_range[1]*x+intercept_confidence_range[1], ':')
        '''
        
        #print len(dataset.speed_at_deceleration)
        
        slope = str(fit[0])[0:5]
        intercept = str(fit[1])[0:5]
        string = 'y = ' + slope + 'x + ' + intercept 
        ax.text(0, 0.5, string)
        string = 'Rsq=' + str(Rsq)[0:4]
        ax.text(0,0.4, string)
        
        if plot_sample_trajecs:
            keys_to_plot = ['2_29065', '2_31060', '8_10323', '6_715']
            for key in keys_to_plot:
                trajec = dataset.trajecs[key]
                ftp = np.arange(trajec.frames[0]-25, trajec.frames[-1]).tolist()
                ax.plot( np.log(trajec.angle_subtended_by_post[ftp]), trajec.speed[ftp], color='black', linewidth=2)
                ax.plot( np.log(trajec.angle_at_deceleration), trajec.speed_at_deceleration, '.', color='black', alpha=1)
        
        fix_angle_log_spine(ax, histograms=False)    
        string = 'N=' + str(len(indices))
        ax.text( 10.*np.pi/180., 0.8, string, color='blue')
        if filename is not None:
            fig.savefig(filename, format='pdf')
    
    return x, y, yminus, yplus
    
  
def get_distance_vs_speed_curve(dataset, plot=False, plot_sample_trajecs=False, post_type=['checkered', 'checkered_angled', 'black', 'black_angled'], filename=None):
    angleok = np.where( dataset.distance_at_deceleration < 1 )[0].tolist()
    
    indices = []
    for i, key in enumerate(dataset.trajectories_at_deceleration):
        if i in angleok:
            if dataset.trajecs[key].post_type in post_type:
                indices.append(i)
            else:
                pass#print dataset.trajecs[key].post_type
    
    fit, variance, intercept_confidence_interval, slope_confidence_interval, Rsq  = floris.linear_fit_type2(dataset.distance_at_deceleration[indices], dataset.speed_at_deceleration[indices], full_output=True)
    print fit
    print 'variance: ', variance, np.sqrt(variance)
    x = np.linspace(np.min(dataset.distance_at_deceleration), np.max(dataset.distance_at_deceleration), 100)
    y = np.polyval(fit, x)
    yplus = np.polyval(fit, x+np.sqrt(variance))+np.sqrt(variance)
    yminus = np.polyval(fit, x-np.sqrt(variance))-np.sqrt(variance)
        
    
    data, bins = np.histogram( dataset.distance_at_deceleration[indices], bins=16, normed=True)
    xvals = np.diff(bins) + bins[0:-1]
    
    butter_b, butter_a = signal.butter(3, 0.3)
    data_filtered = signal.filtfilt(butter_b, butter_a, data)
    
    if plot:
        fig = plt.figure()
        fig.set_facecolor('white')
        ax = fig.add_subplot(111)
        
        ax.plot( (dataset.distance_at_deceleration[indices]), dataset.speed_at_deceleration[indices], '.', color='blue', alpha=0.8)
        ax.plot(x,y, color='blue', linewidth=1)
        
        
        ax.fill_between(x, yplus, yminus, color='blue', linewidth=0, alpha=0.2)
        #ax.plot(x+np.sqrt(variance), y+np.sqrt(variance), ':')
        
        '''
        ax.plot(x,slope_confidence_range[0]*x+intercept_confidence_range[0], ':')
        ax.plot(x,slope_confidence_range[0]*x+intercept_confidence_range[1], ':')
        ax.plot(x,slope_confidence_range[1]*x+intercept_confidence_range[0], ':')
        ax.plot(x,slope_confidence_range[1]*x+intercept_confidence_range[1], ':')
        '''
        
        #print len(dataset.speed_at_deceleration)
        
        slope = str(fit[0])[0:5]
        intercept = str(fit[1])[0:5]
        string = 'y = ' + slope + 'x + ' + intercept 
        ax.text(0, 0.5, string)
        string = 'Rsq=' + str(Rsq)[0:4]
        ax.text(0,0.4, string)
        
        if plot_sample_trajecs:
            keys_to_plot = ['2_29065', '2_31060', '8_10323', '6_715']
            for key in keys_to_plot:
                trajec = dataset.trajecs[key]
                ftp = np.arange(trajec.frames[0]-25, trajec.frames[-1]).tolist()
                ax.plot( (trajec.distance_subtended_by_post[ftp]), trajec.speed[ftp], color='black', linewidth=2)
                ax.plot( (trajec.distance_at_deceleration), trajec.speed_at_deceleration, '.', color='black', alpha=1)
        
        string = 'N=' + str(len(indices))
        ax.text( 10.*np.pi/180., 0.8, string, color='blue')
        if filename is not None:
            fig.savefig(filename, format='pdf')
    
    return x, y, yminus, yplus
    
    
    
def time_to_impact_vs_rrev(dataset, keys=None, plot=True, post_type=['checkered', 'checkered_angled', 'black', 'black_angled'], filename=None, xaxis='rrev', parameter='deceleration'):
    
    if keys is None:
        keys = dataset.trajecs.keys()
    
    indices = []
    if parameter == 'deceleration':
        trajecs_of_interest = dataset.trajectories_at_deceleration
    if parameter == 'saccade':
        trajecs_of_interest = dataset.trajectories_at_saccade
        
    for i, key in enumerate(trajecs_of_interest):
        if dataset.trajecs[key].post_type in post_type:
            if key in keys:
                if parameter == 'deceleration':
                    rrev = dataset.expansion_at_deceleration[i] / dataset.angle_at_deceleration[i]
                elif parameter == 'saccade':
                    rrev = dataset.expansion_at_saccade[i] / dataset.angle_at_saccade[i]
                if rrev > 0:
                    indices.append(i)
        else:
            print dataset.trajecs[key].post_type
    
    
    if plot:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        if xaxis == 'rrev':
            ax.set_autoscale_on(False)
            ax.set_xlim(-100,600)
            ax.set_ylim(-100,400)
            ax.set_aspect('equal')
        elif xaxis =='l_over_v':
            ax.set_autoscale_on(False)
            ax.set_xlim(-20,50)
            ax.set_ylim(-100,400)
            pass
        colormap_norm = matplotlib.colors.Normalize(0, .8, clip=True)
        cmap = plt.get_cmap('jet')
        
        if parameter == 'deceleration':
            l_over_v = 1000.*(dataset.distance_at_deceleration[indices]*np.tan(dataset.angle_at_deceleration[indices]/2.)) / dataset.speed_at_deceleration[indices]
            rrev_inv = 1000.*(dataset.angle_at_deceleration[indices]) / dataset.expansion_at_deceleration[indices]
            angle = 1000.*(dataset.angle_at_deceleration[indices])
            expansion = 1000*dataset.expansion_at_deceleration[indices]
        elif parameter == 'saccade':
            l_over_v = 1000.*(dataset.distance_at_saccade[indices]*np.tan(dataset.angle_at_saccade[indices]/2.)) / dataset.speed_at_deceleration[indices]
            rrev_inv = 1000.*(dataset.angle_at_saccade[indices]) / dataset.expansion_at_saccade[indices]
        
        
        if xaxis == 'rrev':
            x_param = expansion
        elif xaxis == 'l_over_v':
            x_param = l_over_v 
        
        if parameter == 'deceleration':
            tti = 1000*dataset.time_to_impact_at_deceleration[indices]
            speed = dataset.speed_at_deceleration[indices]
            keys = dataset.trajectories_at_deceleration[indices]
        elif parameter == 'saccade':
            tti = 1000*dataset.time_to_impact_at_saccade[indices]
            speed = dataset.speed_at_saccade[indices]
            keys = dataset.trajectories_at_saccade[indices]
        
        speed_threshold = 0.2
        if xaxis == 'rrev':
            indices_for_fit = np.where( (speed>speed_threshold)*(tti>0)*(tti<350)*(x_param<5000) )[0].tolist()
        elif xaxis == 'l_over_v':
            indices_for_fit = np.where( (speed>speed_threshold)*(tti>0)*(tti<350) )[0].tolist()
        #fit = np.polyfit(x_param[indices_for_fit], tti[indices_for_fit], 1)
        fit, variance, intercept_confidence_interval, slope_confidence_interval, Rsq = floris.linear_fit_type2(x_param[indices_for_fit], tti[indices_for_fit], full_output=True)
        x = np.linspace(-50, np.max(x_param[indices_for_fit]), 100)
        y = np.polyval(fit, x)
        
        if xaxis == 'rrev':
            indices_for_fit_for_low_vel = np.where( (speed<speed_threshold)*(tti>0)*(tti<350)*(x_param<4000) )[0].tolist()
        elif xaxis == 'l_over_v':
            indices_for_fit_for_low_vel = np.where( (speed<speed_threshold)*(tti>0)*(tti<200)*(x_param<50)*(x_param>0) )[0].tolist()
        fit_low_vel = np.polyfit(x_param[indices_for_fit_for_low_vel], tti[indices_for_fit_for_low_vel], 1)
        x_low_vel = np.linspace(-50, np.max(x_param[indices_for_fit_for_low_vel]), 100)
        y_low_vel = np.polyval(fit_low_vel, x_low_vel)
        
        print 'fit: ', fit
        
        
        if plot:
            for i in range(len(x_param)):
                c = cmap(colormap_norm(speed[i]))
                
                if speed[i] < speed_threshold:
                    '''
                    if xaxis == 'rrev':
                        pt = patches.Circle( (x_param[i], tti[i]), radius=5, facecolor=c, edgecolor='none', alpha=0.5)
                        ax.add_artist(pt)
                    else:
                    '''
                    ax.plot( x_param[i], tti[i], '.', markerfacecolor=c, markeredgecolor=c, alpha=0.5)
                else:
                    '''
                    if xaxis == 'rrev':
                        pt = patches.Circle( (x_param[i], tti[i]), radius=10, facecolor=c, edgecolor='none', alpha=0.5)
                        ax.add_artist(pt)
                    else:
                    '''
                    ax.plot( x_param[i], tti[i], '.', color=c, alpha=0.5)
            
            # high vel
            #if xaxis == 'rrev':
            ax.plot( x,y, color='blue')
            ax.vlines(0,-100,50,linestyle='dotted', linewidth=0.5)
            ax.hlines(fit[1],-100,25,linestyle='dotted', linewidth=0.5)
            string = 'delay: ' + str(fit[1])[0:3] + 'ms'
            ax.text(25,fit[1],string,horizontalalignment='left', verticalalignment='center')
            
            # low vel
            '''
            ax.plot( x_low_vel,y_low_vel, color='blue', linewidth=0.5)
            ax.vlines(0,-100,25,linestyle='dotted', linewidth=0.5)
            ax.hlines(fit_low_vel[1],-100,50,linestyle='dotted', linewidth=0.5)
            string = 'delay: ' + '+' + str(fit_low_vel[1])[0:3] + 'ms'
            ax.text(0,50,string,horizontalalignment='center', verticalalignment='center')
            '''
            
        adjust_spines(ax,['left', 'bottom'])
        
        if xaxis == 'rrev':
            ax.set_xlabel('angle / expansion, msec')
            xticks = np.arange(ax.get_xlim()[0],ax.get_xlim()[1],100).tolist()
            ax.set_xticks(xticks)
            yticks = np.arange(-100,401,100).tolist()
            ax.set_yticks(yticks)
        elif xaxis == 'l_over_v':
            ax.set_xlabel('half angle / speed, msec')
            
        if parameter == 'deceleration':
            ax.set_ylabel('time to impact at deceleration, msec')
        elif parameter == 'saccade':
            ax.set_ylabel('time to impact at saccade, msec')
        
        if filename is not None:
            #fig.subplots_adjust(bottom=0.2, top=0.8, right=0.8, left=0.2)
            fig.savefig(filename, format='pdf')
    
        return x, y
        
def time_to_impact_vs_rrev_landing(dataset, keys=None, plot=True, post_type=['checkered', 'checkered_angled', 'black', 'black_angled'], filename=None):
    
    if keys is None:
        keys = dataset.trajecs.keys()
    
    indices = []
    trajecs_of_interest = dataset.trajectories_at_deceleration
        
    for i, key in enumerate(trajecs_of_interest):
        if dataset.trajecs[key].post_type in post_type:
            if key in keys:
                rrev = dataset.expansion_at_deceleration[i] / dataset.angle_at_deceleration[i]
                if rrev > 0:
                    indices.append(i)
        else:
            print dataset.trajecs[key].post_type
    
    
    if plot:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_autoscale_on(False)
        ax.set_xlim(-50,50)
        ax.set_ylim(-50,100)
        #ax.set_aspect('equal')
        colormap_norm = matplotlib.colors.Normalize(0, .8, clip=True)
        cmap = plt.get_cmap('jet')
        
        rrev_inv = 1000.*(dataset.angle_at_deceleration[indices]) / dataset.expansion_at_deceleration[indices]
        l_over_v = 1000.*(dataset.distance_at_deceleration[indices]*np.tan(dataset.angle_at_deceleration[indices]/2.)) / dataset.speed_at_deceleration[indices]
        angle = 1000.*(dataset.angle_at_deceleration[indices])
        expansion = 1000*dataset.expansion_at_deceleration[indices]
        x_param = l_over_v
        
        
        tti = 1000*dataset.time_to_impact_at_deceleration[indices]
        speed = dataset.speed_at_deceleration[indices]
        keys = dataset.trajectories_at_deceleration[indices]
        
        speed_threshold = 0.2
        indices_for_fit = np.where( (speed>speed_threshold)*(tti>0)*(tti<350)*(x_param<500) )[0].tolist()
        fit = np.polyfit(x_param[indices_for_fit], tti[indices_for_fit], 1)
        x = np.linspace(-50, np.max(x_param[indices_for_fit]), 100)
        y = np.polyval(fit, x)
        
        indices_for_fit_for_low_vel = np.where( (speed<speed_threshold)*(tti>0)*(tti<350)*(x_param<4000) )[0].tolist()
        fit_low_vel = np.polyfit(x_param[indices_for_fit_for_low_vel], tti[indices_for_fit_for_low_vel], 1)
        x_low_vel = np.linspace(-50, np.max(x_param[indices_for_fit_for_low_vel]), 100)
        y_low_vel = np.polyval(fit_low_vel, x_low_vel)
        
        if plot:
            for i in range(len(x_param)):
                c = cmap(colormap_norm(speed[i]))
                
                if speed[i] < speed_threshold:
                    '''
                    pt = patches.Circle( (x_param[i], tti[i]), radius=5, facecolor=c, edgecolor='none', alpha=0.5)
                    ax.add_artist(pt)
                    '''
                    ax.plot( x_param[i], tti[i], '.', markerfacecolor=c, markeredgecolor=c, alpha=0.5)
                else:
                    '''
                    pt = patches.Circle( (x_param[i], tti[i]), radius=10, facecolor=c, edgecolor='none', alpha=0.5)
                    ax.add_artist(pt)
                    '''
                    ax.plot( x_param[i], tti[i], '.', color=c, alpha=0.5)
            
            # high vel
            #if xaxis == 'rrev':
            ax.plot( x,y, color='blue')
            ax.vlines(0,-100,50,linestyle='dotted', linewidth=0.5)
            ax.hlines(fit[1],-100,25,linestyle='dotted', linewidth=0.5)
            string = 'delay: ' + str(fit[1])[0:3] + 'ms'
            ax.text(25,fit[1],string,horizontalalignment='left', verticalalignment='center')
            
            # low vel
            '''
            ax.plot( x_low_vel,y_low_vel, color='blue', linewidth=0.5)
            ax.vlines(0,-100,25,linestyle='dotted', linewidth=0.5)
            ax.hlines(fit_low_vel[1],-100,50,linestyle='dotted', linewidth=0.5)
            string = 'delay: ' + '+' + str(fit_low_vel[1])[0:3] + 'ms'
            ax.text(0,50,string,horizontalalignment='center', verticalalignment='center')
            '''
            
        ## plot a constant velocity (fake) trajectory
        vels = [0.05, 0.2, 0.4, 0.6, 0.8]
        for vel in vels:
            fps = 1000.
            dist = np.arange(.3, 0.00001, -1/fps*vel)
            radius = 0.009565
            angle = 2*np.arcsin(radius/(dist+radius))
            expansion = sa1.diffa(angle)*fps
            rrev_inv = angle / expansion
            l_over_v = dist*np.tan(angle/2.) / vel
            tti = dist / vel
            ax.plot(1000*l_over_v, 1000*tti, '--', color=cmap(colormap_norm(vel))) 
        
        adjust_spines(ax,['left', 'bottom'])
        
        ax.set_xlabel('angle / expansion, msec')
        xticks = np.arange(ax.get_xlim()[0],ax.get_xlim()[1],100).tolist()
        #ax.set_xticks(xticks)
        yticks = np.arange(-100,401,100).tolist()
        #ax.set_yticks(yticks)
            
        ax.set_ylabel('time to impact at deceleration, msec')
        
        #fig.subplots_adjust(bottom=0.2, top=0.8, right=0.8, left=0.2)
        fig.savefig('tti_vs_rrev_inv_landing_deceleration.pdf', format='pdf')
    
        return x, y
    
def time_to_impact_vs_rrev_flyby(dataset, keys=None, plot=True, post_type=['checkered', 'checkered_angled', 'black', 'black_angled'], filename=None):
    
    if keys is None:
        classified_keys = get_classified_keys(dataset)
        keys = classified_keys['straight']
    
    # correct for radial variations and round object visual distortion (sin vs tan dealio)
    radius=dataset.radius_at_saccade
    phi = np.pi/2. - dataset.angle_at_saccade/2.
    l = np.sin(phi)*radius
    speed = dataset.speed_at_saccade
    
    indices = []
    trajecs_of_interest = dataset.trajectories_at_saccade
        
    for i, key in enumerate(trajecs_of_interest):
        if dataset.trajecs[key].post_type in post_type:
            if key in keys:
                if l[i] > 0:
                    indices.append(i)
        else:
            pass#print dataset.trajecs[key].post_type
    
    

    if plot:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        colormap_norm = matplotlib.colors.Normalize(0, 0.02, clip=True)
        cmap = plt.get_cmap('jet')
        #ax.set_ylim(-10,400)
        #ax.set_xlim(-10,200)
        #ax.set_autoscale_on(False)
        
        # get average speed before saccade:
        
        l_over_v = 1000.*l / speed
        x_param = l_over_v
        
        
        tti = 1000*dataset.time_to_impact_at_saccade
        
        indices_to_fit = np.where( (tti[indices]<500)*(x_param[indices]<200) )[0].tolist()
        
        key_list = dataset.trajectories_at_saccade[indices].tolist()
        trajec_list = [dataset.trajecs[ key ] for key in key_list]
        post_types = [trajec.post_type for trajec in trajec_list]
        indices_for_black = []
        indices_for_checkered = []
        for i, p in enumerate(post_types):
            if 'black' in p and indices[i] in indices_to_fit:
                indices_for_black.append( indices[i] )
            elif 'checkered' in p and indices[i] in indices_to_fit:
                indices_for_checkered.append( indices[i] )
        
        fit_black = np.polyfit(x_param[indices_for_black], tti[indices_for_black], 1)
        x_black = np.linspace(0, 150, 100)
        y_black = np.polyval(fit_black, x_black)
        print 'black: ', fit_black
        print 2*np.arctan2(1,fit_black[0])*180/np.pi, 2*np.arcsin(1/fit_black[0])*180/np.pi
        ax.plot( x_param[indices_for_black], tti[indices_for_black], '.', color='black', alpha=0.5)
        ax.plot( x_black,y_black, '-', color='black')
        
        fit_checkered = np.polyfit(x_param[indices_for_checkered], tti[indices_for_checkered], 1)
        x_checkered = np.linspace(0, 150, 100)
        y_checkered = np.polyval(fit_checkered, x_checkered)
        print 'checkered: ', fit_checkered
        print 2*np.arctan2(1,fit_checkered[0])*180/np.pi, 2*np.arcsin(1/fit_checkered[0])*180/np.pi
        
        for i in indices_for_checkered:
            c = cmap(colormap_norm(speed[i]))
            ax.plot( x_param[i], tti[i], '.', markerfacecolor='white', markeredgecolor=c, alpha=0.5)
        ax.plot( x_checkered,y_checkered, '-', color='black')
        
        
        #ax.vlines(0,-100,50,linestyle='dotted', linewidth=0.5)
        #ax.hlines(fit[1],-100,25,linestyle='dotted', linewidth=0.5)
        #string = 'delay: ' + str(fit[1])[0:3] + 'ms'
        #ax.text(25,fit[1],string,horizontalalignment='left', verticalalignment='center')
        # low vel
        '''
        ax.plot( x_low_vel,y_low_vel, color='blue', linewidth=0.5)
        ax.vlines(0,-100,25,linestyle='dotted', linewidth=0.5)
        ax.hlines(fit_low_vel[1],-100,50,linestyle='dotted', linewidth=0.5)
        string = 'delay: ' + '+' + str(fit_low_vel[1])[0:3] + 'ms'
        ax.text(0,50,string,horizontalalignment='center', verticalalignment='center')
        '''
            
        adjust_spines(ax,['left', 'bottom'])
        
        ax.set_xlabel('l/v (radius over speed), msec')
        #xticks = np.arange(ax.get_xlim()[0],ax.get_xlim()[1],100).tolist()
        #ax.set_xticks(xticks)
        #yticks = np.arange(-100,401,100).tolist()
        #ax.set_yticks(yticks)
            
        ax.set_ylabel('time to impact at saccade, msec')
        
        #fig.subplots_adjust(bottom=0.2, top=0.8, right=0.8, left=0.2)
        fig.savefig('tti_vs_rrev_inv_flyby_saccade.pdf', format='pdf')
    
    
def angle_vs_expansion(dataset): 

    fig = plt.figure()
    ax = fig.add_subplot(111)
    colormap_norm = matplotlib.colors.Normalize(0, .8, clip=True)
    cmap = plt.get_cmap('jet')
        
    for i, key in enumerate(dataset.trajectories_at_deceleration): 
        c = cmap(colormap_norm(dataset.speed_at_deceleration[i]))
        ax.plot( dataset.angle_at_deceleration[i], dataset.expansion_at_deceleration[i], '.', color=c, alpha=0.5)
    fig.savefig('angle_vs_expansion_landing_deceleration.pdf', format='pdf')
    
    
    
def fix_angle_log_spine(ax, set_y=True, histograms=True):

    if histograms:
        yticklabels = ['0.0', '0.2', '0.4', '0.6', '0.8', '1.0']
        probax = ax.twinx()
        adjust_spines(probax, ['right'])
        probax.set_yticklabels(yticklabels, fontsize=9)
        probax.set_ylabel('probability, normalized', fontsize=9)
    deg_ticks = np.array([5, 10, 30, 60, 90, 180])
    deg_tick_strings = [str(d) for d in deg_ticks]
    rad_ticks = deg_ticks*np.pi/180.
    rad_ticks_log = np.log(rad_ticks)
    
    dist_tick_strings = ['(21)', '(10)', '(2.7)', '(0.9)', '(0.4)', '(0)']
    x_tick_strings = []
    for i, d in enumerate(dist_tick_strings):
        x_tick_strings.append( deg_tick_strings[i] + '\n' + dist_tick_strings[i] )
    
    ## plot paramters    
    if set_y:
        ax.set_ylim([0,1.])
    ax.set_xlim(rad_ticks_log[0], rad_ticks_log[-1])
    
    if histograms:
        adjust_spines(ax, ['left', 'bottom', 'right'], color={'bottom':'black'})
    else:
        adjust_spines(ax, ['left', 'bottom'], color={'bottom':'black'})
    ax.set_xlabel('retinal size, deg\n(distance, cm)', labelpad=10)
    ax.set_ylabel('speed, m/s')
    if set_y:
        yticklabels = ['0.0', '0.2', '0.4', '0.6', '0.8', '1.0']
        ax.set_yticklabels(yticklabels)
    
    ax.set_xticks(rad_ticks_log.tolist())
    ax.set_xticklabels(x_tick_strings) 
    
    

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
                pass#spine.set_smart_bounds(True)
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
        
    for line in ax.get_xticklines() + ax.get_yticklines():
        #line.set_markersize(6)
        line.set_markeredgewidth(1)

    
    
def calc_func(dataset, func):
    
    for k, trajec in dataset.trajecs.items():
        func(trajec)
        
def deceleration_plot(dataset, keys=None):    
    
    if keys is None:
        keys = dataset.trajecs.keys() 

    fig = plt.figure()
    ax = fig.add_subplot(111)

    for key in keys:
        trajec = dataset.trajecs[key]
        frames = np.arange(trajec.frame_at_deceleration, trajec.frame_of_landing).tolist()
        ax.plot(np.log(trajec.angle_subtended_by_post[frames]), trajec.speed[frames], color='blue', alpha=0.3)
        
    ax.plot(np.log(dataset.angle_at_deceleration), dataset.speed_at_deceleration, '.', color='blue')
    fix_angle_log_spine(ax, histograms=False)
    fig.savefig('deceleration_plot.pdf', format='pdf')
    
def landing_saccades_plot(dataset, keys=None):    
    
    if keys is None:
        keys = dataset.trajecs.keys() 

    fig = plt.figure()
    ax = fig.add_subplot(111)

    for key in keys:
        trajec = dataset.trajecs[key]
        #frames = np.arange(trajec.frame_at_deceleration, trajec.frame_of_landing).tolist()
        #ax.plot(np.log(trajec.angle_subtended_by_post[frames]), trajec.speed[frames], color='blue', alpha=0.3)
        
        if len(trajec.saccades) == 1:
            s = trajec.saccades[-1]
            ax.plot(  trajec.positions[s,0], trajec.positions[s,1], '.', color='green')

    fig.savefig('landing_saccades_plot.pdf', format='pdf')
    
    
def landing_spagetti_plots(dataset, keys=None):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    
    
    if keys is None:
        classified_keys = get_classified_keys(dataset)
        keys = classified_keys['straight']
        #keys = dataset.trajecs.keys()
        
    for key in keys:
        trajec = dataset.trajecs[key]
        ftp = np.arange(trajec.frames[0]-25, trajec.frames[-1]).tolist()
        ax.plot( np.log(trajec.angle_subtended_by_post[ftp]), trajec.speed[ftp], color='black', linewidth=1, alpha=0.05)

    keys_to_plot = ['2_29065', '2_31060', '8_10323', '6_715']
    for key in keys_to_plot:
        trajec = dataset.trajecs[key]
        ftp = np.arange(trajec.frames[0]-25, trajec.frames[-1]).tolist()
        ax.plot( np.log(trajec.angle_subtended_by_post[ftp]), trajec.speed[ftp], color='black', linewidth=1)
        ax.plot( np.log(trajec.angle_at_deceleration), trajec.speed_at_deceleration, '.', color='blue', alpha=1)
    fix_angle_log_spine(ax, histograms=False)
    fig.savefig('landing_spagetti.pdf', format='pdf')
    
def landing_spagetti_plots_nolog(dataset, keys=None):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    colormap_norm = matplotlib.colors.Normalize(0, .8, clip=True)
    cmap = plt.get_cmap('jet')
        
    if keys is None:
        classified_keys = get_classified_keys(dataset)
        keys = classified_keys['straight']
        #keys = dataset.trajecs.keys()
        
    for key in keys:
        trajec = dataset.trajecs[key]
        ftp = np.arange(trajec.frames[0]-25, trajec.frames[-1]).tolist()
        #ax.plot( trajec.angle_subtended_by_post[ftp], trajec.speed[ftp], color='black', linewidth=1, alpha=0.05)
        if trajec.frame_of_landing - trajec.frame_at_deceleration > 5:
            #dvda = floris.diffa( trajec.speed[trajec.frame_at_deceleration:trajec.frame_of_landing] ) / floris.diffa( (trajec.angle_subtended_by_post[trajec.frame_at_deceleration:trajec.frame_of_landing]) )
            c = cmap(colormap_norm(trajec.speed_at_deceleration))
            #ax.plot(trajec.angle_subtended_by_post[trajec.frame_at_deceleration:trajec.frame_of_landing], dvda, color=c, alpha=0.2)
            dvdd = floris.diffa( trajec.speed[trajec.frame_at_deceleration:trajec.frame_of_landing] ) / floris.diffa( (trajec.dist_to_stim_r_normed[trajec.frame_at_deceleration:trajec.frame_of_landing]) )
            ax.plot(trajec.dist_to_stim_r_normed[trajec.frame_at_deceleration:trajec.frame_of_landing], dvdd, color=c, alpha=0.2)
    '''
    keys_to_plot = ['2_29065', '2_31060', '8_10323', '6_715']
    for key in keys_to_plot:
        trajec = dataset.trajecs[key]
        ftp = np.arange(trajec.frames[0]-25, trajec.frames[-1]).tolist()
        ax.plot( trajec.angle_subtended_by_post[ftp], trajec.speed[ftp], color='black', linewidth=1)
        ax.plot( trajec.angle_at_deceleration, trajec.speed_at_deceleration, '.', color='blue', alpha=1)
    '''
    #ax.set_ylim(-1,0)
    ax.set_ylim(0,50)
    fig.savefig('landing_spagetti_nolog.pdf', format='pdf')
    
def flyby_spagetti_plots(dataset, keys=None):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    if keys is None:
        classified_keys = get_classified_keys(dataset)
        keys = classified_keys['straight']
    
    plotted_trajecs = []
    for key in keys:
        trajec = dataset.trajecs[key]
        if len(trajec.saccades) > 0:
            ftp = np.arange(0, trajec.frames[-1]).tolist()
            ax.plot( np.log(trajec.angle_subtended_by_post[ftp]), trajec.speed[ftp], color='black', linewidth=1, alpha=0.05)
            plotted_trajecs.append(key)
            
    k = [40, 56]
    ks = [plotted_trajecs[i] for i in k]
    for key in ks:
        trajec = dataset.trajecs[key]
        ax.plot( np.log(trajec.angle_subtended_by_post[trajec.frames_of_flyby]), trajec.speed[trajec.frames_of_flyby], color='black', linewidth=1)
        if len(trajec.saccades) > 0:
            saccade_frames = get_saccade_range(trajec, trajec.saccades[-1])
            ax.plot( np.log(trajec.angle_subtended_by_post[saccade_frames]), trajec.speed[saccade_frames], color='green', linewidth=1)  
            ax.plot( np.log(trajec.angle_subtended_by_post[saccade_frames[0]]), trajec.speed[saccade_frames[0]], '.', color='green')
            
            
        if trajec.frame_at_deceleration is not None:
            ax.plot( np.log(trajec.angle_subtended_by_post[trajec.frame_at_deceleration]), trajec.speed[trajec.frame_at_deceleration], '.', color='blue')   
        
    fix_angle_log_spine(ax, histograms=False)
    fig.savefig('flyby_spagetti.pdf', format='pdf')
    
def saccade_deceleration_relation(dataset, keys=None):
    
    if keys is None:
        classified_keys = get_classified_keys(dataset)
        keys = classified_keys['straight']
    
    angle_at_saccade = []
    angle_at_deceleration = []
    speed_at_saccade = []
    speed_at_deceleration = []
    time = []
    for key in keys:
        trajec = dataset.trajecs[key]
        if trajec.frame_at_deceleration is not None:
            if len(trajec.saccades) > 0:
                #s = trajec.saccades[-1]
                saccade_frames = get_saccade_range(trajec, trajec.saccades[-1])
                s = saccade_frames[0]
                d = trajec.frame_at_deceleration
                
                angle_at_saccade.append( trajec.angle_subtended_by_post[s] )
                angle_at_deceleration.append( trajec.angle_subtended_by_post[d] )
                speed_at_saccade.append( trajec.speed[s] )
                speed_at_deceleration.append( trajec.speed[d] )
                time.append( trajec.epoch_time[s] - trajec.epoch_time[d] )
                
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    ax.plot( speed_at_saccade, speed_at_deceleration, '.')
    ax.set_xlabel('speed at saccade')
    ax.set_ylabel('speed at deceleration')
    adjust_spines(ax, ['left', 'bottom'])
    fig.savefig('flyby_speed_at_saccade_vs_deceleration.pdf', format='pdf')
    
    
    # time plot
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    data, bins = np.histogram( time, bins=16, normed=True)
    xvals = np.diff(bins) + bins[0:-1]
    butter_b, butter_a = signal.butter(3, 0.3)
    data_filtered = signal.filtfilt(butter_b, butter_a, data)
    
    axhist = ax.twinx()
    axhist.hist( time, bins=bins, normed=True, facecolor='blue', alpha=0.1, edgecolor='blue')
    axhist.plot( xvals, data_filtered, color='blue' )
    
    ax.plot( time, speed_at_deceleration, '.', color='blue')
    adjust_spines(axhist, ['right'])
    adjust_spines(ax, ['left', 'bottom', 'right'])
    axhist.set_ylabel('probability, normalized')
    ax.set_ylabel('speed at deceleration, m/s')
    ax.set_xlabel('time btwn deceleration and saccade, sec')
    fig.savefig('flyby_time_vs_deceleration.pdf', format='pdf')
    
    
        
def make_figs(dataset_landing, dataset_flyby, movie_dataset=None):
    
    if movie_dataset is not None:
        summary_figure(dataset_landing, dataset_flyby, movie_dataset)
    
    landing_spagetti_plots(dataset_landing)
    flyby_spagetti_plots(dataset_flyby, keys=None)
    
    get_angle_vs_speed_curve(dataset_landing, plot=True, post_type=['checkered', 'checkered_angled'], filename='angle_vs_speed_curve_checkered.pdf')
    get_angle_vs_speed_curve(dataset_landing, plot=True, post_type=['black', 'black_angled'], filename='angle_vs_speed_curve_black.pdf')
    get_angle_vs_speed_curve(dataset_landing, plot=True, filename='angle_vs_speed_curve.pdf', plot_sample_trajecs=False)
    
    saccade_angle_histogram(dataset_flyby, plot=True, post_type=['black','black_angled'], filename='saccade_angle_histogram_black.pdf')
    saccade_angle_histogram(dataset_flyby, plot=True, post_type=['checkered','checkered_angled'], filename='saccade_angle_histogram_checkered.pdf')
    saccade_angle_histogram(dataset_flyby, plot=True, plot_sample_trajecs=False, filename='saccade_angle_histogram.pdf')
    
    saccade_angle_histogram(dataset_landing, plot=True, plot_sample_trajecs=False, filename='saccade_angle_histogram_landing.pdf')
    
    if movie_dataset is not None:
        leg_extension_angle_histogram(movie_dataset, plot=True, behavior='landing')
    deceleration_angle_histogram_flyby(dataset_flyby, plot=True, saccades=True)
    deceleration_angle_histogram_flyby(dataset_flyby, plot=True, saccades=False)
    deceleration_angle_histogram_flyby(dataset_flyby, plot=True, saccades=None)
    
    saccade_deceleration_relation(dataset_flyby, keys=None)
    print 'done'
    return
    
def make_deceleration_angle_histograms_flyby(dataset):
    
    deceleration_angle_histogram_flyby(dataset, keys=None, plot=True, saccades=None, filename='flyby_deceleration_angle_histogram_10deg.pdf', angle_to_post_range=[0,10*np.pi/180.])
    
    deceleration_angle_histogram_flyby(dataset, keys=None, plot=True, saccades=None, filename='flyby_deceleration_angle_histogram_20deg.pdf', angle_to_post_range=[10*np.pi/180.,20*np.pi/180.])
    
    deceleration_angle_histogram_flyby(dataset, keys=None, plot=True, saccades=None, filename='flyby_deceleration_angle_histogram_30deg.pdf', angle_to_post_range=[20*np.pi/180.,30*np.pi/180.])
    
    
def make_landing_angle_vs_speed_figs(dataset_landing):
    classified_keys = get_classified_keys(dataset_landing)
    for key_class in classified_keys.keys():
        filename = 'landing_angle_vs_speed_' + key_class + '.pdf'
        get_angle_vs_speed_curve(dataset_landing, plot=True, filename=filename, keys=classified_keys[key_class])
def make_flyby_angle_vs_speed_figs(dataset_flyby):
    classified_keys = get_classified_keys(dataset_flyby)
    for key_class in classified_keys.keys():
        filename = 'flyby_angle_vs_speed_' + key_class + '.pdf'
        deceleration_angle_histogram_flyby(dataset_flyby, plot=True, saccades=True, filename=filename, keys=classified_keys[key_class])
    
    
def get_saccade_range(trajec, saccade_frame=None, plot=False):
    
    try:
        
        first_frame = np.max([0, saccade_frame-15])
        last_frame = np.min([saccade_frame+15, len(trajec.speed)-1])
        frames = np.arange(first_frame, last_frame).tolist()
        angular_vel = np.abs( sa1.diffa( trajec.heading[frames] ) )*100.*180/np.pi
        saccade_frames_raw = angular_vel > 350
        
        saccade_frames_blob = nim.find_blob_nearest_to_point(saccade_frames_raw, saccade_frame-first_frame)
        saccade_frames = np.where(saccade_frames_blob == 1)[0].tolist()
        
        if plot:
            print saccade_frames
            
            fig = plt.figure()
            ax = fig.add_subplot(111)
            
            ax.plot( frames, np.abs( sa1.diffa( trajec.heading[frames] ) ), color='blue' )
            ax.plot( np.array(saccade_frames)+frames[0], np.abs( sa1.diffa( trajec.heading[np.array(saccade_frames)+frames[0]] ) ), color='green' )
            fig.savefig('test.pdf', format='pdf')
            
        return (np.array(saccade_frames)+frames[0]).tolist()
        
    except:
        print 'excepted'
        return [saccade_frame, saccade_frame+1]
    
    
def histogram_saccade_length(dataset, keys=None):
    if keys is None:
        classified_keys = get_classified_keys(dataset)
        #keys = classified_keys['straight']
        keys = dataset.trajecs.keys()
        
    saccade_lengths = []
    for key in keys:
        trajec = dataset.trajecs[key]
        
        if len(trajec.saccades) > 0:
            saccade_frames = get_saccade_range(trajec, trajec.saccades[-1])
            saccade_lengths.append( float(len(saccade_frames))/100. )
    
    print 'saccade mean length: ', np.mean(saccade_lengths)
    
    fig = plt.figure()
    ax = fig.add_subplot(111) 
    ax.hist(1000*np.array(saccade_lengths), bins=30, edgecolor='none', facecolor='green')
    adjust_spines(ax, ['left', 'bottom'])
    ax.set_xlabel('saccade duration, ms')
    ax.set_ylabel('occurences')
    fig.savefig('saccade_lengths.pdf', format='pdf')
                
                
    
def get_max_deceleration(dataset, post_type=['checkered', 'checkered_angled', 'black', 'black_angled'], filename=None, keys=None):
    angleok = np.where( (dataset.angle_at_deceleration < 2)*(dataset.angle_at_deceleration > .01) )[0].tolist()
    if keys is None:
        classified_keys = get_classified_keys(dataset)
        keys = classified_keys['straight']
    
    indices = []
    for i, key in enumerate(dataset.trajectories_at_deceleration):
        if key in keys:
            if i in angleok:
                if dataset.trajecs[key].post_type in post_type:
                    indices.append(i)
                else:
                    pass
    
    decelerations = []
            
    for i in dataset.trajectories_at_deceleration[indices]:
        trajec = dataset.trajecs[i]
        if trajec.frame_of_landing - trajec.frame_at_deceleration > 5:
            decelerations.append( np.max( np.abs(floris.diffa( trajec.speed[trajec.frame_at_deceleration:trajec.frame_of_landing] ) / floris.diffa( np.log(trajec.angle_subtended_by_post[trajec.frame_at_deceleration:trajec.frame_of_landing]) ))) )
            
    for i, d in enumerate(decelerations):
        if np.isnan(d) or np.isinf(d):
            print d
            del(decelerations[i])
                
    fig = plt.figure()
    ax = fig.add_subplot(111)
    bins = np.linspace(0,1,100)
    ax.hist(np.array(decelerations), bins=bins)
    fig.savefig('decelerations.pdf', format='pdf')
    
    return decelerations
    
        
            
    
    
    
def remove_duplicate_trajectories(dataset, delete=False):
    
    keys = dataset.trajecs.keys()
    deleted_keys = []
    for i, key in enumerate(keys[0:-2]):
        if key not in deleted_keys:
            master_trajec = dataset.trajecs[key]
            tm = master_trajec.epoch_time[0]
            
            for j in range(i+1,len(keys)-1):
                k = keys[j]
                if k not in deleted_keys:
                    test_trajec = dataset.trajecs[k]
                    terr = np.abs(tm - test_trajec.epoch_time[0])
                    if terr < 0.001:
                        print key, k
                        if delete:
                            del(dataset.trajecs[k])
                            deleted_keys.append(k)
            
def fix_sawyers_callibration(dataset):
    keys = dataset.trajecs.keys()
    for key in keys:
        trajec = dataset.trajecs[key]
        trajec.positions[:,2] -= 0.15
    
    
def print_test(dataset):
    keys = dataset.trajecs.keys()
    for key in keys:
        trajec = dataset.trajecs[key]
        if np.max(trajec.positions[:,2]) > 0:
            print key
    
