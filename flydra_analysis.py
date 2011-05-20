import numpy as np

import pickle
import copy
import sys
sys.path.insert(0, '/usr/local/lib/python2.6/dist-packages')
sys.path.append('/home/floris/src/pymovie2')

import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf as pdf

import analysis_plot as ap
import sa1_analysis as sa1 
import flydra_floris_analysis as ffa

REQUIRED_LENGTH = 30
REQUIRED_DIST = 0.08


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
    
    
def find_critical_points(arr, window=5):
    arr_diff = diff_windowed(arr, window)
    return np.where( diff_windowed( np.sign(arr_diff), 1 ) != 0 )[0].tolist()

def cumsum(array):
    cumarr = copy.copy(array)
    for i in range(1, len(cumarr)):
        cumarr[i] += cumarr[i-1]
    return cumarr

def load(filename, prep=True):
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
    
def make_behavior_dataset(dataset, filename='landing_dataset_8cm', behavior='landing'):
    new_dataset = ffa.Dataset(like=dataset)
    if type(behavior) is not list:
        behavior = [behavior]
    for k,trajec in dataset.trajecs.items():
        if trajec.behavior in behavior:
            
            calc_frame_of_landing (trajec, threshold = 0.0005)
            normalize_dist_to_stim_r(trajec)
            
            if trajec.frame_of_landing > REQUIRED_LENGTH:
                if trajec.dist_to_stim_r_normed[0] >= REQUIRED_DIST:
                
                    classify(trajec, dfar=REQUIRED_DIST, dnear=0.005)
                    new_dataset.trajecs.setdefault(k, trajec)
                    
    save(new_dataset, filename)

    return new_dataset
    
def add_datasets(dataset_list):
    
    new_dataset = dataset_list[0]
    
    for d, dataset in enumerate(dataset_list[1:]):
        for k,trajec in dataset.trajecs.items():
            if k not in new_dataset.trajecs.keys():
                new_dataset.trajecs.setdefault(k, trajec)
            else:
                print 'duplicate id: ', k
                
    return new_dataset
    

def prep_dataset(dataset):
    
    for k,trajec in dataset.trajecs.items():
        classify(trajec, dfar=REQUIRED_DIST, dnear=0.005)
        calc_heading(trajec)
        calc_heading_cumsum(trajec)
        trajec.frames = np.arange(get_frame_at_distance(trajec, REQUIRED_DIST), trajec.frame_of_landing).tolist()
        calc_time_to_impact(trajec)
        trajec.key = k
        calc_saccades(trajec)
###
def normalize_dist_to_stim_r(trajec):
    trajec.dist_to_stim_r_normed = trajec.dist_to_stim_r - trajec.dist_to_stim_r[trajec.frame_of_landing]
def get_frame_at_distance(trajec, distance):
    normalize_dist_to_stim_r(trajec)
    frames = np.arange(0, trajec.frame_of_landing).tolist()
    dist_to_post = trajec.dist_to_stim_r_normed[frames]
    dist_crossovers = np.where( sa1.diffa(np.sign(dist_to_post - distance)) != 0 )[0]
    frame = dist_crossovers[-1]
    return frame
def get_speed_at_distance(trajec, distance):
    frame = get_frame_at_distance(trajec, distance)
    speed = trajec.speed[frame]
    return speed
def classify(trajec, dfar=REQUIRED_DIST, dnear=0.005):
    speed_hi_threshold = 0.19
    speed_lo_threshold = 0.15

    if trajec.behavior == 'landing' and trajec.dist_to_stim_r_normed[0] >= REQUIRED_DIST:
        trajec.speed_far = get_speed_at_distance(trajec, dfar)
        trajec.speed_near = get_speed_at_distance(trajec, dnear)
        
        if trajec.speed_near > speed_hi_threshold:
            trajec.classification = 'fast'
        elif trajec.speed_far < speed_lo_threshold:
            trajec.classification = 'slow'
        else:
            trajec.classification = 'mid'
def calc_frame_of_landing (trajec, threshold = 0.0005):
        # search forward in time until we get below threshold, then check to make sure fly stays below threshold for three more frames
        # if it's a flyby, find nearest point to post?
        if trajec.behavior == 'flyby':
            trajec.frame_of_landing = np.argmin(trajec.dist_to_stim)
            return trajec.frame_of_landing
        if trajec.behavior == 'landing':
            diffdist = np.abs(sa1.diffa(trajec.dist_to_stim_r))
            print 'calculating for landing'
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
                        print 'frame of landing: ', frame_of_landing
                        return frame_of_landing
                    elif diffdist[i] < threshold:
                        counter = counter + 1
                    else:
                        counter = 0
                        frame_of_landing = 0
                        
            trajec.frame_of_landing = -1
            print 'frame of landing: ', frame_of_landing
            return -1
        else:
            trajec.frame_of_landing = -1
            return -1
            
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
            print change_diff
            if np.abs(change_diff > .2):
                    trajec.heading_cumsum_pos_pts_of_inflection.append(pt)
                    
    
    if plot:
        plt.plot(trajec.heading_cumsum)
        plt.plot(trajec.heading_cumsum_pos_pts_of_inflection, trajec.heading_cumsum[trajec.heading_cumsum_pos_pts_of_inflection], 'o')
        plt.plot(trajec.heading_cumsum_diff)
        
        
def calc_time_to_impact(trajec):
    trajec.time_to_impact = trajec.dist_to_stim_r_normed / trajec.speed
    
    
def calc_saccades(trajec, magnitude=0.1):
    raw_saccades = find_critical_points(trajec.heading_smooth_diff)
    
    windowed_heading_diff = diff_windowed(trajec.heading, 5)
    
    trajec.saccades = []
    for saccade in raw_saccades:
        print np.abs(windowed_heading_diff[saccade])
        if np.abs(windowed_heading_diff[saccade]) > magnitude:
            trajec.saccades.append(saccade)
    

################################################################################################
## PLOTTING
################################################################################################

###
def pdf_flydra_trajecs_of_classifications(dataset, filename='classified_flydra_trajectories'):
    classified_keys = get_classified_keys(dataset)
    
    page = -1
    plt.close('all')
    pp =  pdf.PdfPages(filename)
    
    for classification, keys in classified_keys.items():
        cl = ap.xy_trajectories(dataset, trajectory=keys, show_saccades=True, trajec_alpha=0.1)
        title = str(classification) + ' course'
        cl.ax0.set_title(title)
        cl.ax0.figure.set_size_inches(2*10,1*10)
        cl.ax0.figure.set_dpi(72)
        pp.savefig()
        plt.close('all')
    
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
    
    
    
    
    
    
    
    
    
    
    
    
    
