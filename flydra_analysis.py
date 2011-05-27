import numpy as np

import time
import datetime

import pickle
import copy
import sys
sys.path.insert(0, '/usr/local/lib/python2.6/dist-packages')
sys.path.append('/home/floris/src/pymovie2')

import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf as pdf
import matplotlib
import matplotlib.colorbar
from mpl_toolkits.axes_grid1 import make_axes_locatable

import analysis_plot as ap
import sa1_analysis as sa1 
import numpyimgproc as nim
import colorline
import flydra_floris_analysis as ffa

REQUIRED_LENGTH = 30
REQUIRED_DIST = 0.08

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
        trajec.key = k
        if trajec.behavior in behavior:
            
            calc_frame_of_landing (trajec, threshold = 0.0005)
            normalize_dist_to_stim_r(trajec)
            trajec.frames = np.arange(get_frame_at_distance(trajec, REQUIRED_DIST), trajec.frame_of_landing).tolist()
            
            if trajec.frame_of_landing > REQUIRED_LENGTH:
                if trajec.dist_to_stim_r_normed[0] >= REQUIRED_DIST:
                    if np.max(trajec.positions[trajec.frames,2]) < 0:
                    
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
    

def prep_dataset(dataset, distance=REQUIRED_DIST, do_classification=True):
    
    for k,trajec in dataset.trajecs.items():
        prep_trajectory(trajec, distance)
        
def prep_trajectory(trajec, distance=REQUIRED_DIST, do_classification=True):
    calc_frame_of_landing(trajec)
    calc_heading(trajec)
    #calc_heading_cumsum(trajec)
    trajec.frames = np.arange(get_frame_at_distance(trajec, distance), trajec.frame_of_landing).tolist()
    calc_time_to_impact(trajec)
    calc_saccades(trajec)
    calc_wv_ratio_cumsum(trajec)
    calc_dist_travelled(trajec)
    calc_dist_travelled_in_window(trajec, window=0.1, window_units='sec')
    if do_classification:
        classify(trajec, dfar=distance, dnear=0.005)
    sa1.calc_post_dynamics_for_flydra_trajectory(trajec)
    trajec.expansion = sa1.diffa(trajec.angle_subtended_by_post)*trajec.fps
    trajec.rrev = sa1.diffa(trajec.angle_subtended_by_post)*trajec.fps / trajec.angle_subtended_by_post
    calc_deceleration_initiation(trajec, plot=False)
###
def normalize_dist_to_stim_r(trajec):
    if trajec.behavior == 'landing':
        trajec.dist_to_stim_r_normed = trajec.dist_to_stim_r - trajec.dist_to_stim_r[trajec.frame_of_landing]
    else:
        trajec.dist_to_stim_r_normed = trajec.dist_to_stim_r
def get_frame_at_distance(trajec, distance, singleframe=True):
    normalize_dist_to_stim_r(trajec)
    frames = np.arange(0, trajec.frame_of_landing).tolist()
    dist_to_post = trajec.dist_to_stim_r_normed[frames]
    dist_crossovers = np.where( sa1.diffa(np.sign(dist_to_post - distance)) != 0 )[0]
    if singleframe:
        return dist_crossovers[-1]
    else:
        return dist_crossovers

def get_speed_at_distance(trajec, distance, singleframe=False):
    frames = get_frame_at_distance(trajec, distance, singleframe=singleframe)
    speed = np.max(trajec.speed[frames])
    return speed
def classify(trajec, dfar=REQUIRED_DIST, dnear=0.005):
    speed_hi_threshold = 0.18
    speed_lo_threshold = 0.18

    if trajec.behavior == 'landing' and trajec.dist_to_stim_r_normed[0] >= REQUIRED_DIST:
        trajec.speed_far = get_speed_at_distance(trajec, dfar, singleframe=True)
        trajec.speed_near = get_speed_at_distance(trajec, dnear, singleframe=True)
        
        rapid_change = ['16_11941', '21_39869', '24_68045', '10_81287', '10_81263', '6_658']
        meander = ['19_19700', '2_29065', '10_77258', '24_68566', '19_20720', '24_68072', '8_13120', '1_25985', '1_25719', '8_10323', '7_5444', '21_27634', '21_40053', '10_77509', '19_19986']
        
        '''
        if trajec.key in rapid_change:
            trajec.classification = 'rapid_change'
        elif trajec.key in meander:
            trajec.classification = 'meander'
        '''
        
        change_in_heading = np.linalg.norm( trajec.heading_smooth[trajec.frames[0]] - trajec.heading_smooth[trajec.frame_of_landing] )
        crash_test = np.min(trajec.smooth_accel[trajec.frames[0:-5]])
        projected_min_dist_to_post = dist_pt_to_line( (0,0), trajec.positions[trajec.frames[0],0:2], trajec.positions[trajec.frames[1],0:2] )
        
        if trajec.speed_far < speed_lo_threshold:
            trajec.classification = 'slow'
        
        elif crash_test >= 0:
            trajec.classification = 'crash_course'        
            
        elif len(trajec.saccades) == 0 and trajec.tortuosity < 1.1:
            trajec.classification = 'straight'
            
        elif len(trajec.saccades) == 1 and projected_min_dist_to_post > 0.03:
            trajec.classification = 'single_saccade'
        
        else:
            trajec.classification = 'mid'
            
    else:
        
        trajec.classification = 'flyby'
            
            
            
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
    vt = sa1.norm_array(trajec.positions[trajec.frames, 0:2] - trajec.positions[ trajec.frames[0]-1:trajec.frames[-1], 0:2])  
    trajec.dist_travelled = cumsum(vt)
    vt -= vt[0]
    euclidean_dist_travelled = np.linalg.norm( trajec.positions[trajec.frames[0], 0:2] - trajec.positions[trajec.frames[-1], 0:2] )
    trajec.tortuosity = trajec.dist_travelled[-1] / euclidean_dist_travelled
    trajec.mean_speed = cumsum(trajec.speed[trajec.frames]) / cumsum(np.ones_like(trajec.speed[trajec.frames]))
    
def calc_saccades(trajec, magnitude=0.15):
    raw_saccades = find_extrema_points(trajec.heading_smooth_diff)
    
    windowed_heading_diff = diff_windowed(trajec.heading, 5)
    
    trajec.saccades = []
    for saccade in raw_saccades:
        if np.abs(windowed_heading_diff[saccade]) > magnitude:
            if np.abs(saccade-trajec.frame_of_landing) > 4:
                if saccade < trajec.frame_of_landing and saccade > trajec.frames[0]:
                    trajec.saccades.append(saccade)
            
            
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
def pdf_flydra_trajecs_of_classifications(dataset, filename='classified_flydra_trajectories'):
    classified_keys = get_classified_keys(dataset)
    
    page = -1
    plt.close('all')
    pp =  pdf.PdfPages(filename)
    
    for classification, keys in classified_keys.items():
        numfliestoplot = 5
        firstfly = 0
        while firstfly < len(keys):
            
            cl = ap.xy_trajectories(dataset, trajectory=keys, show_saccades=True, trajec_alpha=0.8, firstfly=firstfly, numfliestoplot=numfliestoplot, print_obj_ids=True)
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
    
    post_types_black = ['black' for i in range(15)]
    post_types_checkered = ['checkered' for i in range(15)]
    post_types_none = ['none' for i in range(2)]
    post_types = post_types_black + post_types_checkered + post_types_none
    
    for f in dataset.filename:
        stringtime = f[4:19]
        epochtime = convert_stringtime_to_epochtime(stringtime)
        epochtimes.append(epochtime)
    epochtimes = np.array(epochtimes)

    for key, trajec in dataset.trajecs.items():
        time_diff = trajec.epoch_time[0]*np.ones_like(epochtimes) - epochtimes
        time_diff[time_diff<0] = np.inf
        f = np.argmin(time_diff)
        trajec.post_type = post_types[ f ]
    

    
def print_posttype_for_classification(dataset, c):
    classified_keys = get_classified_keys(dataset)
    keys = classified_keys[c]
    for key in keys:
        print dataset.trajecs[key].post_type
        
        
        
        
        
def plot_accel_color(dataset, keys):
    
    figure = None
    ax0_size=[0.1,0.1,0.7,0.7]
    norm=[0,0.6]
    trajec_alpha = 0.9
    #cl = colorline.Colorline(xlim=[0,0.5], ylim =[-.006,.003], norm=norm, colormap = 'jet', figure=figure, hide_colorbar=False, ax0_size=ax0_size)   
    cl = colorline.Colorline(xlim=[0,0.5], ylim =[-10,10], norm=norm, colormap = 'jet', figure=figure, hide_colorbar=False, ax0_size=ax0_size)   
    
    
    # point of deceleration initiation figure
    fig = plt.figure()
    ax = fig.add_subplot(111)
    colormap_norm = matplotlib.colors.Normalize(0, .6, clip=True)
    cmap = plt.get_cmap('jet')
    
    # point of leg extension figure
    figleg = plt.figure()
    axleg = figleg.add_subplot(111)
    
    ## find standard deviation:
    steady_state_acceleration = np.array([])
    for key in keys:
        trajec = dataset.trajecs[key] 
        a = trajec.accel_1d[0:trajec.frame_of_landing]
        t = trajec.time_to_impact[0:trajec.frame_of_landing]
        steady_state_acceleration = np.hstack( (steady_state_acceleration, a[np.where( (t<0.4)*(t>0.2) )[0].tolist()]) )
        std = np.std(steady_state_acceleration)
        
    landing_initialization = []
    
    
    for key in keys:
        trajec = dataset.trajecs[key] 
        
        x = trajec.dist_to_stim_r_normed[0:trajec.frame_of_landing]
        a = trajec.accel_1d[0:trajec.frame_of_landing]
        t = trajec.time_to_impact[0:trajec.frame_of_landing]
        s = trajec.speed[0:trajec.frame_of_landing]
        exp = sa1.diffa(trajec.angle_subtended_by_post)[0:trajec.frame_of_landing] * trajec.fps
        
        #da = sa1.diffa(trajec.smooth_accel[0:trajec.frame_of_landing])
        steady_state_acceleration = a[np.where( (t<0.4)*(t>0.2) )[0].tolist()]
        std = np.std(steady_state_acceleration)
        mean = np.mean(steady_state_acceleration)
    
        try:
            a_past_2std_pts = np.where( (a < mean-2*std)*(sa1.diffa(a)<0) )[0].tolist()
            pt_of_deceleration = a_past_2std_pts[np.where(nim.find_blobs(sa1.diffa(a_past_2std_pts)==1)[-1]==1)[0][0]]
            
            # interpolate to get a better time for the actual 2*std crossover
        
            time_to_impact_at_deceleration = np.interp(mean-2*std, np.abs(a[pt_of_deceleration-1:pt_of_deceleration+1]), t[pt_of_deceleration-1:pt_of_deceleration+1])
            speed_at_deceleration = np.interp(mean-2*std, np.abs(a[pt_of_deceleration-1:pt_of_deceleration+1]), s[pt_of_deceleration-1:pt_of_deceleration+1])
            dist_at_deceleration = np.interp(mean-2*std, np.abs(a[pt_of_deceleration-1:pt_of_deceleration+1]), x[pt_of_deceleration-1:pt_of_deceleration+1])
            expansion_at_deceleration = np.interp(mean-2*std, np.abs(a[pt_of_deceleration-1:pt_of_deceleration+1]), exp[pt_of_deceleration-1:pt_of_deceleration+1])
            angle_at_deceleration = np.interp(mean-2*std, np.abs(a[pt_of_deceleration-1:pt_of_deceleration+1]), trajec.angle_subtended_by_post[pt_of_deceleration-1:pt_of_deceleration+1])
            
            
            trajec.epoch_time_at_deceleration = np.interp(mean-2*std, np.abs(a[pt_of_deceleration-1:pt_of_deceleration+1]), trajec.epoch_time[pt_of_deceleration-1:pt_of_deceleration+1])
            trajec.frame_at_deceleration = pt_of_deceleration
            
            trajec.rrev = sa1.diffa(trajec.angle_subtended_by_post)*trajec.fps / trajec.angle_subtended_by_post
            trajec.time_aligned_to_deccel = trajec.fly_time[0:trajec.frame_of_landing]-trajec.fly_time[trajec.frame_at_deceleration]
            
            #trajec.time_at_deceleration = time_at_deceleration
            c = cmap(colormap_norm(speed_at_deceleration))
            '''
            if trajec.classification == 'straight':
                c = 'red'
            elif trajec.classification == 'slow':
                c = 'blue'
            else:
                c = 'green'
            '''
            
            #other pts:
            '''
            normalizing_factor = trajec.angle_subtended_by_post[pt_of_deceleration]
            ax.plot(trajec.speed[0:pt_of_deceleration+1], sa1.diffa(trajec.angle_subtended_by_post)[0:pt_of_deceleration+1] * trajec.fps / trajec.angle_subtended_by_post[0:pt_of_deceleration+1], color='blue')
            '''
            #normalizing_factor = 50*(speed_at_deceleration-0.7)**2+1
            #normalizing_factor = trajec.angle_subtended_by_post[pt_of_deceleration]
            fp = 0
            x = angle_at_deceleration #trajec.angle_subtended_by_post[trajec.frame_at_deceleration-fp]
            y = expansion_at_deceleration #trajec.rrev[trajec.frame_at_deceleration-fp]
            ax.plot(x*180./np.pi, y*180./np.pi, 'o', color=c)
            
            
            
            
            #landing_initialization.append(expansion_at_deceleration / normalizing_factor)
            #print speed_at_deceleration, 
            #ax.text(speed_at_deceleration, dist_at_deceleration, trajec.key)
            #cl.ax0.plot(time_at_deceleration, -2*std, 'o', color=c)
            #time_of_day = time.localtime(trajec.epoch_time[0]).tm_hour*np.ones_like(x)
            
            y = trajec.smooth_accel[0:trajec.frame_of_landing]
            
            cl.colorline(t, a, s,linewidth=1, norm=norm, alpha=trajec_alpha)
        except:
            trajec.frame_at_deceleration = None
            pass
        ax.set_title('point of deceleration past 2std')
        ax.set_xlabel('retinal size, deg')
        ax.set_ylabel('retinal expansion velocity, deg/s')
        
        cl.ax0.set_title('straight line trajectories n=118')
        cl.ax0.set_xlabel('time to impact, sec')
        cl.ax0.set_ylabel('acceleration, m/s^2')
        cl.ax1.set_ylabel('speed, m/s')
        
    return landing_initialization
    
    
def calc_deceleration_initiation(trajec, plot=False):
    trajec.frame_at_deceleration = None

    x = trajec.dist_to_stim_r_normed[0:trajec.frame_of_landing]
    a = trajec.accel_1d[0:trajec.frame_of_landing]
    t = trajec.time_to_impact[0:trajec.frame_of_landing]
    s = trajec.speed[0:trajec.frame_of_landing]
    trajec.rrev = sa1.diffa(trajec.angle_subtended_by_post)*trajec.fps / trajec.angle_subtended_by_post
    
    steady_state_acceleration = a[np.where( (t<0.4)*(t>0.2) )[0].tolist()]
    std = np.std(steady_state_acceleration)
    mean = np.mean(steady_state_acceleration)
    
    # find where acceleration starts to go down outside of normal flight
    a_past_2std_pts = np.where( (a < mean-2*std)*(sa1.diffa(a)<0) )[0].tolist()
    try:
        pt_of_deceleration = a_past_2std_pts[np.where(nim.find_blobs(sa1.diffa(a_past_2std_pts)==1)[-1]==1)[0][0]]
        
        # search backwards in frames/time to find first frame of deceleration
        f = pt_of_deceleration
        accel = a[f]
        while accel < 0:
            f -= 1
            accel = a[f]
        
        #f = pt_of_deceleration
        trajec.frame_of_deceleration = f-1
        # interpolate to get best estimate of time
        trajec.time_of_deceleration = np.interp(0, a[f-1:f+1], trajec.epoch_time[f-1:f+1]) 
        if plot:
            plt.plot(trajec.epoch_time[0:trajec.frame_of_landing], a)
            plt.plot(trajec.epoch_time[0:trajec.frame_of_landing], trajec.speed[0:trajec.frame_of_landing])
            plt.plot(trajec.time_of_deceleration, 0, 'o', color='red')
            plt.plot(trajec.epoch_time[pt_of_deceleration], a[pt_of_deceleration], 'o', color='blue')
            plt.show()
        #print trajec.key, ' processed', trajec.frame_of_deceleration
    except:
        pt_of_deceleration = None
        trajec.frame_of_deceleration = None
        trajec.time_of_deceleration = None
        #print trajec.key, ' unprocessed'
        
def plot_deceleration_vs_time_to_impact(dataset, keys, show_legs=True):
    figure = None
    ax0_size=[0.1,0.1,0.7,0.7]
    norm=[0,0.6]
    trajec_alpha = 0.9
    
    cl = colorline.Colorline(xlim=[0,0.5], ylim =[-10,10], norm=norm, colormap = 'jet', figure=figure, hide_colorbar=False, ax0_size=ax0_size)
    
    for key in keys:
        trajec = dataset.trajecs[key]
        calc_deceleration_initiation(trajec)
        if trajec.time_of_deceleration is not None:
            a = trajec.accel_1d[0:trajec.frame_of_landing]
            t = trajec.epoch_time[0:trajec.frame_of_landing] - trajec.time_of_deceleration
            s = trajec.speed[0:trajec.frame_of_landing]
            cl.colorline(t, a, s,linewidth=1, norm=norm, alpha=trajec_alpha)
            
            if show_legs:
                acc_at_legext = np.interp(trajec.legextension_time, trajec.epoch_time, trajec.accel_1d)
                timpact_at_legext = np.interp(trajec.legextension_time, trajec.epoch_time[0:trajec.frame_of_landing], t)
                cl.ax0.plot( timpact_at_legext, acc_at_legext, 'o', color='red')
    cl.ax0.set_xlabel('time, sec')
    cl.ax0.set_ylabel('acceleration, m/s^2')
            
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
        if trajec.frame_of_deceleration is not None:
            frame_offset = int(time_offset / trajec.fps)

            c = cmap(colormap_norm(trajec.speed[trajec.frame_of_deceleration+frame_offset]))
            
            frames = np.arange(trajec.frame_of_deceleration-1+frame_offset,trajec.frame_of_deceleration+1+frame_offset).tolist()
            angle_at_deceleration = np.interp(trajec.time_of_deceleration+time_offset, trajec.epoch_time[frames], trajec.angle_subtended_by_post[frames])
            expansion_at_deceleration = np.interp(trajec.time_of_deceleration+time_offset, trajec.epoch_time[frames], trajec.expansion[frames])
            
            speed_at_deceleration = np.interp(trajec.time_of_deceleration+time_offset, trajec.epoch_time[frames], trajec.speed[frames])
            dist_at_deceleration = np.interp(trajec.time_of_deceleration+time_offset, trajec.epoch_time[frames], trajec.dist_to_stim_r[frames])
            
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
        if trajec.frame_of_deceleration is not None:
            frame_offset = int(time_offset / trajec.fps)
            frame_of_interest = trajec.saccades[0] 
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
    
    
    
        
def get_rrev_and_speed(dataset, keys, time_offset=0, plot=False):
    expansion = []
    rrev = []
    speed = []
    angle = []
    for key in keys:
        trajec = dataset.trajecs[key]
        if trajec.frame_of_deceleration is not None:
            frame_offset = int(time_offset / trajec.fps)
            frames = np.arange(trajec.frame_of_deceleration-1+frame_offset,trajec.frame_of_deceleration+1+frame_offset).tolist()
            angle_at_deceleration = np.interp(trajec.time_of_deceleration+time_offset, trajec.epoch_time[frames], trajec.angle_subtended_by_post[frames])
            expansion_at_deceleration = np.interp(trajec.time_of_deceleration+time_offset, trajec.epoch_time[frames], trajec.expansion[frames])
            speed_at_deceleration = np.interp(trajec.time_of_deceleration+time_offset, trajec.epoch_time[frames], trajec.speed[frames])
            rrev_at_deceleration = expansion_at_deceleration / angle_at_deceleration

            distance_at_deceleration = np.interp(trajec.time_of_deceleration+time_offset, trajec.epoch_time[frames], trajec.dist_to_stim_r_normed[frames])
            v_over_d_at_deceleration = speed_at_deceleration / distance_at_deceleration
            
            expansion.append(expansion_at_deceleration)
            rrev.append(rrev_at_deceleration)
            #rrev.append(v_over_d_at_deceleration)
            speed.append(speed_at_deceleration)
            angle.append(angle_at_deceleration)
    if plot:
        colormap_norm = matplotlib.colors.Normalize(0, .6, clip=True)
        cmap = plt.get_cmap('jet')
        fig = plt.figure()
        ax = fig.add_subplot(111)

        for i in range(len(rrev)):
            c = cmap(colormap_norm(angle[i]))
            ax.plot(speed[i], rrev[i], 'o', color=c)
        ax.set_xlabel('velocity, m/s')
        ax.set_ylabel('rrev, 1/s')
        plt.show()      
        
    return np.array(expansion), np.array(rrev), np.array(speed), np.array(angle)
    
    
def plot_exp_vs_angle_for_speed_range(dataset, keys, velrange=[0.4,0.8], time_offset=0, plot=False):
    expansion, rrev, speed, angle = get_rrev_and_speed(dataset, keys, time_offset=0, plot=False)
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

    if plot:
        fig = plt.figure()
        ax = fig.add_subplot(111)
    
    min_number_data_pts = 12
    expansion, rrev, speed, angle = get_rrev_and_speed(dataset, keys, time_offset=0)
    vmax_bin = 0.5
    
    fitted_rrev_array = []
    fitted_speed_array = []
    fitted_residuals = []
    
    v = 0
    while v < vmax_bin:
        vmin = v
        vmax = v
        indices = []
        while len(indices) < min_number_data_pts and vmax < vmax_bin:
            vmax += 0.01
            
            if vmax >= vmax_bin:
                vmax = np.max(speed)
            
            try:
                indices = np.where( (speed>vmin)*(speed<vmax) )[0].tolist()
            except:
                indices = []
        if len(indices) > 2:
            result = np.polyfit(angle[indices], expansion[indices], 1, rcond=None, full=True)
            fit = result[0]
            residuals = result[1][0]
            r = np.sqrt(residuals / len(indices))
            rrev_at_speed = fit[0]
            
            fitted_rrev_array.append(rrev_at_speed)
            fitted_speed_array.append(np.mean(speed[indices]))
            fitted_residuals.append(r)
            
            if plot:
                ax.vlines(np.mean(speed[indices]), rrev_at_speed-r, rrev_at_speed+r)
                ax.plot(np.mean(speed[indices]), rrev_at_speed, 'o', color='black')
                            
        speed_sorted = np.argsort(speed[indices])
        difference_between_two_slowest_speeds_in_region = speed[indices][speed_sorted[1]] - speed[indices][speed_sorted[0]]
        print difference_between_two_slowest_speeds_in_region
        v += difference_between_two_slowest_speeds_in_region
    
    if plot:
        ax.set_xlabel('speed, m/s')
        ax.set_ylabel('rrev, 1/s')
        plt.show()
    
    return np.array(fitted_rrev_array), np.array(fitted_speed_array), np.array(fitted_residuals)
    
def fit_rrev_vs_speed_curve(dataset, keys, plot=False):
        
    fitted_rrev_array, fitted_speed_array, fitted_residuals = get_rrev_vs_speed(dataset, keys, time_offset=0, plot=False)
    
    indices = np.where( (fitted_speed_array<0.49)*(fitted_speed_array>0.19) )[0].tolist()
    
    result = np.polyfit(fitted_speed_array[indices], fitted_rrev_array[indices], 1, rcond=None, full=True)
    coeffs = result[0]
    
    if plot:
        fig = plt.figure()
        ax = fig.add_subplot(111)
    for i in range(len(fitted_speed_array)):
        r = fitted_residuals[i]
        
        if i in indices:
            c = 'black'
        else:
            c = 'gray'
        
        if plot:
            ax.vlines(np.mean(fitted_speed_array[i]), fitted_rrev_array[i]-r, fitted_rrev_array[i]+r, color=c)
            ax.plot(np.mean(fitted_speed_array[i]), fitted_rrev_array[i], 'o', color=c)
            
            x = np.linspace(np.min(fitted_speed_array), np.max(fitted_speed_array), 300)
            y = coeffs[0]*x**1 + coeffs[1]
            ax.plot(x,y, color='blue')
            
            ax.set_title('RREV vs. speed')
            ax.set_xlabel('speed, m/s')
            ax.set_ylabel('RREV, 1/s')
        
    return coeffs
    
    

def calc_coeff_of_var_rrev(dataset, keys, variable='deceleration'):
    fig2 = plt.figure()
    ax2 = fig2.add_subplot(111)
    
    time = np.linspace(-.5, .15, 100)
    
    for t in time:
        
        rrev_list = []
        for key in keys:
            trajec = dataset.trajecs[key] 
            if trajec.frame_of_deceleration is not None:
                if variable == 'deceleration':
                    time_at_poi = trajec.time_of_deceleration
                elif variable == 'saccade':
                    try:
                        time_at_poi = trajec.epoch_time[trajec.saccades[0]]
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
            
        c = cmap(norm(trajec.speed[trajec.saccades[0]-5]))
        plt.plot(trajec.angle_to_post[trajec.saccades[0]-5]*180/np.pi, trajec.dist_to_stim_r_normed[trajec.saccades[0]-5], 'o', color=c)
        

def simulate_rrev_for_fake_trajectory(trajec_example, pos0, vel0, time, behavior='landing'):
    
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
    
    return trajec
    
    
def plot_simulated_rrevs(dataset, keys):
    trajec_example = dataset.trajecs[dataset.trajecs.keys()[0]]
    
    figure=None
    fig = plt.figure(figure)
    ax = fig.add_subplot(111)
    colormap_norm = matplotlib.colors.Normalize(0, .3, clip=True)
    cmap = plt.get_cmap('jet')
    arrow = dict(arrowstyle="->")
    
    min_angle_to_see_post = 20*np.pi/180.
    rrev_threshold = 6.5
    
    angle_subtended_norm = (5, 40)
    cl = colorline.Colorline(xlim=[0,1], ylim =[0,50], norm=angle_subtended_norm, colormap = 'jet', figure=figure, hide_colorbar=True, ax0=ax)   
    
    coeffs = fit_rrev_vs_speed_curve(dataset, keys, plot=False)
    
    ## little hack test
    #rrev_for_speed_interp = [10.5, 12, 11.12, 13.6, 14.2]
    #speed = [.28, .32, .38, .41, .46]
    speed = np.linspace(0.25, 0.5, 5)
    
    delays = [.06]
    
    for delay in delays:
        rrev_for_speed_arr = []
        t_at_delay_arr = []
        for i, s in enumerate(speed):
            # center, slow
            pos0 = [-.2, 0, 0]
            vel0 = [s, 0, 0]
            time = 2
            trajec = simulate_rrev_for_fake_trajectory(trajec_example, pos0, vel0, time)
            frames = np.arange(0,trajec.frame_of_landing-1).tolist()
            
            rrev_for_speed = s*coeffs[0] + coeffs[1]
            #rrev_for_speed = rrev_for_speed_interp[i]
            
            t = np.interp(rrev_for_speed, trajec.rrev[frames], trajec.epoch_time[frames])
            ax.plot(t, rrev_for_speed, 'o', color='black')
            
            tdelayed = t-delay
            tminangle = np.interp(min_angle_to_see_post, trajec.angle_subtended_by_post[frames], trajec.epoch_time[frames])
            tthreshold = np.interp(rrev_threshold, trajec.rrev[frames], trajec.epoch_time[frames])
            
            ttrigger = np.max([tminangle, tthreshold])+delay
            
            rrev_at_time = np.interp(ttrigger, trajec.epoch_time[frames], trajec.rrev[frames])
            rrev_for_speed_arr.append(rrev_at_time)
            t_at_delay_arr.append(ttrigger)
            
            if delay==delays[0]:
                cl.colorline(trajec.epoch_time[frames], trajec.rrev[frames], trajec.angle_subtended_by_post[frames]*180./np.pi,linewidth=1, norm=angle_subtended_norm, alpha=1)
                
                
            if delay==delays[0] and s==speed[-1]:   
                
                #ax.plot(trajec.epoch_time[frames], trajec.rrev[frames], color='gray')
                string = 'RREV at initiation of deceleration'
                string_position = (t-0.05, rrev_for_speed+1)
                ax.annotate(string, (t, rrev_for_speed),
                    xytext=string_position,
                    arrowprops=arrow,
                    horizontalalignment='right', verticalalignment='top')
                    
        #c = cmap(colormap_norm(delay))
        c = 'red'
        ax.plot(t_at_delay_arr, rrev_for_speed_arr, '--', color=c, linewidth=2)
        ax.plot(t_at_delay_arr, rrev_for_speed_arr, 'o', color=c, markeredgecolor=c)
        
    # plot angle grid lines
    angles = [min_angle_to_see_post]
    for min_angle_to_see_post in angles:
        rrev_for_angles = []
        t_for_angles = []
        for i, s in enumerate(speed):
            # center, slow
            pos0 = [-.2, 0, 0]
            vel0 = [s, 0, 0]
            time = 2
            trajec = simulate_rrev_for_fake_trajectory(trajec_example, pos0, vel0, time)
            frames = np.arange(0,trajec.frame_of_landing-1).tolist()
        
            rrev_for_angles.append( np.interp(min_angle_to_see_post, trajec.angle_subtended_by_post[frames], trajec.rrev[frames]) )
            t_for_angles.append( np.interp(min_angle_to_see_post, trajec.angle_subtended_by_post[frames], trajec.epoch_time[frames]) )
    
        c = 'black' #cmap(angles_norm(min_angle_to_see_post*180/np.pi))
        #ax.plot(t_for_angles, rrev_for_angles, '-', color=c)
    ax.hlines(rrev_threshold, 0, 1, color='black')
    ax.text(.1, rrev_threshold, 'RREV threshold')
    
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
                
def plot_movie_dataset_leg_extension(movie_dataset):
    
    straight_ids = ['20101110_C001H001S0038', '20101110_C001H001S0032', '20101110_C001H001S0008', '20101113_C001H001S0020', '20101101_C001H001S0009', '20101030_C001H001S0035', '20101111_C001H001S0035', '20101101_C001H001S0001', '20101101_C001H001S0009', '20101110_C001H001S0035', '20101111_C001H001S0034', '20101110_C001H001S0004', '20101111_C001H001S0001', '20101101_C001H001S0020', '20101101_C001H001S0002', '20101110_C001H001S0039', '20101111_C001H001S0054', '20101111_C001H001S0058', '20101110_C001H001S0027', '20101110_C001H001S0014', '20101111_C001H001S0005']
    
    #fig = plt.figure()
    for movieid in straight_ids:
        movie = movie_dataset.movies[movieid]
        if movie.behavior == 'landing':        
            if movie.legextensionrange is not None:
                try:
                    legextensionframe = movie.legextensionrange[0] - movie.firstframe_ofinterest
                    #legextension_time = sa1.frame_to_timestamp(movie, movie.legextensionrange)[0]
                    angles_subtended = movie.scaled.angle_subtended_by_post[legextensionframe - 1: legextensionframe+1]
                    exp_at_legextension = (angles_subtended[1]-angles_subtended[0])*movie.framerate
                    angle_at_legextension = angles_subtended[1]
                    plt.plot( angle_at_legextension*180/np.pi, (exp_at_legextension / angle_at_legextension)*180/np.pi, 'o', color='black' )
                    #plt.text( angle_at_legextension, exp_at_legextension / angle_at_legextension, movieid)
                except:
                    pass
                
            
                
                
                
                
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
            
    
