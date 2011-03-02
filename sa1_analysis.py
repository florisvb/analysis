import numpy as np
import cPickle as pickle
import time
import datetime
import scipy.optimize

#########################################################################################
def get_movie_dict(movie_info_filename):

    infile = open(movie_info_filename,"r")
    line_number = 0

    movie_dict = {}

    for line in infile.readlines():
        line_number += 1
        entry = map(str, line.split(' : '))
        try:
            entry[-1] = entry[-1].rstrip('\r\n')
        except:
            pass
        
        if len(entry) == 2:
            movie_dict.setdefault(entry[0], entry[1])
        
        if line_number > 18:
            break
            
    movie_date = map(int, movie_dict['Date'].split('/'))
    movie_time = map(int, movie_dict['Time'].split(':'))
    
    date = datetime.date(movie_date[0], movie_date[1], movie_date[2])
    t = date.timetuple()
    timetuple = (t[0], t[1], t[2], movie_time[0], movie_time[1], t[5], t[6], t[7], t[8])
    movie_epoch_time = time.mktime(timetuple)
    movie_dict.setdefault('EpochTime', movie_epoch_time)

    infile.close()
    return movie_dict 
    
def get_timestamp(movie_dict, frame, trigger_stamp):
    
    #trigger_stamp = movie_dict['EpochTime']
    start_frame = int(movie_dict['Start Frame'])
    trigger_frame = int(movie_dict['Correct Trigger Frame'])
    fps = float(movie_dict['Record Rate(fps)'])
    
    frame = (start_frame + frame) - trigger_frame
    dt = float(frame) * 1/float(fps)
    timestamp = trigger_stamp + dt
    
    return timestamp
    
def get_timestamps(movie_dict, trigger_stamp):
    nframes = int(movie_dict['Total Frame'])
    timestamps = []
    for frame in range(nframes):
        timestamp = get_timestamp(movie_dict, frame, trigger_stamp)
        timestamps.append(timestamp)
    movie_dict.setdefault('Timestamps', np.array(timestamps))        
    
def load_movie_info(movie_list):
    movie_list_fd = open(movie_list, 'r')
    movies = {}
    for line in movie_list_fd.readlines():
        entry = line.split(',')
        
        # clean entry
        for i, enum in enumerate(entry):
            entry[i] = entry[i].rstrip()
            entry[i] = entry[i].lstrip()
            try:
                entry[i] = entry[i].rstrip('\n')
            except:
                pass
            entry[i] = entry[i].rstrip()
            entry[i] = entry[i].lstrip()
        #
        #try:
        movie_path = entry[0]
        print
        print movie_path
        print '***************'
        if movie_path[-1] != '/':
            movie_path = movie_path + '/'
        print entry
        movie_id = entry[1]
        movie_info_filename = movie_path + movie_id + '.cih'
        print movie_info_filename
        behavior = entry[2]
        awesome = int(entry[3])
        post_type = entry[4]
        extras = entry[5]
        extras = extras.split('%')
        try:
            movie_dict = get_movie_dict(movie_info_filename)
            movie_dict.setdefault('Behavior', behavior)
            movie_dict.setdefault('Awesome', awesome)
            movie_dict.setdefault('PostType', post_type)
            movie_dict.setdefault('Extras', extras)
            movie_dict.setdefault('Path', movie_path)
            
            movies.setdefault(movie_id, movie_dict)
        except:
            print 'pass'
            pass
    movie_list_fd.close()
    return movies
    
def get_flydra_frame(npmovie, delay=0, find_delay=False):
    # algorithm works, but does not work on the data.. there is some delay in the two timestamps
    sa1_epochtime = npmovie.epochtime+delay
    
    flydra_frames = np.zeros(len(npmovie.uframes))
    errors = np.zeros_like(flydra_frames)
    
    if find_delay is True:
        start = 0
        stop = -1
    else:
        start = 0
        stop = -1       
    
    for f, uframe in enumerate(npmovie.uframes[start:stop]):
        if uframe is not None:
            sa1_time = f*1.0/npmovie.fps + sa1_epochtime
            
            err = np.zeros_like(npmovie.trajec.epoch_time)
            for i, time in enumerate(npmovie.trajec.epoch_time):
                err[i] = np.abs(time-sa1_time)
    
            best_match = np.argmin(err) 
            #print sa1_time, npmovie.trajec.epoch_time[best_match]
            flydra_frames[f] = best_match
            errors[f] = err[best_match]
            if find_delay is False:
                try:
                    delay = npmovie.flydradelay
                except:
                    pass
                uframe.flydraframe = best_match
                try:
                    uframe.timestamp_adjusted = uframe.timestamp+delay
                except:
                    pass
                try:
                    npmovie.timestamps_adjusted = npmovie.timestamps+delay+sa1_epochtime
                except:
                    pass
    if find_delay is True:
        e = np.sum(errors)
        print e
        return e
    else:
        return flydra_frames, errors
        
def fmin_wrapper_get_flydra_frame(delay, npmovie):
    e = get_flydra_frame(npmovie, delay=delay, find_delay=True)
    return e
    
def find_sa1_timestamp_delay(npmovie, guess=30):
    try:
        delay = guess
        results = scipy.optimize.fmin(fmin_wrapper_get_flydra_frame, delay, args=(npmovie,), full_output=True)
        delay = results[0][0]
        npmovie.flydradelay = delay 
        
        get_flydra_frame(npmovie, delay=delay)
        return delay
    except:
        return None
        
def get_flydra_trajectory(npmovie, dataset):

    t = time.localtime(npmovie.epochtime)
    try:
        obj_id = int(npmovie.objid)
    except: 
        npmovie.dataset_id = None
        npmovie.trajec = None
        return None, None
        
    epochtime = npmovie.epochtime
    
    print 'dataset should include: ', time.strftime('%Y%m%d %H%M%S', t)
    print 'obj id: ', obj_id
    
    try:
        for k, trajectory in dataset.trajecs.items():
            t = k.lstrip('1234567890')
            t = t.lstrip('_')
            o = int(t)
            if obj_id == o:
                # check timestamp:
                time_err = np.abs(trajectory.epoch_time[0] - epochtime)
                print time_err
                if time_err < 100:
                    npmovie.dataset_id = k
                    npmovie.trajec = trajectory
                    return trajectory, time_err
            else:
                continue    
    except:
        npmovie.dataset_id = None
        npmovie.trajec = None
        return None, None   
    npmovie.dataset_id = None 
    npmovie.trajec = None
    return None, None
        
        
    
def get_movie_obj_id(movie, obj_id_list):
    movie_time = movie['EpochTime']
    errors = np.zeros_like(obj_id_list[:,1])
    for i in range(len(obj_id_list)):
        errors[i] =  np.abs(obj_id_list[i,1] - movie_time)
        #print i, obj_id_list[i,1]#, movie_time, obj_id_list[i,1]-movie_time
    #print
    #print min(errors)
    if min(errors) < 200:
        obj_id = int(obj_id_list[ np.argmin(errors), 0])
        trigger_stamp = float(obj_id_list[ np.argmin(errors), 1])
    else:
        obj_id = None
        trigger_stamp = 0
    print obj_id, trigger_stamp
    movie.setdefault('Obj_ID', obj_id)
    movie.setdefault('Trigger Stamp', trigger_stamp)
    return obj_id, trigger_stamp
    
    
def get_awesome_movies(movies, behavior='landing', extras=None, awesome=1):
    awesome_list = []
    
    if type(behavior) is not list:
        behavior = [behavior]
    if type(awesome) is not list:
        awesome = [awesome]
        
    if extras is None:
        for k, movie in movies.items():
            if movie['Behavior'] in behavior:
                if movie['Awesome'] in awesome:
                     awesome_list.append((k, movie['Obj_ID'], movie['Date']))
    
    else:
        for k, movie in movies.items():
            if extras in movie['Extras']:
                if movie['Awesome'] in awesome:
                     awesome_list.append((k, movie['Obj_ID'], movie['Date']))
    
    return awesome_list
#########################################################################################


def sa1_analysis():
    if 0:
        sa1_obj_id_files = [    '/home/floris/data/windtunnel/SA1/checkerboard/SA1_20101023', 
                                '/home/floris/data/windtunnel/SA1/checkerboard/SA1_20101024',
                                '/home/floris/data/windtunnel/SA1/black/SA1_20101025',
                                '/home/floris/data/windtunnel/SA1/black/SA1_20101026',
                                '/home/floris/data/windtunnel/SA1/black/SA1_20101026_a',
                                '/home/floris/data/windtunnel/SA1/black/SA1_20101027',
                                '/home/floris/data/windtunnel/SA1/black/SA1_20101028',
                                '/home/floris/data/windtunnel/SA1/checker_angle/SA1_20101029',
                                '/home/floris/data/windtunnel/SA1/checker_angle/SA1_20101031',
                                '/home/floris/data/windtunnel/SA1/checker_angle/SA1_20101101',  
                                '/home/floris/data/windtunnel/SA1/checker_angle/SA1_20101109',  
                                '/home/floris/data/windtunnel/SA1/checker_angle/SA1_20101110',    
                                '/home/floris/data/windtunnel/SA1/black_angle/SA1_20101111', 
                                '/home/floris/data/windtunnel/SA1/black_angle/SA1_20101113',                   
                                ]
    if 1:
        sa1_obj_id_files = [    '/home/floris/Documents/data/sa1_movie_data/SA1_20101023', 
                                '/home/floris/Documents/data/sa1_movie_data/SA1_20101024',
                                '/home/floris/Documents/data/sa1_movie_data/SA1_20101025',
                                '/home/floris/Documents/data/sa1_movie_data/SA1_20101026',
                                '/home/floris/Documents/data/sa1_movie_data/SA1_20101026_a',
                                '/home/floris/Documents/data/sa1_movie_data/SA1_20101027',
                                '/home/floris/Documents/data/sa1_movie_data/SA1_20101028',
                                '/home/floris/Documents/data/sa1_movie_data/SA1_20101029',
                                '/home/floris/Documents/data/sa1_movie_data/SA1_20101031',
                                '/home/floris/Documents/data/sa1_movie_data/SA1_20101101',  
                                '/home/floris/Documents/data/sa1_movie_data/SA1_20101109',  
                                '/home/floris/Documents/data/sa1_movie_data/SA1_20101110',    
                                '/home/floris/Documents/data/sa1_movie_data/SA1_20101111', 
                                '/home/floris/Documents/data/sa1_movie_data/SA1_20101113',                   
                                ]
                            
    movie_list = '/home/floris/Documents/data/sa1_classification.txt'


    obj_id_list = None
    for filename in sa1_obj_id_files:
        ## load SA1 object ID data ##
        infile = open(filename, 'r')
        unpickler = pickle.Unpickler(infile)
        # this bullshit needed to deal with the pickle file, which has multiple dumps.. only the last load() call has the whole file
        new_obj_id_list = unpickler.load()
        if 1:
            while 1:
                try:
                    new_obj_id_list = unpickler.load()
                except:
                    break
        infile.close()
        
        if obj_id_list == None:
            obj_id_list = new_obj_id_list
        else:
            obj_id_list = np.vstack((obj_id_list, new_obj_id_list))
            
    ## load SA1 movie info ##
    if 1:
        movies = load_movie_info(movie_list)

        for k, movie in movies.items():
            obj_id, trigger_stamp = get_movie_obj_id(movie, obj_id_list)
            get_timestamps(movie, trigger_stamp)
            #get_flydra_trajectory(movie, dataset)
            

    return obj_id_list, movies


