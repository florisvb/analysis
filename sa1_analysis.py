import numpy as np
import cPickle as pickle
import time
import datetime


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
        
        if line_number > 3:
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
    
def load_movie_info(movie_path, movie_list):
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
        try:
            movie_id = entry[0]
            movie_info_filename = movie_path + movie_id + '.cih'
            behavior = entry[1]
            awesome = int(entry[2])
            post_type = entry[3]
            movie_dict = get_movie_dict(movie_info_filename)
            movie_dict.setdefault('Behavior', behavior)
            movie_dict.setdefault('Awesome', awesome)
            movie_dict.setdefault('PostType', post_type)
            
            movies.setdefault(movie_id, movie_dict)
        except:
            pass
    movie_list_fd.close()
    return movies
    
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
    else:
        obj_id = None
    print obj_id
    return obj_id
    
    
def get_awesome_movies(movies, behavior='landing'):
    awesome_list = []
    for k, movie in movies.items():
        if movie['Behavior'] == behavior:
            if movie['Awesome'] == 1:
                 awesome_list.append((k, movie['Obj_ID'], movie['Date']))
                 
    return awesome_list
#########################################################################################


if __name__ == '__main__':

    sa1_obj_id_files = [    '/home/floris/data/windtunnel/SA1/checkerboard/SA1_20101023', 
                            '/home/floris/data/windtunnel/SA1/checkerboard/SA1_20101024',
                            '/home/floris/data/windtunnel/SA1/black/SA1_20101025']
    movie_path = '/media/SA1_videos/sa1_videos/'
    movie_list = '/media/SA1_videos/sa1_classification.txt'


    obj_id_list = None
    for filename in sa1_obj_id_files:
        ## load SA1 object ID data ##
        infile = open(filename, 'r')
        unpickler = pickle.Unpickler(infile)
        # this bullshit needed to deal with the pickle file, which has multiple dumps.. only the last load() call has the whole file
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
        movies = load_movie_info(movie_path, movie_list)

        for k, movie in movies.items():
            obj_id = get_movie_obj_id(movie, obj_id_list)
            movie.setdefault('Obj_ID', obj_id)

            




