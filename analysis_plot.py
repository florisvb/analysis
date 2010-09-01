# various handy plotting functions for flydra_floris_analysis
from matplotlib.pyplot import figure, show
from matplotlib.patches import Ellipse, Rectangle
import numpy as np
import colorline
import matplotlib.pyplot as pyplot
import matplotlib.pyplot as plt
import scipy.linalg
import matplotlib.patches
import matplotlib.backends.backend_pdf as pdf
import flydra_floris_analysis as analysis
import pickle
import floris
from matplotlib.patches import Patch


def load(filename, prep=True):
    fname = (filename)
    fd = open( fname, mode='r')
    print 'loading data... '
    dataset = pickle.load(fd)
    
    if prep is True:
        dataset.prep_data()
        
    return dataset
    
    
def save(dataset, filename):
    print 'saving data to file: ', filename
    fname = (filename)  
    fd = open( fname, mode='w' )
    pickle.dump(dataset, fd)
    return 1
    
    
    
def load_raw(filename, dataset=None, center = np.array([0,0,0]), radius = .01913/2., fps = None, experiment = None, gender = None, post_height = 0.3):


    ## experiment specific files ##

    if experiment is 'maimonstraw2006':
        # for straw and maimon data: 
        center = np.array([.4562, .1951, .2798])
        radius = 0.00635
        fps = 100


    if dataset is None:

        flyzone = np.array([0.05, 0.05, 0.05])
        flight_buffer = 0.008

        dataset = analysis.Dataset()
        post = analysis.Post(center, flyzone, radius, flight_buffer)
        dataset.set_stimulus(post)
        
    dataset.load_data(filename = filename, fps = fps)

    return dataset
    
    
def print_vals(dataset, trajectory, feature):
    
    if type(trajectory) is not list:
            trajectory = [trajectory]
    
    for k in trajectory:
        print k, ', ', dataset.trajecs[k].__dict__[feature]
    
######################  count flies  ###########################

def countflies (dataset, behavior = None):

    if behavior is not None:
    
        if type(behavior) is not list:
            behavior = [behavior]
    
        fly = 0
        for k,v in dataset.trajecs.items():
            if dataset.trajecs[k].behavior in behavior:
                fly = fly+1
        return fly
        
    if behavior is None:
        behaviors = ['landing', 'flyby', 'takeoff', 'unclassified']
        for b in behaviors:
            fly = 0
            for k,v in dataset.trajecs.items():
                if dataset.trajecs[k].behavior == b:
                    fly = fly+1
            print b+': ', fly
        return
        
def print_trajectory_numbers(dataset, behavior='landing'):
    if type(behavior) is not list:
        behavior = [behavior]
        
    for k,v in dataset.trajecs.items():
        if dataset.trajecs[k].behavior in behavior:
            print k
            

def get_trajec(dataset, behavior='landing'):
    
    for k,v in dataset.trajecs.items():
        if dataset.trajecs[k].behavior == behavior:
            return k
            
######################  Simple Plot  ################################

def simple_plot (dataset, attribute1, attribute2, behavior = 'landing', figure = None, a2=None):

    if type(behavior) is not list:
        behavior = [behavior]
    
    if figure is None:
        fig = plt.figure()
    else:
        fig = plt.figure(figure)

    for k,v in dataset.trajecs.items():
        if dataset.trajecs[k].behavior in behavior:
            if a2 is None:
                plt.plot( dataset.trajecs[k].__dict__[attribute1], dataset.trajecs[k].__dict__[attribute2] )
            if a2 is not None:
                plt.plot( dataset.trajecs[k].__dict__[attribute1], dataset.trajecs[k].__dict__[attribute2][:,a2] )
                    
    title_str = attribute1 + ' vs. ' + attribute2
    pyplot.title(title_str)
    pyplot.xlabel(attribute1)
    pyplot.ylabel(attribute2)
    
    
######################  Simple Scatter  ################################

def simple_scatter (dataset, attribute1, attribute2, figure = None, a1=None, a2=None, distance=0.06):
    
    if figure is None:
        fig = plt.figure()
    else:
        fig = plt.figure(figure)
        
    behavior = ['landing', 'flyby']

    for k,v in dataset.trajecs.items():
        if dataset.trajecs[k].behavior in behavior:

            distance_index = np.argmin(np.abs(dataset.trajecs[k].bins - distance))
            
            if dataset.trajecs[k].dp_residency[distance_index] < 1000:
        
                plot_style = 'g*'
                if dataset.trajecs[k].behavior == 'landing':
                    plot_style = 'ro'
                if dataset.trajecs[k].behavior == 'flyby':
                    plot_style = 'bo'
            
                x = dataset.trajecs[k].__dict__[attribute1][distance_index,a1]
                y = dataset.trajecs[k].__dict__[attribute2][distance_index,a2]
                
                if type(x) is list:
                    x = x[0]
                if type(y) is list:
                    y = y[0]
                
                plt.plot(x, y, plot_style )
                    
                    
    title_str = attribute1 + ' vs. ' + attribute2
    pyplot.title(title_str)
    pyplot.xlabel(attribute1)
    pyplot.ylabel(attribute2)


######################  Dist vs Speed plot  ###########################

def dist_vs_speed_color ( dataset, behavior = 'landing', figure = None, xlim=(-0.03, 0), ylim=(-.01,.12), norm=(-.04,.03), colormap = 'jet'):


        
    if figure is None:
        cl = colorline.Colorline(xlim=xlim, ylim=ylim, norm=norm, colormap=colormap)
    else:
        cl = colorline.Colorline(figure=figure, xlim=xlim, ylim=ylim, norm=norm, colormap=colormap)
        
    longest_trajec_val = 0
    for k,v in dataset.trajecs.items():
        if dataset.trajecs[k].behavior is behavior:
            if dataset.trajecs[k].frame_of_landing > longest_trajec_val:
                longest_trajec_val = dataset.trajecs[k].frame_of_landing
                longest_trajec_id = k
    # now we know the longest trajectory and its id number
    # find flytime of last landing, we'll add the difference between that and each trajectory to each trajectory to line up landings:
    time_of_last_landing = float(dataset.trajecs[longest_trajec_id].fly_time[longest_trajec_val])

        
    for k,v in dataset.trajecs.items():
        if dataset.trajecs[k].behavior is behavior:
        
            print k
            
            # adjust time
            time_offset = time_of_last_landing - dataset.trajecs[k].fly_time[dataset.trajecs[k].frame_of_landing]
            adjusted_time = dataset.trajecs[k].fly_time + time_offset
            
            # adjust trajectory
            # post isn't perfectly known, so fix that with an adjusted dist_to_stim: 
            # find steady state error, and subtract that
            adjusted_trajectory = dataset.trajecs[k].dist_to_stim - dataset.trajecs[k].dist_to_stim[dataset.trajecs[k].frame_of_landing]
            
            fr_land = dataset.trajecs[k].frame_of_landing
            
            
            cl.colorline(-1*adjusted_trajectory[0:fr_land], dataset.trajecs[k].speed[0:fr_land], dataset.trajecs[k].positions[0:fr_land,2])
            
            
    cl.ax0.set_title('Landings: speed vs. dist to post')
    cl.ax0.set_xlabel('distance from post surface, meters')
    cl.ax0.set_ylabel('speed (3d), meters/sec')
    cl.ax1.set_ylabel('vertical distance from top of post, meters')
    
    return cl
    
    
    
######################  deviation from beeline  ###########################
    
    
def dev_from_beeline (dataset, behavior = 'landing'):

    pyplot.figure()

    for k,v in dataset.trajecs.items():
        if dataset.trajecs[k].behavior is behavior:
            dataset.trajecs[k].calc_dev_from_beeline()
            
    for k,v in dataset.trajecs.items():
        if dataset.trajecs[k].behavior is behavior:
        
            beeline_sum = []
            dist = []
        
            for i in range(dataset.trajecs[k].length):
                beeline_sum.append(  dataset.trajecs[k].dev_from_beeline[i:dataset.trajecs[k].frame_of_landing].sum() / dataset.trajecs[k].dist_to_stim[i]*1000)
                dist.append(  dataset.trajecs[k].dist_to_stim[i]  )
                
                        
            dist_np = np.array(dist)*-1
            beeline_sum_np = np.array(beeline_sum)
            pyplot.plot(dist_np,beeline_sum_np)
                
    
        
    pyplot.title('deviation from projected beeline to post')
    pyplot.xlabel('distance to post, meters')
    pyplot.ylabel('total deviation from projected beeline normalized by dist to post, meters')
                
                
                
                
                
######################  Dist to post edge vs Time plot  ###########################

def dist_to_post_edge ( dataset, behavior = 'landing', figure = None):


        
        
    longest_trajec_val = 0
    for k,v in dataset.trajecs.items():
        if dataset.trajecs[k].behavior is behavior:
            if dataset.trajecs[k].frame_of_landing > longest_trajec_val:
                longest_trajec_val = dataset.trajecs[k].frame_of_landing
                longest_trajec_id = k
    # now we know the longest trajectory and its id number
    # find flytime of last landing, we'll add the difference between that and each trajectory to each trajectory to line up landings:
    time_of_last_landing = float(dataset.trajecs[longest_trajec_id].fly_time[longest_trajec_val])

        
    for k,v in dataset.trajecs.items():
        if dataset.trajecs[k].behavior is behavior:
        
           
            
            # adjust time
            time_offset = time_of_last_landing - dataset.trajecs[k].fly_time[dataset.trajecs[k].frame_of_landing]
            adjusted_time = dataset.trajecs[k].fly_time + time_offset
            
            # adjust trajectory
            # post isn't perfectly known, so fix that with an adjusted dist_to_stim: 
            # find steady state error, and subtract that
            adjusted_trajectory = dataset.trajecs[k].dist_to_post_edge
            
            fr_land = dataset.trajecs[k].frame_of_landing
            
            
            pyplot.plot(adjusted_time[0:fr_land], adjusted_trajectory[0:fr_land])
            
            
    pyplot.title('distance to post vertical edge')
    pyplot.xlabel('time to landing, seconds')
    pyplot.ylabel('dist to post vertical edge, meters')
            
    pyplot.show()
                


######################  Dist to post top vs Time plot  ###########################
                
def dist_to_post_top ( dataset, behavior = 'landing', figure = None):


        
        
    longest_trajec_val = 0
    for k,v in dataset.trajecs.items():
        if dataset.trajecs[k].behavior is behavior:
            if dataset.trajecs[k].frame_of_landing > longest_trajec_val:
                longest_trajec_val = dataset.trajecs[k].frame_of_landing
                longest_trajec_id = k
    # now we know the longest trajectory and its id number
    # find flytime of last landing, we'll add the difference between that and each trajectory to each trajectory to line up landings:
    time_of_last_landing = float(dataset.trajecs[longest_trajec_id].fly_time[longest_trajec_val])

        
    for k,v in dataset.trajecs.items():
        if dataset.trajecs[k].behavior is behavior:
        
           
            
            # adjust time
            time_offset = time_of_last_landing - dataset.trajecs[k].fly_time[dataset.trajecs[k].frame_of_landing]
            adjusted_time = dataset.trajecs[k].fly_time + time_offset
            
            # adjust trajectory
            # post isn't perfectly known, so fix that with an adjusted dist_to_stim: 
            # find steady state error, and subtract that
            adjusted_trajectory = dataset.trajecs[k].dist_to_post_top
            
            fr_land = dataset.trajecs[k].frame_of_landing
            
            
            pyplot.plot(adjusted_time[0:fr_land], adjusted_trajectory[0:fr_land])
            
            
    pyplot.title('distance to post vertical top')
    pyplot.xlabel('time to landing, seconds')
    pyplot.ylabel('dist to post vertical top, meters')
            
    pyplot.show()
                
                
                
                
######################  Dist to nearest edge vs Time plot  ###########################
                
def dist_to_post_nearest_edge ( dataset, behavior = 'landing', figure = None):


        
        
    longest_trajec_val = 0
    for k,v in dataset.trajecs.items():
        if dataset.trajecs[k].behavior is behavior:
            if dataset.trajecs[k].frame_of_landing > longest_trajec_val:
                longest_trajec_val = dataset.trajecs[k].frame_of_landing
                longest_trajec_id = k
    # now we know the longest trajectory and its id number
    # find flytime of last landing, we'll add the difference between that and each trajectory to each trajectory to line up landings:
    time_of_last_landing = float(dataset.trajecs[longest_trajec_id].fly_time[longest_trajec_val])

        
    for k,v in dataset.trajecs.items():
        if dataset.trajecs[k].behavior is behavior:
        
           
            
            # adjust time
            time_offset = time_of_last_landing - dataset.trajecs[k].fly_time[dataset.trajecs[k].frame_of_landing]
            adjusted_time = dataset.trajecs[k].fly_time + time_offset
            
            # adjust trajectory
            # post isn't perfectly known, so fix that with an adjusted dist_to_stim: 
            # find steady state error, and subtract that
            adjusted_trajectory = np.zeros(dataset.trajecs[k].dist_to_post_top.shape)
            for i in range(dataset.trajecs[k].length):  
                adjusted_trajectory[i] = np.min([dataset.trajecs[k].dist_to_post_top[i], dataset.trajecs[k].dist_to_post_edge[i]])
            
            fr_land = dataset.trajecs[k].frame_of_landing
            
            
            pyplot.plot(adjusted_time[0:fr_land], adjusted_trajectory[0:fr_land])
            
            
    pyplot.title('distance to post vertical top')
    pyplot.xlabel('time to landing, seconds')
    pyplot.ylabel('dist to post vertical top, meters')
            
    pyplot.show()
                
                
######################  XY plot  ###########################
                
def xy_trajectories(dataset, behavior = 'landing', figure=None, firstfly=0, numfliestoplot=None, trajectory=None, rotate=False, flip=False, lim=(-.06, .06), zslice=(-.2,.2), colorcode='s', norm=None, show_sacades=False):

    if norm is None:
        if colorcode == 'z':
            norm = (-.2, 0)
        if colorcode == 's':
            norm = (0.02, .3)
        if colorcode == 'r':
            norm = (-.2, .2)

    radius = dataset.stimulus.radius
    
    
    def rotatexy(x,y,R):
        v = np.dot(R, np.array([[x],[y]]))
        x = v[0]
        y = v[1]
        return x,y
    
    
    if type(behavior) is not list:
        behavior = [behavior]
    if trajectory is not None:
        if type(trajectory) is not list:
            trajectory = [trajectory]

    cl = colorline.Colorline(xlim=lim, ylim =lim, norm=norm, colormap = 'jet', figure=figure)
    el = Ellipse((0,0), radius*2, radius*2, facecolor='grey', alpha=0.4)
    cl.ax0.add_artist(el) 

    if trajectory is None:
        trajectory = []
        fly = -1
        for k,v in dataset.trajecs.items():
            if dataset.trajecs[k].behavior in behavior:
            
                if numfliestoplot is not None:
                    fly = fly+1
                    if fly <= firstfly:
                        continue
                    if fly > firstfly+numfliestoplot:
                        break
            
                trajectory.append(k)
                
                
                
    if trajectory is not None:
        for k in trajectory:
            if dataset.trajecs[k].behavior in behavior:
            
                if 1:
                    x = dataset.trajecs[k].initial_state['position'][0]
                    y = dataset.trajecs[k].initial_state['position'][1]
                    theta = -1*np.arctan2(y,x)
                    R = np.array([  [np.cos(theta), -np.sin(theta)],
                                    [np.sin(theta), np.cos(theta)]])
                if 1: 
                    # rotate velocity, pick y component, check sign
                    vtheta = dataset.trajecs[k].initial_state['velocity'][1]
                    F = np.sign(vtheta)
        
                fr_land = -1
                if dataset.trajecs[k].behavior is 'landing':
                    fr_land = dataset.trajecs[k].frame_of_landing
                    x = dataset.trajecs[k].positions[fr_land,0]
                    y = dataset.trajecs[k].positions[fr_land,1]
                    if rotate:
                        x,y = rotatexy(x,y,R)
                    if flip:
                        y = y*F
                    el = Ellipse((x, y), .001, .001, facecolor='black', alpha=1)
                    cl.ax0.add_artist(el)
                    
                x = dataset.trajecs[k].positions[0,0]
                y = dataset.trajecs[k].positions[0,1]
                
                if rotate:
                    x,y = rotatexy(x,y,R)
                    if flip:
                        y = y*F
                    
                el_start = Ellipse((x,y), .001, .001, facecolor='green', edgecolor='green', alpha=1)
                cl.ax0.add_artist(el_start)
                
                if show_sacades is True:
                    sacades = calc_sacades(dataset, k)
                    for sacade in sacades:
                        el_sacade = Ellipse((sacade[0],sacade[1]), .001, .001, facecolor='red', edgecolor='red', alpha=1)
                        cl.ax0.add_artist(el_sacade)
                    
                x = dataset.trajecs[k].positions[0:fr_land,0]-dataset.stimulus.center[0]
                y = dataset.trajecs[k].positions[0:fr_land,1]-dataset.stimulus.center[1]
                s = dataset.trajecs[k].speed[0:fr_land]
                z = dataset.trajecs[k].positions[0:fr_land,2]-dataset.stimulus.center[2]
                
                if rotate:
                    for i in range(len(x)):
                        x[i], y[i] = rotatexy(x[i],y[i],R)
                        if flip:
                            y[i] = y[i]*F
                        
                initial_z = dataset.trajecs[k].initial_state['position'][2]
                if  initial_z < zslice[1] and initial_z > zslice[0]:
                    if colorcode == 'z':
                        c = z
                    if colorcode == 's':
                        c = s
                    if colorcode == 'r':
                        c = dataset.trajecs[k].polar_vel[0:fr_land, 0]
                    cl.colorline(x, y, c,linewidth=1)
                    
        title = 'x-y trajectories for select trajectories'
        cl.ax0.set_title(title)
        cl.ax0.set_xlabel('x dimension, meters')
        cl.ax0.set_ylabel('y dimension, meters')
        if colorcode == 'z':
            clabel = 'z dimension, meters'
        if colorcode == 's':
            clabel = 'speed, meters/sec'
        if colorcode == 'r':
            clabel = 'radial vel, m/s'
        cl.ax1.set_ylabel(clabel)
                
    pyplot.show()

    return cl

######################  RZ plot  ###########################
                
def rz_trajectories(dataset, behavior = 'landing', figure = None, color = 'speed', firstfly=0, numfliestoplot=None, trajectory=None):

    


    if type(behavior) is not list:
        behavior = [behavior]
    if trajectory is not None:
        if type(trajectory) is not list:
            trajectory = [trajectory]
        
        
    radius = dataset.stimulus.radius
    try:
        post_height = dataset.stimulus.post_height
    except:
        post_height = 0.3
        
    #cl = colorline.Colorline(xlim=(0, .17), ylim = (-.15, .02), norm=(0, .2), colormap = 'jet', figure=figure)
    cl = colorline.Colorline(xlim=(0, .2), ylim = (-.15, .05), norm=(0, .2), colormap = 'jet', figure=figure)
    rec = Rectangle( (0,dataset.stimulus.center[2]), radius, -1*post_height, facecolor='grey', alpha=0.4)
    cl.ax0.add_artist(rec) 

    

    
    if trajectory is None:
        fly = -1
        for k,v in dataset.trajecs.items():
            if dataset.trajecs[k].behavior in behavior:
            
                if numfliestoplot is not None:
                    fly = fly+1
                    if fly <= firstfly:
                        continue
                    if fly > firstfly+numfliestoplot:
                        break
            
                fr_land = -1
                if dataset.trajecs[k].behavior is 'landing':
                    fr_land = dataset.trajecs[k].frame_of_landing
                    el = Ellipse((dataset.trajecs[k].dist_to_stim_r[fr_land]+radius, dataset.trajecs[k].positions[fr_land, 2]), .001, .001, facecolor='black', alpha=1)
                    cl.ax0.add_artist(el)
                
                if color is 'speed':
                    colorcode = dataset.trajecs[k].speed[0:fr_land]
                elif color is 'error':
                    colorcode = dataset.trajecs[k].covariance_pos_sum[0:fr_land]
                cl.colorline(dataset.trajecs[k].dist_to_stim_r[0:fr_land]+radius, dataset.trajecs[k].positions[0:fr_land, 2], colorcode,linewidth=1)
                
        if len(behavior) is 1:
            behaviortext = behavior[0]
        else:
            behaviortext = str(behavior)
        title = 'z-polar trajectories for ' + behaviortext + ' behavior'
        cl.ax0.set_title(title)
        cl.ax0.set_xlabel('radial position, meters')
        cl.ax0.set_ylabel('z position, meters')
        cl.ax1.set_ylabel('3D speed, meters/sec')        
        
        
    if trajectory is not None:
    
        for k in trajectory:
      
            if dataset.trajecs[k].behavior in behavior:
        
                fr_land = -1
                if dataset.trajecs[k].behavior == 'landing':
                    fr_land = dataset.trajecs[k].frame_of_landing
                    el = Ellipse((dataset.trajecs[k].dist_to_stim_r[fr_land]+radius, dataset.trajecs[k].dist_to_stim_z[fr_land]), .001, .001, facecolor='black', alpha=1)
                    cl.ax0.add_artist(el)
                    
                el_start = Ellipse((dataset.trajecs[k].dist_to_stim_r[0]+radius, dataset.trajecs[k].dist_to_stim_z[0]), .001, .001, facecolor='green', edgecolor='green', alpha=1)
                cl.ax0.add_artist(el_start)
                
                if color is 'speed':
                    colorcode = dataset.trajecs[k].speed[0:fr_land]
                elif color is 'error':
                    colorcode = dataset.trajecs[k].covariance_pos_sum[0:fr_land]
                cl.colorline(dataset.trajecs[k].dist_to_stim_r[0:fr_land]+radius, dataset.trajecs[k].dist_to_stim_z[0:fr_land], colorcode,linewidth=1)
                    
        
        title = 'z-polar trajectories for select trajectories' 
        cl.ax0.set_title(title)
        cl.ax0.set_xlabel('radial position, meters')
        cl.ax0.set_ylabel('z position, meters')
        cl.ax1.set_ylabel('3D speed, meters/sec') 

    
    pyplot.show()
    
    
######################  probability plot  ###########################

# probability that at a given distance from the post, a fly's trajectory will end in a landing
                
def prob_to_land(dataset, bins = None, figure = None):

    if figure is not None:
        print 'FIGURE NUMBER not supported yet'



    # histogram for landing and flybys
    # landing
    landing_hist = np.zeros(nbins-1)
    for k,v in dataset.trajecs.items():
        if dataset.trajecs[k].behavior is 'landing':
            n, hist, patches = pyplot.hist(dataset.trajecs[k].dist_to_stim, bins, normed = False, visible = False)  
            binary_hist = [floris.binarize(x) for x in n] 
            landing_hist = landing_hist + binary_hist
               
    # flyby
    flyby_hist = np.zeros(nbins-1)
    for k,v in dataset.trajecs.items():
        if dataset.trajecs[k].behavior is 'flyby':
            n, hist, patches = pyplot.hist(dataset.trajecs[k].dist_to_stim, bins, normed = False, visible = False)  
            binary_hist = [floris.binarize(x) for x in n] 
            flyby_hist = flyby_hist + binary_hist
        
            #for d in dataset.trajecs[k].dist_to_stim:
            #    dist_flyby.append( d )
    
    prob_landing = landing_hist / (landing_hist + flyby_hist)
    prob_landing = prob_landing / max(prob_landing)
    
    plt.figure(3)
    plt.plot(bins[0:len(bins)-1], prob_landing)

    pyplot.title('probability of landing')
    pyplot.xlabel('distance from post, meters')
    pyplot.ylabel('p(landing)')

    return landing_hist, flyby_hist
    
    
    
######################  residency plot  ###########################

def residency_color(dataset, behavior = 'landing', avg = True, figure = None, normalized = False, firstfly=0, numfliestoplot=None):

    if figure is not None:
        print 'FIGURE NUMBER not supported yet'


    max_dist_to_post = binrange[1]
    min_dist_to_post = binrange[0]


    
    
    cl = colorline.Colorline(xlim=(binrange[0],binrange[1]), ylim = (-1, 100), norm=(0, .2), colormap = 'jet', figure=figure)
    
    nflies = countflies(dataset, behavior = behavior)
    avg_residency = np.zeros(len(bins))

    fly = -1
    for k,v in dataset.trajecs.items():
        if dataset.trajecs[k].behavior is behavior:
        
            
            if numfliestoplot is not None:
                fly = fly+1
                if fly <= firstfly:
                    continue
                if fly > firstfly+numfliestoplot:
                    break
        
           
            
            if normalized is False:
                avg_residency = avg_residency + dataset.trajecs[k].dp_residency
            if normalized is True:
                avg_residency = avg_residency + dataset.trajecs[k].dp_residency_normalized
                
            cl.colorline(bins, dataset.trajecs[k].dp_residency, dataset.trajecs[k].dp_speed,linewidth=1)
            
            
    #avg_residency = (avg_residency / nflies)
    #avg_residency = avg_residency / max(avg_residency)
    if 0: #avg is True:
        plt.plot(bins[0:len(avg_residency)-1], avg_residency[0:len(avg_residency)-1], 'black', linewidth='2')
        title_str = 'residency vs. dist to post for ' + behavior + ' behavior'
        plt.title(title_str)
        xlabel_str = 'distance to post, meters (bin size: ' + str(binsize) + ')'
        plt.xlabel(xlabel_str)
        plt.ylabel('avg residency, normalized')

    
    

            


######################  avg diff residency plot  ###########################

def avg_diff_residency(dataset, behavior = 'landing'):

    if figure is not None:
        print 'FIGURE NUMBER not supported yet'


    bins = np.linspace(binrange[0],binrange[1],binrange[1]/binsize)
    #min_dist_to_post_index = int(binsize/min_dist_to_post)
    
    
    nflies = countflies(dataset, behavior = behavior)
    avg_residency = np.zeros(len(bins))

    start = 0.035
    end = 0.04

    for k,v in dataset.trajecs.items():
        if dataset.trajecs[k].behavior is behavior:
            dataset.trajecs[k].calc_avg_diff_residency( int(start/binsize), int(end/binsize))
            plt.plot(1,dataset.trajecs[k].avg_diff_residency,'*')
            print dataset.trajecs[k].avg_diff_residency
            
    plt.show()



######################  rdot vs r (from raw data)  ################################

def rdot_vs_r (dataset, behavior = 'flyby', figure = None, firstfly=0, numfliestoplot=None):

    if type(behavior) is not list:
        behavior = [behavior]
    
    if figure is None:
        fig = plt.figure()
    else:
        fig = plt.figure(figure)
        
    fly = -1
    for k,v in dataset.trajecs.items():
        if dataset.trajecs[k].behavior in behavior:
            
            if numfliestoplot is not None:
                fly = fly+1
                if fly <= firstfly:
                    continue
                if fly > firstfly+numfliestoplot:
                    break
            
            plt.plot( dataset.trajecs[k].dist_to_stim, dataset.trajecs[k].polar_vel[:,0])
                
    title_str = 'rdot vs r'
    pyplot.title(title_str)
    pyplot.xlabel('r, m')
    pyplot.ylabel('rdot, m/s')


######################  rdot vs r (binned by distance to allow averaging)  ################################

def rdot_vs_r_avg (dataset, behavior = 'flyby', figure = None):

    if type(behavior) is not list:
        behavior = [behavior]
    
    if figure is None:
        fig = plt.figure()
    else:
        fig = plt.figure(figure)
        
    
        
    nflies = countflies(dataset, behavior = behavior)
    total_residency = np.zeros(len(bins))
    total_polar_vel = np.zeros(len(bins))
    for k,v in dataset.trajecs.items():
        if dataset.trajecs[k].behavior in behavior:
            
            
            total_polar_vel = total_polar_vel + dataset.trajecs[k].dp_polar_vel
            binarized_residency = [floris.binarize(x) for x in dataset.trajecs[k].dp_residency]
            total_residency = total_residency + binarized_residency
            
            plt.plot( dataset.trajecs[k].bins, dataset.trajecs[k].dp_polar_vel, 'gray', linewidth='1' )
            
    total_residency_floats = [float(x) for x in total_residency]
    avg_polar_vel = np.zeros(len(bins))
    for i in range(len(bins)):
        if total_residency_floats[i] != 0:
            avg_polar_vel[i] = total_polar_vel[i] / total_residency_floats[i]

    plt.plot( bins, avg_polar_vel, 'black', linewidth='4' )
            
    plt.show()
        
                
    title_str = 'rdot vs r'
    pyplot.title(title_str)
    x_str = 'r, m (binsize: ' + str(binsize) + ')'
    pyplot.xlabel(x_str)
    pyplot.ylabel('rdot, m/s')



######################  speed vs r (binned by distance to allow averaging)  ################################

def speed_vs_r (dataset, behavior = 'flyby', figure = None, firstfly=0, numfliestoplot=None):

    if type(behavior) is not list:
        behavior = [behavior]
    
    if figure is None:
        fig = plt.figure()
    else:
        fig = plt.figure(figure)
        
    
    fly = -1
    for k,v in dataset.trajecs.items():
        if dataset.trajecs[k].behavior in behavior:
        
        
            if numfliestoplot is not None:
                fly = fly+1
                if fly <= firstfly:
                    continue
                if fly > firstfly+numfliestoplot:
                    break
        
            
            
            # find region where residency is not zero, and only plot that section
            
            s = dataset.trajecs[k].dp_speed
            
            s_mean = np.mean( s[np.nonzero( s[4:len(s)] )[0]+4] )
            
            plt.plot(s_mean,'*')
            
    
    plt.show()
        
                
    title_str = 'speed vs r'
    pyplot.title(title_str)
    x_str = 'r, m (binsize: ' + str(binsize) + ')'
    pyplot.xlabel(x_str)
    pyplot.ylabel('speed, m/s')



#############################   find similar trajectory   ##########################################
# find a trajectory of behavior 'behavior' that most closely matches the master_trajectory
# uses the initial state from the trajectory class for comparisons

def find_similar_trajec(dataset, master_trajectory, behavior = 'flyby', ntrajecs=5):
    if type(behavior) is not list:
        behavior = [behavior]

    master = dataset.trajecs[master_trajectory].initial_state


    best_match = ''
    best_err = 100
    err_arr = []
    trajec_arr = []
    for k,v in dataset.trajecs.items():
        if dataset.trajecs[k].behavior in behavior:
        
            match = dataset.trajecs[k].initial_state
            err = np.sum(   np.abs(master['velocity'] - match['velocity']) + np.abs(master['height']-match['height']) + np.abs(master['angle']-match['angle']))
            
            trajec_arr.append(k)
            err_arr.append(err)
            
            if err < best_err:
                best_err = err
                best_match = k
                #print 'behavior: ', behavior, ' individual: ', k, ' error: ', err
                
    err_arr_np = np.array(err_arr)      
    err_arr_sorting = [i for i in err_arr_np.argsort()]
    sorted_trajecs = [trajec_arr[i] for i in err_arr_sorting]
     
    return best_match, best_err, sorted_trajecs[0:ntrajecs]
            

def plot_similar_trajecs(dataset, master_trajectory=None, behavior='flyby', rotate=True, flip=True, ntrajecs=5):

    if master_trajectory is None:
        master_trajectory = get_trajec(dataset, behavior='landing')

    best_match, best_err, similar_trajecs = find_similar_trajec(dataset, master_trajectory, behavior=behavior, ntrajecs=ntrajecs)
    similar_trajecs.append(master_trajectory)
    
    
    cl = xy_trajectories(dataset, behavior=['landing', 'flyby'], trajectory=similar_trajecs, rotate=rotate, flip=flip, lim=(-.15, .15))
    x = 0.06+dataset.stimulus.radius
    y = 0
    el = Ellipse((x,y), .005, .005, facecolor='red', alpha=1)
    cl.ax0.add_artist(el) 
    #rz_trajectories(dataset, behavior=['landing', 'flyby'], trajectory=similar_trajecs)
    
def pdf_similar_trajecs(dataset):
    
    # Initialize:
    pp =  pdf.PdfPages('similartrajecs.pdf')

    # As many times as you like, create a figure fig, then either:
    f = 0
    for k,v in dataset.trajecs.items():
        if dataset.trajecs[k].behavior is 'landing':
            f += 1
            fig = plt.figure(f)
            plot_similar_trajecs(dataset, master_trajectory=k, flip=True)
            pp.savefig(fig)
            plt.close(fig)

    # Once you are done, remember to close the object:
    pp.close()



def similar_trajec_divergence(dataset):

    for k,v in dataset.trajecs.items():
        if dataset.trajecs[k].behavior is 'landing':
            master = k
            match, error, matchlist = find_similar_trajec(dataset, master)
            
            if error > 0.05:
                continue
                
            
            
            # for each attribute, plot error between match and master vs. bins
            err_polar_vel_r = np.abs(dataset.trajecs[master].dp_polar_vel[:,0] - dataset.trajecs[match].dp_polar_vel[:,0])
            #err_polar_vel_r = np.nan_to_num(err_polar_vel_r)
            if max(np.abs(err_polar_vel_r)) > 10:
                continue
                
            
            fig=plt.figure(1)
            plt.plot(dataset.trajecs[master].bins,err_polar_vel_r)
            
            
            
            
            
            
            # for each attribute, plot error between match and master vs. bins
            err_polar_vel_theta = np.abs(dataset.trajecs[master].dp_polar_vel[:,1] - dataset.trajecs[match].dp_polar_vel[:,1])
            #err_polar_vel_r = np.nan_to_num(err_polar_vel_r)
            if max(np.abs(err_polar_vel_theta)) > 10:
                continue
                
            
            fig=plt.figure(2)
            plt.plot(dataset.trajecs[master].bins,err_polar_vel_theta)
            
            
            
            
            # for each attribute, plot error between match and master vs. bins
            err_polar_vel_z = np.abs(dataset.trajecs[master].dp_polar_vel[:,2] - dataset.trajecs[match].dp_polar_vel[:,2])
            #err_polar_vel_r = np.nan_to_num(err_polar_vel_r)
            if max(np.abs(err_polar_vel_z)) > 10:
                continue
                
            
            fig=plt.figure(3)
            plt.plot(dataset.trajecs[master].bins,err_polar_vel_z)
            
            


            # for each attribute, plot error between match and master vs. bins
            err_z = np.abs(dataset.trajecs[master].dp_dist_to_post_z[:] - dataset.trajecs[match].dp_dist_to_post_z[:])
            #err_polar_vel_r = np.nan_to_num(err_polar_vel_r)
            if max(np.abs(err_z)) > 10:
                continue
                
            
            fig=plt.figure(4)
            plt.plot(dataset.trajecs[master].bins,err_z)




# for various distances, plot N vs. speed.... histograms

def feature_hist(dataset, feature_string, behavior = ['landing', 'flyby'], distance = 0.04, binsize = 0.01, binrange = [0.03,0.2], figure = None, feature_element = None):
    # feature string should be something like: "dp_speed"

    if type(behavior) is not list:
        behavior = [behavior]
        
    if type(distance) is not list:
        distance = [distance]
    
    if figure is None:
        fig = plt.figure()
    else:
        fig = plt.figure(figure)
        
    bins = np.linspace(binrange[0], binrange[1], int((binrange[1]-binrange[0])/binsize) )
    
    for d in distance:
        
        # for each behavior, create a dataset of the desired features, at the given distance
        for b in behavior:
            feature_array = []
            for k,v in dataset.trajecs.items():
                if dataset.trajecs[k].behavior is b:
                
                    distance_index = np.argmin(np.abs(dataset.trajecs[k].bins - d))
                    
                    if dataset.trajecs[k].dp_residency[distance_index] < 1000:
                        
                        try:
                            if feature_element is None:
                                feature_array.append( dataset.trajecs[k].__dict__[feature_string][distance_index] )
                            else:
                                feature_array.append( dataset.trajecs[k].__dict__[feature_string][distance_index, feature_element] )
                        except:
                            continue
            
            hist, bin_edges = np.histogram(feature_array, bins)
            normalized_hist = hist / float(max(hist))
            
            bin_centers = np.zeros(len(bin_edges)-1)
            for i in range(len(bin_edges)-1):
                bin_centers[i] = bin_edges[i] + (bin_edges[i+1]-bin_edges[i]) / 2.0
            
            line_style = '-'
            if b is 'landing':
                line_style = '-r'
            elif b is 'flyby':
                line_style = '--b'
            plt.plot(bin_centers, normalized_hist, line_style )

                    
    title_str = feature_string
    pyplot.title(title_str)
    pyplot.xlabel(feature_string)
    pyplot.ylabel('N')
    
    legend_str = []
    for b in behavior:
        legend_str.append(str(b))
    plt.legend(legend_str)







def predict_behavior(dataset, distance=0.06, method = 'all'):

    if method == 'all':
        method = ['polar vel', 'orientation']
    if type(method) is not list:
        method = [method]

    # [i,j] = guess, real
    # [0,0] = landing, landing
    # [0,1] = landing, flyby
    # [1,0] = flyby, landing
    # [1,1] = flyby, flyby
    classifier_matrix = np.zeros([2,2])
    unclassified_trajecs = []
    misclassified_trajecs = []
    correctlyclassified_trajecs = []

    

    polar_threshold = -0.055
    speed_threshold = 0.12
    orientation_threshold_upper = 1.3
    orientation_threshold_lower = 0.05
    angular_threshold = 0.1
    
    m = .1
    m_pos = m
    m_neg = -1*m
    
    
    
    for k,v in dataset.trajecs.items():
        
        distance_index = np.argmin(np.abs(dataset.trajecs[k].bins - distance))
        
        if dataset.trajecs[k].dp_residency[distance_index] < 1000:
        
        
            ## classifier ##
            classifier = 0    
            # landing:
            
            if 'polar vel' in method:
                #print 'using polar velocity classifier'
                if dataset.trajecs[k].dp_polar_vel[distance_index,0] > polar_threshold:
                    classifier = classifier + 1
                if dataset.trajecs[k].dp_polar_vel[distance_index,0] > np.abs(polar_threshold):
                    classifier = classifier - 1
                    
                if dataset.trajecs[k].dp_speed[distance_index] > (speed_threshold):
                    classifier = classifier - 1
                    
            if 'angular vel' in method:
                #print 'using polar velocity classifier'
                if np.abs(dataset.trajecs[k].dp_polar_vel[distance_index,1]) > angular_threshold:
                    classifier = classifier - 1
                
             
             
            if 'orientation'  in method:  
                #print 'using orientation classifier' 
                if dataset.trajecs[k].dp_angle_to_post[distance_index] > orientation_threshold_upper:
                    classifier = 0
                if dataset.trajecs[k].dp_angle_to_post[distance_index] < orientation_threshold_lower:
                    classifier = 0
                    
                    
            if 'speed'  in method:  
                #print 'using speed classifier' 
                if dataset.trajecs[k].dp_speed[distance_index] < speed_threshold:
                    classifier = classifier + 1
                    
            if 'angle and vel' in method:  
                if dataset.trajecs[k].dp_polar_vel[distance_index,1] > dataset.trajecs[k].dp_angle_to_post[distance_index]*m_pos:
                    classifier = classifier - 1
                if dataset.trajecs[k].dp_polar_vel[distance_index,1] < dataset.trajecs[k].dp_angle_to_post[distance_index]*m_neg:
                    classifier = classifier - 1
             
                    
            
                    
            classifier = floris.binarize(classifier / float(len(method)), threshold = 0.3)  
                    
            #################        
                    
            
            if classifier == 1: #landing       
                if dataset.trajecs[k].behavior is 'flyby':
                    misclassified_trajecs.append(k)
                    classifier_matrix[0,1] = classifier_matrix[0,1]+1
                if dataset.trajecs[k].behavior is 'landing':
                    correctlyclassified_trajecs.append(k)
                    classifier_matrix[0,0] = classifier_matrix[0,0]+1
                    
            else: # flyby
                if dataset.trajecs[k].behavior is 'flyby':
                    correctlyclassified_trajecs.append(k)
                    classifier_matrix[1,1] = classifier_matrix[1,1]+1
                if dataset.trajecs[k].behavior is 'landing':
                    misclassified_trajecs.append(k)
                    classifier_matrix[1,0] = classifier_matrix[1,0]+1

        else:
            unclassified_trajecs.append(k)

                


        
        
        
        
        
    percent_correct = np.trace(classifier_matrix) / float(np.sum(classifier_matrix))
    percent_flyby_correct = classifier_matrix[1,1] / float(classifier_matrix[1,1]+classifier_matrix[0,1])
    percent_landing_correct = classifier_matrix[0,0] / float(classifier_matrix[0,0]+classifier_matrix[1,0])
    
    tmp = { 'classifier_matrix': classifier_matrix, 
            'percent_correct': percent_correct, 
            'percent_flyby_correct': percent_flyby_correct, 
            'percent_landing_correct': percent_landing_correct, 
            'correctlyclassified_trajecs': correctlyclassified_trajecs, 
            'misclassified_trajecs': misclassified_trajecs, 
            'unclassified_trajecs': unclassified_trajecs,
            }
            
    print classifier_matrix
    print 'percent_flyby_correct: ', percent_flyby_correct
    print 'percent_landing_correct: ', percent_landing_correct
    
    return tmp
        












#############################    plot scatter of timestamps   ###############################33


def plot_timestamps (dataset, trajectory=None, behavior='landing', yval = 1, figure=None):

    if trajectory is not None:
        if type(trajectory) is not list:
                trajectory = [trajectory]
            
    if figure is None:
        fig = plt.figure()
    else:
        fig = plt.figure(figure)
    
    
    time_arr_landing = np.zeros(24)
    time_arr_flyby = np.zeros(24)
    
    
    if trajectory is None:
        for k,v in dataset.trajecs.items():
            t = int(dataset.trajecs[k].timestamp[9:11])
            if dataset.trajecs[k].behavior == 'landing':
                time_arr_landing[t] = time_arr_landing[t]+1
            if dataset.trajecs[k].behavior == 'flyby':
                time_arr_flyby[t] = time_arr_flyby[t]+1
    
    if trajectory is not None:
        for k in trajectory:
            t = int(dataset.trajecs[k].timestamp[9:11])
            if dataset.trajecs[k].behavior == 'landing':
                time_arr_landing[t] = time_arr_landing[t]+1
            if dataset.trajecs[k].behavior == 'flyby':
                time_arr_flyby[t] = time_arr_flyby[t]+1
            
    bins = np.linspace(0,23,24)
    plt.plot(bins, time_arr_landing / float(max(time_arr_landing)), 'ro')
    plt.plot(bins, time_arr_flyby / float(max(time_arr_flyby)), 'bo')
                    
                    
    title_str = 'timestamps histogram'
    pyplot.title(title_str)
    pyplot.xlabel('time (hour)')
    pyplot.ylabel('n')

        
        
#############################    heatmap   ###############################

def heatmap (dataset, trajectory=None, behavior='landing', zslice=[-.15,0], figure=None, normalize=True, plane='xy'):
    
    if behavior == 'all':
        behavior = ['landing', 'flyby', 'takeoff', 'unclassified']
    elif type(behavior) is not list:
        behavior = [behavior]
        
    radius = dataset.stimulus.radius
    
    xyrange = [-.1, .1]
    zrange = [-.15, .15]
    rrange = [0, 0.3]
    bins = np.linspace(xyrange[0], xyrange[1], 50)
    bins_x = bins
    bins_y = bins
    bins_z = np.linspace(zrange[0], zrange[1], 50)
    bins_r = np.linspace(rrange[0], rrange[1], 100)
    
    residency = np.zeros([len(bins_x),len(bins_y),len(bins_z)])
    rz_residency = np.zeros([len(bins_r), len(bins_z)])
    
    t = 1
    n = 0
    if trajectory is None:
        for k,v in dataset.trajecs.items():
            if dataset.trajecs[k].behavior in behavior:    
            
                n += 1       
                
                bin_assignment_x = np.searchsorted(bins_x, dataset.trajecs[k].positions[:,0])-1
                bin_assignment_y = np.searchsorted(bins_y, dataset.trajecs[k].positions[:,1])-1
                bin_assignment_z = np.searchsorted(bins_z, dataset.trajecs[k].positions[:,2])-1
                radial_pos = (dataset.trajecs[k].positions[:,0]**2 + dataset.trajecs[k].positions[:,1]**2)**(0.5)
                bin_assignment_r = np.searchsorted(bins_r, radial_pos)-1
                
                residency_of_trajectory = np.zeros_like(residency)
                rz_residency_of_trajectory = np.zeros_like(rz_residency)
                for i in range(dataset.trajecs[k].length):
                    x = bin_assignment_x[i]
                    y = bin_assignment_y[i]
                    z = bin_assignment_z[i]
                    r = bin_assignment_r[i]
                    residency_of_trajectory[x,y,z] += 1
                    rz_residency_of_trajectory[r,z] += 1
         
                # print just one trajectory for debugging
                if t == 0:
                    t = 1
                    print bins_z
                    print bin_assignment_z
                    print dataset.trajecs[k].positions[:,2]
                    print rz_residency_of_trajectory.sum(axis=0)
                    
                # plot residency, not dwell time
                if normalize: 
                    residency_of_trajectory = floris.binarize(residency_of_trajectory)
                    rz_residency_of_trajectory = floris.binarize(rz_residency_of_trajectory)
                residency += residency_of_trajectory
                rz_residency += rz_residency_of_trajectory
                
    
    
    residency = residency[0:len(bins_x)-1, 0:len(bins_y)-1, 0:len(bins_z)-1]
    rz_residency = rz_residency[0:len(bins_r)-1, 0:len(bins_z)-1]
    
    #print rz_residency.sum(axis=0)
    
    # if XY:
    if plane == 'xy':
        xy_residency = np.zeros([residency.shape[0], residency.shape[1]])
        zslicebins = np.searchsorted(bins_z,zslice)-1
        for i in range(zslicebins[0],zslicebins[1]):
            xy_residency += residency[:,:,i] 

        ax = pyplot.subplot(111)
        im = pyplot.imshow(xy_residency.T, cmap=pyplot.cm.jet, extent=(xyrange[0], xyrange[1], xyrange[0], xyrange[1]), origin='lower' )
        #im.set_interpolation('nearest')
        #im.set_interpolation('bicubic')
        im.set_interpolation('bilinear')
        
        el = Ellipse( (0,0), radius*2, radius*2, facecolor='white', edgecolor='none', alpha=0.8)
        ax.add_artist(el)

        pyplot.xlabel('x position, m')
        pyplot.ylabel('y position, m')
        
    # if RZ:
    if plane == 'rz':

        ax = pyplot.subplot(111)
        im = pyplot.imshow(rz_residency.T, cmap=pyplot.cm.jet, extent=(rrange[0], rrange[1], zrange[0], zrange[1]), origin='lower' )
        #im.set_interpolation('nearest')
        im.set_interpolation('bicubic')
        #im.set_interpolation('bilinear')
        
        rec = Rectangle( (0,0), radius, -1*dataset.stimulus.height, facecolor='white', edgecolor='none', alpha=0.8)
        ax.add_artist(rec)

        pyplot.xlabel('r position, m')
        pyplot.ylabel('z position, m')
        
    s = ', behavior: ' + behavior[0]
        
    if normalize:
        s = 'flux probablity' + s
        pyplot.title(s)
    else:
        s = 'dwell time probability' + s
        pyplot.title(s)
    show()

    print 'n: ', n

        
    return rz_residency
    
        
def state_heatmap (dataset, trajectory=None, behavior='landing', zslice=[-.15,0], state=None, figure=None, normalize=True, plane='xy'):
    
    if behavior == 'all':
        behavior = ['landing', 'flyby', 'takeoff', 'unclassified']
    elif type(behavior) is not list:
        behavior = [behavior]
        
    radius = dataset.stimulus.radius
    
    xyrange = [-.1, .1]
    zrange = [-.15, .15]
    rrange = [0, 0.3]
    bins = np.linspace(xyrange[0], xyrange[1], 50)
    bins_x = bins
    bins_y = bins
    bins_z = np.linspace(zrange[0], zrange[1], 50)
    bins_r = np.linspace(rrange[0], rrange[1], 100)
    
    residency = np.zeros([len(bins_x),len(bins_y),len(bins_z)])
    state_heatmap = np.zeros_like(residency)
    rz_residency = np.zeros([len(bins_r), len(bins_z)])
    rz_state_heatmap = np.zeros_like(rz_residency)
    
    t = 1
    n = 0
    if trajectory is None:
        for k,v in dataset.trajecs.items():
            if dataset.trajecs[k].behavior in behavior:    
            
                n += 1       
                
                bin_assignment_x = np.searchsorted(bins_x, dataset.trajecs[k].positions[:,0])-1
                bin_assignment_y = np.searchsorted(bins_y, dataset.trajecs[k].positions[:,1])-1
                bin_assignment_z = np.searchsorted(bins_z, dataset.trajecs[k].positions[:,2])-1
                radial_pos = (dataset.trajecs[k].positions[:,0]**2 + dataset.trajecs[k].positions[:,1]**2)**(0.5)
                bin_assignment_r = np.searchsorted(bins_r, radial_pos)-1
                
                residency_of_trajectory = np.zeros_like(residency)
                state_heatmap_of_trajectory = np.zeros_like(residency)
                rz_residency_of_trajectory = np.zeros_like(rz_residency)
                rz_state_heatmap_of_trajectory = np.zeros_like(rz_residency)
                for i in range(dataset.trajecs[k].length):
                    x = bin_assignment_x[i]
                    y = bin_assignment_y[i]
                    z = bin_assignment_z[i]
                    r = bin_assignment_r[i]
                    residency_of_trajectory[x,y,z] += 1
                    state_heatmap_of_trajectory[x,y,z] += dataset.trajecs[k].__dict__[state][i]
                    rz_residency_of_trajectory[r,z] += 1
                    rz_state_heatmap_of_trajectory[r,z] += dataset.trajecs[k].__dict__[state][i]
         
                # print just one trajectory for debugging
                if t == 0:
                    t = 1
                    print bins_z
                    print bin_assignment_z
                    print dataset.trajecs[k].positions[:,2]
                    print rz_residency_of_trajectory.sum(axis=0)
                    
                residency += residency_of_trajectory
                state_heatmap += state_heatmap_of_trajectory
                rz_residency += rz_residency_of_trajectory
                rz_state_heatmap += rz_state_heatmap_of_trajectory
    
    # normalize for residency
    if normalize: 
        residency += (residency==0) # set all zeros to one to prevent divide by zero (state value will be zero, so answer will still be zero)
        state_heatmap = state_heatmap / residency

        rz_residency += (rz_residency==0) # set all zeros to one to prevent divide by zero (state value will be zero, so answer will still be zero)
        rz_state_heatmap = rz_state_heatmap / rz_residency
            
    state_heatmap = state_heatmap[0:len(bins_x)-1, 0:len(bins_y)-1, 0:len(bins_z)-1]
    rz_state_heatmap = rz_state_heatmap[0:len(bins_r)-1, 0:len(bins_z)-1]
    
    #print rz_residency.sum(axis=0)
    
    # if XY:
    if plane == 'xy':
        xy_state_heatmap = np.zeros([state_heatmap.shape[0], state_heatmap.shape[1]])
        zslicebins = np.searchsorted(bins_z,zslice)-1
        for i in range(zslicebins[0],zslicebins[1]):
            xy_state_heatmap += state_heatmap[:,:,i] 

        ax = pyplot.subplot(111)
        im = pyplot.imshow(xy_state_heatmap.T, cmap=pyplot.cm.jet, extent=(xyrange[0], xyrange[1], xyrange[0], xyrange[1]), origin='lower' )
        #im.set_interpolation('nearest')
        #im.set_interpolation('bicubic')
        im.set_interpolation('bilinear')
        
        el = Ellipse( (0,0), radius*2, radius*2, facecolor='white', edgecolor='none', alpha=0.8)
        ax.add_artist(el)

        pyplot.xlabel('x position, m')
        pyplot.ylabel('y position, m')
        
    # if RZ:
    if plane == 'rz':

        ax = pyplot.subplot(111)
        im = pyplot.imshow(rz_state_heatmap.T, cmap=pyplot.cm.jet, extent=(rrange[0], rrange[1], zrange[0], zrange[1]), origin='lower' )
        #im.set_interpolation('nearest')
        im.set_interpolation('bicubic')
        #im.set_interpolation('bilinear')
        
        rec = Rectangle( (0,0), radius, -1*dataset.stimulus.height, facecolor='white', edgecolor='none', alpha=0.8)
        ax.add_artist(rec)

        pyplot.xlabel('r position, m')
        pyplot.ylabel('z position, m')
        
    
    show()

    print 'n: ', n

        
    return rz_state_heatmap


def parameterize_by_dist_to_post(dataset, behavior=['landing', 'flyby']):
    for k,v in dataset.trajecs.items():
        if dataset.trajecs[k].behavior in behavior:    
            dataset.trajecs[k].calc_polar_speed()
            dataset.trajecs[k].calc_angle_to_post()
            dataset.trajecs[k].calc_parameterize_by_dist()
            dataset.trajecs[k].calc_initial_state(0.06)

    
def avg_landing_dist(dataset):


    nflies = 0
    tot_landing_dist = 0
    for k,v in dataset.trajecs.items():
        if dataset.trajecs[k].behavior is 'landing':
        
            nflies += 1
            print dataset.trajecs[k].dist_to_stim[-1]
            tot_landing_dist += dataset.trajecs[k].dist_to_stim[-1]
            
    avg_landing_dist = tot_landing_dist / float(nflies)
    return avg_landing_dist
    
    
    
def pdf_heatmaps(dataset):
    
    pp =  pdf.PdfPages('heatmaps.pdf')

    # As many times as you like, create a figure fig, then either:
    f = 0
    fig = plt.figure(f)
    heatmap(dataset, behavior='flyby', plane='xy')
    pp.savefig(fig)
    plt.close(fig)
    
    f += 1
    fig = plt.figure(f)
    heatmap(dataset, behavior='flyby', plane='rz')
    pp.savefig(fig)     
    plt.close(fig)
    
    f += 1
    fig = plt.figure(f)
    heatmap(dataset, behavior='landing', plane='xy')
    pp.savefig(fig)     
    plt.close(fig)
    
    f += 1
    fig = plt.figure(f)
    heatmap(dataset, behavior='landing', plane='rz')
    pp.savefig(fig)     
    plt.close(fig)

    # Once you are done, remember to close the object:
    pp.close()
    
def pdf_sacademaps(dataset):
    pp =  pdf.PdfPages('sacademaps.pdf')
    
    f = 0
    fig = plt.figure(f)
    sacade_heatmap(dataset, behavior='flyby',plane='xy',zslice=(-.15,.15))
    pp.savefig(fig)
    plt.close(fig)
    
    f += 1
    fig = plt.figure(f)
    sacade_heatmap(dataset, behavior='landing',plane='xy',zslice=(-.15,.15))
    pp.savefig(fig)
    plt.close(fig)
    
    f += 1
    fig = plt.figure(f)
    heatmap(dataset, behavior='flyby', plane='xy', zslice=(-.15,.15))
    pp.savefig(fig)     
    plt.close(fig)
    
    f += 1
    fig = plt.figure(f)
    heatmap(dataset, behavior='landing', plane='xy', zslice=(-.15,.15))
    pp.savefig(fig)     
    plt.close(fig)
    
    pp.close()
    
    
def calc_sacades(dataset, trajec):
    # find u-turns where the fly turns sharply away from the post
    # look for places in trajectory where polar radial velocity changes from negative to positive
    trajectory = dataset.trajecs[trajec]
    window = 5 
    magnitude = 0.1
    sacades = []
    for i in range(trajectory.length-window):
        if i < window:
            continue
        old_vel = trajectory.polar_vel[i-1,0]
        new_vel = trajectory.polar_vel[i,0]
        if old_vel < 0:
            if new_vel > 0:
                # check magnitude:
                    if trajectory.polar_vel[i-window,0] < -1*magnitude:
                                
                        x = trajectory.positions[i,0]
                        y = trajectory.positions[i,1]
                        z = trajectory.positions[i,2]
                        sacades.append( (x,y,z) )
                
    return sacades
    
    
def sacade_heatmap (dataset, behavior='landing', zslice=[-.15,0], plane='xy'):
    
    if behavior == 'all':
        behavior = ['landing', 'flyby', 'takeoff', 'unclassified']
    elif type(behavior) is not list:
        behavior = [behavior]
        
    radius = dataset.stimulus.radius
    
    xyrange = [-.1, .1]
    zrange = [-.15, .15]
    rrange = [0, 0.3]
    bins = np.linspace(xyrange[0], xyrange[1], 50)
    bins_x = bins
    bins_y = bins
    bins_z = np.linspace(zrange[0], zrange[1], 50)
    bins_r = np.linspace(rrange[0], rrange[1], 100)
    
    sacade_map = np.zeros([len(bins_x),len(bins_y),len(bins_z)])
    sacade_map_rz = np.zeros([len(bins_r), len(bins_z)])
    n = 0
    for k,v in dataset.trajecs.items():
        if dataset.trajecs[k].behavior in behavior:    
        
            n += 1       
            
            sacades = calc_sacades(dataset, k)
            for sacade in sacades:        
                                    
                x = np.searchsorted(bins_x, sacade[0])-1
                y = np.searchsorted(bins_y, sacade[1])-1
                z = np.searchsorted(bins_z, sacade[2])-1
                radial_pos = (sacade[0]**2 + sacade[1]**2)**(0.5)
                r = np.searchsorted(bins_r, radial_pos)-1
            
                sacade_map[x,y,z] += 1
                sacade_map_rz[r,z] += 1
    
    # if XY:
    if plane == 'xy':
        sacade_map_xy = np.zeros([sacade_map.shape[0], sacade_map.shape[1]])
        zslicebins = np.searchsorted(bins_z,zslice)-1
        for i in range(zslicebins[0],zslicebins[1]):
            sacade_map_xy += sacade_map[:,:,i] 

        ax = pyplot.subplot(111)
        im = pyplot.imshow(sacade_map_xy.T, cmap=pyplot.cm.jet, extent=(xyrange[0], xyrange[1], xyrange[0], xyrange[1]), origin='lower' )
        #im.set_interpolation('nearest')
        im.set_interpolation('bicubic')
        #im.set_interpolation('bilinear')
        
        el = Ellipse( (0,0), radius*2, radius*2, facecolor='white', edgecolor='none', alpha=0.8)
        ax.add_artist(el)

        pyplot.xlabel('x position, m')
        pyplot.ylabel('y position, m')
        
    # if RZ:
    if plane == 'rz':

        ax = pyplot.subplot(111)
        im = pyplot.imshow(sacade_map_rz.T, cmap=pyplot.cm.jet, extent=(rrange[0], rrange[1], zrange[0], zrange[1]), origin='lower' )
        #im.set_interpolation('nearest')
        im.set_interpolation('bicubic')
        #im.set_interpolation('bilinear')
        
        rec = Rectangle( (0,0), radius, -1*dataset.stimulus.height, facecolor='white', edgecolor='none', alpha=0.8)
        ax.add_artist(rec)

        pyplot.xlabel('r position, m')
        pyplot.ylabel('z position, m')
        
    title = 'sacade map, behavior: ' + behavior[0]
    pyplot.title(title)
        
    
    show()

    print 'n: ', n

        
    return sacade_map

    
    
    
    
    
    
    
    
    
    
    

