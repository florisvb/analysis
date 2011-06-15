import flydra.a2.core_analysis as core_analysis
import flydra.analysis.result_utils as result_utils
import numpy as np
import scipy.linalg
import time
import enthought.mayavi.mlab as mlab
import matplotlib.pyplot as pyplot
import matplotlib.patches
import floris

def load_trajectory(source=None, obj=None, num=1):

    #source = '/home/floris/data/windtunnel/SA1/checkerboard/DATA20101023_135300.h5'
    #source_type = None
    #obj = 6061

    # num refers to the prefix number - relevant if using a dataset class

    if source.__class__ is Dataset:
        source_type = 'dataset'
    else:
        kalman_smoothing = True
        objs = None
        
        dataset = Dataset()
        dataset.load_data(filename = source, kalman_smoothing = kalman_smoothing, objs=obj)
        source_type = 'h5_file'
            
    if source_type is None:
        raise ValueError('please enter a valid source_type: .h5 file or dataset class')

    obj = str(num) + '_' + str(obj)
    trajectory = dataset.trajecs[obj]
    
    return trajectory

class Dataset:

    def __init__(self, like=None):
    
        if like is None:
            # trajectory related initializations
            self.trajecs = {}
            self.stimulus = None
            self.n_artificial_trajecs = 0
            
            self.datasets = []
            self.filename = []
        else:
            self.trajecs = {}
            self.stimulus = like.stimulus
            self.n_artificial_trajecs = 0
            self.datasets = like.datasets
            self.filename = like.filename
        
    def set_stimulus (self, stimulus):
        self.stimulus = stimulus

        
    def load_data (self, filename, calibration_file = None, objs = None, obj_filelist = None, kalman_smoothing = True, fps = None, dynamic_model = None, gender = None, post_type=None):
    
        self.datasets.append(len(self.datasets)+1)
        self.filename.append(filename)
        
        if calibration_file is not None:
            print "Calibration Files not yet supported!!!"
            
        # set up analyzer
        ca = core_analysis.get_global_CachingAnalyzer()
        (obj_ids, use_obj_ids, is_mat_file, data_file, extra) = ca.initial_file_load(filename)
        
        # data set defaults
        if fps is None:
            fps = result_utils.get_fps(data_file)
        if dynamic_model is None:
            try:
                dyn_model = extra['dynamic_model_name']
            except:
                dyn_model = 'EKF mamarama, units: mm'
        if dynamic_model is not None:
            dyn_model = dynamic_model
        
        # if kalman smoothing is on, then we cannot use the EKF model - remove that from the model name
        print '**Kalman Smoothing is: ', kalman_smoothing, ' **'
        if kalman_smoothing is True:
            dyn_model = dyn_model[4:]
        print 'using dynamic model: ', dyn_model
    
        if objs is None and obj_filelist is None:
            print "running through all object id's, this might take a while..."
            obj_only = use_obj_ids # this is all the unique object id's 
        if obj_filelist is not None:
            tmp = np.loadtxt(obj_filelist,delimiter=',')
            obj_only = np.array(tmp[:,0], dtype='int')
        elif objs is not None:
            if type(objs) is not list:
                objs = [objs] 
            obj_only = np.array(objs)
            
            
        print 'loading data.... '
        for obj_id in obj_only:
        
            try: 
                # kalman rows = [
                kalman_rows =  ca.load_data( obj_id, data_file,
                                     dynamic_model_name = dyn_model,
                                     use_kalman_smoothing= kalman_smoothing,
                                     frames_per_second= fps)
                                      
                                    
            except:
                print 'object id failed to load (probably no data): ', obj_id
                continue
            
            #print kalman_rows[0]
          
            # couple object ID dictionary with trajectory objects
            traj_id = (str(self.datasets[-1])+'_'+str(obj_id))
            self.trajecs.setdefault(traj_id, Trajectory(kalman_rows, extra, stimulus = self.stimulus, fps = fps, post_type=post_type) )
            #print 'loaded trajectory for object id: ', obj_id
        
            #self.trajecs.setdefault(obj, 
            
    def cull_trajectories(self, behavior = 'unknown'):
    
        if type(behavior) is not list:
            behavior = [behavior]

        for k,v in self.trajecs.items():
            if self.trajecs[k].behavior in behavior:
                print k
                del(self.trajecs[k])
                
                
            
    def create_trajectory (self, nframes, h5filebase, start, finish):
    
        # need to make some kalman rows!
        # kalman rows = [obj_id, frame, time of acquisition, positions[3], velocities[3], other crap]
        self.n_artificial_trajecs = self.n_artificial_trajecs + 1
        
        # get extra, pretty much just for the time model
        ca = core_analysis.get_global_CachingAnalyzer()
        (obj_ids, use_obj_ids, is_mat_file, data_file, extra) = ca.initial_file_load(h5filebase)
        
        
        start = np.array(start)
        finish = np.array(finish)
        
        mean_vel = (finish - start)
        
        delta_t = extra['time_model'].framestamp2timestamp(2) - extra['time_model'].framestamp2timestamp(1)

        
        kalman_rows = np.zeros([nframes, 9])
        for n in range(nframes):
            
            kalman_rows[n,0] = 1
            kalman_rows[n,1] = n
            kalman_rows[n,2] = 0
            
            
            kalman_rows[n,6:9] = mean_vel
            kalman_rows[n,3:6] = start + mean_vel*delta_t
            
        traj_id = (str('artificial'+'_'+str(self.n_artificial_trajecs)))
        self.trajecs.setdefault(traj_id, Trajectory(kalman_rows, extra, stimulus = self.stimulus) )
            
            
    def plot3dtraj (self, behavior = 'landing', fig = 1, stim = 1):
    
        # initialize 3d plot:
        
        #f = mlab.figure(fig)
        
        # plot the stimulus
        if stim:

            radius = self.stimulus.radius
            post_height = self.stimulus.center[2]
            post_bottom = -.1

            step = 0.0005
            x, y = np.mgrid[-1*radius-step:radius+step:step, -1*radius-step:radius+step:step]
            z = x*0

            for i in range(x.shape[0]):
                for j in range(y.shape[1]):
                    if (x[i,j]**2 + y[i,j]**2) < radius**2:
                        z[i,j] = post_height    
                    else:
                        z[i,j] = post_bottom

            s = mlab.surf(x, y, z, color=(0.01,0.01,0.01) )
        
        
        if 1:
            # find max speed for normalizing colors
            max_speed = 0
            for k,v in self.trajecs.items():
                if self.trajecs[k].behavior is behavior: 
                    max_speed_tmp = self.trajecs[k].speed.max()
                    if max_speed_tmp > max_speed:
                        max_speed = max_speed_tmp

            # plot trajectories of chosen behavior
            for k,v in self.trajecs.items():
                if self.trajecs[k].behavior is behavior: 
                    colors = self.trajecs[k].speed / max_speed
                    x = self.trajecs[k].positions[:,0]
                    y = self.trajecs[k].positions[:,1]
                    z = self.trajecs[k].positions[:,2]
                    mlab.points3d(x, y, z, colors, colormap="jet", scale_mode='none')
                    
        print "REMEMBER TO RUN IPYTHON USING --WTHREAD"
                
    def prep_data (self, behavior = 'landing'):
    
        for k,v in self.trajecs.items():
            self.trajecs[k].calc_polar()
            self.trajecs[k].calc_dist_to_stim()
            self.trajecs[k].calc_dist_to_stim_r()
            #self.trajecs[k].calc_dist_to_stim_z()
            #self.trajecs[k].calc_dist_to_post_top()
            #self.trajecs[k].calc_dist_to_post_edge()
            self.trajecs[k].classify(method = 'velocity')
            #self.trajecs[k].calc_polar_vel()
            #self.trajecs[k].calc_polar_speed()
            #self.trajecs[k].calc_angle_to_post()
            
            
            '''
            if self.trajecs[k].behavior is 'landing':
                #self.trajecs[k].calc_orientation_to_post_2d()
                self.trajecs[k].find_frame_of_landing()
                self.trajecs[k].classify_landingtype()           
            if self.trajecs[k].behavior is 'takeoff':
                self.trajecs[k].find_frame_of_takeoff()
            if self.trajecs[k].behavior is 'flyby':
                self.trajecs[k].find_frame_of_landing()
                '''
            if self.trajecs[k].behavior is not 'unclassified':
                self.trajecs[k].calc_accel_1d()
                
                
            
            '''
            try:
                self.trajecs[k].calc_parameterize_by_dist()
                self.trajecs[k].calc_initial_state(0.06)
            except:
                continue
            #self.trajecs[k].calc_covariance_pos()
                '''
                
    def plot_raw_data (self, figure = None, behavior = None):
    
        if figure is not None:
            pyplot.figure(figure)
        else:
            pyplot.figure()
                
                
        if behavior is None:
            behavior = ['takeoff', 'landing', 'flyby', 'unknown']
            
        if type(behavior) is not list:
            behavior = [behavior]
            
        print 'plotting trajectories for the following behaviors: ', behavior
                
                
        for k,v in self.trajecs.items():
        
            if self.trajecs[k].behavior is 'takeoff' and self.trajecs[k].behavior in behavior:
                pyplot.plot(-1*self.trajecs[k].dist_to_stim, self.trajecs[k].speed, color='blue')
            if self.trajecs[k].behavior is 'landing' and self.trajecs[k].behavior in behavior:
                pyplot.plot(-1*self.trajecs[k].dist_to_stim, self.trajecs[k].speed, color='red')
            if self.trajecs[k].behavior is 'flyby' and self.trajecs[k].behavior in behavior:
                pyplot.plot(-1*self.trajecs[k].dist_to_stim, self.trajecs[k].speed, color='green')
            if self.trajecs[k].behavior is 'unknown' and self.trajecs[k].behavior in behavior:
                pyplot.plot(-1*self.trajecs[k].dist_to_stim, self.trajecs[k].speed, color='black')
        
        pyplot.title('all trajectories: speed vs. dist to post')
        pyplot.xlabel('distance to post surface, meters')
        pyplot.ylabel('speed (3d), meters/sec')


        blue = matplotlib.patches.Patch(color='blue')
        red = matplotlib.patches.Patch(color='red')
        green = matplotlib.patches.Patch(color='green')
        black = matplotlib.patches.Patch(color='black')


        pyplot.figlegend( (blue, red, green, black), ('takeoff', 'landing', 'flyby', 'unknown' ), 'upper right')
        
             
            
class Post:


    def __init__(self, center = None, flyzone = None, radius = None, flight_buffer = None, post_height = 0.3):
    
    
        # default post:
        if center is None:
            print 'using default post values'
            center = np.array([0,0,0])
            flyzone = np.array([0.025, 0.025, 0.015])

            radius = 0.03177
            flight_buffer = 0.015
    
    
        # landing post stuff
        self.center = center #np.array([0,0,0])
        self.flyzone = flyzone #np.array([0,0,0])
        self.radius = radius #np.array([0,0,0])
        self.buffer = flight_buffer #0
        self.height = post_height
        
    def dist_to_post (self, point):
    
        dist_z = point[2] - self.center[2]
        dist_r = scipy.linalg.norm(self.center[0:2] - point[0:2])
    
        # fly is above post: distance to the top plane
        if dist_z >= 0 and dist_r <= self.radius:
            dist = dist_z
        elif dist_z >= 0 and dist_r > self.radius:
            dist = np.sqrt(dist_z**2 + (dist_r-self.radius)**2)
        elif dist_z < 0:
            dist = dist_r-self.radius
    
        return dist
        
    def dist_to_post_z (self, point):
    
        dist_z = point[2] - self.center[2]
        return dist_z
        
    def dist_to_post_r (self, point):
    
        dist_r = scipy.linalg.norm(self.center[0:2] - point[0:2])-self.radius
        return dist_r
    
   
    


class Trajectory:

    def __init__(self, kalman_rows, extra, stimulus = None, fps = None, post_type=None):
        """
        kalman rows =   [0] = obj_id
                        [1] = frame
                        [2] = timestamp
                        [3:6] = positions
                        [6:9] = velocities
                        [9:12] = P00-P02
                        [12:18] = P11,P12,P22,P33,P44,P55
                        [18:21] = rawdir_pos
                        [21:24] = dir_pos
                        
                        dtype=[('obj_id', '<u4'), ('frame', '<i8'), ('timestamp', '<f8'), ('x', '<f8'), ('y', '<f8'), ('z', '<f8'), ('xvel', '<f8'), ('yvel', '<f8'), ('zvel', '<f8'), ('P00', '<f8'), ('P01', '<f8'), ('P02', '<f8'), ('P11', '<f8'), ('P12', '<f8'), ('P22', '<f8'), ('P33', '<f8'), ('P44', '<f8'), ('P55', '<f8'), ('rawdir_x', '<f4'), ('rawdir_y', '<f4'), ('rawdir_z', '<f4'), ('dir_x', '<f4'), ('dir_y', '<f4'), ('dir_z', '<f4')])
                        
                        
        Covariance Matrix (P):
        
          [ xx  xy  xz
            xy  yy  yz
            xz  yz  zz ]
            
            full covariance for trajectories as velc as well
                        
                        
        """
        self.post_type = post_type
        self.obj_id = kalman_rows[0][0]
        self.kalman_rows = kalman_rows
        
        try:
            self.timestamp = time.strftime( '%Y%m%d_%H%M%S', time.localtime(extra['time_model'].framestamp2timestamp(kalman_rows[0][1])) )
        except:
            self.timestamp = None
        
        print len(kalman_rows)
        
        # trajectory attributes
        self.frame_numbers = np.zeros([len(kalman_rows)], dtype=int)
        self.epoch_time = np.zeros([len(kalman_rows)])
        self.fly_time = np.zeros([len(kalman_rows)])
        self.positions = np.zeros([len(kalman_rows), 3])
        self.velocities = np.zeros([len(kalman_rows), 3])
        self.speed = np.zeros([len(kalman_rows)])
        self.length = len(self.speed)
        
        if fps is not None:
            self.fps = float(fps)
        
        
        
        
        
        # stimulus related statistics
        self.behavior = 'unclassified'
        self.landingtype = 'unclassified'
        
        
        # stimulus attributes
        self.stimulus = stimulus

        
        for i in range(len(kalman_rows)):
   
            self.frame_numbers[i] = int(kalman_rows[i][1])
            
            try:
                self.epoch_time[i] = extra['time_model'].framestamp2timestamp(kalman_rows[i][1])
                self.fly_time[i] = self.epoch_time[i] - self.epoch_time[0]
            except:
                if i > 0:
                    self.fly_time[i] = self.fly_time[i-1]+1/fps 
            
            
            for j in range(3):
                self.positions[i][j] = kalman_rows[i][j+3]
                self.velocities[i][j] = kalman_rows[i][j+6]
            self.speed[i] = np.sqrt(kalman_rows[i][6]**2+kalman_rows[i][7]**2+kalman_rows[i][8]**2)
            
            
        # average time between frames (should be constant!)
        self.dt = np.mean(np.diff(self.fly_time))
        
        
    def calc_covariance_pos(self):
        
        self.covariance_pos_matrix = np.ones([self.length, 3, 3])
        for i in range(self.length):
            
            P00 = self.kalman_rows[i][9]
            P01 = self.kalman_rows[i][10]
            P02 = self.kalman_rows[i][11]
            P11 = self.kalman_rows[i][12]
            P12 = self.kalman_rows[i][13]
            P22 = self.kalman_rows[i][14]
        
            self.covariance_pos_matrix[i] = np.array([  [P00, P01, P02],
                                                        [P01, P11, P12],
                                                        [P02, P12, P22] ])
                                        
        self.covariance_pos_sum = [np.sum(self.covariance_pos_matrix[i]) for i in range(self.length)]
            
    def calc_dist_to_stim (self):
        
        self.dist_to_stim = np.zeros([self.length])
        for i in range(self.length):
            self.dist_to_stim[i] = self.stimulus.dist_to_post(self.positions[i,:])
            
    def calc_dist_to_stim_z (self):
        self.dist_to_stim_z = np.zeros([self.length])
        for i in range(self.length):
            self.dist_to_stim_z[i] = self.stimulus.dist_to_post_z(self.positions[i,:])
      
    def calc_dist_to_stim_r (self):
        self.dist_to_stim_r = np.zeros([self.length])
        for i in range(self.length):
            self.dist_to_stim_r[i] = self.stimulus.dist_to_post_r(self.positions[i,:])    
            
    def calc_pathlength_prior_landing (self):
    
        self.landing_path = np.zeros([self.length])
                
        for i in range(self.frame_of_landing+1):
           
            self.landing_path[self.frame_of_landing - i] = scipy.linalg.norm(self.positions[self.frame_of_landing - i] - self.positions[self.frame_of_landing - (i+1)]) + self.landing_path[self.frame_of_landing - i + 1]    
            
            print self.frame_of_landing - i, scipy.linalg.norm(self.positions[self.frame_of_landing - i] - self.positions[self.frame_of_landing - (i+1)]), self.landing_path[self.frame_of_landing - i]
            
            
    def calc_dist_to_post_edge (self):
    
        self.dist_to_post_edge = ( (self.dist_to_stim_r-self.stimulus.radius)**2 + self.stimulus.radius**2)**(0.5)
        
    def calc_dist_to_post_top (self):
        
        self.dist_to_post_top = ( (self.dist_to_stim_r-self.stimulus.radius)**2 + self.dist_to_stim_z**2)**(0.5)
        
    
    def calc_parameterize_by_dist (self, binsize = 0.005, binrange = 0.1, bins = None, stop_at_closest_dist=False, smoothing = False, monotonic = False):
        # note: there's a problem, we're parameterizing by distance as if the fly always flies towards the post, so if it makes a u-turn, those values will average out. 
    
        # create a list of the distances
        if bins is None:
            self.bins = np.linspace(0,binrange,binrange/binsize)
        else:
            self.bins = bins
        
        self.dp_residency = np.zeros(self.bins.shape)
        self.dp_speed = np.zeros(self.bins.shape)
        self.dp_accel_1d = np.zeros(self.bins.shape)
        self.dp_polar_vel = np.zeros([self.bins.shape[0],3])
        self.dp_dist_to_post_z = np.zeros(self.bins.shape)
        self.dp_angle_to_post = np.zeros(self.bins.shape)
        self.dp_positions = np.zeros([self.bins.shape[0], 3])
        self.dp_z = np.zeros(self.bins.shape)
        # self.dp_accel = np.zeros(self.bins.shape)
        # self.dp_angvel = np.zeros(self.bins.shape)
        
        # searchsorted: returns bin assignment
        bin_assignment = np.searchsorted(self.bins, self.dist_to_stim)
        bin_assignment = bin_assignment-1

        for i in range(self.length):
            if stop_at_closest_dist is True:
                if i>self.frame_of_landing:
                    break
                
            if monotonic is True:
                if self.dp_residency[bin_assignment[i]] != 0:
                    continue
            self.dp_residency[bin_assignment[i]] += 1
            self.dp_speed[bin_assignment[i]] = self.speed[i]+self.dp_speed[bin_assignment[i]]
            self.dp_accel_1d[bin_assignment[i]] = self.accel_1d[i]+self.dp_accel_1d[bin_assignment[i]]
            self.dp_polar_vel[bin_assignment[i]] = self.polar_vel[i]+self.dp_polar_vel[bin_assignment[i]]
            self.dp_dist_to_post_z[bin_assignment[i]] = self.positions[i,2]+self.dp_dist_to_post_z[bin_assignment[i]]
            self.dp_angle_to_post[bin_assignment[i]] = self.angle_to_post[i]+self.dp_angle_to_post[bin_assignment[i]]
            self.dp_positions[bin_assignment[i]] = self.positions[i]+self.dp_positions[bin_assignment[i]]
            self.dp_z[bin_assignment[i]] = self.positions[i,2]*self.polar_vel[i,2]*(-1.)+self.dp_z[bin_assignment[i]]
            
            
        for i in range(len(self.dp_residency)):
            if self.dp_residency[i] == 0:
                self.dp_residency[i] = 1000
            
        if monotonic is not True:
               
            self.dp_z = self.dp_z / self.dp_residency
            self.dp_z = np.nan_to_num(self.dp_z)
                
            self.dp_angle_to_post = self.dp_angle_to_post / self.dp_residency
            self.dp_angle_to_post = np.nan_to_num(self.dp_angle_to_post)
            
            for i in range(3):
                self.dp_positions[:,i] = self.dp_positions[:,i] / self.dp_residency
            self.dp_positions = np.nan_to_num(self.dp_positions)
                
            self.dp_speed = self.dp_speed / self.dp_residency
            self.dp_speed = np.nan_to_num(self.dp_speed)
            if smoothing is True:
                for i in range(1,len(self.dp_speed)-1):
                    if self.dp_speed[i] == 0:
                        self.dp_speed[i] = (self.dp_speed[i-1] + self.dp_speed[i+1]) / 2
            
            self.dp_accel_1d = self.dp_accel_1d / self.dp_residency
            self.dp_accel_1d = np.nan_to_num(self.dp_accel_1d)
            if smoothing is True:
                for i in range(1,len(self.dp_accel_1d)-1):
                    if self.dp_accel_1d[i] == 0:
                        self.dp_accel_1d[i] = (self.dp_accel_1d[i-1] + self.dp_accel_1d[i+1]) / 2
                        
            for i in range(3):
                self.dp_polar_vel[:,i] = self.dp_polar_vel[:,i] / self.dp_residency
            self.dp_polar_vel = np.nan_to_num(self.dp_polar_vel)
            if smoothing is True:
                for i in range(1,len(self.dp_polar_vel)-1):
                    for j in range(3):
                        if self.dp_polar_vel[i,j] == 0:
                            self.dp_polar_vel[i,j] = (self.dp_polar_vel[i-1,j] + self.dp_polar_vel[i+1,j]) / 2

            self.dp_dist_to_post_z = self.dp_dist_to_post_z / self.dp_residency
            self.dp_dist_to_post_z = np.nan_to_num(self.dp_dist_to_post_z)
            if smoothing is True:
                for i in range(1,len(self.dp_dist_to_post_z)-1):
                    if self.dp_dist_to_post_z[i] == 0:
                        self.dp_dist_to_post_z[i] = (self.dp_dist_to_post_z[i-1] + self.dp_dist_to_post_z[i+1]) / 2

                        
            self.dp_residency_normalized = self.dp_residency / max(self.dp_residency)
            #self.dp_residency_normalized = np.nan_to_num(self.dp_residency_normalized)
            
                    
            #print self.dp_speed
        
    
    def calc_diff_residency(self):
    
        tmp = np.diff(self.dp_residency)
        self.dp_diff_residency = np.zeros(len(tmp)+1)
        self.dp_diff_residency[1:-1] = tmp[0:-1]
        
    def calc_avg_diff_residency(self,start,end):
        self.calc_diff_residency()
        self.avg_diff_residency = np.sum(self.dp_diff_residency[start:end]) / len(self.dp_diff_residency[start:end])
        
        
    ################ initial state #####################
    
    def calc_initial_state(self, dist_to_post_r):
    
        # first need the following values in dist_to_post parameters:
        # polar velocity (r,theta,z), height
    
        # find the index of the bins for dist_to_post_r
        initial_index = np.argmin(np.abs(self.bins-dist_to_post_r))
    
        self.initial_state = {  'velocity': self.dp_polar_vel[initial_index,:],
                                'height': self.dp_dist_to_post_z[initial_index],
                                'position': self.dp_positions[initial_index,:],
                                'angle': self.dp_angle_to_post[initial_index],
                                'dist': dist_to_post_r,
                             }
        
                
            
    def find_frame_of_landing (self, threshold = 0.003):
    
        # search forward in time until we get below threshold, then check to make sure fly stays below threshold for three more frames
        
        # if it's a flyby, find nearest point to post?
        if self.behavior == 'flyby':
            self.frame_of_landing = np.argmin(self.dist_to_stim)
            return self.frame_of_landing

        
        if self.behavior == 'landing':
                
            frame_of_landing = 0
            counter = 0
            for i in range(self.length):
                if frame_of_landing == 0:
                    if self.speed[i] < threshold and self.dist_to_stim[i] < 0.01:
                        frame_of_landing = i
                if frame_of_landing > 0:
                    if counter >= 3:
                        self.frame_of_landing = frame_of_landing
                        self.time_of_landing = self.fly_time[frame_of_landing]
                        print 'frame of landing: ', frame_of_landing
                        return frame_of_landing
                    elif self.speed[i] < threshold:
                        counter = counter + 1
                    else:
                        counter = 0
                        frame_of_landing = 0
                        
            self.frame_of_landing = -1
            print 'frame of landing: ', frame_of_landing
            return -1
                
                
    def find_frame_of_takeoff (self, threshold = 0.035):
    
    
        frame_of_takeoff = 0
        counter = 0
        for i in range(self.length):
            if frame_of_takeoff == 0:
                if self.speed[i] > threshold and self.dist_to_stim[i] < 0.01:
                    frame_of_takeoff = i
            if frame_of_takeoff > 0:
                if counter >= 3:
                    self.frame_of_takeoff = frame_of_takeoff
                    self.time_of_takeoff = self.fly_time[frame_of_takeoff]
                    return frame_of_takeoff
                elif self.speed[i] > threshold:
                    counter = counter + 1
                else:
                    counter = 0
                    frame_of_takeoff = 0
                    
        self.frame_of_takeoff = -1
        return -1
        
                
    def calc_orientation (self):
    
        self.orientation = np.ones([self.length, 3])
        
        for i in range(self.length):
            self.orientation[i,:] = self.velocities[i,:] / scipy.linalg.norm(self.velocities[i,:])
            
    def calc_orientation_to_post_2d (self):
        
        self.calc_orientation()
        self.orientation_to_post_2d = np.ones([self.length])
        
        for i in range(self.length-1):
            orientation_2d = self.orientation[i,0:2] / scipy.linalg.norm(self.orientation[i,0:2])
            vector_to_post = (self.stimulus.center[0:2] - self.positions[i,0:2]) / scipy.linalg.norm(self.stimulus.center[0:2] - self.positions[i,0:2]) 
            self.orientation_to_post_2d[i] = np.dot(orientation_2d, vector_to_post)
            if self.speed[i] < 0.005:
                self.orientation_to_post_2d[i] = 1
                
    def calc_angle_to_post (self):
    
        self.calc_orientation_to_post_2d()
        self.angle_to_post = np.arccos(self.orientation_to_post_2d)
        
    def calc_polar_vel (self):
    
        # r,theta,z
        self.polar_vel = np.zeros([self.length, 3])
    
        for i in range(self.length):
            # unit vectors
            r = (self.positions[i,0:2] - self.stimulus.center[0:2]) / scipy.linalg.norm(self.positions[i,0:2] - self.stimulus.center[0:2])
            theta = np.array([-1*r[1], r[0]])
            
            self.polar_vel[i,0] = np.dot(self.velocities[i,0:2],r)
            self.polar_vel[i,1] = np.dot(self.velocities[i,0:2],theta)
            self.polar_vel[i,2] = self.velocities[i,2]
            
    def calc_polar_speed(self):
        
        self.polar_speed = np.zeros(self.length)
    
        for i in range(self.length):
            self.polar_speed[i] = np.linalg.norm(self.polar_vel[i,1])
        
    
    def calc_polar (self):
    
        # r,theta,z
        self.polar = np.zeros([self.length, 3])
            
        self.polar[:,0] = np.sqrt( self.positions[:,0]**2 + self.positions[:,1]**2 )
        
        for i in range(self.length): 
            y = self.positions[i,1]
            if np.abs(y) > 0.0005:
                self.polar[i,1] = np.arcsin ( y / self.polar[i,0] )
            else:
                self.polar[i,1] = 0
        
        self.polar[:,2] = self.positions[:,2] 
        
    def calc_accel_1d (self):
    
        self.accel_1d = np.zeros(self.length)
        
        for i in range(self.length-1):
            self.accel_1d[i] = (self.speed[i+1] - self.speed[i]) / (self.dt)
    


    def classify (self, method = None):
    
        if method is None:
                method = 'velocity'

        if self.length >= 15:

            
             
            if method is 'velocity':
            
                # place constraints on initial, final velocity, initial, final position, trajectory length?
            
            
                # find initial and final average velocities:
                
                n_frames_to_avg = 1
                # NOTE: first two frames seem to have the incorrect velocity, ie. zero
                initial_offset = 3
                
                initial_mean_vel = 0
                for i in range(n_frames_to_avg):
                    initial_mean_vel = initial_mean_vel + self.speed[i+initial_offset]
                initial_mean_vel = initial_mean_vel/n_frames_to_avg
                
                
                final_mean_vel = 0
                for i in range(n_frames_to_avg):
                    final_mean_vel = final_mean_vel + self.speed[i+len(self.speed)-3]
                final_mean_vel = final_mean_vel/n_frames_to_avg 
                
                
                flying_vel = 0.01
                notmoving_vel = 0.006
                
                # place some kind of volume restriction - distance from center axis of post:
                
                
                #distance to post:               
                
                
                
                initial_dist = self.dist_to_stim[0]
                final_dist = self.dist_to_stim[-1]
                
                min_dist_to_post = self.dist_to_stim.min()             
                    
                # start out as unclassified
                self.behavior = 'unknown'
                    
                if initial_dist > 0.04 and final_dist < 0.01 and final_mean_vel < notmoving_vel and initial_mean_vel > notmoving_vel:
                    self.behavior = 'landing'
                if initial_mean_vel <= notmoving_vel and final_mean_vel >= flying_vel and initial_dist-final_dist < -0.01 and initial_dist < 0.03:
                    self.behavior = 'takeoff'
                if min_dist_to_post < 0.03 and initial_dist > 0.04 and final_dist > 0.04:
                    self.behavior = 'flyby'                    
                #if initial_mean_vel <= notmoving_vel and final_mean_vel <= notmoving_vel:
                #    self.behavior = 'walking'                                       
            
            if method is 'volume':
                min_dist_to_post = min(self.dist_to_stim)
                frame_of_min_dist_to_post = np.argmin(self.dist_to_stim)
                landing_threshold = 0.005
                landing_speed = 0.005
                #if min_dist_to_post < landing_threshold:
                    
                    
                                     
    # classify: landing, take off, walking, fly by

    # for initial position: is inside or outside fly zone?

    
    def classify_landingtype (self):
    
        if self.behavior is not 'landing':
            print 'needs to be a landing behavior to classify type!'
    
        if self.behavior is 'landing':
            
            if self.dist_to_stim_z[self.frame_of_landing] >= -0.001:
                self.landingtype = 'top'
            else:
                self.landingtype = 'side'
            


    def calc_dev_from_beeline (self):
    
        # calculate deviation from straightline - first calculate straight path: line in 3D space
        # line in 3D space: paramterized, but t-scale doesn't matter
        
        Li = self.positions[0,:]
        Lf = self.positions[self.frame_of_landing,:]
        t_scale = Lf-Li
        self.dev_from_beeline = np.zeros(self.length)

        for i in range(self.length):
            Q = self.positions[i,:]
            t =  ((Li - Q) / t_scale.sum()).sum()
            d = scipy.linalg.norm( t_scale*t + Q - Li ) #( (t_scale[0]*t+Q[0]-Li[0])**2 + (t_scale[1]*t+Q[1]-Li[1])**2 + (t_scale[2]*t+Q[2]-Li[2])**2 )**(1/2)         
            self.dev_from_beeline[i] = d
            
                                              




