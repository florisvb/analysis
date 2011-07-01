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
import floris




def saccade_angle(dataset, keys=None):
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    fig_attractiveness = plt.figure()
    ax_attractiveness = fig_attractiveness.add_subplot(111)
    
    colormap_norm = matplotlib.colors.Normalize(0, 0.1, clip=True)
    colormap_angle_norm = matplotlib.colors.Normalize(0, np.pi, clip=True)
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
        if len(trajec.saccades) > 2:
            keys_with_saccades.append(key)
        else:
            keys_without_saccades.append(key)
            
    speed_during_saccade_array = []
    angle_of_saccade_array = []
            
    for key in keys_with_saccades:
        trajec = dataset.trajecs[key]
        #sa1.calc_post_dynamics_for_flydra_trajectory(trajec)
        
        #for s in trajec.saccades:
        if 1:
            s = trajec.saccades[-1]
            sac_range = fa.get_saccade_range(trajec, s)
            
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
            
            if signed_angleofsaccade < 0:
                saccade_shifted = signed_angleofsaccade + np.pi*2
            else:
                saccade_shifted = signed_angleofsaccade
            
            '''
            while angle_of_saccade > np.pi:
                angle_of_saccade -= np.pi
            while angle_of_saccade < -1*np.pi:
                angle_of_saccade += np.pi
            '''
            #angle_after = trajec.angle_to_post[sac_range[0]] + angle_of_saccade
            
            dist_at_saccade = trajec.dist_to_stim_r_normed[s]
            c = cmap(colormap_norm(0.1 - dist_at_saccade))
            #c = cmap(colormap_norm(trajec.speed[s]))
            #c = cmap(colormap_norm( np.log(trajec.angle_subtended_by_post[s])))
            
            speed_during_saccade_array.append(trajec.speed[s])
            angle_of_saccade_array.append(np.abs(signed_angleofsaccade))
            
            sac = patches.Circle( (angle_prior, signed_angleofsaccade), radius=0.05, facecolor=c, edgecolor='none', alpha=0.8, zorder=100)
            ax.add_artist(sac)
            
            #attractiveness = floris.dist_point_to_line((angle_prior, signed_angleofsaccade), [-np.pi,-np.pi], [np.pi,np.pi], sign=False)
            attractiveness = np.abs(angle_prior - signed_angleofsaccade) / np.pi
            
            c = cmap(colormap_angle_norm( np.pi-np.abs(angle_prior) ))
            #pt = patches.Circle( ( np.log(trajec.angle_subtended_by_post[f0]), attractiveness), radius=0.05, facecolor='green', edgecolor='none', alpha=0.2, zorder=100)
            #ax_attractiveness.add_artist(pt)
            
            
            ax_attractiveness.plot( np.log(trajec.angle_subtended_by_post[f0]), attractiveness, '.', color='green', alpha=0.2, markeredgecolor='white') 
    
    # regular ax
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
    
    fig.subplots_adjust(bottom=0.2, top=0.9, right=0.9, left=0.2)
    ax.set_aspect('equal')
    fig.savefig('saccades.pdf', format='pdf')
        
    
    # attractiveness ax
    ax_attractiveness.set_ylim(0, 1)
    ax_attractiveness.set_xlim( np.log(5*np.pi/180.), np.log(np.pi))
    ax_attractiveness.set_autoscale_on(False)
    
    
    
    fa.fix_angle_log_spine(ax_attractiveness, set_y=False, histograms=False)
    ax_attractiveness.set_yticks([0, 0.5, 1])
    #ax_attractiveness.set_xlabel('retinal size, deg')
    ax_attractiveness.set_ylabel('attractiveness')
    
    fig_attractiveness.savefig('attractiveness.pdf', format='pdf')


    # histogram stuff
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
    
    

def colorbar():
    fig = plt.figure()
    cax = fig.add_subplot(111)
    
    cticks = np.linspace(0, 0.1, 5, endpoint=True)
    cb = matplotlib.colorbar.ColorbarBase(cax, cmap=plt.get_cmap('jet'), norm=plt.Normalize(0, 0.1), orientation='horizontal', boundaries=None, ticks=cticks, drawedges=False)
    
    fig.savefig('colorbar_for_saccades.pdf', format='pdf')
    
    
    
    
    
    
def get_decision_range(dataset):
    if keys is None:
        classified_keys = fa.get_classified_keys(dataset)
        keys = classified_keys['straight']
        keys = dataset.trajecs.keys()
        
    keys_with_saccades = []
    keys_without_saccades = []
    for key in keys:
        trajec = dataset.trajecs[key]    
        if len(trajec.saccades) > 1:
            keys_with_saccades.append(key)
        else:
            keys_without_saccades.append(key)
            
    speed_during_saccade_array = []
    angle_of_saccade_array = []
            
    for key in keys_with_saccades:
        trajec = dataset.trajecs[key]
        #sa1.calc_post_dynamics_for_flydra_trajectory(trajec)
        
        #for s in trajec.all_saccades:
        s0 = trajec.saccades[-2]
        s1 = trajec.saccades[-1]
        
        
        
    
    
    
        
        
        
def free_flight_fixation(dataset, keys=None):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
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
    
    print keys_with_saccades
    
    
    for key in keys:
        
        trajec = dataset.trajecs[key]
        sa1.calc_post_dynamics_for_flydra_trajectory(trajec)
        '''
        post_angle = []
        saccade_angle = []
        
        for s in trajec.saccades:
            
            sac_range = fa.get_saccade_range(trajec, s)
            
            post_angle.append( trajec.angle_to_post[sac_range[0]]*-1 )
            
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
            
            saccade_angle.append(signed_angleofsaccade)
        '''
        edge = trajec.angle_to_post - np.sign(trajec.angle_to_post)*trajec.angle_subtended_by_post/2.
        last_frame = np.argmin(trajec.dist_to_stim_r_normed)
        frames = np.arange(0, last_frame).tolist()
        ax.plot(edge[frames]*180/np.pi, trajec.dist_to_stim_r_normed[frames], '-', color='black', alpha=0.1)
    
    ax.set_ylim(0, 0.4)
    ax.set_xlim(-180, 180)
    ax.set_autoscale_on(False)
    
    fa.adjust_spines(ax, ['left', 'bottom'])
    xticks = [-180, -90, 0, 90, 180]
    ax.set_xticks(xticks)

    ax.set_xlabel('retinal position of post center, deg')
    ax.set_ylabel('distance to post, m')
    
    fig.savefig('fixation.pdf', format='pdf')
