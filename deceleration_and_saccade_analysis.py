
import matplotlib.pyplot as plt
import floris_plot_lib as fpl
import numpy as np
import saccade_analysis as sac

def get_frame_of_decel_and_last_saccade(trajec):
    
    frame_of_last_saccade_initiation = None
    
    if trajec.behavior == 'landing':
        for sac in trajec.sac_ranges:
            if sac[-1] < trajec.frame_of_landing:
                frame_of_last_saccade_initiation = sac[0]
            else:
                break
    if trajec.behavior == 'flyby':
        for sac in trajec.sac_ranges:
            if sac[0] < trajec.frame_nearest_to_post:
                frame_of_last_saccade_initiation = sac[0]
            else:
                break
    
    
    return trajec.frame_at_deceleration, frame_of_last_saccade_initiation
    
def get_saccade_range_from_first_frame(trajec, s):
    
    for sac_range in trajec.sac_ranges:
        if s in sac_range:
            return sac_range
            
def print_frame_diff_btwn_decel_and_saccade(dataset):
    
    neg = 0
    pos = 0
    none = 0
    
    fig_neg = plt.figure()
    ax_neg = fig_neg.add_subplot(111)
    
    fig_pos = plt.figure()
    ax_pos = fig_pos.add_subplot(111)
    
    fig_none = plt.figure()
    ax_none = fig_none.add_subplot(111)
    
    fig_saccades = plt.figure()
    ax_saccades = fig_saccades.add_subplot(111)
    
    fig_deceleration = plt.figure()
    ax_deceleration = fig_deceleration.add_subplot(111)
        
    for k, trajec in dataset.trajecs.items():
        fd, fs = get_frame_of_decel_and_last_saccade(trajec)
        pos_frame_minus_landing = []
        
        if trajec.classification == 'straight':
            none += 1
            ax_none.plot(np.log(trajec.angle_subtended_by_post[fd]), trajec.speed[fd], '.', color='purple')
            ax_deceleration.plot(np.log(trajec.angle_subtended_by_post[fd]), trajec.speed[fd], '.', color='purple')
        else:
            try:
                difference = fd - fs
                if difference < 0:
                    neg += 1
                    ax_neg.plot(np.log(trajec.angle_subtended_by_post[fd]), trajec.speed[fd], '.', color='purple')
                    ax_neg.plot(np.log(trajec.angle_subtended_by_post[fs]), trajec.speed[fs], '.', color='green')
                    ax_neg.plot([np.log(trajec.angle_subtended_by_post[fd]), np.log(trajec.angle_subtended_by_post[fs])], [trajec.speed[fd], trajec.speed[fs]], '-', color='black')
                    ax_deceleration.plot(np.log(trajec.angle_subtended_by_post[fd]), trajec.speed[fd], '.', color='red')
                    
                    sac_range = get_saccade_range_from_first_frame(trajec, fs)
                    angle = sac.get_angle_of_saccade(trajec, sac_range)
                    print trajec.angle_to_post, angle
                    ax_saccades.plot(trajec.angle_to_post[fs], angle, '.')
                    print trajec.angle_to_post, angle
                    
                elif difference > 0:
                    pos += 1
                    pos_frame_minus_landing.append(trajec.frame_of_landing - fs)
                    ax_pos.plot(np.log(trajec.angle_subtended_by_post[fd]), trajec.speed[fd], '.', color='purple')
                    ax_pos.plot(np.log(trajec.angle_subtended_by_post[fs]), trajec.speed[fs], '.', color='green')
                    #ax_pos.plot([np.log(trajec.angle_subtended_by_post[fd]), np.log(trajec.angle_subtended_by_post[fs])], [trajec.speed[fd], trajec.speed[fs]], '-', color='black')
                    ax_pos.plot(np.log(trajec.angle_subtended_by_post[0:trajec.frame_of_landing]), trajec.speed[0:trajec.frame_of_landing], '-', color='black', linewidth=0.5)
                    ax_deceleration.plot(np.log(trajec.angle_subtended_by_post[fd]), trajec.speed[fd], '.', color='green')
            except:
                none += 1
                ax_none.plot(np.log(trajec.angle_subtended_by_post[fd]), trajec.speed[fd], '.', color='purple')
                ax_deceleration.plot(np.log(trajec.angle_subtended_by_post[fd]), trajec.speed[fd], '.', color='purple')
            
            
    yticks = [0, 0.2, 0.4, 0.6, 0.8, 1.0]
    xticks = np.log(np.array([10, 20, 30, 60, 90, 180])*np.pi/180.).tolist()
    
    set_log_angle_ticks(ax_neg)
    set_log_angle_ticks(ax_pos)
    set_log_angle_ticks(ax_none)
    set_log_angle_ticks(ax_deceleration)
            
    fpl.adjust_spines(ax_saccades, ['left', 'bottom'], yticks=[-np.pi, -np.pi/2., 0, np.pi/2., np.pi], xticks=[-np.pi, -np.pi/2., 0, np.pi/2., np.pi], smart_bounds=True)
            
    fig_neg.savefig('decel_and_saccade_neg.pdf', format='pdf')
    fig_pos.savefig('decel_and_saccade_pos.pdf', format='pdf')    
    fig_none.savefig('decel_and_saccade_none.pdf', format='pdf')   
    fig_deceleration.savefig('decel_for_all.pdf', format='pdf')    
    fig_saccades.savefig('saccade_for_saccade_after_deceleration.pdf', format='pdf')
    
    
    print 'neg: ', neg, 'pos: ', pos, 'none: ', none
    print pos_frame_minus_landing
    
    
def set_log_angle_ticks(ax):
    yticks = [0, 0.2, 0.4, 0.6, 0.8, 1.0]
    xticks = np.log(np.array([5, 10, 20, 30, 60, 90, 180])*np.pi/180.).tolist()
    fpl.adjust_spines(ax, ['left', 'bottom'], yticks=yticks, xticks=xticks, smart_bounds=True)
    
    xticklabels = [str(i) for i in [5, 10, 20, 30, 60, 90, 180]]
    ax.set_xticklabels(xticklabels)
    
    ax.set_xlabel('Retinal size, deg')
    ax.set_ylabel('Speed, m/s')
    
    
def get_sac_ranges_in_frame_range(trajec):
    sac_ranges_in_range = []
    sac_ranges = trajec.sac_ranges
    for sac_range in sac_ranges:
        if sac_range[-1] in trajec.frames:
            sac_ranges_in_range.append(sac_range)
    return sac_ranges_in_range
    
def plot_deceleration_not_near_saccades(dataset, restrict_angle=False, show_decel_if_before_saccade=False):
    
    speed = []
    angle_subtended = []
    angle_to_post = []
    
    def plot(trajec):
        fd = trajec.frame_at_deceleration
        angleok = True
        if restrict_angle:
            if np.abs(trajec.angle_to_post[fd]) < 20.*np.pi/180.:
                angleok=True
            else:
                angleok=False 
        if angleok:
            trajec.classification = 'no_saccade_after_deceleration'
            speed.append(trajec.speed[fd])
            angle_subtended.append(trajec.angle_subtended_by_post[fd])
            angle_to_post.append(trajec.angle_to_post[fd])
    
    for k, trajec in dataset.trajecs.items():
        trajec.classification = 'saccade_after_deceleration'
        fd = trajec.frame_at_deceleration
        
        sac_ranges = get_sac_ranges_in_frame_range(trajec)
        
        if len(sac_ranges) > 0:
            sac_frames = []
            for sac_range in sac_ranges:
                sac_frames.extend(sac_range)
            decel_sac_diff = np.array(sac_frames) - fd
        
            # if no saccades after deceleration:
            if np.max(decel_sac_diff) <= 0:
                plot(trajec)
            
            else:
                indices = np.where(decel_sac_diff>0)[0].tolist()
                nearest_sac_to_fd_after_fd = np.min(decel_sac_diff[indices])
                if nearest_sac_to_fd_after_fd > 75:
                    plot(trajec)
        
        else:
            plot(trajec)
                
                
            
    speed = np.array(speed)
    angle_subtended = np.array(angle_subtended)
    angle_to_post = np.array(angle_to_post)
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    print speed.shape, angle_subtended.shape
    
    ax.scatter(np.log(angle_subtended), speed, color='purple', linewidth=0) #, color=angle_to_post, colornorm=(-90*np.pi/180., 90*np.pi/180.))
    
    
    
    # add points for trajecs with saccades after deceleration:
    if show_decel_if_before_saccade:
        keys = get_keys_for_landings_with_saccade_after_deceleration(dataset)
        speed = []
        angle_subtended = []
        angle_to_post = []
        for key in keys:
            trajec = dataset.trajecs[key]
            fd = trajec.frame_at_deceleration
            speed.append(trajec.speed[fd])
            angle_subtended.append(trajec.angle_subtended_by_post[fd])
            angle_to_post.append(trajec.angle_to_post[fd])
        speed = np.array(speed)
        angle_subtended = np.array(angle_subtended)
        angle_to_post = np.array(angle_to_post)
        ax.scatter(np.log(angle_subtended), speed, color='black', linewidth=0) #, color=angle_to_post, colornorm=(-90*np.pi/180., 90*np.pi/180.))
    
    set_log_angle_ticks(ax)
    
    fig.savefig('deceleration_for_flybys_not_near_saccades.pdf', format='pdf')
    
    
def get_keys_for_landings_with_saccade_after_deceleration(dataset_landing):
    keys = []
    for key, trajec in dataset_landing.trajecs.items():
        if trajec.classification == 'saccade_after_deceleration':
            keys.append(key)
    return keys
    
    
## deprecated
def plot_deceleration_for_flybys(dataset):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    fig_saccade_flyby = plt.figure()
    ax_saccade_flyby = fig_saccade_flyby.add_subplot(111)
    
    for k, trajec in dataset.trajecs.items():
        fd, fs = get_frame_of_decel_and_last_saccade(trajec)
        if np.abs(trajec.angle_to_post[fd]) < 180.*np.pi/180.:
            
            if 0:#fs is not None:
                
                if trajec.angle_subtended_by_post[fs] > 20*np.pi/180.:
                
                    if fs > fd:
                        ax.plot(np.log(trajec.angle_subtended_by_post[fd]), trajec.speed[fd], '.', color='green')
                        #ax.plot(np.log(trajec.angle_subtended_by_post[fs]), trajec.speed[fs], '.', color='green')
                        #ax.plot([np.log(trajec.angle_subtended_by_post[fd]), np.log(trajec.angle_subtended_by_post[fs])], [trajec.speed[fd], trajec.speed[fs]], '-', color='black')
                        sac_range = get_saccade_range_from_first_frame(trajec, fs)
                        angle = sac.get_angle_of_saccade(trajec, sac_range)
                        ax_saccade_flyby.plot(trajec.angle_to_post[fs], angle, '.', color='green')
                    else:
                        print trajec.key # saccade before deceleration
                        ax.plot(np.log(trajec.angle_subtended_by_post[fd]), trajec.speed[fd], '.', color='red')
                        sac_range = get_saccade_range_from_first_frame(trajec, fs)
                        angle = sac.get_angle_of_saccade(trajec, sac_range)
                        ax_saccade_flyby.plot(trajec.angle_to_post[fs], angle, '.', color='red')
                    
            if fs is not None:
                if np.abs(fs-fd) > 10 and np.abs(trajec.angle_to_post[fd]) < 60.*np.pi/180.:
                    ax.plot(np.log(trajec.angle_subtended_by_post[fd]), trajec.speed[fd], '.', color='purple')
                    
            
            
        
        # try looking for trajectories that slowed down after saccading?
    
    set_log_angle_ticks(ax)
    fig.savefig('deceleration_for_flybys.pdf', format='pdf') 
    
    fpl.adjust_spines(ax_saccade_flyby, ['left', 'bottom'], yticks=[-np.pi, -np.pi/2., 0, np.pi/2., np.pi], xticks=[-np.pi, -np.pi/2., 0, np.pi/2., np.pi], smart_bounds=True)
    
    ax_saccade_flyby.set_xlim([-np.pi, np.pi])
    ax_saccade_flyby.set_ylim([-np.pi, np.pi])
    
    fig_saccade_flyby.savefig('last_saccade_of_flyby.pdf', format='pdf')
    
    
    
    
# heatmap of straightlanding trajectories

def heatmap_of_straight_landing_trajecs(dataset_landing, dataset_flyby):

    speeds = []
    angles_subtended = []
    for k, trajec in dataset_landing.trajecs.items():
        if trajec.classification == 'no_saccade_after_deceleration':
            speeds.extend(trajec.speed[0:trajec.frame_of_landing].tolist())
            angles_subtended.extend(trajec.angle_subtended_by_post[0:trajec.frame_of_landing].tolist())
        
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    fpl.histogram2d(ax, np.log(np.array(angles_subtended)), np.array(speeds), bins=50, normed=False, histrange=None, weights=None, logcolorscale=True, colormap='jet', interpolation='bicubic')
    
    fig_saccades_flyby = plt.figure()
    ax_saccades_flyby = fig_saccades_flyby.add_subplot(111)
    
    
    # now look at speed and retinal size for flyby trajectories at last saccade if headed towards post:
    for k, trajec in dataset_flyby.trajecs.items():
        fd, fs = get_frame_of_decel_and_last_saccade(trajec)
        for sac_range in trajec.sac_ranges:
            if fs in sac_range:
                if trajec.angle_subtended_by_post[sac_range[0]] > 0.*np.pi/180.:
                    
                #if np.abs(trajec.angle_to_post[sac_range[0]]) < 180.*np.pi/180.:
                    ax.plot(np.log(trajec.angle_subtended_by_post[sac_range[0]]), trajec.speed[sac_range[0]], '.', markersize=2, color='white', markeredgecolor='black', linewidth=0.5)

                    angle = sac.get_angle_of_saccade(trajec, sac_range)
                    ax_saccades_flyby.plot(-1*trajec.angle_to_post[sac_range[0]], angle, '.', color='red', markersize=2)

    set_log_angle_ticks(ax)
    ax.set_aspect('auto')
        
    fig.savefig('landing_deceleration_heatmap.pdf')
        
    fpl.adjust_spines(ax_saccades_flyby, ['left', 'bottom'], yticks=[-np.pi, -np.pi/2., 0, np.pi/2., np.pi], xticks=[-np.pi, -np.pi/2., 0, np.pi/2., np.pi], smart_bounds=True)
    ax_saccades_flyby.set_xlim([-np.pi, np.pi])
    ax_saccades_flyby.set_ylim([-np.pi, np.pi])
    
    deg_ticks = ['-180', '-90', '0', '90', '180']
    ax_saccades_flyby.set_xticklabels(deg_ticks)
    ax_saccades_flyby.set_yticklabels(deg_ticks)
    
    ax_saccades_flyby.set_xlabel('Angle to post, deg')
    ax_saccades_flyby.set_ylabel('Turn angle, deg')
    
    fig_saccades_flyby.savefig('saccades_flyby.pdf', format='pdf')
        


def heatmap_of_flyby_trajecs(dataset_flyby):

    speeds = []
    angles_subtended = []
    for k, trajec in dataset_flyby.trajecs.items():
        speeds.extend(trajec.speed[0:trajec.frame_of_landing].tolist())
        angles_subtended.extend(trajec.angle_subtended_by_post[0:trajec.frame_of_landing].tolist())
        
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    fpl.histogram2d(ax, np.log(np.array(angles_subtended)), np.array(speeds), bins=100, normed=False, histrange=None, weights=None, logcolorscale=True, colormap='jet', interpolation='bicubic')
    
    set_log_angle_ticks(ax)
    ax.set_aspect('auto')
        
    fig.savefig('flyby_deceleration_heatmap.pdf')
        
