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

import copy

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib

import scipy.signal as signal

import colorline
import flydra_analysis as fa
import sa1_analysis as sa1
import numpy as np
import floris



def get_angle_of_saccade(trajec, s):
    sac_range = fa.get_saccade_range(trajec, s)
    if sac_range is None:
        return None
    
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
    
    
    
def calc_saccade_stats(dataset, keys=None):
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
            
    for key in keys:
        trajec = dataset.trajecs[key]
        #sa1.calc_post_dynamics_for_flydra_trajectory(trajec)
        
        if len(trajec.all_saccades) > 0:
            try:
                frames = trajec.frames_below_post
            except:
                frames = np.arange(0, trajec.frame_nearest_to_post).tolist()
        
            for s in trajec.all_saccades:
                if s in frames:
                    sac_range = fa.get_saccade_range(trajec, s)
                    if sac_range is None:
                        continue
                    if sac_range[0] <= trajec.frame_nearest_to_post:
                        angle_prior = trajec.angle_to_post[sac_range[0]]*-1
                        signed_angleofsaccade = get_angle_of_saccade(trajec, s)
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
                            
    dataset.speed_at_saccade_array = np.array(speed_at_saccade_array)
    dataset.angle_of_saccade_array = np.array(angle_of_saccade_array)
    dataset.flipped_angle_of_saccade_array = np.array(flipped_angle_of_saccade_array)
    dataset.angle_prior_array = np.array(angle_prior_array)
    dataset.flipped_angle_prior_array = np.array(flipped_angle_prior_array)
    dataset.angle_subtended_array = np.array(angle_subtended_array)
    dataset.attractiveness_array = np.array(attractiveness_array)
    dataset.dist_at_saccade_array = np.array(dist_at_saccade_array)
    dataset.post_types_saccade_array = post_types_saccade_array
    
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
                alpha = colormap_norm_dist(0.2 - dataset.dist_at_saccade_array[i])
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
    
    
def plot_saccade_angle(dataset, keys=None, angle_subtended_range=[0,180], retinal_position_range=[-160, 160], show_evasive_thresholds=False, show_tanh_evasive_thresholds=False, edge='far', colorscale=False, speedscale=False, regress=None, plot_regress=False, plot=True, plot_color_regress=False, mirror=False, postcolorscale=False, return_left_right=False, recalc_saccade_stats=False):
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
                    alpha = colormap_norm_dist(0.2 - dataset.dist_at_saccade_array[i])
                elif postcolorscale is True:
                    if 'checkered' in dataset.post_types_saccade_array[i]:
                        c = 'red'
                        alpha = 0.8
                    elif 'black' in dataset.post_types_saccade_array[i]:
                        c = 'blue'
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
    
    bins, hists, curves = floris.histogram(ax, [np.log(angle_subtended[evasive]), np.log(angle_subtended[nonevasive])], bins=bins, colors=['orange', 'blue'], return_vals=True)
    
    #ax.hist( np.log(angle_subtended[evasive]), bins=bins, facecolor='orange', edgecolor='none', alpha=0.5, normed=False)
    #ax.hist( np.log(angle_subtended[nonevasive]), bins=bins, facecolor='blue', edgecolor='none', alpha=0.5, normed=False)
    
    evasive_hist = hists[0]
    nonevasive_hist = hists[1]
    
    ymax = np.max([np.max(evasive_hist), np.max(nonevasive_hist)])
    
    deg_ticks = np.array([5, 10, 30, 60, 90, 180])
    deg_tick_strings = [str(d) for d in deg_ticks]
    rad_ticks = deg_ticks*np.pi/180.
    rad_ticks_log = np.log(rad_ticks)
    
    dist_tick_strings = ['(21)', '(10)', '(2.7)', '(0.9)', '(0.4)', '(0)']
    x_tick_strings = []
    for i, d in enumerate(dist_tick_strings):
        x_tick_strings.append( deg_tick_strings[i] + '\n' + dist_tick_strings[i] )
    
    ax.set_xlim(rad_ticks_log[0], rad_ticks_log[-1])
    ax.set_ylim(0,ymax)
    
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
    
    plot_saccade_angle(dataset_landing, show_evasive_thresholds=False, edge='far', retinal_position_range=[-160,160], regress=dataset_landing, plot_regress=False, colorscale=True, speedscale=True, plot_color_regress=True)
    plot_saccade_angle(dataset_flyby, show_evasive_thresholds=False, edge='far', retinal_position_range=[-160,160], regress=dataset_landing, plot_regress=False, colorscale=True, speedscale=True, plot_color_regress=True)
    plot_saccade_angle(dataset_nopost, show_evasive_thresholds=False, edge='far', retinal_position_range=[-160,160], regress=dataset_landing, plot_regress=False, colorscale=True, speedscale=True, plot_color_regress=True)    
    
def save_post_type_far_edge_all_behaviors(dataset_landing, dataset_flyby, dataset_nopost):
    
    plot_saccade_angle(dataset_landing, show_evasive_thresholds=False, edge='far', retinal_position_range=[-160,160], regress=dataset_landing, plot_regress=False, colorscale=False, speedscale=False, plot_color_regress=True, postcolorscale=True)
    
    plot_saccade_angle(dataset_flyby, show_evasive_thresholds=False, edge='far', retinal_position_range=[-160,160], regress=dataset_landing, plot_regress=False, colorscale=False, speedscale=False, plot_color_regress=True, postcolorscale=True)
    
    plot_saccade_angle(dataset_nopost, show_evasive_thresholds=False, edge='far', retinal_position_range=[-160,160], regress=dataset_landing, plot_regress=False, colorscale=False, speedscale=False, plot_color_regress=True, postcolorscale=True)    
    
def save_flyby_evasive_fit(dataset_flyby):

    plot_saccade_angle(dataset_flyby, show_evasive_thresholds=False, edge='far', retinal_position_range=[-160,160], angle_subtended_range=[20,180], regress=dataset_flyby, plot_regress=True, colorscale=True, speedscale=True, plot_color_regress=False, mirror=True)


def save_left_right_turn_histogram(dataset_flyby, mirror=False):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    left, right = plot_saccade_angle(dataset_flyby, show_evasive_thresholds=False, edge='far', retinal_position_range=[-160,160], angle_subtended_range=[20,180], regress=dataset_flyby, plot_regress=True, colorscale=True, speedscale=True, plot_color_regress=False, mirror=True, return_left_right=True)
    
    print '*** left ***'
    print left
    
    bins = np.linspace(-160*np.pi/180., 160*np.pi/180., 14, endpoint=True)
    bins, hists, curves = floris.histogram(ax, [left, right], bins=bins, colors=['crimson', 'green'], return_vals=True, show_smoothed=False)
    
    left_hist = np.array(hists[0], dtype=float)
    right_hist = np.array(hists[1], dtype=float)
    
    if mirror:
        percent_turn_left = (left_hist + right_hist[::-1])  / (left_hist+right_hist+right_hist[::-1]+left_hist[::-1])
        percent_turn_right = (right_hist + left_hist[::-1])  / (left_hist+right_hist+right_hist[::-1]+left_hist[::-1])
    else:
        percent_turn_left = left_hist / (left_hist+right_hist)
        percent_turn_right = right_hist / (left_hist+right_hist)
    
    ymax = np.max([np.max(left_hist), np.max(right_hist)])
    ax.set_ylim(0, 60)
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
    percax.plot(bincenters, percent_turn_left, color='crimson')
    percax.plot(bincenters, percent_turn_right, color='green')
    
    percax.set_ylim(0, 1)
    percax.set_xlim(-np.pi, np.pi)
    percax.set_autoscale_on(False)
    rad_ticks = [-np.pi, -np.pi/2., 0, np.pi/2., np.pi]
    deg_tick_strings = ['-180', '-90', '0', '90', '180']
    fa.adjust_spines(percax, ['left', 'bottom'])
    percax.set_xticks(rad_ticks)
    percax.set_xticklabels(deg_tick_strings) 
    percax.set_xlabel('angle to far edge')
    percax.set_ylabel('probability of turn direction')
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
    
    
    
    
    
    
        
    
    
    
    
