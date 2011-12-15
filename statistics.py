import floris
import saccade_analysis as sac
import numpy as np
import scipy.stats
import flydra_analysis as fa

import matplotlib.pyplot as plt

def ANCOVA_for_deceleration(dataset):
    
    black_decel_speed = []
    black_decel_angle = []
    checkered_decel_speed = []
    checkered_decel_angle = []
    
    keys = dataset.trajecs.keys()

    for key in keys:
        trajec = dataset.trajecs[key]
        if 'black' in trajec.post_type:
            black_decel_speed.append(trajec.speed_at_deceleration)
            black_decel_angle.append(trajec.angle_at_deceleration)
        elif 'checkered' in trajec.post_type:
            checkered_decel_speed.append(trajec.speed_at_deceleration)
            checkered_decel_angle.append(trajec.angle_at_deceleration)

    black_decel_speed = np.array(black_decel_speed) 
    black_decel_angle = np.array(black_decel_angle)
    checkered_decel_speed = np.array(checkered_decel_speed)
    checkered_decel_angle = np.array(checkered_decel_angle)
    groups = [[np.log(black_decel_angle), black_decel_speed], [np.log(checkered_decel_angle), checkered_decel_speed]]

    floris.ANCOVA(groups)

def t_test_for_evasive_saccades(dataset):
    sac.calc_saccade_stats(dataset)
    
    black_saccade_angle = []
    checkered_saccade_angle = []
    
    keys = dataset.keys_saccade_array
    
    evasive_threshold = 0.836248003234

    for i, key in enumerate(keys):
        trajec = dataset.trajecs[key]
        a = dataset.angle_of_saccade_array[i]
        ap = dataset.angle_prior_array[i]
        
        if np.abs(a-ap) > evasive_threshold:
            if 'black' in trajec.post_type:
                black_saccade_angle.append(dataset.angle_subtended_array[i])
            elif 'checkered' in trajec.post_type:
                checkered_saccade_angle.append(dataset.angle_subtended_array[i])

    black_saccade_angle = np.array(black_saccade_angle)
    checkered_saccade_angle = np.array(checkered_saccade_angle)
    
    if 1:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        bins = np.linspace( np.log(5*np.pi/180.), np.log(180*np.pi/180.), 20)
        bins, hists, curves = floris.histogram(ax, [np.log(black_saccade_angle), np.log(checkered_saccade_angle)], bins=bins, colors=['black', 'teal'], return_vals=True, normed_occurences=True, curve_line_alpha=0)

        deg_ticks = np.array([5, 10, 30, 60, 90, 180])
        deg_tick_strings = [str(d) for d in deg_ticks]
        rad_ticks = deg_ticks*np.pi/180.
        rad_ticks_log = np.log(rad_ticks)
        
        dist_tick_strings = ['(21)', '(10)', '(2.7)', '(0.9)', '(0.4)', '(0)']
        x_tick_strings = []
        for i, d in enumerate(dist_tick_strings):
            x_tick_strings.append( deg_tick_strings[i] + '\n' + dist_tick_strings[i] )
        
        ax.set_xlim(rad_ticks_log[0], rad_ticks_log[-1])
        ax.set_ylim(0,1)
        
        fa.adjust_spines(ax,['left', 'bottom'])
        ax.set_xlabel('Retinal size, deg\n(distance, cm)', labelpad=10)
        ax.set_ylabel('Occurances')
        
        ax.set_xticks(rad_ticks_log.tolist())
        ax.set_xticklabels(x_tick_strings) 
        
        
        fig.subplots_adjust(bottom=0.25, top=0.9, right=0.9, left=0.2)
        filename = 'saccade_histogram_black_vs_checkered' + '.pdf'
        fig.savefig(filename, format='pdf')

    print 't test: ', scipy.stats.ttest_ind(black_saccade_angle, checkered_saccade_angle)
    print 'mann whitney u: ', scipy.stats.mannwhitneyu(black_saccade_angle, checkered_saccade_angle)
    print 'black: ', np.mean(black_saccade_angle), 'n: ', len(black_saccade_angle)
    print 'checkered: ', np.mean(checkered_saccade_angle), 'n: ', len(checkered_saccade_angle)
    print 'total n - 2: ', len(black_saccade_angle) + len(checkered_saccade_angle) - 2
    
def ANCOVA_for_attractive_saccades(dataset):
    sac.calc_saccade_stats(dataset)
    
    black_saccade_angle = []
    black_angle_prior = []
    checkered_saccade_angle = []
    checkered_angle_prior = []
    
    keys = dataset.keys_saccade_array
    
    evasive_threshold = 0.836248003234

    for i, key in enumerate(keys):
        trajec = dataset.trajecs[key]
        a = dataset.angle_of_saccade_array[i]
        ap = dataset.angle_prior_array[i]
        
        if 1:#np.abs(a-ap) < evasive_threshold:
        
            if 'black' in trajec.post_type:
                black_saccade_angle.append(a)
                black_angle_prior.append(ap)
                
            elif 'checkered' in trajec.post_type:
                checkered_saccade_angle.append(a)
                checkered_angle_prior.append(ap)

    black_saccade_angle = np.array(black_saccade_angle)
    black_angle_prior = np.array(black_angle_prior)
    checkered_saccade_angle = np.array(checkered_saccade_angle) 
    checkered_angle_prior = np.array(checkered_angle_prior) 

    groups = [[black_angle_prior, black_saccade_angle], [checkered_angle_prior, checkered_saccade_angle]]

    floris.ANCOVA(groups)    
    
    
def t_test_for_leg_ext(movie_dataset, behavior='landing'):
    
    black_legext_angle = []
    checkered_legext_angle = []
    
    keys = movie_dataset.movies.keys()
    
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
                
                if 'crash' not in movie.subbehavior and tmp and legextensionframe < movie.landingframe_relative: #tmp: #movie.trajec.classification == 'straight':
                    if 'black' in movie.posttype:
                        black_legext_angle.append(movie.scaled.angle_subtended_by_post[legextensionframe][0])
                    elif 'checkered' in movie.posttype:
                        checkered_legext_angle.append(movie.scaled.angle_subtended_by_post[legextensionframe][0])
                    n += 1
        
    black_legext_angle = np.array(black_legext_angle)
    checkered_legext_angle = np.array(checkered_legext_angle)
    
    if 1:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        bins = np.linspace( np.log(5*np.pi/180.), np.log(180*np.pi/180.), 16)
        bins, hists, curves = floris.histogram(ax, [np.log(black_legext_angle), np.log(checkered_legext_angle)], bins=bins, colors=['black', 'teal'], return_vals=True, normed_occurences=False, curve_line_alpha=0, show_smoothed=True)

        deg_ticks = np.array([5, 10, 30, 60, 90, 180])
        deg_tick_strings = [str(d) for d in deg_ticks]
        rad_ticks = deg_ticks*np.pi/180.
        rad_ticks_log = np.log(rad_ticks)
        
        dist_tick_strings = ['(21)', '(10)', '(2.7)', '(0.9)', '(0.4)', '(0)']
        x_tick_strings = []
        for i, d in enumerate(dist_tick_strings):
            x_tick_strings.append( deg_tick_strings[i] + '\n' + dist_tick_strings[i] )
        
        ax.set_xlim(rad_ticks_log[0], rad_ticks_log[-1])
        ax.set_ylim(0,10)
        
        fa.adjust_spines(ax,['left', 'bottom'])
        ax.set_xlabel('Retinal size, deg\n(distance, cm)', labelpad=10)
        ax.set_ylabel('Occurances')
        
        ax.set_xticks(rad_ticks_log.tolist())
        ax.set_xticklabels(x_tick_strings) 
        
        
        fig.subplots_adjust(bottom=0.25, top=0.9, right=0.9, left=0.2)
        filename = 'legext_histogram_black_vs_checkered' + '.pdf'
        fig.savefig(filename, format='pdf')
        
    print 'T-test', scipy.stats.ttest_ind(black_legext_angle, checkered_legext_angle)
    print 'Mann Whitney U', scipy.stats.mannwhitneyu(black_legext_angle, checkered_legext_angle)
    print 'black: ', np.mean(black_legext_angle), 'n: ', len(black_legext_angle)
    print 'checkered: ', np.mean(checkered_legext_angle), 'n: ', len(checkered_legext_angle)
    print 'total n - 2: ', len(black_legext_angle) + len(checkered_legext_angle) - 2
    
    
    
    
def get_speeds_at_dist(dataset, distance, post_type=['black', 'black_angled', 'checkered', 'checkered_angled', 'none']):    
    speeds = []    
    for k, trajec in dataset.trajecs.items():
        if trajec.post_type in post_type:
            s = fa.get_speed_at_distance(trajec, distance, return_none_if_no_frame=True)
            if s is not None:
                speeds.append( s )
    return speeds

def plot_speeds_at_dist(dataset_landing, dataset_flyby, dataset_nopost, distance=None):
    
    if distance is None:
        distance = np.linspace(0.15, 0.01, 15, endpoint=True)
    
    mean_landing_speeds_black = []
    std_landing_speeds_black = []
    mean_landing_speeds_checkered = []
    std_landing_speeds_checkered = []
    mean_flyby_speeds_black = []
    mean_flyby_speeds_checkered = []
    mean_nopost_speeds = []
    
    histfig = plt.figure()
    histax = histfig.add_subplot(111) 
    
    for d in distance:
        landing_speeds_black = get_speeds_at_dist(dataset_landing, d, post_type=['black', 'black_angled'])
        mean_landing_speeds_black.append(np.mean(landing_speeds_black))
        std_landing_speeds_black.append(np.std(landing_speeds_black))
        landing_speeds_checkered = get_speeds_at_dist(dataset_landing, d, post_type=['checkered', 'checkered_angled'])
        mean_landing_speeds_checkered.append(np.mean(landing_speeds_checkered))
        std_landing_speeds_checkered.append(np.std(landing_speeds_checkered))
        
        print 'landing: ', scipy.stats.mannwhitneyu(landing_speeds_black, landing_speeds_checkered), np.mean(landing_speeds_black), np.mean(landing_speeds_checkered)
        
        bins = np.linspace(0.1, 0.5, 15)
        bins, hists, hist_std, curves = floris.histogram(histax, [landing_speeds_black, landing_speeds_checkered], bins=bins, colors=['black', 'teal'], bin_width_ratio=0.9, edgecolor='none', bar_alpha=0.8, curve_fill_alpha=0.2, curve_line_alpha=0, return_vals=True, normed_occurences=False, bootstrap_std=True)
        
        flyby_speeds_black = get_speeds_at_dist(dataset_flyby, d, post_type=['black', 'black_angled'])
        mean_flyby_speeds_black.append(np.mean(flyby_speeds_black))
        flyby_speeds_checkered = get_speeds_at_dist(dataset_flyby, d, post_type=['checkered', 'checkered_angled'])
        mean_flyby_speeds_checkered.append(np.mean(flyby_speeds_checkered))
        
        print 'flyby: ', scipy.stats.mannwhitneyu(flyby_speeds_black, flyby_speeds_checkered), np.mean(flyby_speeds_black), np.mean(flyby_speeds_checkered)
        
        nopost_speeds = get_speeds_at_dist(dataset_nopost, d)
        mean_nopost_speeds.append(np.mean(nopost_speeds))
        
        print 'flyby vs landing (black): ', scipy.stats.mannwhitneyu(flyby_speeds_black, landing_speeds_black)
        print 'flyby vs landing (checkered): ', scipy.stats.mannwhitneyu(flyby_speeds_checkered, landing_speeds_checkered)
        print 'flyby (black) vs landing (checkered): ', scipy.stats.mannwhitneyu(flyby_speeds_black, landing_speeds_checkered)
        print 'flyby (checkered) vs landing (black): ', scipy.stats.mannwhitneyu(flyby_speeds_checkered, landing_speeds_black)
        
        print 'flyby vs nopost (black): ', scipy.stats.mannwhitneyu(flyby_speeds_black, nopost_speeds)
        print 'flyby vs nopost (checkered): ', scipy.stats.mannwhitneyu(flyby_speeds_checkered, nopost_speeds)
        print 'landing vs nopost (black): ', scipy.stats.mannwhitneyu(landing_speeds_black, nopost_speeds)
        print 'landing vs nopost (checkered): ', scipy.stats.mannwhitneyu(landing_speeds_checkered, nopost_speeds)
        
        print
        print 'mean flight speed black landing: ', np.mean(landing_speeds_black), np.std(landing_speeds_black)
        print 'n: ', len(landing_speeds_black)
        all_others = []
        all_others.extend(landing_speeds_checkered)
        all_others.extend(flyby_speeds_black)
        all_others.extend(flyby_speeds_checkered)
        all_others.extend(nopost_speeds)
        print 'mean flight all others: ', np.mean(all_others), np.std(all_others)
        print 'n: ', len(all_others)

    mean_landing_speeds_black = np.array(mean_landing_speeds_black)
    std_landing_speeds_black = np.array(std_landing_speeds_black) 
    mean_landing_speeds_checkered = np.array(mean_landing_speeds_checkered) 
    std_landing_speeds_checkered = np.array(std_landing_speeds_checkered)
    mean_flyby_speeds_black = np.array(mean_flyby_speeds_black) 
    mean_flyby_speeds_checkered = np.array(mean_flyby_speeds_checkered) 
    mean_nopost_speeds = np.array(mean_nopost_speeds) 

        
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    ax.plot(distance, mean_landing_speeds_black, 'black')
    ax.fill_between(distance, mean_landing_speeds_black+std_landing_speeds_black, mean_landing_speeds_black-std_landing_speeds_black, alpha=0.2, color='black', edgecolor='none')
    ax.fill_between(distance, mean_landing_speeds_checkered+std_landing_speeds_checkered, mean_landing_speeds_checkered-std_landing_speeds_checkered, alpha=0.2, color='teal', edgecolor='none')
    
    ax.plot(distance, mean_landing_speeds_checkered, '-', color='teal')
    ax.plot(distance, mean_flyby_speeds_black, 'red')
    ax.plot(distance, mean_flyby_speeds_checkered, ':', color='red')
    ax.plot(distance, mean_nopost_speeds, 'blue')
    
    ax.set_xlabel('distance, m')
    ax.set_ylabel('speed, m/s')
    
    fig.savefig('flight_speeds', format='pdf')
    histfig.savefig('flight_speeds_hist', format='pdf')
    
    
    
def t_test_speed_at_deceleration(dataset):
    
    fa.calc_stats_at_deceleration(dataset, keys=None, time_offset=0, return_vals=False)
    angleok = np.where( (dataset.angle_at_deceleration < 2)*(dataset.angle_at_deceleration > .01) )[0].tolist()
    
    
    black_post = []
    checkered_post = []
    for i, p in enumerate(dataset.post_type_at_deceleration):
        if i in angleok:
            if 'black' in p:
                black_post.append(i)
            if 'checkered' in p:
                checkered_post.append(i)
            
    print '*'*40
    print len(black_post), len(checkered_post)
    
    
    print scipy.stats.ttest_ind(dataset.speed_at_deceleration[black_post], dataset.speed_at_deceleration[checkered_post])
    print 'black: ', np.mean(dataset.speed_at_deceleration[black_post]), 'n: ', len(dataset.speed_at_deceleration[black_post])
    print 'checkered: ', np.mean(dataset.speed_at_deceleration[checkered_post]), 'n: ', len(dataset.speed_at_deceleration[checkered_post])
    print 'total n - 2: ', len(black_post) + len(checkered_post) - 2
    
    
    
    
    
    

