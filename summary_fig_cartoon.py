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





def landing(dataset_landing, movie_dataset, speed=0.2):

    behavior = 'landing'
    
    
    fps = 1000.
    dt = 1/fps
    r = 0.009565
    radius = r
    pos0 = [-0.2, 0]
    vel = [speed, 0]
    dvda = -0.2
    
    nf = 5000
    positions = np.zeros([nf, 2])
    positions[0] = pos0
    velocities = np.zeros([nf, 2])
    velocities[0] = vel
    speed = np.zeros([nf])
    speed[0] = np.linalg.norm(velocities[0])
    distance = np.zeros([nf])
    angle_subtended_by_post = np.zeros([nf])
    leg_ext = np.zeros([nf])
    frames = [0]
    frame_at_deceleration = None
    deceleration_initiated = False
    
    for f in range(1,nf): 
        if np.linalg.norm(positions[f-1])-radius <= 0.0001:
            landed = True
        else:
            landed = False
            
            
        if not landed:
            frames.append(f)
            positions[f] = positions[f-1] + velocities[f-1]*dt
            distance[f] = np.linalg.norm(positions[f]) - radius
            angle_subtended_by_post[f] = 2*np.arcsin( radius / (distance[f]+radius) )
            
            if f>5:
                #velocities[f,0] = -.21*np.log(angle_subtended_by_post[f])+.2
                da = np.log(angle_subtended_by_post[f])-np.log(angle_subtended_by_post[f-1])
                
                a = angle_subtended_by_post
                af = np.min([a[f], 3])
                exp0 = (a[f]-a[f-1])/dt #/ (-2.*np.tan(a[f]/2.))
                exp1 = (a[f-1]-a[f-2])/dt #/ (-2.*np.tan(a[f-1]/2.))
                
                m = -0.21/radius
                b = 0.159/radius
                expthreshold = (m*np.log(af)+b)*(2*np.tan(af/2.)*np.sin(af/2.))
        
                exp0 -= expthreshold
                exp1 -= expthreshold
                
                exp0 = np.max([exp0, 0])
                exp1 = np.max([exp1, 0])
                
                #c = -1*exp0 / 3500.
                
                dda = (exp1-exp0)/dt
                c = dda / 150000.
                print dda, velocities[f-1,0]
                
                c = np.min([c,0])
                v = np.max([speed[f-1] + c, 0.0])
                
                velocities[f,0] = v
            else:
                velocities[f] = velocities[f-1]

            speed[f] = np.linalg.norm(velocities[f])
            if speed[f] > -0.21*np.log(angle_subtended_by_post[f])+0.159:
                deceleration_initiated = True
                if frame_at_deceleration is None:
                    frame_at_deceleration = f
            else:
                deceleration_initiated = False
                
            if angle_subtended_by_post[f]*180/np.pi > 70 or np.isnan(angle_subtended_by_post[f]):
                leg_ext[f] = 1

    
    
            
            
    fig2 = plt.figure()  
    ax2 = fig2.add_subplot(111)
    
    fit, Rsq, x, y, yminus, yplus = fa.get_angle_vs_speed_curve(dataset_landing, plot=False)
    ax2.plot( x, y, color='blue', alpha=0.3)
    ax2.fill_between(x, yplus, yminus, color='blue', linewidth=0, alpha=0.2)
    
    angle_at_leg_extension, bins, data_filtered, xvals = fa.leg_extension_angle_histogram(movie_dataset, plot=False)
    ax2.plot(xvals, data_filtered, color='red', alpha=0.3)
    ax2.fill_between(xvals, data_filtered, np.zeros_like(xvals), color='red', linewidth=0, alpha=0.2)
    
    ax2.plot(np.log(angle_subtended_by_post), speed, color='black')
    
    fa.fix_angle_log_spine(ax2, histograms=False) 
    
    fig2.subplots_adjust(bottom=0.3, top=0.8, right=0.9, left=0.25)
    
    filename = 'landing_cartoon_plot.pdf'
    fig2.savefig(filename, format='pdf')
            
            





def flyby(dataset_flyby, speed=0.5):

    behavior = 'flyby'
    
    
    fps = 100.
    dt = 1/fps
    r = 0.009565
    radius = r
    pos0 = [-0.2, 0]
    vel = [speed, 0]
    dvda = -0.2
    
    nf = 200
    positions = np.zeros([nf, 2])
    positions[0] = pos0
    velocities = np.zeros([nf, 2])
    velocities[0] = vel
    speed = np.zeros([nf])
    speed[0] = np.linalg.norm(velocities[0])
    distance = np.zeros([nf])
    angle_subtended_by_post = np.zeros([nf])
    frames = [0]
    frame_at_deceleration = None
    frame_at_saccade = None
    deceleration_initiated = False
    saccade_time = 0
    saccading = False
    
    for f in range(1,nf): 

        frames.append(f)
        positions[f] = positions[f-1] + velocities[f-1]*dt
        distance[f] = np.linalg.norm(positions[f]) - radius
        angle_subtended_by_post[f] = 2*np.arcsin( radius / (distance[f]+radius) )
        
        if saccading:
            saccade_rate = 400*np.pi/180.
            s = np.linalg.norm(velocities[f-1])
            velocities[f,1] = np.sin(saccade_rate)*s
            velocities[f,0] = np.cos(saccade_rate)*s
        elif deceleration_initiated:
            worldangle = np.arctan2(velocities[f-1,1], velocities[f-1,0])
            da = np.log(angle_subtended_by_post[f])-np.log(angle_subtended_by_post[f-1])
            velocities[f,1] = velocities[f-1,1] + dvda*da*np.sin(worldangle)
            velocities[f,0] = velocities[f-1,0] + dvda*da*np.cos(worldangle)
        else:
            velocities[f] = velocities[f-1]

        speed[f] = np.linalg.norm(velocities[f])
        if speed[f] > -0.1*np.log(angle_subtended_by_post[f])+0.212 and saccade_time == 0:
            deceleration_initiated = True
            print 'decelerating triggered'
            if frame_at_deceleration is None:
                frame_at_deceleration = f
        else:
            pass
            #deceleration_initiated = False
            
        if angle_subtended_by_post[f]*180/np.pi > 30 or np.isnan(angle_subtended_by_post[f]):
            if frame_at_saccade is None:
                frame_at_saccade = f
                saccading = True
            saccade_time += dt
            if saccade_time > 0.17:
                saccading = False
            
        if saccade_time > 0:
            deceleration_initiated = False
            
             

            
                
                
    # plot
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    flycon_flying = plt.imread('flying_flycon.png')
    
    for f in frames:
        '''
        if not leg_ext[f]:
            color = 'black'
        else:
            color = 'red'            
        if f == frame_at_deceleration:
            color = 'blue'
        '''
        color='black'
                 
        x = positions[f,0]
        y = positions[f,1]
        w = flycon_flying.shape[0]/(7*1000.*.8)
        h = flycon_flying.shape[1]/(5.*1000.*.8)
        
        if f == frame_at_deceleration:
            #ax.imshow(flycon_flying, extent=(x-2*w, x, y-h, y+h), aspect='equal')
            
            # show angle subtended
            d = np.linalg.norm([x-.001,y])
            half_angle_to_post = np.arcsin( radius / d )
            world_angle = np.arctan2(y,x-.001)
            
            a = half_angle_to_post - world_angle
            visual_intercept_1 = [0+np.cos(np.pi/2.-a)*radius, 0+np.sin(np.pi/2.-a)*radius]
            
            a = half_angle_to_post + world_angle
            visual_intercept_2 = [0+np.cos(np.pi/2.-a)*radius, 0-np.sin(np.pi/2.-a)*radius]
            
            xy = np.vstack( (visual_intercept_1, visual_intercept_2, [x-.001,y]) )
            triangle = patches.Polygon( xy, facecolor='blue', edgecolor='none', zorder=-10, alpha=0.2 )
            ax.add_artist(triangle)
            
            '''
            arc = patches.Arc( (x,y), .05, .05, 180, (world_angle - half_angle_to_post)*180/np.pi, (half_angle_to_post + world_angle)*180/np.pi, edgecolor='blue', linewidth=1)
            ax.add_artist(arc)
            '''
            
        elif f == frame_at_saccade:
            #ax.imshow(flycon_flying, extent=(x-2*w, x, y-h, y+h), aspect='equal')
            
            # show angle subtended
            d = np.linalg.norm([x-.001,y])
            half_angle_to_post = np.arcsin( radius / d )
            world_angle = np.arctan2(y,x-.001)
            
            a = half_angle_to_post - world_angle
            visual_intercept_1 = [0+np.cos(np.pi/2.-a)*radius, 0+np.sin(np.pi/2.-a)*radius]
            
            a = half_angle_to_post + world_angle
            visual_intercept_2 = [0+np.cos(np.pi/2.-a)*radius, 0-np.sin(np.pi/2.-a)*radius]
            
            xy = np.vstack( (visual_intercept_1, visual_intercept_2, [x-.001,y]) )
            triangle = patches.Polygon( xy, facecolor='green', edgecolor='none', zorder=-10, alpha=0.2 )
            ax.add_artist(triangle)
        
            
        else:
            pt = patches.Circle( (positions[f,0],positions[f,1]), radius=0.0005, facecolor=color, edgecolor='none')
            ax.add_artist(pt)  
            
    
         
    # post      
    post = patches.Circle( (0,0), radius=radius, facecolor='black', edgecolor='none')
    ax.add_artist(post)
            
    
    
    
            
    ax.set_aspect('equal')
    ax.set_xlim([-0.2, 0.01])
    ax.set_ylim([-0.08, 0.08])
    ax.set_axis_off()
    fig.set_size_inches(6.5,6.5)
    fig.subplots_adjust(bottom=0., top=1, right=0.95, left=0.05)
    filename = 'flyby_cartoon.pdf'
    fig.savefig(filename, format='pdf')
            
            
            
            
            
    fig2 = plt.figure()  
    ax2 = fig2.add_subplot(111)
    
    angle, bins, data_filtered, xvals, curve, yplus, yminus = fa.deceleration_angle_histogram_flyby(dataset_flyby, keys=None, plot=False, saccades=True)
    ax2.plot( curve[:,0], curve[:,1], color='blue', alpha=0.3)
    ax2.fill_between(curve[:,0], yplus, yminus, color='blue', linewidth=0, alpha=0.2)
    
    angle, bins, data_filtered, xvals = fa.saccade_angle_histogram(dataset_flyby, keys=None, plot=False)
    ax2.plot(xvals, data_filtered, color='green', alpha=0.3)
    ax2.fill_between(xvals, data_filtered, np.zeros_like(xvals), color='green', linewidth=0, alpha=0.2)
    
    ax2.plot(np.log(angle_subtended_by_post), speed, color='black')
    ax2.plot(np.log(angle_subtended_by_post[frame_at_saccade]), speed[frame_at_saccade], '.', color='green')
    
    fa.fix_angle_log_spine(ax2, histograms=False) 
    
    fig2.subplots_adjust(bottom=0.3, top=0.8, right=0.9, left=0.25)
    
    filename = 'flyby_cartoon_plot.pdf'
    fig2.savefig(filename, format='pdf') 
            
            


