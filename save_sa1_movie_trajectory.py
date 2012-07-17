def plot_trajectory(movieinfo, figure=None, show_post=True, show_flydra=True):
    fig = plt.figure(figure)
    ax0 = fig.add_axes([.1,.1,.8,.8])
    post_pos = movieinfo.post_pos
    frames = sa1.get_frames_until_landing(movieinfo)
    
    ax0.plot(movieinfo.scaled.positions[frames,0], movieinfo.scaled.positions[frames,1])
    if show_flydra:
        ax0.plot(movieinfo.trajec.positions[:,0], movieinfo.trajec.positions[:,1], '.', markersize=2)
    
    ax0.set_xlim([-.1, .1])
    ax0.set_ylim([-.1, .1])
    ax0.set_aspect('equal')
    
    strobe = sa1.strobe_from_movieinfo(movieinfo, interval=200)
    ax0.imshow(strobe.T, plt.get_cmap('gray'), origin='lower', extent = [-1*movieinfo.post_pos[0]*movieinfo.scale, (1024-movieinfo.post_pos[0])*movieinfo.scale, -1*movieinfo.post_pos[1]*movieinfo.scale, (1024-movieinfo.post_pos[1])*movieinfo.scale])
    
    # plot body orientation vector
    interval = 50
    i = 0
    while i < frames[-1]:
        ax0.plot(movieinfo.scaled.head_pos[i,0],movieinfo.scaled.head_pos[i,1], '.', color='black', zorder=10, markersize=2)
        center = movieinfo.scaled.positions[i]
        long_axis = movieinfo.smooth_ori[i]
        
        factor = .001
                
        dx = long_axis[0]*factor
        dy = long_axis[1]*factor
        
        arrow = Arrow(center[0], center[1], dx, dy, width=.0001, color='red')
        ax0.add_artist(arrow)
        
        i += interval        
        
    if show_post:
        circle = patches.Circle( (0,0), radius=movieinfo.post_radius, facecolor='none', edgecolor='green')
        ax0.add_artist(circle)
    
    title = movieinfo.id + ' ' + movieinfo.behavior + ' ' + movieinfo.subbehavior[0]
    ax0.set_title(title)
    
    
    
    return ax0
