# this is an example file that shows how to generate a plot in matplotlib, and save the file to a layered scalable vector graphics pdf that can be read and edited in Adobe Illustrator. usage: either run the file as python code, or import this file and run function example()

# Known issues:
# so far, I have not found a way to actually make the font carry through properly.. if you find a way, let me know! in the end, this is not a huge deal, as the editable text can be set to the appropriate font in illustrator rather easily.

# written by Floris van Breugel, 6/2011

#############################
# default figure parameters #
fig_width = 4 # width in inches
fig_height = 4  # height in inches
fontsize = 9
fig_size =  (fig_width, fig_height)
#############################

# set defaults
# note: this needs to be performed BEFORE importing pyplot or pylab. If you reload... reload pyplot/pylab after as well. 
# all of these paramters, except 'backend', 'ps.usedistiller', 'text.usetex', can be adjusted elsewhere.. changing them here just makes it convenient if you are doing many figures with the same settings. 
from matplotlib import rcParams
params = {'backend': 'Agg', # must be Agg, PDF, or PS. Only Agg allows for gui-viewing plots in ipython
          'ps.usedistiller': 'xpdf', # may need to apt-get install xpdf
          'ps.fonttype' : 3, # 3: approximate font, 42: exact font
          'pdf.fonttype' : 3,
          'font.family' : 'sans-serif',
          'font.serif' : 'Times, Palatino, New Century Schoolbook, Bookman, Computer Modern Roman',
          'font.sans-serif' : 'Helvetica, Avant Garde, Computer Modern Sans serif', # uses the first font in the series
          'font.cursive' : 'Zapf Chancery',
          'font.monospace' : 'Courier, Computer Modern Typewriter',
          'font.size' : fontsize,
          'text.fontsize': fontsize,
          'axes.labelsize': fontsize,
          'axes.linewidth': 1.0,
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
          'figure.subplot.right': 0.9,
          'figure.subplot.bottom': 0.2,
          'figure.subplot.top': 0.9,
          'figure.subplot.wspace': 0.0,
          'figure.subplot.hspace': 0.0,
          'lines.linewidth': 2.0,
          'text.usetex': True, # THIS IS THE KEY. You will need Latex, Ghostscript, etc. installed..
          }
rcParams.update(params) 
# a more lightweight option is to use something like:
# >>from matplotlib import rc
# >>rc('text', usetex=True)
# but you will need to have xpdf and an appropriate backend set 

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches

def example():
    
    # example sin function for plotting
    x = np.linspace(0, np.pi, 100)
    y = 0.9*np.sin(x)
    
    # set up a figure and axes instance to use for plotting
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    # set limits, and aspect of axis
    ax.set_xlim(0,3)
    ax.set_ylim(0,1)
    ax.set_aspect('equal') # make the x and y axis with the same scale
    
    ######################
    #    basic plot      #
    ax.plot(x,y, color='black')
    ######################
    
    ######################
    # example annotation #
    string = 'look, the top\nof the graph'
    string_position = (np.pi/2., 0.5)
    arrow_position = (np.pi/2,0.9)
    arrow = {'facecolor':'red', 'arrowstyle':"->", 'edgecolor':'red'}
    ax.annotate(string, arrow_position,
        xytext=string_position,
        arrowprops=arrow,
        horizontalalignment='center', verticalalignment='center', color='red')
    ######################
    
    ######################
    # example shape      #
    pts = np.vstack( ([.2,.5], [.5,.8], [2,.1]) )
    triangle = patches.Polygon( pts, facecolor='purple', edgecolor='none', zorder=-10, alpha=0.2 )
    circle = patches.Circle( (2.5,.2), 0.1, facecolor='green', edgecolor='green', alpha=1)
    ax.add_artist(triangle)
    ax.add_artist(circle)
    ######################
    
    #######################
    # example text        #
    string = 'hi, i am a\npurple triangle'
    ax.text(0.9, 0.2, string, horizontalalignment='center', verticalalignment='center', color='purple')
    #######################
    
    
    #######################
    # fix and label stuff #
    
    # note, this call to adjust_spines needs to happen before you label the axes or edit the ticklabels
    # the color and spine_locations dicts allow you to adjust the color and location of individual spines. if you leave them blank they default to black and 10. 
    adjust_spines(ax, ['left', 'bottom'], color={'left':'purple'}, spine_locations={'bottom':20}) 
    
    # x labels and ticks
    ax.set_xlabel('time, s')

    # y labels and ticks
    yticks = ax.get_yticks()
    # make custom ticks for y axis
    nticks = 3
    newyticks = np.linspace(yticks[0], yticks[-1], nticks, endpoint=True)
    ax.set_yticks(newyticks)
    ax.set_yticklabels([str(s) for s in newyticks], color='purple')
    ax.set_ylabel('sin(t)', color='purple')
    #######################
    
    
    #######################
    # save                #
    fig.savefig('example_figure.pdf', format='pdf') 
    # use pdf, as eps has issues with text output, and cannot handle opacity (alpha) settings. these pdfs are layered and perfectly editable in adobe illustrator.. from there you can output an eps file.
    #######################
    
    return fig
    

# function to fix axes spines so they are 'Dickinson Style'
def adjust_spines(ax,spines, color={}, spine_locations={}, smart_bounds=False):
    if type(spines) is not list:
        spines = [spines]
    spine_locations_dict = {'top': 10, 'right': 10, 'left': 10, 'bottom': 10}
    for key in spine_locations.keys():
        spine_locations_dict[key] = spine_locations[key]
        
    if 'none' in spines:
        for loc, spine in ax.spines.iteritems():
            spine.set_color('none') # don't draw spine
        ax.yaxis.set_ticks([])
        ax.xaxis.set_ticks([])
        return
    
    for loc, spine in ax.spines.iteritems():
        if loc in spines:
            spine.set_position(('outward',spine_locations_dict[loc])) # outward by x points
            
            if loc in color.keys():
                c = color[loc]
            else:
                c = 'black'
            
            spine.set_color(c)
            #spine.set_smart_bounds(True)
            if loc == 'bottom' and smart_bounds:
                spine.set_smart_bounds(True)
        else:
            spine.set_color('none') # don't draw spine

    # turn off ticks where there is no spine
    if 'left' in spines:
        ax.yaxis.set_ticks_position('left')
    elif 'right' in spines:
        ax.yaxis.set_ticks_position('right')
    else:
        # no yaxis ticks
        ax.yaxis.set_ticks([])

    if 'bottom' in spines:
        ax.xaxis.set_ticks_position('bottom')
    else:
        # no xaxis ticks
        ax.xaxis.set_ticks([])    
        
    for line in ax.get_xticklines() + ax.get_yticklines():
        #line.set_markersize(6)
        line.set_markeredgewidth(1)
        
        
# ---------------------------------------------------------------------------------
if __name__ == '__main__':
    example()




