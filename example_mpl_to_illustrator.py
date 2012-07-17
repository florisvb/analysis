# this is an example file that shows how to generate a plot in matplotlib, and save the file to a layered scalable vector graphics pdf that can be read and edited in Adobe Illustrator. usage: either run the file as python code, or import this file and run function example()

# Known issues:
# so far, I have not found a way to actually make the font carry through properly.. if you find a way, let me know! in the end, this is not a huge deal, as the editable text can be set to the appropriate font in illustrator rather easily.

# written by Floris van Breugel, 6/2011

# set defaults
# note: this needs to be performed BEFORE importing pyplot or pylab. If you reload... reload pyplot/pylab after as well. 
# all of these paramters, except 'backend', 'ps.usedistiller', 'text.usetex', can be adjusted elsewhere.. changing them here just makes it convenient if you are doing many figures with the same settings. 
from matplotlib import rcParams
params = {'backend': 'Agg', # must be Agg, PDF, or PS. Only Agg allows for gui-viewing plots in ipython
          'ps.usedistiller': 'xpdf', # may need to apt-get install xpdf
          'text.usetex': True, # THIS IS THE KEY. You will need Latex, Ghostscript, etc. installed..
          }
rcParams.update(params) 

import numpy as np
import matplotlib.pyplot as plt

def example():
    
    # example sin function for plotting
    x = np.linspace(0, np.pi, 100)
    y = 0.9*np.sin(x)
    
    # set up a figure and axes instance to use for plotting
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    ax.plot(x,y, color='black')
    ax.set_xlabel('time, s')
    ax.set_ylabel('sin(t)')
    
    fig.savefig('example_figure.pdf', format='pdf') 
    # use pdf, as eps has issues with text output, and cannot handle opacity (alpha) settings. these pdfs are layered and perfectly editable in adobe illustrator.. from there you can output an eps file.
    
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




