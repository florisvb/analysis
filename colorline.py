import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import matplotlib.colorbar



class Colorline ():

    def __init__ (self, figure = None, colormap = 'copper', ylim = (-1,1), xlim = (-1,1), norm = None, hide_colorbar=False, ax0=None, ax0_size=[-.15,0.1,0.7,0.7]):
        
        # "figure" should be a matplotlib figure
        self.hide_colorbar = hide_colorbar
        
        # initialize figure and axes
        if figure is None:
            self.fig = plt.figure()
        else:
            self.fig = plt.figure(figure)
        if ax0 is None:
            self.ax0 = self.fig.add_axes(ax0_size)
        else:
            self.ax0 = self.fig.add_axes(ax0)
        
        if self.hide_colorbar is False:
            self.ax1 = self.fig.add_axes([0.85,0.1,0.05,0.7])

        
        self.ax0.set_ylim(ylim)
        self.ax0.set_xlim(xlim)
        
        
        # plot paramters
        self.norm = norm
        if norm is not None:
            self.norm = plt.Normalize(norm[0], norm[1])
        self.cmap = plt.get_cmap(colormap)
        self.cb = None
        
          


    def colorline(self, x,y,z,linewidth=3, colormap = None, norm=None ):
        
        if colormap is None:
            cmap = self.cmap
        else:
            cmap = plt.get_cmap(colormap)
        
        if norm is None:
            self.norm = plt.Normalize(z.min(), z.max())

        if self.hide_colorbar is False:
            if self.cb is None:
                self.cb = matplotlib.colorbar.ColorbarBase(self.ax1, cmap=cmap, norm=self.norm, orientation='vertical', boundaries=None)
            
            
        
        # Create a set of line segments so that we can color them individually
        # This creates the points as a N x 1 x 2 array so that we can stack points
        # together easily to get the segments. The segments array for line collection
        # needs to be numlines x points per line x 2 (x and y)
        
        
        points = np.array([x, y]).T.reshape(-1, 1, 2)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)
        
        # Create the line collection object, setting the colormapping parameters.
        # Have to set the actual values used for colormapping separately.
        lc = LineCollection(segments, cmap=cmap,norm=self.norm )
        lc.set_array(z)
        lc.set_linewidth(linewidth)
        
        self.ax0.add_collection(lc)
        
        
        
        
        
        
        
    def show(self):
        plt.show()

        

# ---------------------------------------------------------------------------------
if __name__ == '__main__':

    def tent(x):
        """
        A simple tent map
        """
        if x < 0.5:
            return x
        else:
            return -1.0*x + 1
    
    pi = np.pi
    
    t = np.linspace(0, 1, 200)
    y = np.sin(2*pi*t)
    z = np.array([tent(x) for x in t]) 
    
    cl = Colorline()
    cl.colorline(t,y,z)

