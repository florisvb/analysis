import matplotlib
import matplotlib.pyplot as plt
import numpy as np




class Colorgrid:
    
    def __init__(self, X, figure = None, colormap = 'jet', ylim = (-1,1), xlim = (-1,1), norm = None, origin='lower', interpolation='nearest', color_label='colorscale'):
        
        # "figure" should be a matplotlib figure
        
        # initialize figure and axes
        if figure is None:
            self.fig = plt.figure()
        else:
            self.fig = plt.figure(figure)
            
        self.ax0 = self.fig.add_axes([0.1,0.1, 0.75, 0.75])
        self.ax1 = self.fig.add_axes([0.85,0.1,0.05,0.75])
        
        self.cmap = plt.get_cmap(colormap)
        
        # make the heatmap
        im = self.ax0.imshow(   X, 
                                cmap=self.cmap, 
                                extent=(xlim[0], xlim[1], ylim[0], ylim[1]), 
                                origin=origin, 
                                interpolation=interpolation)
        #im.set_interpolation(interpolation)
        
        # make the colorbar
        norm = plt.Normalize(X.min(), X.max())
        cb = matplotlib.colorbar.ColorbarBase(self.ax1, cmap=self.cmap, norm=norm, orientation='vertical', boundaries=None)
        self.ax1.set_ylabel(color_label)
