import numpy as np
import matplotlib.pyplot as pyplot
import matplotlib.pylab as pylab
import matplotlib.mlab as mlab
import matplotlib.cm as colormap
from load_pickled_data import load
import analysis_plot
import matplotlib.colorbar
import flymodel_getfz

try:
    dataset_type = type(dataset)
except:
    print 'loading data!....'
    dataset = load('dataset_tmp4')


# generate flight envelope for flying flies
accels = []
for k,v in dataset.trajecs.items():
    if dataset.trajecs[k].behavior is 'flyby':
        for a in dataset.trajecs[k].accel_1d:
            accels.append( a )

pyplot.figure(10)
n, bins, patches = pyplot.hist(accels, 100, facecolor='green', normed = True, alpha=1, histtype = 'stepfilled')


# calc accel required to stop given velocity and distance to stop
def accel (x,v):
    a = v**2 / (2*x)
    return a
    
# find probability of that acceleration actually occuring in a flying fly
def find_alpha (arr, n, bins):
    
    alpha = np.zeros(arr.shape)
    
    for i in range(arr.shape[0]):
        for j in range(arr.shape[1]):
            index = np.argmin( np.abs(bins - arr[i,j]) )
            a = n[index]
            alpha[i,j] = a
    
    
    return alpha
    
mesh_x,mesh_y = pylab.meshgrid( np.linspace(-.01,-.0001,300), np.linspace(.0001,.1,300) )
c = accel(mesh_x, mesh_y)







alpha = find_alpha(c, n, bins)
alpha = alpha - alpha.max()/2
alpha = np.sin(np.pi/2 * (alpha / alpha.max()) )

cm = colormap.ScalarMappable(cmap='gray')
rgb = cm.to_rgba(-1*np.exp(c))
rgba = rgb
rgba[:,:,3] = alpha





################    plot   ###############

fig = pyplot.figure(1)

cl = analysis_plot.dist_vs_speed_color(dataset, figure=1, room_bottom = True)
ax2 = fig.add_axes([0.12,0.05,0.70,0.05])
cl.ax0.pcolor(mesh_x, mesh_y, alpha, cmap='gray')
cb = matplotlib.colorbar.ColorbarBase(ax2, cmap='gray', norm=pyplot.Normalize(0,1), orientation='horizontal')
ax2.set_xlabel('probability of being on this decceleration isoline in normal flight (n = 166)')
    
    
################   flat wing deceleration  #############3

dx = 0.001
vi = 1
v = [vi]
x = [0]
vf = 1
n_iterations = 0
while vf > 0.001:
    n_iterations = n_iterations + 1
    a = 2+ 2.652*vi 
    vf = -1*((a)*2*dx)**(0.5)+vi # max deceleration of fly if just holding wings out flat
    vi = vf
    if vf > 0.001:
        v.append(vf)
        x.append(x[-1]+dx)

xf = x[-1]
print x
x = np.array(x)
v = np.array(v)
x = x-xf
    
cl.ax0.plot(x,v)


















