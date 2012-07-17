import numpy as np
import matplotlib.pyplot as plt


def test_l_over_v():
    fig = plt.figure()
    ax = fig.add_subplot(111)

    vels = [2, 4, 6, 8, 10]
    for vel in vels:
        fps = 1000.
        l = 140
        
        dist = np.arange(1000, 0.00001, -1/fps*vel)
        #radius = 0.009565
        #angle = 2*np.arcsin(radius/(dist+radius))
        #expansion = sa1.diffa(angle)*fps
        #rrev_inv = angle / expansion
        l_over_v = l / vel
        tti = dist / vel
        
        ax.plot(l_over_v, tti, '--')
