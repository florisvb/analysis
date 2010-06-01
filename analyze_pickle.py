import flydra_floris_analysis as analysis
import numpy as np
import matplotlib.pyplot as pyplot
import scipy.linalg
import time
import pickle

center = np.array([0,0,0])
flyzone = np.array([0.05, 0.05, 0.05])
radius = 0.031/2.0
flight_buffer = 0.008


dataset = analysis.Dataset()
post = analysis.Post(center, flyzone, radius, flight_buffer)
dataset.set_stimulus(post)
dataset.load_data(filename ='/home/floris/rekalmanized_2_264.h5' )




dataset.load_data(filename = '/home/floris/rekalmanized_1_26630.h5'  )
dataset.load_data(filename =  '/home/floris/rekalmanized_1_17818.h5' )
dataset.load_data(filename =  '/home/floris/rekalmanized_1_709.h5' )
dataset.load_data(filename = '/home/floris/rekalmanized_1_22988.h5'  )
dataset.load_data(filename = '/home/floris/rekalmanized_3_412.h5'  )
dataset.load_data(filename = '/home/floris/rekalmanized_3_17454.h5'  )
dataset.load_data(filename =  '/home/floris/rekalmanized_3_15227.h5' )
dataset.load_data(filename = '/home/floris/rekalmanized_3_15860.h5'  )
dataset.load_data(filename = '/home/floris/rekalmanized_3_3671.h5'  )

dataset.load_data(filename = '/home/floris/DATA20100319_173739.h5'  )

dataset.load_data(filename = '/home/floris/DATA20100319_173739.h5'  )
dataset.load_data(filename = '/home/floris/DATA20100328_181250.h5'  )
dataset.load_data(filename = '/home/floris/DATA20100329_170127.h5' )

if 0:
   

    # 60 fps
    # crap > dataset.load_data(filename  = '/home/floris/DATA20100330_181406.h5' )
    dataset.load_data(filename  =  '/home/floris/DATA20100331_181103.h5')
    dataset.load_data(filename  =  '/home/floris/DATA20100401_182548.h5')
    dataset.load_data(filename  = '/home/floris/DATA20100402_190041.h5' )
    dataset.load_data(filename  =  '/home/floris/DATA20100403_114457.h5')
    dataset.load_data(filename  =  '/home/floris/DATA20100404_120907.h5')
    dataset.load_data(filename  =  '/home/floris/DATA20100404_162638.h5' )











