import sys
import numpy as np

sys.path.append("/usr/share/pyshared/flydra/a2")
sys.path.append("/home/floris/src/analysis")

import flydra_floris_analysis as ffa


def match_obj_to_kalmanized_obj(h5source, original_obj_id, kalmanizedsource=None):

    original_trajectory = ffa.load_trajectory(source=h5source, obj=original_obj_id)

    if kalmanizedsource == None:
        kalmanizedsource = h5source
        kalmanizedsource = kalmanizedsource.rstrip('.h5')
        kalmanizedsource = kalmanizedsource + '.kalmanized.h5'
        
    if kalmanizedsource.__class__ is not ffa.Dataset:
        dataset = ffa.Dataset()
        dataset.load_data(kalmanizedsource,kalman_smoothing = True)
    else:
        dataset = kalmanizedsource
        
    best_match = find_match(original_trajectory, dataset)
    return best_match


def find_match(original_trajectory, dataset, volume_sloppiness = 0.01):
    original_time = original_trajectory.epoch_time[int(original_trajectory.length/2.)]
    original_position = original_trajectory.positions[int(original_trajectory.length/2.),:]

    best_time_diff = 1000
    best_match = None
    best_volume_error = None
    for k, trajectory in dataset.trajecs.items():
        for i in range(trajectory.length):
            time_diff = np.abs( original_time - trajectory.epoch_time[i] )
            if time_diff < best_time_diff:
                # check to make sure the object is in the same area, roughly:
                volume_error = np.linalg.norm(original_position - trajectory.positions[i,:])
                if volume_error < volume_sloppiness:            
                    best_time_diff = time_diff
                    best_volume_error = volume_error
                    best_match = k
                #else:
                    #print 'good match found: ', k, ' but volume error exceeds ', volume_sloppiness
    print best_match, best_time_diff, best_volume_error
    return best_match
    
#########################################################################################


if __name__ == '__main__':
        
    h5source = '/home/floris/data/windtunnel/SA1/checkerboard/DATA20101024_154426.h5'
    original_obj_id = 16967
    kalmanizedsource = None
    
    kalmanized_id = match_obj_to_kalmanized_obj(h5source, original_obj_id)
