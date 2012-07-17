from play_array_list import *

def save_imgs_landing(movie_dataset):

    movie = movie_dataset.movies['20101101_C001H001S0019']

    im = movie.frames[150].uimg
    plt.imsave('flight', im, cmap=plt.get_cmap('gray'))
    
    legextensionframe = movie.legextensionrange[0] - movie.firstframe_ofinterest
    im = movie.frames[legextensionframe].uimg
    plt.imsave('leg extension', im, cmap=plt.get_cmap('gray'))
    
    im = movie.frames[movie.landingframe_relative].uimg
    plt.imsave('landing', im, cmap=plt.get_cmap('gray'))
    
    im = movie.frames[-1].uimg
    plt.imsave('sitting', im, cmap=plt.get_cmap('gray'))
    
def save_imgs_crash(movie_dataset):

    movie = movie_dataset.movies['20101111_C001H001S0032']

    im = movie.frames[580].uimg
    plt.imsave('crash_flight', im, cmap=plt.get_cmap('gray'))
    
    im = movie.frames[780].uimg
    plt.imsave('crash_collision', im, cmap=plt.get_cmap('gray'))
    
    im = movie.frames[900].uimg
    plt.imsave('crash_legextension', im, cmap=plt.get_cmap('gray'))
    
    im = movie.frames[-1].uimg
    plt.imsave('crash_sitting', im, cmap=plt.get_cmap('gray'))
    
def save_imgs_wingcrash(movie_dataset):
    
    movie = movie_dataset.movies['20101111_C001H001S0040']

    im = movie.frames[-148+646].uimg
    plt.imsave('wingcrash_flight', im, cmap=plt.get_cmap('gray'))
    
    im = movie.frames[-89+646].uimg
    plt.imsave('wingcrash_collision', im, cmap=plt.get_cmap('gray'))
    
    im = movie.frames[55+646].uimg
    plt.imsave('wingcrash_legext', im, cmap=plt.get_cmap('gray'))
    
    im = movie.frames[-1].uimg
    plt.imsave('sitting', im, cmap=plt.get_cmap('gray'))
