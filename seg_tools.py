import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from scipy import signal


def init():
    plt.close('all')


def draw_uint8(I):
    fig = plt.figure();
    im = plt.imshow(I, cmap = "gray", vmin = -0.5, vmax = 255.5)
    plt.colorbar(im);
    plt.show()
    
    
def draw_float(I):
    fig = plt.figure();
    im = plt.imshow(I, cmap = "gray", vmin = np.min(I), vmax = np.max(I))
    plt.colorbar(im);
    plt.show()


def draw_labels(L):
    nb = np.max(L)
    colors = np.random.random_sample((nb+1, 3))
    colors = np.concatenate((colors, np.ones((nb+1, 1))), axis=1)
    my_map = ListedColormap(colors)
    fig = plt.figure()
    im = plt.imshow(L, cmap = my_map, vmin = -0.5, vmax = nb+0.5)
    plt.colorbar(im) 
    plt.show()
    
    
def gauss1d_smooth(sigma):
    r = np.ceil(3*sigma)
    v = np.arange(-r, r+1)/sigma
    G = np.exp(-v*v/2)
    G = G/np.sum(G)
    G = np.matrix(G)
    return G


def gauss1d_deriv1(sigma):
    r = np.ceil(3*sigma)
    v = np.arange(-r, r+1)/sigma
    G = (-v)*np.exp(-v*v/2)
    x = np.arange(-r, r+1)
    G = G/np.sum(-x*G)
    G = np.matrix(G)
    return G

    
def norm_gradient(I, sigma):
    G1 = gauss1d_smooth(sigma)
    G2 = gauss1d_deriv1(sigma)
    Ix = signal.convolve2d(signal.convolve2d(I, G2, mode='same'), G1.T, mode='same') 
    Iy = signal.convolve2d(signal.convolve2d(I, G1, mode='same'), G2.T, mode='same')
    return np.sqrt(Ix*Ix+Iy*Iy)