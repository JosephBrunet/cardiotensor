"""supporting function ported over for the structure-tensor library"""

import numpy as np
from scipy.ndimage import filters
import matplotlib.pyplot as plt
#test for GPU
flag_GPU = 1
try:
    import cupy as cp
except:
    flag_GPU = 0

#import 3D structure tensor and utils (helper functions)
if not(flag_GPU):
    from structure_tensor import eig_special_3d, structure_tensor_3d    #CPU version
else:
    from structure_tensor.cp import eig_special_3d, structure_tensor_3d #GPU version  


def tensor_vector_distance(S, u):
    """ Caclulating pairwise distance between tensors and vectors
    Arguments:
        S: an array with shape (6,N) containing tensor
        v: an array with shape (M,3) containing vectors
    Returns:
        v: an array with shape (N,M) containing pairwise distances
    Author: vand@dtu.dk, 2019
    """
    if flag_GPU:
        if type(S) == cp.core.core.ndarray:
            S = S.get()
        if type(u) == cp.core.core.ndarray:
            u = u.get()        
    dist = np.dot(np.moveaxis(S[0:3], 0, -1), u**2) + 2*np.dot(np.moveaxis(S[3:], 0, -1), u[[0,0,1]]*u[[1,2,2]])
    return dist




def arrow_navigation(event,z,Z):
    if event.key == "up":
        z = min(z+1,Z-1)
    elif event.key == 'down':
        z = max(z-1,0)
    elif event.key == 'right':
        z = min(z+10,Z-1)
    elif event.key == 'left':
        z = max(z-10,0)
    elif event.key == 'pagedown':
        z = min(z+50,Z+1)
    elif event.key == 'pageup':
        z = max(z-50,0)
    return z




def show_vol(V,cmap='gray'): 
    """
    Shows volumetric data and colored orientation for interactive inspection.
    @author: vand at dtu dot dk
    """
    if flag_GPU:
        if type(V) == cp.core.core.ndarray:
            V = V.get()
    def update_drawing():
        ax.images[0].set_array(V[z])
        ax.set_title(f'slice z={z}')
        fig.canvas.draw()
    def key_press(event):
        nonlocal z
        z = arrow_navigation(event,z,Z)
        update_drawing()
    Z = V.shape[0]
    z = (Z-1)//2
    fig, ax = plt.subplots()
    vmin = np.min(V)
    vmax = np.max(V)
    ax.imshow(V[z], cmap=cmap, vmin=vmin, vmax=vmax)
    ax.set_title(f'slice z={z}')
    fig.canvas.mpl_connect('key_press_event', key_press)



def show_vol_flow(V, fxy, s=5, double_arrow = False): 
    """
    Shows volumetric data and xy optical flow for interactive inspection.
    Arguments:
         V: volume
         fxy: flow in x and y direction
         s: spacing of quiver arrows
    @author: vand at dtu dot dk
    """
    if flag_GPU:
        print('in')
        if type(V) == cp.core.core.ndarray:
            V = V.get()
        if type(fxy) == cp.core.core.ndarray:
            fxy = fxy.get()
    def update_drawing():
        ax.images[0].set_array(V[:,:,z])
        ax.collections[0].U = fxy[0,s//2::s,s//2::s,z].ravel()
        ax.collections[0].V = fxy[1,s//2::s,s//2::s,z].ravel()
        if double_arrow:
            ax.collections[1].U = -fxy[0,s//2::s,s//2::s,z].ravel()
            ax.collections[1].V = -fxy[1,s//2::s,s//2::s,z].ravel()
        ax.set_title(f'slice z={z}')
        fig.canvas.draw()
    def key_press(event):
        nonlocal z
        z = arrow_navigation(event,z,Z)
        update_drawing()
    Z = V.shape[2]
    z = (Z-1)//2
    xmesh, ymesh = np.meshgrid(np.arange(V.shape[0]), np.arange(V.shape[1]), indexing='ij')
            # TODO: figure out exactly why this ij later needs 'xy' 
    fig, ax = plt.subplots()
    ax.imshow(V[:,:,z],cmap='gray')
    ax.quiver(ymesh[s//2::s,s//2::s], xmesh[s//2::s,s//2::s],
              fxy[0,s//2::s,s//2::s,z], fxy[1,s//2::s,s//2::s,z],
              color='r', angles='xy')
    if double_arrow:
        ax.quiver(ymesh[s//2::s,s//2::s], xmesh[s//2::s,s//2::s],
          -fxy[0,s//2::s,s//2::s,z], -fxy[1,s//2::s,s//2::s,z],
          color='r', angles='xy')
    ax.set_title(f'slice z={z}')
    fig.canvas.mpl_connect('key_press_event', key_press)
    
    
    
    
def fan_coloring(vec):
    """
    Fan-based colors for orientations in xy plane
    Arguments:
        vec: an array with shape (3,N) containing orientations
    Returns:
        rgba: an array with shape (4,N) containing rgba colors
     @author:vand@dtu.dk
    """
    if flag_GPU:
        if type(vec) == cp.core.core.ndarray:
            vec = vec.get()
    h = (vec[2]**2)[:,:,:,np.newaxis] # no discontinuity and less gray
    s = np.mod(np.arctan(vec[0]/vec[1]),np.pi) # hue angle
    hue = plt.cm.hsv(s/np.pi)
    rgba = hue*(1-h) + 0.5*h
    rgba[:,3] = 1 # fixing alpha value
    return rgba
    



def show_vol_orientation(V, vec, 
                         coloring = lambda v : np.c_[abs(v).T,np.ones((v.shape[1],1))], 
                         blending = lambda g,c : 0.5*(g+c)): 
    """
    Shows volumetric data and colored orientation for interactive inspection.
    @author: vand at dtu dot dk
    """
    print(V.shape)
    if flag_GPU:
        if type(V) == cp.core.core.ndarray:
            V = V.get()
        if type(vec) == cp.core.core.ndarray:
            vec = vec.get()
    rgba = coloring(vec)
    def update_drawing():
        ax.images[0].set_array(blending(plt.cm.gray(V[:,:,z]), rgba[:,:,z]))
        ax.set_title(f'slice z={z}')
        fig.canvas.draw()
    def key_press(event):
        nonlocal z
        z = arrow_navigation(event,z,Z)
        update_drawing()
    Z = V.shape[2]
    z = (Z-1)//2
    fig, ax = plt.subplots()
    ax.imshow(blending(plt.cm.gray(V[:,:,z]), rgba[:,:,z]))
    ax.set_title(f'slice z={z}')
    fig.canvas.mpl_connect('key_press_event', key_press)


def cart2sph(x,y,z):
    azimuth = np.arctan2(y,x)
    elevation = np.arctan2(z,np.sqrt(x**2 + y**2))
    #r = np.sqrt(x**2 + y**2 + z**2)
    return azimuth, elevation

def histogramSphere(eigVec, nBin):
    
    # Convert eigenvectors from xyz to directions on sphere (azimuth and elevation)
    sphDir = np.empty([2,eigVec.shape[1]], dtype='float')
    sphDir[0,:], sphDir[1,:] = cart2sph(eigVec[0,:],eigVec[1,:],eigVec[2,:])
    
    # Define uv-histogram (edges)
    cAz  = np.linspace(-np.pi,np.pi,nBin[0]+1)
    cEle = np.linspace(-np.pi/2,np.pi/2,nBin[1]+1)
    
    # Define bin center:
    binC_az = (cAz[:-1] + cAz[1:]) / 2
    binC_ele = (cEle[:-1] + cEle[1:]) / 2
    
    # Area of bins (on the sphere):
    binArea = np.outer((cAz[:-1] - cAz[1:]), \
                  np.sign(np.cos(cEle[:-1])) * np.sin(cEle[:-1]) - \
                  np.sign(np.cos(cEle[1:])) * np.sin(cEle[1:]) )
    
    # Count stats: 
    binCount = np.histogram2d(sphDir[0,:], sphDir[1,:], [cAz, cEle], density=None)[0]
    
    # Normalization:
    binVal = np.empty(binCount.shape,dtype=float)
    totalCount = np.sum(binCount)
    binIdx = np.logical_and(binCount > 1, binArea > 0.0001/np.prod(nBin))
    
    # only 'pdf' for now:
    binVal[binIdx] = binCount[binIdx] / (totalCount * binArea[binIdx]) # area weighting
    
    return binVal, binC_az, binC_ele





