"""
Utility functions for operations involving masks and mask searches
==================================================================
"""
import sys
import warnings
import numpy as np
from scipy.signal import convolve
import geocontour.check as gcc
import geocontour.grid as gcg

def bbox(latitudes,longitudes,boundary):
    """
    Check input dimensions (lat/lon) against input boundary and return
    min/max indicies of bounding box

    Mainly used internally for mask search functions

    Parameters
    ----------
    latitudes : ndarray
        1D Nx1 array of latitude points (degrees)
    longitudes : ndarray
        1D Nx1 array of longitude points (degrees)
    boundary : ndarray
        2D Nx2 array of latitude/longitude points (degrees) with the
        last point equal to the first

    Returns
    -------
    boxlatmin : int
        the minimum bounding box latitude index
    boxlatmax : int
        the maximum bounding box latitude index
    boxlonmin : int
        the minimum bounding box longitude index
    boxlonmax : int
        the maximum bounding box longitude index
    """
    gcc.cdim(latitudes)
    gcc.cdim(longitudes)
    gcc.cboundary(boundary)
    latdir=gcg.clatdir(latitudes)
    loninput=gcg.clonrng(longitudes)
    boundloninput=gcg.clonrng(boundary[:,1])
    if (loninput=='neg' and boundloninput=='pos') or (loninput=='pos' and boundloninput=='neg'):
        sys.exit('ERROR - Longitude input range is '+loninput+' and boundary longitude range is '+boundloninput)
    boundlatmin=boundary.min(axis=0)[0]
    boundlatmax=boundary.max(axis=0)[0]
    boundlonmin=boundary.min(axis=0)[1]
    boundlonmax=boundary.max(axis=0)[1]
    if latdir=='inc':
        if boundlatmin<=latitudes.min():
            warnings.warn('WARNING - boundary latitude minimum on or outside of minimum latitude')
            boxlatmin=0
        else:
            boxlatmin=((latitudes-boundlatmin)<0).nonzero()[0][-1]
        if boundlatmax>=latitudes.max():
            warnings.warn('WARNING - boundary latitude maximum on or outside of maximum latitude')
            boxlatmax=len(latitudes)-1
        else:
            boxlatmax=((latitudes-boundlatmax)>0).nonzero()[0][0]
    elif latdir=='dec':
        if boundlatmin<=latitudes.min():
            warnings.warn('WARNING - boundary latitude minimum on or outside of minimum latitude')
            boxlatmin=len(latitudes)-1
        else:
            boxlatmin=((latitudes-boundlatmin)<0).nonzero()[0][0]
        if boundlatmax>=latitudes.max():
            warnings.warn('WARNING - boundary latitude maximum on or outside of maximum latitude')
            boxlatmax=0
        else:
            boxlatmax=((latitudes-boundlatmax)>0).nonzero()[0][-1]
    if boundlonmin<=longitudes.min():
        warnings.warn('WARNING - boundary longitude minimum on or outside of minimum longitude')
        boxlonmin=0
    else:
        boxlonmin=((longitudes-boundlonmin)<0).nonzero()[0][-1]
    if boundlonmax>=longitudes.max():
        warnings.warn('WARNING - boundary longitude maximum on or outside of maximum longitude')
        boxlonmax=len(longitudes)-1
    else:
        boxlonmax=((longitudes-boundlonmax)>0).nonzero()[0][0]
    return boxlatmin, boxlatmax, boxlonmin, boxlonmax

def edge(mask,latitudes=None,longitudes=None,connectivity=8):
    """
    Find a mask's edge cells only, and optionally an array of the edge
    cells 
    
    Parameters
    ----------
    mask : ndarray
        2D MxN bool array where M=len(`latitudes`) and
        N=len(`longitudes`) 
    latitudes : ndarray, optional
        1D Nx1 array of latitude points (degrees)
    longitudes : ndarray, optional
        1D Nx1 array of longitude points (degrees)
    connectivity : {4, 8}, default=8
        connectivity parameter for finding edge cells

            ``4``
                test only lateral neighbors
            ``8``
                test lateral and diagonal neighbors

    Returns
    -------
    edgemask : ndarray
        2D bool array of the same dimensions as mask input
    edgecells : ndarray, optional
        2D Nx2 bool array of latitude/longitude points (degrees) of edge
        cells (unordered), where N is number of edge cells 

    See Also
    --------
    vertex
    conn

    Notes
    -----
    `edgecells` only returned if optional parameters `latitudes` and
    `longitudes` are provided
    """
    if latitudes is not None and longitudes is not None:
        gcc.cmask(mask,latitudes,longitudes)
    else:
        gcc.cmask(mask)
    kernel=np.ones((3,3),dtype='int')
    kernel[1,1]=0
    if connectivity==8:
        pass
    elif connectivity==4:
        kernel[0,0]=0
        kernel[2,0]=0
        kernel[0,2]=0
        kernel[2,2]=0
    else:
        sys.exit('ERROR - connectivity='+str(connectivity)+' is not a valid selection, valid selections are 4/8')
    neighborcount=convolve(mask.astype('int'),kernel,mode='same')
    edgemask=(neighborcount<connectivity)*mask
    if latitudes is not None and longitudes is not None:
        edgecells=np.vstack((latitudes[edgemask.nonzero()[0]],longitudes[edgemask.nonzero()[1]])).T
        return edgemask, edgecells
    else:
        return edgemask

def vertex(mask,latitudes,longitudes):
    """
    Find the vertex points of all cells in a mask, and the vertex points
    of only the mask edge

    Parameters
    ----------
    mask : ndarray
        2D MxN bool array where M=len(`latitudes`) and
        N=len(`longitudes`) 
    latitudes : ndarray, optional
        1D Nx1 array of latitude points (degrees)
    longitudes : ndarray, optional
        1D Nx1 array of longitude points (degrees)

    Returns
    -------
    vertexpoints : ndarray 
        2D Nx2 array of latitude/longitude points (degrees) of all
        vertices of mask cells
    edgevertexpoints : ndarray
        2D Nx2 array of latitude/longitude points (degrees) of all
        vertices of cells at the mask edge (8-connected)

    See Also
    --------
    edge
    conn
    """
    gcc.cmask(mask,latitudes,longitudes)
    vertexmask=np.full((tuple(np.array(mask.shape)+1)),0)
    maskint=mask.astype('int')
    vertexmask[:-1,:-1]+=maskint
    vertexmask[:-1,1:]+=maskint
    vertexmask[1:,:-1]+=maskint
    vertexmask[1:,1:]+=maskint
    latspc=gcg.spacing(latitudes)
    lonspc=gcg.spacing(longitudes)
    latdir=gcg.clatdir(latitudes)
    if latdir=='inc':
        vertexlatitudes=np.append(latitudes-latspc/2,latitudes[-1]+latspc/2)
    elif latdir=='dec':
        vertexlatitudes=np.append(latitudes+latspc/2,latitudes[-1]-latspc/2)
    vertexlongitudes=np.append(longitudes-lonspc/2,longitudes[-1]+lonspc/2)
    vertexpoints=np.vstack((vertexlatitudes[vertexmask.nonzero()[0]],vertexlongitudes[vertexmask.nonzero()[1]])).T
    edgevertexmask=((vertexmask<4)*(vertexmask>0))
    edgevertexpoints=np.vstack((vertexlatitudes[edgevertexmask.nonzero()[0]],vertexlongitudes[edgevertexmask.nonzero()[1]])).T
    return vertexpoints, edgevertexpoints

def neighbors(cell,connectivity=8,direction='cw'):
    """
    Find the neighbors of a cell, with selected connectivity and
    direction

    Parameters
    ----------
    cell : ndarray 
        1x2 numpy array describing the indices of the cell
    connectivity : {4, 8}, default=8
        connectivity parameter for finding neighbor cells

            ``4``
                test only lateral neighbors
            ``8``
                test lateral and diagonal neighbors
    direction : {'cw', 'ccw'}, default='cw'
        select the direction of the returned neighbors

    Returns
    -------
    neighbors : ndarray
        8x2 or 4x2 array of the neighboring cell indices for the input
        cell 

    See Also
    --------
    conn
    """
    if connectivity==8:
        neighborrange=np.array([[-1,0],[-1,1],[0,1],[1,1],[1,0],[1,-1],[0,-1],[-1,-1]])
    elif connectivity==4:
        neighborrange=np.array([[-1,0],[0,1],[1,0],[0,-1]])
    else:
        sys.exit('ERROR - connectivity='+str(connectivity)+' is not a valid selection, valid selections are 4/8')
    if direction=='cw':
        pass
    elif direction=='ccw':
        neighborrange=neighborrange[::-1,:]
    else:
        sys.exit('ERROR - direction=\''+direction+'\' is not a valid selection, valid selections are \'cw\'/\'ccw\'')
    neighbors=cell+neighborrange
    return neighbors

def conn(mask,checkcells='full',connectivity=8):
    """
    Determine whether a mask or its inverse are connected

    Parameters
    ----------
    mask : ndarray
        2D MxN bool array where M=len(`latitudes`) and
        N=len(`longitudes`) 
    checkcells : {'full', 'empty'}, default='full'
        select the mask cells to test ('empty' would select the mask
        inverse)
    connectivity : {4, 8}, default=8
        connectivity parameter for testing connectivity

            ``4``
                test only lateral neighbors
            ``8``
                test lateral and diagonal neighbors
    
    Returns
    -------
    connected : bool
        descriptor for whether the input `mask` is connected under the input conditions
    """
    gcc.cmask(mask)
    if checkcells=='full':
        maskcells=np.argwhere(mask==True)
    elif checkcells=='empty':
        maskcells=np.argwhere(mask==False)
    else:
        sys.exit('ERROR checkcells=\''+checkcells+'\' is not a valid selection, valid selections are \'full\'/\'empty\'')
    startcell=maskcells[0]
    checked=[]
    tocheck=[]
    tocheck.append(startcell)
    while len(tocheck)>0:
        checked.append(tocheck[0])
        cellneighbors=neighbors(tocheck.pop(0),connectivity=connectivity)
        for k in cellneighbors:
            if (k==maskcells).all(axis=1).any() and not any(np.array_equal(k,x) for x in checked) and not any(np.array_equal(k,x) for x in tocheck):
                tocheck.append(k)
    if len(checked)==len(maskcells):
        connected=True
    else:
        connected=False
    return connected

