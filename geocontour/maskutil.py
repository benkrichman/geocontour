import sys
import warnings
import numpy as np
from scipy.signal import convolve
import geocontour.check as gcc
import geocontour.grid as gcg

def boxset(latitudes,longitudes,boundary):
    """
    Checks input dimensions (lat/lon) against input boundary and returns min/max indicies of bounding box

    Inputs (Required):
        latitudes - An evenly spaced numpy array of latitude points (degrees)
        longitudes - An evenly spaced numpy array of longitude points (degrees)
        boundary - A 2-d Nx2 numpy array of latitude/longitude points (degrees)

    Outputs:
        boxlatmin - The minimum bounding box latitude index
        boxlatmax - The maximum bounding box latitude index
        boxlonmin - The minimum bounding box longitude index
        boxlonmax - The maximum bounding box longitude index
    """
    gcc.cdim(latitudes)
    gcc.cdim(longitudes)
    gcc.cboundary(boundary)
    latdir=gcg.checklatitudedirection(latitudes)
    loninput=gcg.checklongituderange(longitudes)
    boundloninput=gcg.checklongituderange(boundary[:,1])
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

def maskedgecells(mask,latitudes=None,longitudes=None,connectivity=8):
    """
    Returns a mask of only the edge cells, and if latitudes and longitudes are provided also returns an array of the edge cells 
    
    Inputs (Required):
        mask - A 2-d boolean numpy array of dimension MxN where M=len(latitudes) and N=len(longitudes)

    Inputs (Optional):
        latitudes - An evenly spaced numpy array of latitude points (degrees)
        longitudes - An evenly spaced numpy array of longitude points (degrees)
        connectivity (4/8) - The connectivity parameter for testing edge cells, default=8
            4 tests only lateral neighbors while 8 tests lateral and diagonal neighbors

    Outputs:
        edgemask - A 2-d boolean numpy array of the same dimensions as mask input
        edgecells - A 2-d Nx2 numpy array of latitude/longitude points (degrees) of edge cells (unordered), where N is number of edge cells 
            Only output if optional latitude and longitude arguments are provided
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

def maskvertexpoints(mask,latitudes,longitudes):
    """
    Returns the vertex points of all cells in the input mask, and the vertex points of only the mask edge

    Inputs (Required):
        mask - A 2-d boolean numpy array of dimension MxN where M=len(latitudes) and N=len(longitudes)
        latitudes - An evenly spaced numpy array of latitude points (degrees)
        longitudes - An evenly spaced numpy array of longitude points (degrees)

    Outputs:
        vertexpoints - A 2-d Nx2 numpy array of latitude/longitude points (degrees) of all vertices of mask cells
        edgevertexpoints - A 2-d Nx2 numpy array of latitude/longitude points (degrees) of all vertices of cells at the mask edge (8-connected)
    """
    gcc.cmask(mask,latitudes,longitudes)
    vertexmask=np.full((tuple(np.array(mask.shape)+1)),0)
    maskint=mask.astype('int')
    vertexmask[:-1,:-1]+=maskint
    vertexmask[:-1,1:]+=maskint
    vertexmask[1:,:-1]+=maskint
    vertexmask[1:,1:]+=maskint
    latspc=gcg.gridspacing(latitudes)
    lonspc=gcg.gridspacing(longitudes)
    latdir=gcg.checklatitudedirection(latitudes)
    if latdir=='inc':
        vertexlatitudes=np.append(latitudes-latspc/2,latitudes[-1]+latspc/2)
    elif latdir=='dec':
        vertexlatitudes=np.append(latitudes+latspc/2,latitudes[-1]-latspc/2)
    vertexlongitudes=np.append(longitudes-lonspc/2,longitudes[-1]+lonspc/2)
    vertexpoints=np.vstack((vertexlatitudes[vertexmask.nonzero()[0]],vertexlongitudes[vertexmask.nonzero()[1]])).T
    edgevertexmask=((vertexmask<4)*(vertexmask>0))
    edgevertexpoints=np.vstack((vertexlatitudes[edgevertexmask.nonzero()[0]],vertexlongitudes[edgevertexmask.nonzero()[1]])).T
    return vertexpoints, edgevertexpoints

def findneighbors(cell,connectivity=8,direction='cw'):
    """
    Returns the neighbors of a cell, with selected connectivity and direction

    Inputs (Required):
        cell - A 1x2 numpy array describing the index of the cell

    Inputs (Optional):
        connectivity (4/8) - The connectivity parameter for testing edge cells, default=8
            4 tests only lateral neighbors while 8 tests lateral and diagonal neighbors
        direction ('cw'/'ccw') - A string selecting the direction of the returned neighbors, default='cw'

    Outputs:
        neighbors - An 8x2 or 4x2 numpy array of the neighboring cell indices for the input cell 
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

def checkconnectivity(mask,checkcells='full',connectivity=8):
    """
    Returns whether a mask or its inverse are connected

    Inputs (Required):
        mask - A 2-d boolean numpy array of dimension MxN where M=len(latitudes) and N=len(longitudes)

    Inputs (Optional):
        checkcells ('full'/'empty') - A string selecting the mask cells to test (True/False), default='full'
        connectivity (4/8) - The connectivity parameter for testing edge cells, default=8
            4 tests only lateral neighbors while 8 tests lateral and diagonal neighbors
    
    Outputs:
        connected (True/False) - A boolean describing whether the input mask is connected under the input conditions
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
        cellneighbors=findneighbors(tocheck.pop(0),connectivity=connectivity)
        for k in cellneighbors:
            if (k==maskcells).all(axis=1).any() and not any(np.array_equal(k,x) for x in checked) and not any(np.array_equal(k,x) for x in tocheck):
                tocheck.append(k)
    if len(checked)==len(maskcells):
        connected=True
    else:
        connected=False
    return connected

