"""
Functions for construction of geocontours from contours
=======================================================

Notes
-----
Is it nonsensical that a package named geocontour contains a module
named geocontour which is used to build something referred to as a
geocontour? Yes. Is the developer sorry they chose to do it this way?
Also yes.
"""
import sys
import numpy as np
import geocontour.check as gcc
import geocontour.grid as gcg

def build(contour,latitudes,longitudes,connecttype='cell',simplify=False):
    """
    Construct a geocontour from a contour input

    Parameters
    ----------
    contour : ndarray
        2D Nx2 array of ordered latitude/longitude points (degrees)
        describing the edge of a mask
    latitudes : ndarray, optional
        1D Nx1 array of latitude points (degrees)
    longitudes : ndarray, optional
        1D Nx1 array of longitude points (degrees)
    connecttype : {'cell', 'center'}, default='cell'
        method of linking contour cells
        
            ``cell``
                implies drawing a single connection through a cell from
                the preceding to the following cell, resulting in a
                geocontour the same length as the input contour
            ``center``
                implies drawing two lines through a cell, intersecting
                the center of the cell, resulting in a geocontour double
                the length of the input contour
    simplify : bool, default=False
        select whether to remove cells for which outward vectors sum to
        0, and merge cells containing multiple contour segments

    Returns
    -------
    geocontour : ndarray
        3D Nx2x5 array defining a list of N contour cells (column 1),
        their edge points (columns 2,3), segment lengths (column 4), and
        outward unit vectors (column 5)

    Notes
    -----
    - Both `connecttype` methods functionally provide the same output when
      using the geocontour with data
    - Simplification will functionally reduce size and provide the same
      output when using the geocontour with data, but the directional
      information contained in the contour may not be preserved
    """
    gcc.ccontour(contour,latitudes,longitudes)
    contourdiff=np.diff(contour,axis=0)
    if connecttype=='center':
        geocontour=np.full((2*(contour.shape[0]-1),2,5),np.nan)
        geocontour[:,:,0]=np.stack((contour[:-1],contour[1:]),axis=1).reshape(-1,2)
        geocontour[:,:,1]=np.stack((contour[:-1],contour[:-1]+contourdiff/2),axis=1).reshape(-1,2)
        geocontour[:,:,2]=np.stack((contour[1:]-contourdiff/2,contour[1:]),axis=1).reshape(-1,2)
        if simplify:
            deleteindices=np.array([],dtype='int')
            geocontour_unique,countindex=np.unique(geocontour[:,:,0],axis=0,return_counts=True)
            for k in geocontour_unique[countindex>1]:
                indices=np.nonzero(((geocontour[:,:,0]-k)==0).all(axis=1))[0]
                starts=geocontour[indices,:,1]
                ends=geocontour[indices,:,2]
                uniquestart=starts[(starts-ends[:,None]).any(axis=2).all(axis=0)]
                uniqueend=ends[(ends-starts[:,None]).any(axis=2).all(axis=0)]
                if uniquestart.size==0 and uniqueend.size==0:
                    deleteindices=np.append(deleteindices,indices)
                else:
                    startkeep=np.argwhere((geocontour[:,:,1]==uniquestart[0]).all(axis=1))[0]==indices
                    endkeep=np.argwhere((geocontour[:,:,2]==uniqueend[0]).all(axis=1))[0]==indices
                    deleteindices=np.append(deleteindices,indices[~(startkeep+endkeep)])
            geocontour=np.delete(geocontour,deleteindices,axis=0)
    elif connecttype=='cell':
        geocontour=np.full((contour.shape[0]-1,2,5),np.nan)
        geocontour[:,:,0]=contour[:-1]
        geocontour[:,:,1]=geocontour[:,:,0]-np.roll(contourdiff/2,1,axis=0)
        geocontour[:,:,2]=geocontour[:,:,0]+contourdiff/2
        if simplify:
            deleteindices=np.array([],dtype='int')
            geocontour_unique,countindex=np.unique(geocontour[:,:,0],axis=0,return_counts=True)
            for k in geocontour_unique[countindex>1]:
                indices=np.nonzero(((geocontour[:,:,0]-k)==0).all(axis=1))[0]
                starts=geocontour[indices,:,1]
                ends=geocontour[indices,:,2]
                uniquestart=starts[(starts-ends[:,None]).any(axis=2).all(axis=0)]
                uniqueend=ends[(ends-starts[:,None]).any(axis=2).all(axis=0)]
                if uniquestart.size==0 and uniqueend.size==0:
                    deleteindices=np.append(deleteindices,indices)
                elif uniquestart.size==2 and uniqueend.size==2:
                    geocontour[indices[0],:,1]=uniquestart[0]
                    geocontour[indices[0],:,2]=uniqueend[0]
                    deleteindices=np.append(deleteindices,indices[1:])
            deleteindices=np.append(deleteindices,np.nonzero(((geocontour[:,:,1]-geocontour[:,:,2])==0).all(axis=1))[0])
            geocontour=np.delete(geocontour,deleteindices,axis=0)
    else:
        sys.exit('ERROR - unrecognized connecttype input '+connecttype+', valid options are \'cell\'/\'center\'')
    gridnorm=np.sqrt(np.sum((geocontour[:,:,2]-geocontour[:,:,1])**2,axis=1))
    geocontour[:,0,3]=gridnorm
    lonlengths=gcg.lonlen(geocontour[:,0,0])
    latlengths=gcg.latlen(geocontour[:,0,0])
    lengths=np.stack((latlengths,lonlengths),axis=1)
    geocontour[:,1,3]=np.sqrt(np.sum(((geocontour[:,:,2]-geocontour[:,:,1])*lengths)**2,axis=1))
    winding=np.sum((contour[:-1,0]+contour[1:,0])*np.diff(contour[:,1]))
    if winding<0:
        rotation=np.array([-1,1])
    elif winding>0:
        rotation=np.array([1,-1])
    else:
        sys.exit('ERROR - Can\'t determine orientation of contour')
    gridnorm[gridnorm==0]=1
    geocontour[:,:,4]=(geocontour[:,:,2]-geocontour[:,:,1])[:,::-1]*rotation/gridnorm[:,None]
    return geocontour

