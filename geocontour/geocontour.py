import sys
import numpy as np
import geocontour.check as gcc
import geocontour.grid as gcg

def build(contour,latitudes,longitudes,connecttype='cell',simplify=False):
    """
    Inputs (Required):
        contour - A 2-d Nx2 numpy array of ordered latitude/longitude points (degrees) describing the edge of a mask
        latitudes - An evenly spaced numpy array of latitude points (degrees)
        longitudes - An evenly spaced numpy array of longitude points (degrees)

    Inputs (Optional):
        connecttype ('cell'/'center') - method of linking contour cells, default='cell'
            'cell' implies drawing a single connection through a cell from the preceeding to the following cell, resulting in a geocontour the same length as the input contour
            'center' implies drawing two lines through a cell, intersecting the center of the cell, resulting in a geocontour double the length of the input contour
            Both methods functionally provide the same output when using the geocontour with data
        simplify (True/False) - remove cells for which outward vectors sum to 0 and merge cells containing multiple contour segments, default=False

    Outputs:
        geocontour - A 3-d Nx2x5 numpy array defining a list of N contour cells and their edge points, lengths, and outward unit vectors
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
    lonlengths=gcg.longitudelength(geocontour[:,0,0])
    latlengths=gcg.latitudelength(geocontour[:,0,0])
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

