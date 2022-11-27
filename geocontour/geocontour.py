import sys
import numpy

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
    geocontour.check.checkcontour(contour,latitudes,longitudes)
    contourdiff=numpy.diff(contour,axis=0)
    if connecttype=='center':
        geocontour=numpy.full((2*(contour.shape[0]-1),2,5),numpy.nan)
        geocontour[:,:,0]=numpy.stack((contour[:-1],contour[1:]),axis=1).reshape(-1,2)
        geocontour[:,:,1]=numpy.stack((contour[:-1],contour[:-1]+contourdiff/2),axis=1).reshape(-1,2)
        geocontour[:,:,2]=numpy.stack((contour[1:]-contourdiff/2,contour[1:]),axis=1).reshape(-1,2)
        if simplify:
            deleteindices=numpy.array([],dtype='int')
            geocontour_unique,countindex=numpy.unique(geocontour[:,:,0],axis=0,return_counts=True)
            for k in geocontour_unique[countindex>1]:
                indices=numpy.nonzero(((geocontour[:,:,0]-k)==0).all(axis=1))[0]
                starts=geocontour[indices,:,1]
                ends=geocontour[indices,:,2]
                uniquestart=starts[(starts-ends[:,None]).any(axis=2).all(axis=0)]
                uniqueend=ends[(ends-starts[:,None]).any(axis=2).all(axis=0)]
                if uniquestart.size==0 and uniqueend.size==0:
                    deleteindices=numpy.append(deleteindices,indices)
                else:
                    startkeep=numpy.argwhere((geocontour[:,:,1]==uniquestart[0]).all(axis=1))[0]==indices
                    endkeep=numpy.argwhere((geocontour[:,:,2]==uniqueend[0]).all(axis=1))[0]==indices
                    deleteindices=numpy.append(deleteindices,indices[~(startkeep+endkeep)])
            geocontour=numpy.delete(geocontour,deleteindices,axis=0)
    elif connecttype=='cell':
        geocontour=numpy.full((contour.shape[0]-1,2,5),numpy.nan)
        geocontour[:,:,0]=contour[:-1]
        geocontour[:,:,1]=geocontour[:,:,0]-numpy.roll(contourdiff/2,1,axis=0)
        geocontour[:,:,2]=geocontour[:,:,0]+contourdiff/2
        if simplify:
            deleteindices=numpy.array([],dtype='int')
            geocontour_unique,countindex=numpy.unique(geocontour[:,:,0],axis=0,return_counts=True)
            for k in geocontour_unique[countindex>1]:
                indices=numpy.nonzero(((geocontour[:,:,0]-k)==0).all(axis=1))[0]
                starts=geocontour[indices,:,1]
                ends=geocontour[indices,:,2]
                uniquestart=starts[(starts-ends[:,None]).any(axis=2).all(axis=0)]
                uniqueend=ends[(ends-starts[:,None]).any(axis=2).all(axis=0)]
                if uniquestart.size==0 and uniqueend.size==0:
                    deleteindices=numpy.append(deleteindices,indices)
                elif uniquestart.size==2 and uniqueend.size==2:
                    geocontour[indices[0],:,1]=uniquestart[0]
                    geocontour[indices[0],:,2]=uniqueend[0]
                    deleteindices=numpy.append(deleteindices,indices[1:])
            deleteindices=numpy.append(deleteindices,numpy.nonzero(((geocontour[:,:,1]-geocontour[:,:,2])==0).all(axis=1))[0])
            geocontour=numpy.delete(geocontour,deleteindices,axis=0)
    else:
        sys.exit('ERROR - unrecognized connecttype input '+connecttype+', valid options are \'cell\'/\'center\'')
    gridnorm=numpy.sqrt(numpy.sum((geocontour[:,:,2]-geocontour[:,:,1])**2,axis=1))
    geocontour[:,0,3]=gridnorm
    lonlengths=geocontour.grid.longitudelength(geocontour[:,0,0])
    latlengths=geocontour.grid.latitudelength(geocontour[:,0,0])
    lengths=numpy.stack((latlengths,lonlengths),axis=1)
    geocontour[:,1,3]=numpy.sqrt(numpy.sum(((geocontour[:,:,2]-geocontour[:,:,1])*lengths)**2,axis=1))
    winding=numpy.sum((contour[:-1,0]+contour[1:,0])*numpy.diff(contour[:,1]))
    if winding<0:
        rotation=numpy.array([-1,1])
    elif winding>0:
        rotation=numpy.array([1,-1])
    else:
        sys.exit('ERROR - Can\'t determine orientation of contour')
    gridnorm[gridnorm==0]=1
    geocontour[:,:,4]=(geocontour[:,:,2]-geocontour[:,:,1])[:,::-1]*rotation/gridnorm[:,None]
    return geocontour

