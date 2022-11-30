import sys
import numpy as np
import geocontour.check as gcc
import geocontour.grid as gcg

def findstart(mask,searchdir='ru'):
    """
    Returns a starting cell for a contour, given a mask and a search criteria

    Inputs (Required):
        mask - A 2-d boolean numpy array of dimension MxN where M=len(latitudes) and N=len(longitudes)
        searchdir - A 2 character string describing the search directions
            e.g. 'dr' will search rows moving downward and columns moving rightward, in that order
            appropriate searches for cw contours: 'dr','ul','ru','ld' and ccw contours: 'dl','ur','rd','lu'

    Outputs:
        startcell - A 1x2 numpy array describing the index of the start cell
        startorientation - A 1x2 numpy array describing the orientation (entry direction) of the start cell
    """
    gcc.cmask(mask)
    if mask.sum()==0:
        sys.exit('ERROR - Input mask has no full cells')
    if searchdir=='dr':
        sortcode=[0.,1.]
    elif searchdir=='dl':
        sortcode=[0.,-1.]
    elif searchdir=='ur':
        sortcode=[-0.,1.]
    elif searchdir=='ul':
        sortcode=[-0.,-1.]
    elif searchdir=='rd':
        sortcode=[1.,0.]
    elif searchdir=='ru':
        sortcode=[1.,-0.]
    elif searchdir=='ld':
        sortcode=[-1.,0.]
    elif searchdir=='lu':
        sortcode=[-1.,-0.]
    else:
        sys.exit('ERROR - searchdir=\''+searchdir+'\' is not a valid selection, valid selections are \'dr\'\'dl\'\'ur\'\'ul\'\'rd\'\'ru\'\'ld\'\'lu\'')
    s0s=int(np.copysign(1,sortcode[0]))
    s0c=int(sortcode[0])
    s1s=int(np.copysign(1,sortcode[1]))
    s1c=int(sortcode[1])
    startcell=sorted(np.argwhere(mask==True),key=lambda x: (s0s*x[s0c],s1s*x[s1c]))[0]
    startorientation=(abs(np.array(sortcode))*np.copysign(1,sortcode[1])).astype('int')
    return startcell,startorientation

def parsestart(start,buffermask):
    """
    Checks start input for contour tracing
    Mainly used internally for contour trace functions

    Inputs (Required):
        start ('auto'/array) - A string or a 2x2 array describind the start cell and orientation, default='auto'
            'auto' will find a proper starting cell/orientation based on the direction input
            if not 'auto' input should be a 2x2 numpy array or nested list with row 1 describing the index of the start cell and row 2 describing the start orientation (e.g. [0,1] points right, [1,0] points down)
        buffermask - A 2-d boolean numpy array of dimension M+1xN+1 where M=len(latitudes) and N=len(longitudes)
        
    Outputs:
        startcell - A 1x2 numpy array describing the index of the start cell
        startorientation - A 1x2 numpy array describing the orientation (entry direction) of the start cell
    """
    start=np.array(start)
    if start.dtype!='int':
        sys.exit('ERROR - Start input is wrong datatype, provide start input as int only')
    if start.shape!=(2,2):
        sys.exit('ERROR - Start input in wrong format, provide a numpy array of shape (2,2) (or a 2x2 nested list), with row 1 containing start cell indices/coordinates and row 2 containing start orientation unit vector')
    startcell=start[0,:]+1
    startorientation=start[1,:]
    if buffermask[startcell[0],startcell[1]]==False:
        sys.exit('ERROR - Startcell '+str(startcell-1)+' is not in mask, select a start cell in the mask')
    if not ((startorientation==0)+(startorientation==1)+(startorientation==-1)).all() or np.linalg.norm(startorientation)!=1:
        sys.exit('ERROR - Startorientation '+str(startorientation)+' should be a unit vector pointing up/down/left/right')
    return startcell,startorientation

def setstop(stop,startvisits,startcell,startorientation):
    """
    Returns a stopping function for use in contour tracing while loop
    Mainly used internally for contour trace functions

    Inputs (Required):
        stop ('Elisoff'/'Nvisits'/'either') - A string selecting the stopping criterion, default='either'
            'Elisoff' stops when the start cell has been re-visited with the same orientation as started with
            'Nvisits' stops when the start cell has been re-visited N number of times (N set by startvisits)
            'either' stops when either Elisoff or Nvisits has been satisfied
        startvisits - An integer condition for the number of times re-visiting the start cell will trigger an end to the search, default=3
            only applicable when stop='Nvisits' or 'either'
        startcell - A 1x2 numpy array describing the index of the start cell
        startorientation - A 1x2 numpy array describing the orientation (entry direction) of the start cell

    Outputs:
        checkbreak - A stopping function that takes a cell, orientation, and start visit counter as input
    """
    if type(startvisits) is not int:
        sys.exit('ERROR - startvisits is not an integer, must be an integer')
    if stop=='either':
        def checkbreak(cell,orientation,Nvisits):
            if (Nvisits>=startvisits) or ((cell==startcell).all() and (orientation==startorientation).all()):
                return True
            else:
                return False
    elif stop=='Elisoff':
        def checkbreak(cell,orientation,Nvisits):
            if (cell==startcell).all() and (orientation==startorientation).all():
                return True
            else:
                return False
    elif stop=='Nvisits':
        def checkbreak(cell,orientation,Nvisits):
            if Nvisits>=startvisits:
                return True
            else:
                return False
    else:
        sys.exit('ERROR - stop=\''+stop+'\' is not a valid selection, valid selections are \'Elisoff\'/\'Nvisits\'/\'either\'')
    return checkbreak

def clean(contourcells,searchcells,latitudes=None,longitudes=None,closecontour=True,remcontourrepeat=True,remsearchrepeat=False):
    """
    Returns a cleaned contour that will pass checks
    Mainly used internally for contour trace functions

    Inputs (Required):
        contourcells - A python list of 1x2 numpy arrays containing contour indices
        searchcells - A python list of 1x2 numpy arrays containing contour search indices

    Inputs (Optional):
        latitudes - An evenly spaced numpy array of latitude points (degrees)
        longitudes - An evenly spaced numpy array of longitude points (degrees)
            Note: if latitudes/longitudes not provided, output will be in indices of input mask
        remcontourrepeat (True/False) - A boolean for selecting whether to remove consecutive repeating cells in the output contour, default=True
        remsearchrepeat (True/False) - A boolean for selecting whether to remove consecutive repeating cells in the output contoursearch, default=False
        closecontour (True/False) - A boolean for selecting whether to close the contour (first cell = last cell), default=True

    Outputs:
        contour - A 2-d Nx2 numpy array of ordered latitude/longitude points (degrees) describing the contour trace of a mask
        contoursearch - A 2-d Nx2 numpy array of ordered latitude/longitude points (degrees) describing the cells searched during contour tracing
            Note: If latitudes/longitudes not provided, returned contour/contoursearch will be indices of the input mask
    """
    if closecontour:
        contourcells.append(contourcells[0])
        searchcells.append(searchcells[0])
    contour=np.array(contourcells)
    contoursearch=np.array(searchcells)
    if remcontourrepeat:
        contour=contour[(np.diff(contour,axis=0,prepend=np.nan)!=0).any(axis=1)]
    if remsearchrepeat:
        contoursearch=contoursearch[(np.diff(contoursearch,axis=0,prepend=np.nan)!=0).any(axis=1)]
    if latitudes is not None and longitudes is not None:
        latspc=gcg.gridspacing(latitudes)
        lonspc=gcg.gridspacing(longitudes)
        latdir=gcg.checklatitudedirection(latitudes)
        if latdir=='inc':
            latitudes_ext=np.concatenate(([latitudes[0]-latspc],latitudes,[latitudes[-1]+latspc]))
        elif latdir=='dec':
            latitudes_ext=np.concatenate(([latitudes[0]+latspc],latitudes,[latitudes[-1]-latspc]))
        longitudes_ext=np.concatenate(([longitudes[0]-lonspc],longitudes,[longitudes[-1]+lonspc]))
        contour=np.stack((latitudes_ext[contour[:,0]],longitudes_ext[contour[:,1]]),axis=1)
        contoursearch=np.stack((latitudes_ext[contoursearch[:,0]],longitudes_ext[contoursearch[:,1]]),axis=1)
    else:
        contour-=1
        contoursearch-=1
    return contour,contoursearch

