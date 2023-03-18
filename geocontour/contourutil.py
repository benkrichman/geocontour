"""
Utility functions for operations involving contours and contour searches
========================================================================
"""
import sys
import numpy as np
import geocontour.check as gcc
import geocontour.grid as gcg

def findstart(mask,searchdir='ru'):
    """
    Find starting cell for contour, given a mask and a search criteria

    Parameters
    ----------
    mask : ndarray
        2D MxN bool array where M=len(`latitudes`) and
        N=len(`longitudes`)
    searchdir : {'dr', 'ul', 'ru', 'ld', 'dl', 'ur', 'rd', 'lu'}
        2 char string describing the search directions
        [e.g. 'dr' will search rows moving downward and columns moving
        rightward, in that order]

            for a desired clockwise search:
                {'dr', 'ul', 'ru', 'ld'}
            for a desired counterclockwise search:
                {'dl', 'ur', 'rd', 'lu'}

    Returns
    -------
    startcell : ndarray 
        1x2 array containing the indices of the start cell
    startorientation : ndarray
        1x2 array describing the orientation (entry direction) of the
        start cell

    See Also
    --------
    parsestart
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
    Check start input for contour tracing
    
    Mainly used internally for contour trace functions

    Parameters
    ----------
    start : array_like 
        2x2 array describing the start cell and orientation
            
            row 1: indices of the start cell
                e.g. for start cell [2,3] second row (lat) and third
                column (lon)
            row 2: the start orientation
                e.g. orientation [0,1] points right, [1,0] points down
    buffermask : ndarray 
        2D M+1xN+1 bool array where M=len(`latitudes`) and
        N=len(`longitudes`)
        
    Returns
    -------
    startcell : ndarray 
        1x2 array containing the indices of the start cell
    startorientation : ndarray
        1x2 array describing the orientation (entry direction) of the
        start cell

    See Also
    --------
    findstart
    """
    start=np.array(start)
    if start.dtype!='int':
        sys.exit('ERROR - Start input is wrong datatype, provide start input as ints only')
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
    Create a stopping function for use in contour tracing while loop
    
    Mainly used internally for contour trace functions

    Parameters
    ----------
    stop : {'Elisoff', 'Nvisits', 'either'}, default='either'
        selector for the stopping criterion

            ``Elisoff``
                stops when the start cell has been re-visited with the
                same orientation as started with
            ``Nvisits``
                stops when the start cell has been re-visited N number
                of times (N set by `startvisits` parameter)
            ``either``
                stops when either Elisoff or Nvisits has been satisfied
    startvisits : int, default=3
        the number of times re-visiting the start cell will trigger an
        end to the search
    startcell : ndarray 
        1x2 array containing the indices of the start cell
    startorientation : ndarray
        1x2 array describing the orientation (entry direction) of the
        start cell

    Returns
    -------
    checkbreak : function
        stopping function that takes a cell, orientation, and start
        visit counter as input

    See Also
    --------
    parsestart

    Notes
    -----
    `startvisits` parameter will only be utilized when `stop` parameter
    is set to 'Nvisits' or 'either'
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

    Parameters
    ----------
    contourcells : array_like
        list of 1x2 numpy arrays containing contour indices
    searchcells : array_like
        list of 1x2 numpy arrays containing contour search indices
    latitudes : ndarray, optional
        1D Nx1 array of latitude points (degrees)
    longitudes : ndarray, optional
        1D Nx1 array of longitude points (degrees)
    remcontourrepeat : bool, default=True
        select whether to remove consecutive repeating cells in the
        output `contour`
    remsearchrepeat : bool, default=False
        select whether to remove consecutive repeating cells in the
        output `contoursearch`
    closecontour : bool, default=True
        select whether to close the output `contour` (first cell = last
        cell)

    Returns
    -------
    contour : ndarray
        2D Nx2 array of ordered latitude/longitude points (degrees)
        describing the edge of a mask
    contoursearch : ndarray
        2D Nx2 array of ordered latitude/longitude points (degrees)
        describing the cells searched during contour tracing

    Notes
    -----
    If `latitudes`/`longitudes` not provided, returned
    `contour`/`contoursearch` will be indices of the input mask
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
        latspc=gcg.spacing(latitudes)
        lonspc=gcg.spacing(longitudes)
        latdir=gcg.clatdir(latitudes)
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

def fancysearch(contoursearch,contraction=0.2,shift=0.25):
    """
    Create a contour search that is visually more easy to follow via:

        - contraction, which ranges from 0 (no change) to 0.5
          (completely flattened to contour cell edges)
        - shift, which ranges from 0 (no change) to 0.5 (half a grid
          cell)

    Parameters
    ----------
    contoursearch : ndarray
        2D Nx2 array of ordered latitude/longitude points (degrees)
        describing the cells searched during contour tracing
    contraction : float, default=0.2 
        value determining how much `contoursearch` "shrinks" towards
        `contour` cell edges
    shift : float, default=0.25 
        value determining how much `contoursearch` "shifts" to avoid
        doubling back on itself

    Returns
    -------
    fancycontoursearch : ndarray 
        2D Nx2 array of ordered latitude/longitude points describing the
        cells searched during contour tracing, with `contraction` and
        `shift` applied

    See Also
    --------
    geocontour.output.plot
    """
    fdiff=np.diff(contoursearch,axis=0,prepend=[contoursearch[0,:]])
    if contraction>0.5:
        sys.exit('ERROR - contraction cannot be greater than 0.5, valid range is 0 to 0.5')
    elif contraction<0:
        sys.exit('ERROR - contraction cannot be less than 0, valid range is 0 to 0.5')
    else:
        zrange=np.arange(fdiff.shape[0])
        lzero=np.stack((zrange,zrange),axis=1)
        lzero[fdiff==0]=0
        lzero=np.maximum.accumulate(lzero,axis=0)
        filldiff=np.stack((fdiff[lzero[:,0],0],fdiff[lzero[:,1],1]),axis=1)
        fancycontoursearch=contoursearch.astype('float')-filldiff*contraction
    if shift>0.5:
        sys.exit('ERROR - shift cannot be greater than 0.5, valid range is 0 to 0.5')
    elif shift<0:
        sys.exit('ERROR - shift cannot be less than 0, valid range is 0 to 0.5')
    elif shift>0:
        fddby=np.convolve(contoursearch[:,0],[1,0,-1],mode='valid')
        fddbx=np.convolve(contoursearch[:,1],[1,0,-1],mode='valid')
        dbINDS=np.where(((fddby==0) & (fddbx==0)))[0]+1
        dbspans=fdiff[dbINDS+2]*shift
        fancycontoursearch[dbINDS+1]+=dbspans
        fancycontoursearch=np.insert(fancycontoursearch,dbINDS+1,fancycontoursearch[dbINDS]+dbspans,axis=0)
    return fancycontoursearch
