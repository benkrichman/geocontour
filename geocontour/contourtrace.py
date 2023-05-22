"""
Functions for tracing a contour on a lat/lon grid from an input mask
====================================================================
"""
import sys
import warnings
import numpy as np
import geocontour.check as gcc
import geocontour.maskutil as gcmu
import geocontour.contourutil as gccu

def square(mask,latitudes=None,longitudes=None,direction='cw',start='auto',stop='either',startvisits=3,checkconn=False,remcontourrepeat=True,remsearchrepeat=False,closecontour=True):
    """
    Find the contour of a mask using the square tracing algorithm
    
    Also known as:

        - Papert's turtle algorithm
        - simple boundary follower (SBF)

    Parameters
    ----------
    mask : ndarray
        2D MxN bool array where M=len(`latitudes`) and
        N=len(`longitudes`)
    latitudes : ndarray, optional
        1D Nx1 array of latitude points (degrees)
    longitudes : ndarray, optional
        1D Nx1 array of longitude points (degrees)
    direction : {'cw', 'ccw'}, default='cw'
        select the direction of contour tracing
    start : {'auto', array_like}, default='auto' 
        either a selection for automatic start cell assignment or a
        manual assignment via a 2x2 array describing the start cell and
        orientation

            row 1: indices of the start cell
                e.g. for start cell [2,3] second row (lat) and third
                column (lon)
            row 2: the start orientation
                e.g. orientation [0,1] points right, [1,0] points down
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
    checkconn : bool, default=False
        select whether to check connectivity and warn the user of
        potential issues, default=False
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
    
    See Also
    --------
    geocontour.contourutil.findstart
    geocontour.contourutil.setstop
    moore
    moore_imp
    pavlidis
    pavlidis_imp
    MSBF
    ISBF
    TSR

    Notes
    -----
    - If `latitudes`/`longitudes` not provided, returned
      `contour`/`contoursearch` will be indices of the input mask
    - `startvisits` parameter will only be utilized when `stop`
      parameter is set to 'Nvisits' or 'either'

    References
    ----------
    Ghuneim, A.G. (2000). Contour Tracing. McGill University.
    <https://www.imageprocessingplace.com/downloads_V3/root_downloads/tutorials/contour_tracing_Abeer_George_Ghuneim/alg.html>

    Gose, E., Johnsonbaugh, R., & Jost, S. (1996). Pattern recognition
    and image analysis. Prentice Hall PTR.

    Papert, S. (1973). Uses of Technology to Enhance Education (No. 298;
    Artificial Intelligence). Massachusetts Institute of Technology.
    <https://dspace.mit.edu/handle/1721.1/6213>

    Seo, J., Chae, S., Shim, J., Kim, D., Cheong, C., & Han, T.-D.
    (2016). Fast Contour-Tracing Algorithm Based on a Pixel-Following
    Method for Image Sensors. Sensors, 16(3), 353.
    <https://doi.org/10.3390/s16030353>

    Toussaint, G.T. (2010). Grids Connectivity and Contour Tracing
    [Lesson Notes]. McGill University.
    <http://www-cgrl.cs.mcgill.ca/~godfried/teaching/mir-reading-assignments/Chapter-2-Grids-Connectivity-Contour-Tracing.pdf>
    """
    if latitudes is not None and longitudes is not None:
        gcc.cmask(mask,latitudes,longitudes)
    else:
        gcc.cmask(mask)
    buffermask=np.full((tuple(np.array(mask.shape)+2)),False)
    buffermask[1:-1,1:-1]=mask
    if checkconn:
        fullconnectivity=gcmu.conn(buffermask,checkcells='full',connectivity=4)
        emptyconnectivity=gcmu.conn(buffermask,checkcells='empty',connectivity=4)
        if not fullconnectivity or not emptyconnectivity:
            warnings.warn('WARNING - Mask and/or non-mask not 4-connected: square tracing may not extract the full contour')
    if direction=='cw':
        def outsideturn(orientation):
            return orientation[::-1]*np.array([-1,1])
        def insideturn(orientation):
            return orientation[::-1]*np.array([1,-1])
        searchdir='ru'
    elif direction=='ccw':
        def outsideturn(orientation):
            return orientation[::-1]*np.array([1,-1])
        def insideturn(orientation):
            return orientation[::-1]*np.array([-1,1])
        searchdir='rd'
    else:
        sys.exit('ERROR - direction=\''+direction+'\' is not a valid selection, valid selections are \'cw\'/\'ccw\'')
    if type(start) is not str:
        startcell,startorientation=gccu.parsestart(start,buffermask)
    elif start=='auto':
        startcell,startorientation=gccu.findstart(buffermask,searchdir=searchdir)
    else:
        sys.exit('ERROR - start=\''+start+'\' is not a valid selection, valid selections are \'auto\' or a 2x2 array describing start cell index and start orientation')
    checkbreak=gccu.setstop(stop,startvisits,startcell,startorientation)
    searchcells=[]
    contourcells=[]
    orientation=startorientation
    cell=startcell
    Nvisits=0
    breakloop=False
    while not breakloop:
        searchcells.append(cell)
        if buffermask[cell[0],cell[1]]==True:
            contourcells.append(cell)
            orientation=outsideturn(orientation)
        else:
            orientation=insideturn(orientation)
        cell=cell+orientation
        if (cell==startcell).all():
            Nvisits+=1
        breakloop=checkbreak(cell,orientation,Nvisits)
    contour,contoursearch=gccu.clean(contourcells,searchcells,latitudes=latitudes,longitudes=longitudes,closecontour=closecontour,remcontourrepeat=remcontourrepeat,remsearchrepeat=remsearchrepeat)
    return contour,contoursearch

def moore(mask,latitudes=None,longitudes=None,direction='cw',start='auto',stop='either',startvisits=3,remcontourrepeat=True,remsearchrepeat=False,closecontour=True):
    """
    Find the contour of a mask using the Moore neighbor tracing
    algorithm

    Parameters
    ----------
    mask : ndarray
        2D MxN bool array where M=len(`latitudes`) and
        N=len(`longitudes`)
    latitudes : ndarray, optional
        1D Nx1 array of latitude points (degrees)
    longitudes : ndarray, optional
        1D Nx1 array of longitude points (degrees)
    direction : {'cw', 'ccw'}, default='cw'
        select the direction of contour tracing
    start : {'auto', array_like}, default='auto' 
        either a selection for automatic start cell assignment or a
        manual assignment via a 2x2 array describing the start cell and
        orientation

            row 1: indices of the start cell
                e.g. for start cell [2,3] second row (lat) and third
                column (lon)
            row 2: the start orientation
                e.g. orientation [0,1] points right, [1,0] points down
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
    
    See Also
    --------
    geocontour.contourutil.findstart
    geocontour.contourutil.setstop
    square
    moore_imp
    pavlidis
    pavlidis_imp
    MSBF
    ISBF
    TSR

    Notes
    -----
    - If `latitudes`/`longitudes` not provided, returned
      `contour`/`contoursearch` will be indices of the input mask
    - `startvisits` parameter will only be utilized when `stop`
      parameter is set to 'Nvisits' or 'either'

    References
    ----------
    Ghuneim, A.G. (2000). Contour Tracing. McGill University.
    <https://www.imageprocessingplace.com/downloads_V3/root_downloads/tutorials/contour_tracing_Abeer_George_Ghuneim/alg.html>

    Toussaint, G.T. (2010). Grids Connectivity and Contour Tracing
    [Lesson Notes]. McGill University.
    <http://www-cgrl.cs.mcgill.ca/~godfried/teaching/mir-reading-assignments/Chapter-2-Grids-Connectivity-Contour-Tracing.pdf>
    """
    if latitudes is not None and longitudes is not None:
        gcc.cmask(mask,latitudes,longitudes)
    else:
        gcc.cmask(mask)
    buffermask=np.full((tuple(np.array(mask.shape)+2)),False)
    buffermask[1:-1,1:-1]=mask
    if direction=='cw':
        searchdirstrings=['ru','dr','ld','ul']
        def insideturn(startorientation):
            return startorientation[::-1]*np.array([1,-1])
    elif direction=='ccw':
        searchdirstrings=['rd','ur','lu','dl']
        def insideturn(startorientation):
            return startorientation[::-1]*np.array([-1,1])
    else:
        sys.exit('ERROR - direction=\''+direction+'\' is not a valid selection, valid selections are \'cw\'/\'ccw\'')
    if type(start) is not str:
        startcell,startorientation=gccu.parsestart(start,buffermask)
    elif start=='auto':
        if stop=='either' or stop=='Elisoff':
            for ct,searchdir in enumerate(searchdirstrings):
                startcell,startorientation=gccu.findstart(buffermask,searchdir=searchdir)
                I=startcell+insideturn(startorientation)
                RI=I-startorientation
                RRI=RI-startorientation
                if buffermask[I[0],I[1]]==True or buffermask[RI[0],RI[1]]==True:
                    break
                else:
                    if ct<3:
                        continue
                    else:
                        if stop=='Elisoff':
                            sys.exit('ERROR - Suitable start cell for Elisoff stop not found, try reversing direction or manually select start')
                        if stop=='either':
                            warnings.warn('Warning - Suitable start cell for Elisoff stop not found, trace will stop on Nvisits only')
        else:
            startcell,startorientation=gccu.findstart(buffermask,searchdir='ru')
    else:
        sys.exit('ERROR - start=\''+start+'\' is not a valid selection, valid selections are \'auto\' or a 2x2 array describing start cell index and start orientation')
    checkbreak=gccu.setstop(stop,startvisits,startcell,startorientation)
    searchcells=[]
    contourcells=[]
    orientation=startorientation
    cell=startcell
    Nvisits=0
    breakloop=False
    while not breakloop:
        searchcells.append(cell)
        if buffermask[cell[0],cell[1]]==True:
            contourcells.append(cell)
            neighbors=gcmu.neighbors(cell,connectivity=8,direction=direction)
            nextneighborindex=(np.nonzero(((cell-orientation)==neighbors).all(axis=1))[0][0])%8
        else:
            nextneighborindex=(nextneighborindex+1)%8
        orientation=neighbors[nextneighborindex]-cell
        cell=neighbors[nextneighborindex]
        if (cell==startcell).all():
            Nvisits+=1
        breakloop=checkbreak(cell,orientation,Nvisits)
    contour,contoursearch=gccu.clean(contourcells,searchcells,latitudes=latitudes,longitudes=longitudes,closecontour=closecontour,remcontourrepeat=remcontourrepeat,remsearchrepeat=remsearchrepeat)
    return contour,contoursearch

def moore_imp(mask,latitudes=None,longitudes=None,direction='cw',start='auto',stop='either',startvisits=3,remcontourrepeat=True,remsearchrepeat=False,closecontour=True):
    """
    Find the contour of a mask using an improved Moore neighbor tracing
    algorithm that reliably captures inside corners

    Parameters
    ----------
    mask : ndarray
        2D MxN bool array where M=len(`latitudes`) and
        N=len(`longitudes`)
    latitudes : ndarray, optional
        1D Nx1 array of latitude points (degrees)
    longitudes : ndarray, optional
        1D Nx1 array of longitude points (degrees)
    direction : {'cw', 'ccw'}, default='cw'
        select the direction of contour tracing
    start : {'auto', array_like}, default='auto' 
        either a selection for automatic start cell assignment or a
        manual assignment via a 2x2 array describing the start cell and
        orientation

            row 1: indices of the start cell
                e.g. for start cell [2,3] second row (lat) and third
                column (lon)
            row 2: the start orientation
                e.g. orientation [0,1] points right, [1,0] points down
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
    
    See Also
    --------
    geocontour.contourutil.findstart
    geocontour.contourutil.setstop
    square
    moore
    pavlidis
    pavlidis_imp
    MSBF
    ISBF
    TSR

    Notes
    -----
    - If `latitudes`/`longitudes` not provided, returned
      `contour`/`contoursearch` will be indices of the input mask
    - `startvisits` parameter will only be utilized when `stop`
      parameter is set to 'Nvisits' or 'either'
    """
    if latitudes is not None and longitudes is not None:
        gcc.cmask(mask,latitudes,longitudes)
    else:
        gcc.cmask(mask)
    buffermask=np.full((tuple(np.array(mask.shape)+2)),False)
    buffermask[1:-1,1:-1]=mask
    if direction=='cw':
        searchdirstrings=['ru','dr','ld','ul']
        def insideturn(startorientation):
            return startorientation[::-1]*np.array([1,-1])
    elif direction=='ccw':
        searchdirstrings=['rd','ur','lu','dl']
        def insideturn(startorientation):
            return startorientation[::-1]*np.array([-1,1])
    else:
        sys.exit('ERROR - direction=\''+direction+'\' is not a valid selection, valid selections are \'cw\'/\'ccw\'')
    if type(start) is not str:
        startcell,startorientation=gccu.parsestart(start,buffermask)
    elif start=='auto':
        if stop=='either' or stop=='Elisoff':
            for ct,searchdir in enumerate(searchdirstrings):
                startcell,startorientation=gccu.findstart(buffermask,searchdir=searchdir)
                I=startcell+insideturn(startorientation)
                RI=I-startorientation
                RRI=RI-startorientation
                if buffermask[I[0],I[1]]==True or buffermask[RI[0],RI[1]]==True:
                    break
                else:
                    if ct<3:
                        continue
                    else:
                        if stop=='Elisoff':
                            sys.exit('ERROR - Suitable start cell for Elisoff stop not found, try reversing direction or manually select start')
                        if stop=='either':
                            warnings.warn('Warning - Suitable start cell for Elisoff stop not found, trace will stop on Nvisits only')
        else:
            startcell,startorientation=gccu.findstart(buffermask,searchdir='ru')
    else:
        sys.exit('ERROR - start=\''+start+'\' is not a valid selection, valid selections are \'auto\' or a 2x2 array describing start cell index and start orientation')
    checkbreak=gccu.setstop(stop,startvisits,startcell,startorientation)
    cell=startcell
    orientation=startorientation
    searchcells=[]
    contourcells=[]
    searchcells.append(cell)
    contourcells.append(cell)
    neighbors=gcmu.neighbors(cell,connectivity=8,direction=direction)
    nextneighborindex=(np.nonzero(((cell-orientation)==neighbors).all(axis=1))[0][0]+1)%8
    Nvisits=0
    breakloop=False
    while not breakloop:
        orientation=neighbors[nextneighborindex]-cell
        cell=neighbors[nextneighborindex]
        searchcells.append(cell)
        if buffermask[cell[0],cell[1]]==True:
            nextnextneighbor=neighbors[(nextneighborindex+1)%8]
            nextnextinmask=(buffermask[nextnextneighbor[0],nextnextneighbor[1]]==True)
            normdist=np.linalg.norm(contourcells[-1]-cell)
            if nextnextinmask and normdist==np.sqrt(2):
                searchcells.insert(-1,nextnextneighbor)
                contourcells.append(nextnextneighbor)
            contourcells.append(cell)
            neighbors=gcmu.neighbors(cell,connectivity=8,direction=direction)
            nextneighborindex=(np.nonzero(((cell-orientation)==neighbors).all(axis=1))[0][0])%8
        else:
            nextneighborindex=(nextneighborindex+1)%8
        if (cell==startcell).all():
            Nvisits+=1
        breakloop=checkbreak(cell,orientation,Nvisits)
    contour,contoursearch=gccu.clean(contourcells,searchcells,latitudes=latitudes,longitudes=longitudes,closecontour=closecontour,remcontourrepeat=remcontourrepeat,remsearchrepeat=remsearchrepeat)
    return contour,contoursearch

def pavlidis(mask,latitudes=None,longitudes=None,direction='cw',start='auto',stop='Nvisits',startvisits=1,remcontourrepeat=True,remsearchrepeat=False,closecontour=True):
    """
    Find the contour of a mask using the Pavlidis algorithm

    Parameters
    ----------
    mask : ndarray
        2D MxN bool array where M=len(`latitudes`) and
        N=len(`longitudes`)
    latitudes : ndarray, optional
        1D Nx1 array of latitude points (degrees)
    longitudes : ndarray, optional
        1D Nx1 array of longitude points (degrees)
    direction : {'cw', 'ccw'}, default='cw'
        select the direction of contour tracing
    start : {'auto', array_like}, default='auto' 
        either a selection for automatic start cell assignment or a
        manual assignment via a 2x2 array describing the start cell and
        orientation

            row 1: indices of the start cell
                e.g. for start cell [2,3] second row (lat) and third
                column (lon)
            row 2: the start orientation
                e.g. orientation [0,1] points right, [1,0] points down
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
    startvisits : int, default=1
        the number of times re-visiting the start cell will trigger an
        end to the search
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
    
    See Also
    --------
    geocontour.contourutil.findstart
    geocontour.contourutil.setstop
    square
    moore
    moore_imp
    pavlidis_imp
    MSBF
    ISBF
    TSR

    Notes
    -----
    - If `latitudes`/`longitudes` not provided, returned
      `contour`/`contoursearch` will be indices of the input mask
    - `startvisits` parameter will only be utilized when `stop`
      parameter is set to 'Nvisits' or 'either'

    References
    ----------
    Ghuneim, A.G. (2000). Contour Tracing. McGill University.
    <https://www.imageprocessingplace.com/downloads_V3/root_downloads/tutorials/contour_tracing_Abeer_George_Ghuneim/alg.html>

    Pavlidis, T. (1982) Algorithms for Graphics and Image Processing.
    Computer Science Press, New York, NY.
    <https://doi.org/10.1007/978-3-642-93208-3>
    """
    if latitudes is not None and longitudes is not None:
        gcc.cmask(mask,latitudes,longitudes)
    else:
        gcc.cmask(mask)
    buffermask=np.full((tuple(np.array(mask.shape)+2)),False)
    buffermask[1:-1,1:-1]=mask
    if direction=='cw':
        def outsideturn(orientation):
            return orientation[::-1]*np.array([-1,1])
        def insideturn(orientation):
            return orientation[::-1]*np.array([1,-1])
        searchdirstrings=['ru','dr','ld','ul']
    elif direction=='ccw':
        def outsideturn(orientation):
            return orientation[::-1]*np.array([1,-1])
        def insideturn(orientation):
            return orientation[::-1]*np.array([-1,1])
        searchdirstrings=['rd','ur','lu','dl']
    else:
        sys.exit('ERROR - direction=\''+direction+'\' is not a valid selection, valid selections are \'cw\'/\'ccw\'')
    if type(start) is not str:
        startcell,startorientation=gccu.parsestart(start,buffermask)
    elif start=='auto':
        for ct,searchdir in enumerate(searchdirstrings):
            startcell,startorientation=gccu.findstart(buffermask,searchdir=searchdir)
            RI=startcell-startorientation+insideturn(startorientation)
            if buffermask[RI[0],RI[1]]==False:
                break
            else:
                if ct<3:
                    continue
                else:
                    warnings.warn('WARNING - Suitable start cell not found, trace may miss cells')
    else:
        sys.exit('ERROR - start=\''+start+'\' is not a valid selection, valid selections are \'auto\' or a 2x2 array describing start cell index and start orientation')
    if stop=='Elisoff':
        warnings.warn('WARNING - Pavlidis tracing often will not stop on the Elisoff condition, tracing may become stuck in an infinite loop')
    checkbreak=gccu.setstop(stop,startvisits,startcell,startorientation)
    searchcells=[]
    contourcells=[]
    searchcells.append(startcell)
    contourcells.append(startcell)
    orientation=startorientation
    cell=startcell
    rotations=0
    Nvisits=0
    breakloop=False
    while not breakloop:
        FO=cell+orientation+outsideturn(orientation)
        F=cell+orientation
        FI=cell+orientation+insideturn(orientation)
        if buffermask[FO[0],FO[1]]==True:
            searchcells.append(FO)
            contourcells.append(FO)
            cell=FO
            orientation=outsideturn(orientation)
            rotations=0
        elif buffermask[F[0],F[1]]==True:
            searchcells.append(FO)
            searchcells.append(F)
            contourcells.append(F)
            cell=F
            rotations=0
        elif buffermask[FI[0],FI[1]]==True:
            searchcells.append(FO)
            searchcells.append(F)
            searchcells.append(FI)
            contourcells.append(FI)
            cell=FI
            rotations=0
        else:
            searchcells.append(FO)
            searchcells.append(F)
            searchcells.append(FI)
            orientation=insideturn(orientation)
            rotations+=1
        if rotations>3:
            sys.exit('ERROR - Stuck on isolated cell '+str(cell))
        if (cell==startcell).all() and len(contourcells)>1:
            Nvisits+=1
        breakloop=checkbreak(cell,orientation,Nvisits)
    contour,contoursearch=gccu.clean(contourcells,searchcells,latitudes=latitudes,longitudes=longitudes,closecontour=closecontour,remcontourrepeat=remcontourrepeat,remsearchrepeat=remsearchrepeat)
    return contour,contoursearch

def pavlidis_imp(mask,latitudes=None,longitudes=None,direction='cw',start='auto',stop='Nvisits',startvisits=1,remcontourrepeat=True,remsearchrepeat=False,closecontour=True):
    """
    Find the contour of a mask using an improved Pavlidis algorithm that
    reliably captures inside corners

    Parameters
    ----------
    mask : ndarray
        2D MxN bool array where M=len(`latitudes`) and
        N=len(`longitudes`)
    latitudes : ndarray, optional
        1D Nx1 array of latitude points (degrees)
    longitudes : ndarray, optional
        1D Nx1 array of longitude points (degrees)
    direction : {'cw', 'ccw'}, default='cw'
        select the direction of contour tracing
    start : {'auto', array_like}, default='auto' 
        either a selection for automatic start cell assignment or a
        manual assignment via a 2x2 array describing the start cell and
        orientation

            row 1: indices of the start cell
                e.g. for start cell [2,3] second row (lat) and third
                column (lon)
            row 2: the start orientation
                e.g. orientation [0,1] points right, [1,0] points down
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
    startvisits : int, default=1
        the number of times re-visiting the start cell will trigger an
        end to the search
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
    
    See Also
    --------
    geocontour.contourutil.findstart
    geocontour.contourutil.setstop
    square
    moore
    moore_imp
    pavlidis
    MSBF
    ISBF
    TSR

    Notes
    -----
    - If `latitudes`/`longitudes` not provided, returned
      `contour`/`contoursearch` will be indices of the input mask
    - `startvisits` parameter will only be utilized when `stop`
      parameter is set to 'Nvisits' or 'either'
    """
    if latitudes is not None and longitudes is not None:
        gcc.cmask(mask,latitudes,longitudes)
    else:
        gcc.cmask(mask)
    buffermask=np.full((tuple(np.array(mask.shape)+2)),False)
    buffermask[1:-1,1:-1]=mask
    if direction=='cw':
        def outsideturn(orientation):
            return orientation[::-1]*np.array([-1,1])
        def insideturn(orientation):
            return orientation[::-1]*np.array([1,-1])
        searchdirstrings=['ru','dr','ld','ul']
    elif direction=='ccw':
        def outsideturn(orientation):
            return orientation[::-1]*np.array([1,-1])
        def insideturn(orientation):
            return orientation[::-1]*np.array([-1,1])
        searchdirstrings=['rd','ur','lu','dl']
    else:
        sys.exit('ERROR - direction=\''+direction+'\' is not a valid selection, valid selections are \'cw\'/\'ccw\'')
    if type(start) is not str:
        startcell,startorientation=gccu.parsestart(start,buffermask)
    elif start=='auto':
        for ct,searchdir in enumerate(searchdirstrings):
            startcell,startorientation=gccu.findstart(buffermask,searchdir=searchdir)
            RI=startcell-startorientation+insideturn(startorientation)
            if buffermask[RI[0],RI[1]]==False:
                break
            else:
                if ct<3:
                    continue
                else:
                    warnings.warn('WARNING - Suitable start cell not found, trace may miss cells')
    else:
        sys.exit('ERROR - start=\''+start+'\' is not a valid selection, valid selections are \'auto\' or a 2x2 array describing start cell index and start orientation')
    if stop=='Elisoff':
        warnings.warn('WARNING - Pavlidis tracing often will not stop on the Elisoff condition, tracing may become stuck in an infinite loop')
    checkbreak=gccu.setstop(stop,startvisits,startcell,startorientation)
    searchcells=[]
    contourcells=[]
    searchcells.append(startcell)
    contourcells.append(startcell)
    orientation=startorientation
    cell=startcell
    rotations=0
    Nvisits=0
    breakloop=False
    while not breakloop:
        FO=cell+orientation+outsideturn(orientation)
        F=cell+orientation
        FI=cell+orientation+insideturn(orientation)
        I=cell+insideturn(orientation)
        if buffermask[FO[0],FO[1]]==True:
            if buffermask[F[0],F[1]]==True:
                contourcells.append(F)
            searchcells.append(F)
            searchcells.append(FO)
            contourcells.append(FO)
            cell=FO
            rotations=0
            orientation=outsideturn(orientation)
        elif buffermask[F[0],F[1]]==True:
            searchcells.append(FO)
            searchcells.append(F)
            contourcells.append(F)
            cell=F
            rotations=0
        elif buffermask[FI[0],FI[1]]==True:
            if buffermask[I[0],I[1]]==True:
                contourcells.append(I)
            searchcells.append(FO)
            searchcells.append(F)
            searchcells.append(I)
            searchcells.append(FI)
            contourcells.append(FI)
            cell=FI
            rotations=0
        else:
            searchcells.append(FO)
            searchcells.append(F)
            searchcells.append(FI)
            orientation=insideturn(orientation)
            rotations+=1
        if rotations>3:
            sys.exit('ERROR - Stuck on isolated cell '+str(cell))
        if (cell==startcell).all() and len(contourcells)>1:
            Nvisits+=1
        breakloop=checkbreak(cell,orientation,Nvisits)
    contour,contoursearch=gccu.clean(contourcells,searchcells,latitudes=latitudes,longitudes=longitudes,closecontour=closecontour,remcontourrepeat=remcontourrepeat,remsearchrepeat=remsearchrepeat)
    return contour,contoursearch

def MSBF(mask,latitudes=None,longitudes=None,direction='cw',start='auto',stop='either',startvisits=3,remcontourrepeat=True,remsearchrepeat=False,closecontour=True):
    """
    Find the contour of a mask using the modified simple boundary
    following algorithm

    Parameters
    ----------
    mask : ndarray
        2D MxN bool array where M=len(`latitudes`) and
        N=len(`longitudes`)
    latitudes : ndarray, optional
        1D Nx1 array of latitude points (degrees)
    longitudes : ndarray, optional
        1D Nx1 array of longitude points (degrees)
    direction : {'cw', 'ccw'}, default='cw'
        select the direction of contour tracing
    start : {'auto', array_like}, default='auto' 
        either a selection for automatic start cell assignment or a
        manual assignment via a 2x2 array describing the start cell and
        orientation

            row 1: indices of the start cell
                e.g. for start cell [2,3] second row (lat) and third
                column (lon)
            row 2: the start orientation
                e.g. orientation [0,1] points right, [1,0] points down
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
    checkconn : bool, default=False
        select whether to check connectivity and warn the user of
        potential issues, default=False
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
    
    See Also
    --------
    geocontour.contourutil.findstart
    geocontour.contourutil.setstop
    square
    moore
    moore_imp
    pavlidis
    pavlidis_imp
    ISBF
    TSR

    Notes
    -----
    - If `latitudes`/`longitudes` not provided, returned
      `contour`/`contoursearch` will be indices of the input mask
    - `startvisits` parameter will only be utilized when `stop`
      parameter is set to 'Nvisits' or 'either'

    References
    ----------
    Cheong, C.-H., & Han, T.-D. (2006). Improved Simple Boundary
    Following Algorithm. Journal of KIISE: Software and Applications,
    33(4), 427–439.
    <https://koreascience.kr/article/JAKO200622219415761.pdf>

    Cheong, C.-H., Seo, J., & Han, T.-D. (2006). Advanced Contour
    Tracing Algorithms based on Analysis of Tracing Conditions.
    Proceedings of the 33rd KISS Fall Conference, 33, 431–436.
    <https://koreascience.kr/article/CFKO200614539217302.pdf>

    Gose, E., Johnsonbaugh, R., & Jost, S. (1996). Pattern recognition
    and image analysis. Prentice Hall PTR.
    """
    if latitudes is not None and longitudes is not None:
        gcc.cmask(mask,latitudes,longitudes)
    else:
        gcc.cmask(mask)
    buffermask=np.full((tuple(np.array(mask.shape)+2)),False)
    buffermask[1:-1,1:-1]=mask
    if direction=='cw':
        def outsideturn(orientation):
            return orientation[::-1]*np.array([-1,1])
        def insideturn(orientation):
            return orientation[::-1]*np.array([1,-1])
        searchdir='ru'
    elif direction=='ccw':
        def outsideturn(orientation):
            return orientation[::-1]*np.array([1,-1])
        def insideturn(orientation):
            return orientation[::-1]*np.array([-1,1])
        searchdir='rd'
    else:
        sys.exit('ERROR - direction=\''+direction+'\' is not a valid selection, valid selections are \'cw\'/\'ccw\'')
    if type(start) is not str:
        startcell,startorientation=gccu.parsestart(start,buffermask)
    elif start=='auto':
        startcell,startorientation=gccu.findstart(buffermask,searchdir=searchdir)
    else:
        sys.exit('ERROR - start=\''+start+'\' is not a valid selection, valid selections are \'auto\' or a 2x2 array describing start cell index and start orientation')
    checkbreak=gccu.setstop(stop,startvisits,startcell,startorientation)
    orientation=startorientation
    cell=startcell
    searchcells=[]
    contourcells=[]
    contourcells.append(cell)
    Nvisits=0
    breakloop=False
    tag=False
    while not breakloop:
        searchcells.append(cell)
        if buffermask[cell[0],cell[1]]==True:
            contourcells.append(cell)
            O=cell+outsideturn(orientation)
            R=cell-orientation
            RO=O-orientation
            if tag==False and buffermask[O[0],O[1]]==False and buffermask[RO[0],RO[1]]==True and buffermask[R[0],R[1]]==False:
                orientation=-orientation
                cell=RO
                tag=True
            else:
                orientation=outsideturn(orientation)
                cell=O
                tag=False
        else:
            orientation=insideturn(orientation)
            cell=cell+orientation
            tag=False
        if (cell==startcell).all():
            Nvisits+=1
        breakloop=checkbreak(cell,orientation,Nvisits)
    contour,contoursearch=gccu.clean(contourcells,searchcells,latitudes=latitudes,longitudes=longitudes,closecontour=closecontour,remcontourrepeat=remcontourrepeat,remsearchrepeat=remsearchrepeat)
    return contour,contoursearch

def ISBF(mask,latitudes=None,longitudes=None,direction='cw',start='auto',stop='either',startvisits=3,remcontourrepeat=True,remsearchrepeat=False,closecontour=True):
    """
    Find the contour of a mask using the improved simple boundary
    following algorithm

    Parameters
    ----------
    mask : ndarray
        2D MxN bool array where M=len(`latitudes`) and
        N=len(`longitudes`)
    latitudes : ndarray, optional
        1D Nx1 array of latitude points (degrees)
    longitudes : ndarray, optional
        1D Nx1 array of longitude points (degrees)
    direction : {'cw', 'ccw'}, default='cw'
        select the direction of contour tracing
    start : {'auto', array_like}, default='auto' 
        either a selection for automatic start cell assignment or a
        manual assignment via a 2x2 array describing the start cell and
        orientation

            row 1: indices of the start cell
                e.g. for start cell [2,3] second row (lat) and third
                column (lon)
            row 2: the start orientation
                e.g. orientation [0,1] points right, [1,0] points down
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
    checkconn : bool, default=False
        select whether to check connectivity and warn the user of
        potential issues, default=False
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
    
    See Also
    --------
    geocontour.contourutil.findstart
    geocontour.contourutil.setstop
    square
    moore
    moore_imp
    pavlidis
    pavlidis_imp
    MSBF
    TSR

    Notes
    -----
    - If `latitudes`/`longitudes` not provided, returned
      `contour`/`contoursearch` will be indices of the input mask
    - `startvisits` parameter will only be utilized when `stop`
      parameter is set to 'Nvisits' or 'either'

    References
    ----------
    Cheong, C.-H., & Han, T.-D. (2006). Improved Simple Boundary
    Following Algorithm. Journal of KIISE: Software and Applications,
    33(4), 427–439.
    <https://koreascience.kr/article/JAKO200622219415761.pdf>

    Cheong, C.-H., Seo, J., & Han, T.-D. (2006). Advanced Contour
    Tracing Algorithms based on Analysis of Tracing Conditions.
    Proceedings of the 33rd KISS Fall Conference, 33, 431–436.
    <https://koreascience.kr/article/CFKO200614539217302.pdf>

    Seo, J., Chae, S., Shim, J., Kim, D., Cheong, C., & Han, T.-D.
    (2016). Fast Contour-Tracing Algorithm Based on a Pixel-Following
    Method for Image Sensors. Sensors, 16(3), 353.
    <https://doi.org/10.3390/s16030353>
    """
    if latitudes is not None and longitudes is not None:
        gcc.cmask(mask,latitudes,longitudes)
    else:
        gcc.cmask(mask)
    buffermask=np.full((tuple(np.array(mask.shape)+2)),False)
    buffermask[1:-1,1:-1]=mask
    if direction=='cw':
        def outsideturn(orientation):
            return orientation[::-1]*np.array([-1,1])
        def insideturn(orientation):
            return orientation[::-1]*np.array([1,-1])
        searchdir='ru'
    elif direction=='ccw':
        def outsideturn(orientation):
            return orientation[::-1]*np.array([1,-1])
        def insideturn(orientation):
            return orientation[::-1]*np.array([-1,1])
        searchdir='rd'
    else:
        sys.exit('ERROR - direction=\''+direction+'\' is not a valid selection, valid selections are \'cw\'/\'ccw\'')
    if type(start) is not str:
        startcell,startorientation=gccu.parsestart(start,buffermask)
    elif start=='auto':
        startcell,startorientation=gccu.findstart(buffermask,searchdir=searchdir)
    else:
        sys.exit('ERROR - start=\''+start+'\' is not a valid selection, valid selections are \'auto\' or a 2x2 array describing start cell index and start orientation')
    checkbreak=gccu.setstop(stop,startvisits,startcell,startorientation)
    orientation=startorientation
    cell=startcell
    searchcells=[]
    contourcells=[]
    searchcells.append(cell)
    contourcells.append(cell)
    Nvisits=0
    breakloop=False
    #following 2006 paper/conf (2016 paper missing tag component - goes into infinite loop on case 4b)
    tag=False
    while not breakloop:
        O=cell+outsideturn(orientation)
        if buffermask[O[0],O[1]]==True: #Case A
            searchcells.append(O)
            contourcells.append(O)
            cell=O
            orientation=outsideturn(orientation)
            tag=False
        else:
            R=cell-orientation
            RO=O-orientation
            if buffermask[RO[0],RO[1]]==True and buffermask[R[0],R[1]]==False and tag==False: #Case B
                searchcells.append(O)
                searchcells.append(R)
                searchcells.append(RO)
                contourcells.append(RO)
                cell=RO
                orientation=-orientation
                tag=True
            else:
                searchcells.append(RO)
                searchcells.append(O)
                F=cell+orientation
                FO=O+orientation
                if buffermask[FO[0],FO[1]]==True:
                    searchcells.append(F)
                    searchcells.append(FO)
                    if buffermask[F[0],F[1]]==True: #Case D
                        contourcells.append(F)
                        contourcells.append(FO)
                    else: #Case C
                        contourcells.append(FO)
                    cell=FO
                elif buffermask[F[0],F[1]]==True: #Case E
                    searchcells.append(FO)
                    searchcells.append(F)
                    contourcells.append(F)
                    cell=F
                    orientation=insideturn(orientation)
                else: #Case F
                    searchcells.append(FO)
                    searchcells.append(F)
                    searchcells.append(cell)
                    orientation=-orientation
                tag=False
        if (cell==startcell).all():
            Nvisits+=1
        breakloop=checkbreak(cell,orientation,Nvisits)
    contour,contoursearch=gccu.clean(contourcells,searchcells,latitudes=latitudes,longitudes=longitudes,closecontour=closecontour,remcontourrepeat=remcontourrepeat,remsearchrepeat=remsearchrepeat)
    return contour,contoursearch

def TSR(mask,latitudes=None,longitudes=None,direction='cw',start='auto',stop='either',startvisits=3,remcontourrepeat=True,remsearchrepeat=False,closecontour=True):
    """
    Find the contour of a mask using the two-step representative tracing
    algorithm

    Parameters
    ----------
    mask : ndarray
        2D MxN bool array where M=len(`latitudes`) and
        N=len(`longitudes`)
    latitudes : ndarray, optional
        1D Nx1 array of latitude points (degrees)
    longitudes : ndarray, optional
        1D Nx1 array of longitude points (degrees)
    direction : {'cw', 'ccw'}, default='cw'
        select the direction of contour tracing
    start : {'auto', array_like}, default='auto' 
        either a selection for automatic start cell assignment or a
        manual assignment via a 2x2 array describing the start cell and
        orientation

            row 1: indices of the start cell
                e.g. for start cell [2,3] second row (lat) and third
                column (lon)
            row 2: the start orientation
                e.g. orientation [0,1] points right, [1,0] points down
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
    checkconn : bool, default=False
        select whether to check connectivity and warn the user of
        potential issues, default=False
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
    
    See Also
    --------
    geocontour.contourutil.findstart
    geocontour.contourutil.setstop
    square
    moore
    moore_imp
    pavlidis
    pavlidis_imp
    MSBF
    ISBF

    Notes
    -----
    - If `latitudes`/`longitudes` not provided, returned
      `contour`/`contoursearch` will be indices of the input mask
    - `startvisits` parameter will only be utilized when `stop`
      parameter is set to 'Nvisits' or 'either'

    References
    ----------
    Seo, J., Chae, S., Shim, J., Kim, D., Cheong, C., & Han, T.-D.
    (2016). Fast Contour-Tracing Algorithm Based on a Pixel-Following
    Method for Image Sensors. Sensors, 16(3), 353.
    <https://doi.org/10.3390/s16030353>
    """
    if latitudes is not None and longitudes is not None:
        gcc.cmask(mask,latitudes,longitudes)
    else:
        gcc.cmask(mask)
    buffermask=np.full((tuple(np.array(mask.shape)+2)),False)
    buffermask[1:-1,1:-1]=mask
    if direction=='cw':
        def outsideturn(orientation):
            return orientation[::-1]*np.array([-1,1])
        def insideturn(orientation):
            return orientation[::-1]*np.array([1,-1])
        searchdir='ru'
    elif direction=='ccw':
        def outsideturn(orientation):
            return orientation[::-1]*np.array([1,-1])
        def insideturn(orientation):
            return orientation[::-1]*np.array([-1,1])
        searchdir='rd'
    else:
        sys.exit('ERROR - direction=\''+direction+'\' is not a valid selection, valid selections are \'cw\'/\'ccw\'')
    if type(start) is not str:
        startcell,startorientation=gccu.parsestart(start,buffermask)
    elif start=='auto':
        startcell,startorientation=gccu.findstart(buffermask,searchdir=searchdir)
    else:
        sys.exit('ERROR - start=\''+start+'\' is not a valid selection, valid selections are \'auto\' or a 2x2 array describing start cell index and start orientation')
    checkbreak=gccu.setstop(stop,startvisits,startcell,startorientation)
    orientation=startorientation
    cell=startcell
    searchcells=[]
    contourcells=[]
    searchcells.append(cell)
    contourcells.append(cell)
    Nvisits=0
    breakloop=False
    while not breakloop:
        O=cell+outsideturn(orientation)
        RO=O-orientation
        #Stage 1
        if buffermask[RO[0],RO[1]]==True:
            searchcells.append(O)
            searchcells.append(RO)
            if buffermask[O[0],O[1]]==True: #Case 1
                contourcells.append(O)
                contourcells.append(RO)
            else: #Case 2
                contourcells.append(RO)
            cell=RO
            orientation=-orientation
        else:
            searchcells.append(RO)
            searchcells.append(O)
            if buffermask[O[0],O[1]]==True: #Case 3
                contourcells.append(O)
                cell=O
                orientation=outsideturn(orientation)
            else: #Case 4
                pass
        #Stage 2
        F=cell+orientation
        FO=F+outsideturn(orientation)
        if buffermask[FO[0],FO[1]]==True:
            searchcells.append(F)
            searchcells.append(FO)
            if buffermask[F[0],F[1]]==True: #Case 6
                contourcells.append(F)
                contourcells.append(FO)
            else: #Case 5
                contourcells.append(FO)
            cell=FO
        else:
            searchcells.append(FO)
            searchcells.append(F)
            if buffermask[F[0],F[1]]==True: #Case 7
                contourcells.append(F)
                cell=F
                orientation=insideturn(orientation)
            else: #Case 8
                orientation=-orientation
        if (cell==startcell).all():
            Nvisits+=1
        breakloop=checkbreak(cell,orientation,Nvisits)
    contour,contoursearch=gccu.clean(contourcells,searchcells,latitudes=latitudes,longitudes=longitudes,closecontour=closecontour,remcontourrepeat=remcontourrepeat,remsearchrepeat=remsearchrepeat)
    return contour,contoursearch

