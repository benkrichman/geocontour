import sys
import warnings
import numpy
import geocontour.check
import geocontour.maskutil
import geocontour.contourutil

def square(mask,latitudes=None,longitudes=None,direction='cw',start='auto',stop='either',startvisits=3,checkconn=False,remcontourrepeat=True,remsearchrepeat=False,closecontour=True):
    """
    Returns the contour trace of a mask input using the square tracing algorithm
    Source:
        https://www.imageprocessingplace.com/downloads_V3/root_downloads/tutorials/contour_tracing_Abeer_George_Ghuneim/square.html
        http://www-cgrl.cs.mcgill.ca/~godfried/teaching/mir-reading-assignments/Chapter-2-Grids-Connectivity-Contour-Tracing.pdf

    Inputs (Required):
        mask - A 2-d boolean numpy array of dimension MxN where M=len(latitudes) and N=len(longitudes)

    Inputs (Optional):
        latitudes - An evenly spaced numpy array of latitude points (degrees)
        longitudes - An evenly spaced numpy array of longitude points (degrees)
            Note: if latitudes/longitudes not provided, output will be in indices of input mask
        direction ('cw'/'ccw') - A string selecting the direction of the returned neighbors, default='cw'
        start ('auto'/array) - A string or a 2x2 array describind the start cell and orientation, default='auto'
            'auto' will find a proper starting cell/orientation based on the direction input
            if not 'auto' input should be a 2x2 numpy array or nested list with row 1 describing the index of the start cell and row 2 describing the start orientation (e.g. [0,1] points right, [1,0] points down)
        stop ('Elisoff'/'Nvisits'/'either') - A string selecting the stopping criterion, default='either'
            'Elisoff' stops when the start cell has been re-visited with the same orientation as started with
            'Nvisits' stops when the start cell has been re-visited N number of times (N set by startvisits)
            'either' stops when either Elisoff or Nvisits has been satisfied
        startvisits - An integer condition for the number of times re-visiting the start cell will trigger an end to the search, default=3
            only applicable when stop='Nvisits' or 'either'
        checkconn (True/False) - A boolean for selecting whether to check connectivity and warn the user of potential issues, default=False
        remcontourrepeat (True/False) - A boolean for selecting whether to remove consecutive repeating cells in the output contour, default=True
        remsearchrepeat (True/False) - A boolean for selecting whether to remove consecutive repeating cells in the output contoursearch, default=False
        closecontour (True/False) - A boolean for selecting whether to close the contour (first cell = last cell), default=True

    Outputs:
        contour - A 2-d Nx2 numpy array of ordered latitude/longitude points (degrees) describing the contour trace of a mask
        contoursearch - A 2-d Nx2 numpy array of ordered latitude/longitude points (degrees) describing the cells searched during contour tracing
            Note: If latitudes/longitudes not provided, returned contour/contoursearch will be indices of the input mask
    """
    if latitudes is not None and longitudes is not None:
        geocontour.check.checkmask(mask,latitudes,longitudes)
    else:
        geocontour.check.checkmask(mask)
    buffermask=numpy.full((tuple(numpy.array(mask.shape)+2)),False)
    buffermask[1:-1,1:-1]=mask
    if checkconn:
        fullconnectivity=geocontour.maskutils.checkconnectivity(buffermask,checkcells='full',connectivity=4)
        emptyconnectivity=geocontour.maskutils.checkconnectivity(buffermask,checkcells='empty',connectivity=4)
        if not fullconnectivity or not emptyconnectivity:
            warnings.warn('WARNING - Mask and/or non-mask not 4-connected: square tracing may not extract the full contour')
    if direction=='cw':
        def outsideturn(orientation):
            return orientation[::-1]*numpy.array([-1,1])
        def insideturn(orientation):
            return orientation[::-1]*numpy.array([1,-1])
        searchdir='ru'
    elif direction=='ccw':
        def outsideturn(orientation):
            return orientation[::-1]*numpy.array([1,-1])
        def insideturn(orientation):
            return orientation[::-1]*numpy.array([-1,1])
        searchdir='rd'
    else:
        sys.exit('ERROR - direction=\''+direction+'\' is not a valid selection, valid selections are \'cw\'/\'ccw\'')
    if type(start) is not str:
        startcell,startorientation=geocontour.contourutil.parsestart(start,buffermask)
    elif start=='auto':
        startcell,startorientation=geocontour.contourutil.findstart(buffermask,searchdir=searchdir)
    else:
        sys.exit('ERROR - start=\''+start+'\' is not a valid selection, valid selections are \'auto\' or a 2x2 array describing start cell index and start orientation')
    checkbreak=geocontour.contourutil.setstop(stop,startvisits,startcell,startorientation)
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
    contour,contoursearch=geocontour.contourutil.clean(contourcells,searchcells,latitudes=latitudes,longitudes=longitudes,closecontour=closecontour,remcontourrepeat=remcontourrepeat,remsearchrepeat=remsearchrepeat)
    return contour,contoursearch

def moore(mask,latitudes=None,longitudes=None,direction='cw',start='auto',stop='either',startvisits=3,remcontourrepeat=True,remsearchrepeat=False,closecontour=True):
    """
    Returns the contour trace of a mask input using the Moore neighbor tracing algorithm
    Source:
        https://www.imageprocessingplace.com/downloads_V3/root_downloads/tutorials/contour_tracing_Abeer_George_Ghuneim/square.html
        http://www-cgrl.cs.mcgill.ca/~godfried/teaching/mir-reading-assignments/Chapter-2-Grids-Connectivity-Contour-Tracing.pdf

    Inputs (Required):
        mask - A 2-d boolean numpy array of dimension MxN where M=len(latitudes) and N=len(longitudes)

    Inputs (Optional):
        latitudes - An evenly spaced numpy array of latitude points (degrees)
        longitudes - An evenly spaced numpy array of longitude points (degrees)
            Note: if latitudes/longitudes not provided, output will be in indices of input mask
        direction ('cw'/'ccw') - A string selecting the direction of the returned neighbors, default='cw'
        start ('auto'/array) - A string or a 2x2 array describind the start cell and orientation, default='auto'
            'auto' will find a proper starting cell/orientation based on the direction input, stop input, and mask shape, or return a warning if unable
            if not 'auto' input should be a 2x2 numpy array or nested list with row 1 describing the index of the start cell and row 2 describing the start orientation (e.g. [0,1] points right, [1,0] points down)
        stop ('Elisoff'/'Nvisits'/'either') - A string selecting the stopping criterion, default='either'
            'Elisoff' stops when the start cell has been re-visited with the same orientation as started with
            'Nvisits' stops when the start cell has been re-visited N number of times (N set by startvisits)
            'either' stops when either Elisoff or Nvisits has been satisfied
        startvisits - An integer condition for the number of times re-visiting the start cell will trigger an end to the search, default=3
            only applicable when stop='Nvisits' or 'either'
        remcontourrepeat (True/False) - A boolean for selecting whether to remove consecutive repeating cells in the output contour, default=True
        remsearchrepeat (True/False) - A boolean for selecting whether to remove consecutive repeating cells in the output contoursearch, default=False
        closecontour (True/False) - A boolean for selecting whether to close the contour (first cell = last cell), default=True

    Outputs:
        contour - A 2-d Nx2 numpy array of ordered latitude/longitude points (degrees) describing the contour trace of a mask
        contoursearch - A 2-d Nx2 numpy array of ordered latitude/longitude points (degrees) describing the cells searched during contour tracing
            Note: If latitudes/longitudes not provided, returned contour/contoursearch will be indices of the input mask
    """
    if latitudes is not None and longitudes is not None:
        geocontour.check.checkmask(mask,latitudes,longitudes)
    else:
        geocontour.check.checkmask(mask)
    buffermask=numpy.full((tuple(numpy.array(mask.shape)+2)),False)
    buffermask[1:-1,1:-1]=mask
    if direction=='cw':
        searchdirstrings=['ru','dr','ld','ul']
        def insideturn(startorientation):
            return startorientation[::-1]*numpy.array([1,-1])
    elif direction=='ccw':
        searchdirstrings=['rd','ur','lu','dl']
        def insideturn(startorientation):
            return startorientation[::-1]*numpy.array([-1,1])
    else:
        sys.exit('ERROR - direction=\''+direction+'\' is not a valid selection, valid selections are \'cw\'/\'ccw\'')
    if type(start) is not str:
        startcell,startorientation=geocontour.contourutil.parsestart(start,buffermask)
    elif start=='auto':
        if stop=='either' or stop=='Elisoff':
            for ct,searchdir in enumerate(searchdirstrings):
                startcell,startorientation=geocontour.contourutil.findstart(buffermask,searchdir=searchdir)
                I=startcell+insideturn(startorientation)
                RI=I-startorientation
                RRI=RI-startorientation
                if buffermask[I[0],I[1]]==True and buffermask[RI[0],RI[1]]==False:
                    C1=True
                else:
                    C1=False
                if any(RRI<0) or RRI[0]>=buffermask.shape[0] or RRI[1]>=buffermask.shape[1]:
                    C2=False
                else:
                    if buffermask[RI[0],RI[1]]==True and buffermask[RRI[0],RRI[1]]==False:
                        C2=True
                    else:
                        C2=False
                if C1 or C2:
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
            startcell,startorientation=geocontour.contourutil.findstart(buffermask,searchdir='ru')
    else:
        sys.exit('ERROR - start=\''+start+'\' is not a valid selection, valid selections are \'auto\' or a 2x2 array describing start cell index and start orientation')
    checkbreak=geocontour.contourutil.setstop(stop,startvisits,startcell,startorientation)
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
            neighbors=geocontour.maskutils.findneighbors(cell,connectivity=8,direction=direction)
            nextneighborindex=(numpy.nonzero(((cell-orientation)==neighbors).all(axis=1))[0][0]+1)%8
        else:
            nextneighborindex=(nextneighborindex+1)%8
        orientation=neighbors[nextneighborindex]-cell
        cell=neighbors[nextneighborindex]
        if (cell==startcell).all():
            Nvisits+=1
        breakloop=checkbreak(cell,orientation,Nvisits)
    contour,contoursearch=geocontour.contourutil.clean(contourcells,searchcells,latitudes=latitudes,longitudes=longitudes,closecontour=closecontour,remcontourrepeat=remcontourrepeat,remsearchrepeat=remsearchrepeat)
    return contour,contoursearch

def moore_imp(mask,latitudes=None,longitudes=None,direction='cw',start='auto',stop='either',startvisits=3,remcontourrepeat=True,remsearchrepeat=False,closecontour=True):
    """
    Returns the contour trace of a mask input using an improved Moore neighbor tracing algorithm
    Captures inside corners missed by Moore neighbor tracing

    Inputs (Required):
        mask - A 2-d boolean numpy array of dimension MxN where M=len(latitudes) and N=len(longitudes)

    Inputs (Optional):
        latitudes - An evenly spaced numpy array of latitude points (degrees)
        longitudes - An evenly spaced numpy array of longitude points (degrees)
            Note: if latitudes/longitudes not provided, output will be in indices of input mask
        direction ('cw'/'ccw') - A string selecting the direction of the returned neighbors, default='cw'
        start ('auto'/array) - A string or a 2x2 array describind the start cell and orientation, default='auto'
            'auto' will find a proper starting cell/orientation based on the direction input, stop input, and mask shape, or return a warning if unable
            if not 'auto' input should be a 2x2 numpy array or nested list with row 1 describing the index of the start cell and row 2 describing the start orientation (e.g. [0,1] points right, [1,0] points down)
        stop ('Elisoff'/'Nvisits'/'either') - A string selecting the stopping criterion, default='either'
            'Elisoff' stops when the start cell has been re-visited with the same orientation as started with
            'Nvisits' stops when the start cell has been re-visited N number of times (N set by startvisits)
            'either' stops when either Elisoff or Nvisits has been satisfied
        startvisits - An integer condition for the number of times re-visiting the start cell will trigger an end to the search, default=3
            only applicable when stop='Nvisits' or 'either'
        remcontourrepeat (True/False) - A boolean for selecting whether to remove consecutive repeating cells in the output contour, default=True
        remsearchrepeat (True/False) - A boolean for selecting whether to remove consecutive repeating cells in the output contoursearch, default=False
        closecontour (True/False) - A boolean for selecting whether to close the contour (first cell = last cell), default=True

    Outputs:
        contour - A 2-d Nx2 numpy array of ordered latitude/longitude points (degrees) describing the contour trace of a mask
        contoursearch - A 2-d Nx2 numpy array of ordered latitude/longitude points (degrees) describing the cells searched during contour tracing
            Note: If latitudes/longitudes not provided, returned contour/contoursearch will be indices of the input mask
    """
    if latitudes is not None and longitudes is not None:
        geocontour.check.checkmask(mask,latitudes,longitudes)
    else:
        geocontour.check.checkmask(mask)
    buffermask=numpy.full((tuple(numpy.array(mask.shape)+2)),False)
    buffermask[1:-1,1:-1]=mask
    if direction=='cw':
        searchdirstrings=['ru','dr','ld','ul']
        def insideturn(startorientation):
            return startorientation[::-1]*numpy.array([1,-1])
    elif direction=='ccw':
        searchdirstrings=['rd','ur','lu','dl']
        def insideturn(startorientation):
            return startorientation[::-1]*numpy.array([-1,1])
    else:
        sys.exit('ERROR - direction=\''+direction+'\' is not a valid selection, valid selections are \'cw\'/\'ccw\'')
    if type(start) is not str:
        startcell,startorientation=geocontour.contourutil.parsestart(start,buffermask)
    elif start=='auto':
        if stop=='either' or stop=='Elisoff':
            for ct,searchdir in enumerate(searchdirstrings):
                startcell,startorientation=geocontour.contourutil.findstart(buffermask,searchdir=searchdir)
                I=startcell+insideturn(startorientation)
                RI=I-startorientation
                RRI=RI-startorientation
                if buffermask[I[0],I[1]]==True and buffermask[RI[0],RI[1]]==False:
                    C1=True
                else:
                    C1=False
                if any(RRI<0) or RRI[0]>=buffermask.shape[0] or RRI[1]>=buffermask.shape[1]:
                    C2=False
                else:
                    if buffermask[RI[0],RI[1]]==True and buffermask[RRI[0],RRI[1]]==False:
                        C2=True
                    else:
                        C2=False
                if C1 or C2:
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
            startcell,startorientation=geocontour.contourutil.findstart(buffermask,searchdir='ru')
    else:
        sys.exit('ERROR - start=\''+start+'\' is not a valid selection, valid selections are \'auto\' or a 2x2 array describing start cell index and start orientation')
    checkbreak=geocontour.contourutil.setstop(stop,startvisits,startcell,startorientation)
    cell=startcell
    orientation=startorientation
    searchcells=[]
    contourcells=[]
    searchcells.append(cell)
    contourcells.append(cell)
    neighbors=geocontour.maskutils.findneighbors(cell,connectivity=8,direction=direction)
    nextneighborindex=(numpy.nonzero(((cell-orientation)==neighbors).all(axis=1))[0][0]+1)%8
    Nvisits=0
    breakloop=False
    while not breakloop:
        orientation=neighbors[nextneighborindex]-cell
        cell=neighbors[nextneighborindex]
        searchcells.append(cell)
        if buffermask[cell[0],cell[1]]==True:
            nextnextneighbor=neighbors[(nextneighborindex+1)%8]
            nextnextinmask=(buffermask[nextnextneighbor[0],nextnextneighbor[1]]==True)
            normdist=numpy.linalg.norm(contourcells[-1]-cell)
            if nextnextinmask and normdist==numpy.sqrt(2):
                searchcells.insert(-1,nextnextneighbor)
                contourcells.append(nextnextneighbor)
            contourcells.append(cell)
            neighbors=geocontour.maskutils.findneighbors(cell,connectivity=8,direction=direction)
            nextneighborindex=(numpy.nonzero(((cell-orientation)==neighbors).all(axis=1))[0][0]+1)%8
        else:
            nextneighborindex=(nextneighborindex+1)%8
        if (cell==startcell).all():
            Nvisits+=1
        breakloop=checkbreak(cell,orientation,Nvisits)
    contour,contoursearch=geocontour.contourutil.clean(contourcells,searchcells,latitudes=latitudes,longitudes=longitudes,closecontour=closecontour,remcontourrepeat=remcontourrepeat,remsearchrepeat=remsearchrepeat)
    return contour,contoursearch

def pavlidis(mask,latitudes=None,longitudes=None,direction='cw',start='auto',stop='Nvisits',startvisits=1,remcontourrepeat=True,remsearchrepeat=False,closecontour=True):
    """
    Returns the contour trace of a mask input using the Pavlidis tracing algorithm
    Source:
        10.1007/978-3-642-93208-3
        https://www.imageprocessingplace.com/downloads_V3/root_downloads/tutorials/contour_tracing_Abeer_George_Ghuneim/theo.html

    Inputs (Required):
        mask - A 2-d boolean numpy array of dimension MxN where M=len(latitudes) and N=len(longitudes)

    Inputs (Optional):
        latitudes - An evenly spaced numpy array of latitude points (degrees)
        longitudes - An evenly spaced numpy array of longitude points (degrees)
            Note: if latitudes/longitudes not provided, output will be in indices of input mask
        direction ('cw'/'ccw') - A string selecting the direction of the returned neighbors, default='cw'
        start ('auto'/array) - A string or a 2x2 array describind the start cell and orientation, default='auto'
            'auto' will find a proper starting cell/orientation based on the direction input, stop input, and mask shape, or return a warning if unable
            if not 'auto' input should be a 2x2 numpy array or nested list with row 1 describing the index of the start cell and row 2 describing the start orientation (e.g. [0,1] points right, [1,0] points down)
        stop ('Elisoff'/'Nvisits'/'either') - A string selecting the stopping criterion, default='Nvisits'
            Note: The Elisoff stopping criterion is usable in this function, but the Pavlidis tracing algorithm is not designed to utilize it and if used it is likely to produce an infinite loop
            'Elisoff' stops when the start cell has been re-visited with the same orientation as started with
            'Nvisits' stops when the start cell has been re-visited N number of times (N set by startvisits)
            'either' stops when either Elisoff or Nvisits has been satisfied
        startvisits - An integer condition for the number of times re-visiting the start cell will trigger an end to the search, default=1
            only applicable when stop='Nvisits' or 'either'
        remcontourrepeat (True/False) - A boolean for selecting whether to remove consecutive repeating cells in the output contour, default=True
        remsearchrepeat (True/False) - A boolean for selecting whether to remove consecutive repeating cells in the output contoursearch, default=False
        closecontour (True/False) - A boolean for selecting whether to close the contour (first cell = last cell), default=True

    Outputs:
        contour - A 2-d Nx2 numpy array of ordered latitude/longitude points (degrees) describing the contour trace of a mask
        contoursearch - A 2-d Nx2 numpy array of ordered latitude/longitude points (degrees) describing the cells searched during contour tracing
            Note: If latitudes/longitudes not provided, returned contour/contoursearch will be indices of the input mask
    """
    if latitudes is not None and longitudes is not None:
        geocontour.check.checkmask(mask,latitudes,longitudes)
    else:
        geocontour.check.checkmask(mask)
    buffermask=numpy.full((tuple(numpy.array(mask.shape)+2)),False)
    buffermask[1:-1,1:-1]=mask
    if direction=='cw':
        def outsideturn(orientation):
            return orientation[::-1]*numpy.array([-1,1])
        def insideturn(orientation):
            return orientation[::-1]*numpy.array([1,-1])
        searchdirstrings=['ru','dr','ld','ul']
    elif direction=='ccw':
        def outsideturn(orientation):
            return orientation[::-1]*numpy.array([1,-1])
        def insideturn(orientation):
            return orientation[::-1]*numpy.array([-1,1])
        searchdirstrings=['rd','ur','lu','dl']
    else:
        sys.exit('ERROR - direction=\''+direction+'\' is not a valid selection, valid selections are \'cw\'/\'ccw\'')
    if type(start) is not str:
        startcell,startorientation=geocontour.contourutil.parsestart(start,buffermask)
    elif start=='auto':
        for ct,searchdir in enumerate(searchdirstrings):
            startcell,startorientation=geocontour.contourutil.findstart(buffermask,searchdir=searchdir)
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
    checkbreak=geocontour.contourutil.setstop(stop,startvisits,startcell,startorientation)
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
            if (cell==startcell).all():
                Nvisits+=1
            rotations=0
            breakloop=checkbreak(cell,orientation,Nvisits)
            orientation=outsideturn(orientation)
        elif buffermask[F[0],F[1]]==True:
            searchcells.append(FO)
            searchcells.append(F)
            contourcells.append(F)
            cell=F
            if (cell==startcell).all():
                Nvisits+=1
            rotations=0
            breakloop=checkbreak(cell,orientation,Nvisits)
        elif buffermask[FI[0],FI[1]]==True:
            searchcells.append(FO)
            searchcells.append(F)
            searchcells.append(FI)
            contourcells.append(FI)
            cell=FI
            if (cell==startcell).all():
                Nvisits+=1
            rotations=0
            breakloop=checkbreak(cell,orientation,Nvisits)
        else:
            searchcells.append(FO)
            searchcells.append(F)
            searchcells.append(FI)
            orientation=insideturn(orientation)
            rotations+=1
        if rotations>3:
            sys.exit('ERROR - Stuck on isolated cell '+str(cell))
    contour,contoursearch=geocontour.contourutil.clean(contourcells,searchcells,latitudes=latitudes,longitudes=longitudes,closecontour=closecontour,remcontourrepeat=remcontourrepeat,remsearchrepeat=remsearchrepeat)
    return contour,contoursearch

def pavlidis_imp(mask,latitudes=None,longitudes=None,direction='cw',start='auto',stop='Nvisits',startvisits=1,remcontourrepeat=True,remsearchrepeat=False,closecontour=True):
    """
    Returns the contour trace of a mask input using an improved Pavlidis tracing algorithm
    Captures inside corners missed by Pavlidis tracing

    Inputs (Required):
        mask - A 2-d boolean numpy array of dimension MxN where M=len(latitudes) and N=len(longitudes)

    Inputs (Optional):
        latitudes - An evenly spaced numpy array of latitude points (degrees)
        longitudes - An evenly spaced numpy array of longitude points (degrees)
            Note: if latitudes/longitudes not provided, output will be in indices of input mask
        direction ('cw'/'ccw') - A string selecting the direction of the returned neighbors, default='cw'
        start ('auto'/array) - A string or a 2x2 array describind the start cell and orientation, default='auto'
            'auto' will find a proper starting cell/orientation based on the direction input, stop input, and mask shape, or return a warning if unable
            if not 'auto' input should be a 2x2 numpy array or nested list with row 1 describing the index of the start cell and row 2 describing the start orientation (e.g. [0,1] points right, [1,0] points down)
        stop ('Elisoff'/'Nvisits'/'either') - A string selecting the stopping criterion, default='Nvisits'
            Note: The Elisoff stopping criterion is usable in this function, but the Pavlidis tracing algorithm is not designed to utilize it and if used it is likely to produce an infinite loop
            'Elisoff' stops when the start cell has been re-visited with the same orientation as started with
            'Nvisits' stops when the start cell has been re-visited N number of times (N set by startvisits)
            'either' stops when either Elisoff or Nvisits has been satisfied
        startvisits - An integer condition for the number of times re-visiting the start cell will trigger an end to the search, default=1
            only applicable when stop='Nvisits' or 'either'
        remcontourrepeat (True/False) - A boolean for selecting whether to remove consecutive repeating cells in the output contour, default=True
        remsearchrepeat (True/False) - A boolean for selecting whether to remove consecutive repeating cells in the output contoursearch, default=False
        closecontour (True/False) - A boolean for selecting whether to close the contour (first cell = last cell), default=True

    Outputs:
        contour - A 2-d Nx2 numpy array of ordered latitude/longitude points (degrees) describing the contour trace of a mask
        contoursearch - A 2-d Nx2 numpy array of ordered latitude/longitude points (degrees) describing the cells searched during contour tracing
            Note: If latitudes/longitudes not provided, returned contour/contoursearch will be indices of the input mask
    """
    if latitudes is not None and longitudes is not None:
        geocontour.check.checkmask(mask,latitudes,longitudes)
    else:
        geocontour.check.checkmask(mask)
    buffermask=numpy.full((tuple(numpy.array(mask.shape)+2)),False)
    buffermask[1:-1,1:-1]=mask
    if direction=='cw':
        def outsideturn(orientation):
            return orientation[::-1]*numpy.array([-1,1])
        def insideturn(orientation):
            return orientation[::-1]*numpy.array([1,-1])
        searchdirstrings=['ru','dr','ld','ul']
    elif direction=='ccw':
        def outsideturn(orientation):
            return orientation[::-1]*numpy.array([1,-1])
        def insideturn(orientation):
            return orientation[::-1]*numpy.array([-1,1])
        searchdirstrings=['rd','ur','lu','dl']
    else:
        sys.exit('ERROR - direction=\''+direction+'\' is not a valid selection, valid selections are \'cw\'/\'ccw\'')
    if type(start) is not str:
        startcell,startorientation=geocontour.contourutil.parsestart(start,buffermask)
    elif start=='auto':
        for ct,searchdir in enumerate(searchdirstrings):
            startcell,startorientation=geocontour.contourutil.findstart(buffermask,searchdir=searchdir)
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
    checkbreak=geocontour.contourutil.setstop(stop,startvisits,startcell,startorientation)
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
            if (cell==startcell).all():
                Nvisits+=1
            rotations=0
            breakloop=checkbreak(cell,orientation,Nvisits)
            orientation=outsideturn(orientation)
        elif buffermask[F[0],F[1]]==True:
            searchcells.append(FO)
            searchcells.append(F)
            contourcells.append(F)
            cell=F
            if (cell==startcell).all():
                Nvisits+=1
            rotations=0
            breakloop=checkbreak(cell,orientation,Nvisits)
        elif buffermask[FI[0],FI[1]]==True:
            if buffermask[I[0],I[1]]==True:
                contourcells.append(I)
            searchcells.append(FO)
            searchcells.append(F)
            searchcells.append(I)
            searchcells.append(FI)
            contourcells.append(FI)
            cell=FI
            if (cell==startcell).all():
                Nvisits+=1
            rotations=0
            breakloop=checkbreak(cell,orientation,Nvisits)
        else:
            searchcells.append(FO)
            searchcells.append(F)
            searchcells.append(FI)
            orientation=insideturn(orientation)
            rotations+=1
        if rotations>3:
            sys.exit('ERROR - Stuck on isolated cell '+str(cell))
    contour,contoursearch=geocontour.contourutil.clean(contourcells,searchcells,latitudes=latitudes,longitudes=longitudes,closecontour=closecontour,remcontourrepeat=remcontourrepeat,remsearchrepeat=remsearchrepeat)
    return contour,contoursearch

def TSR(mask,latitudes=None,longitudes=None,direction='cw',start='auto',stop='either',startvisits=4,remcontourrepeat=True,remsearchrepeat=False,closecontour=True):
    """
    Returns the contour trace of a mask input using fast representative tracing
    Source: 10.3390/s16030353

    Inputs (Required):
        mask - A 2-d boolean numpy array of dimension MxN where M=len(latitudes) and N=len(longitudes)

    Inputs (Optional):
        latitudes - An evenly spaced numpy array of latitude points (degrees)
        longitudes - An evenly spaced numpy array of longitude points (degrees)
            Note: if latitudes/longitudes not provided, output will be in indices of input mask
        direction ('cw'/'ccw') - A string selecting the direction of the returned neighbors, default='cw'
        start ('auto'/array) - A string or a 2x2 array describind the start cell and orientation, default='auto'
            'auto' will find a proper starting cell/orientation based on the direction input
            if not 'auto' input should be a 2x2 numpy array or nested list with row 1 describing the index of the start cell and row 2 describing the start orientation (e.g. [0,1] points right, [1,0] points down)
        stop ('Elisoff'/'Nvisits'/'either') - A string selecting the stopping criterion, default='either'
            'Elisoff' stops when the start cell has been re-visited with the same orientation as started with
            'Nvisits' stops when the start cell has been re-visited N number of times (N set by startvisits)
            'either' stops when either Elisoff or Nvisits has been satisfied
        startvisits - An integer condition for the number of times re-visiting the start cell will trigger an end to the search, default=4
            only applicable when stop='Nvisits' or 'either'
        checkconn (True/False) - A boolean for selecting whether to check connectivity and warn the user of potential issues, default=False
        remcontourrepeat (True/False) - A boolean for selecting whether to remove consecutive repeating cells in the output contour, default=True
        remsearchrepeat (True/False) - A boolean for selecting whether to remove consecutive repeating cells in the output contoursearch, default=False
        closecontour (True/False) - A boolean for selecting whether to close the contour (first cell = last cell), default=True

    Outputs:
        contour - A 2-d Nx2 numpy array of ordered latitude/longitude points (degrees) describing the contour trace of a mask
        contoursearch - A 2-d Nx2 numpy array of ordered latitude/longitude points (degrees) describing the cells searched during contour tracing
            Note: If latitudes/longitudes not provided, returned contour/contoursearch will be indices of the input mask
    """
    if latitudes is not None and longitudes is not None:
        geocontour.check.checkmask(mask,latitudes,longitudes)
    else:
        geocontour.check.checkmask(mask)
    buffermask=numpy.full((tuple(numpy.array(mask.shape)+2)),False)
    buffermask[1:-1,1:-1]=mask
    if direction=='cw':
        def outsideturn(orientation):
            return orientation[::-1]*numpy.array([-1,1])
        def insideturn(orientation):
            return orientation[::-1]*numpy.array([1,-1])
        searchdir='ru'
    elif direction=='ccw':
        def outsideturn(orientation):
            return orientation[::-1]*numpy.array([1,-1])
        def insideturn(orientation):
            return orientation[::-1]*numpy.array([-1,1])
        searchdir='rd'
    else:
        sys.exit('ERROR - direction=\''+direction+'\' is not a valid selection, valid selections are \'cw\'/\'ccw\'')
    if type(start) is not str:
        startcell,startorientation=geocontour.contourutil.parsestart(start,buffermask)
    elif start=='auto':
        startcell,startorientation=geocontour.contourutil.findstart(buffermask,searchdir=searchdir)
    else:
        sys.exit('ERROR - start=\''+start+'\' is not a valid selection, valid selections are \'auto\' or a 2x2 array describing start cell index and start orientation')
    checkbreak=geocontour.contourutil.setstop(stop,startvisits,startcell,startorientation)
    orientation=startorientation
    cell=startcell
    searchcells=[]
    contourcells=[]
    searchcells.append(cell)
    contourcells.append(cell)
    Nvisits=0
    breakloop=False
    while not breakloop:
        L=cell+outsideturn(orientation)
        RL=cell-orientation+outsideturn(orientation)
        #Stage 1
        if buffermask[RL[0],RL[1]]==True:
            searchcells.append(L)
            searchcells.append(RL)
            if buffermask[L[0],L[1]]==True: #Case 1
                contourcells.append(L)
                contourcells.append(RL)
            else: #Case 2
                contourcells.append(RL)
            cell=RL
            orientation=-orientation
        else:
            searchcells.append(RL)
            searchcells.append(L)
            if buffermask[L[0],L[1]]==True: #Case 3
                contourcells.append(L)
                cell=L
                orientation=outsideturn(orientation)
            else: #Case 4
                pass
        #Stage 2
        F=cell+orientation
        FL=cell+orientation+outsideturn(orientation)
        if buffermask[FL[0],FL[1]]==True:
            searchcells.append(F)
            searchcells.append(FL)
            if buffermask[F[0],F[1]]==True: #Case 6
                contourcells.append(F)
                contourcells.append(FL)
            else: #Case 5
                contourcells.append(FL)
            cell=FL
        else:
            searchcells.append(FL)
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
    contour,contoursearch=geocontour.contourutil.clean(contourcells,searchcells,latitudes=latitudes,longitudes=longitudes,closecontour=closecontour,remcontourrepeat=remcontourrepeat,remsearchrepeat=remsearchrepeat)
    return contour,contoursearch

