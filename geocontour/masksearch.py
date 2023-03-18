"""
Functions for finding a mask on a lat/lon grid from an input boundary
=====================================================================
"""
import warnings
import numpy as np
import shapely as sh
import matplotlib.path as mplp
import geocontour.grid as gcg
import geocontour.maskutil as gcmu

def center(latitudes,longitudes,boundary,precision=1e-5):
    """
    Find a mask on a lat/lon grid from an input boundary
        
    Criteria for inclusion of a cell is whether the center of the cell
    falls within the boundary

    Parameters
    ----------
    latitudes : ndarray
        1D Nx1 array of latitude points (degrees)
    longitudes : ndarray
        1D Nx1 array of longitude points (degrees)
    boundary : ndarray
        2D Nx2 array of latitude/longitude points (degrees) with the
        last point equal to the first
    precision : float, default=1e-5
        value by which `boundary` is expanded, capturing indeterminate
        points that may fall on/near `boundary`

    Returns
    -------
    mask : ndarray
        2D MxN bool array where M=len(`latitudes`) and
        N=len(`longitudes`)

    See Also
    --------
    center2
    nodes
    nodes2
    area

    Notes
    -----
    Regarding `precision`: Shapely can't beat machine precision, and can
    thus give "incorrect" results for very close points or shapes. The
    `precision` input errs on being more inclusive, in particular to
    capture points falling directly on a boundary. A decent rule is to
    set the `precision` value as high as you can without impeding the
    accuracy. For instance, the default of 1e-5 (degrees) translates to
    roughly 1m precision at the equator. The buffer can be negated by
    setting this input very low (to machine precision).
    """
    latdir=gcg.clatdir(latitudes)
    if latdir=='dec':
        latitudes=np.flip(latitudes)
    boxlatmin, boxlatmax, boxlonmin, boxlonmax = gcmu.bbox(latitudes,longitudes,boundary)
    boxmask=np.full((boxlatmax-boxlatmin+1,boxlonmax-boxlonmin+1),False)
    boundpoly=sh.polygons(boundary).buffer(precision)
    ysp, xsp = np.meshgrid(latitudes[boxlatmin:boxlatmax+1],longitudes[boxlonmin:boxlonmax+1], indexing='ij')
    searchpoints=np.hstack((ysp.reshape((-1,1)), xsp.reshape((-1,1))))
    shsearchpoints=sh.points(searchpoints)
    boxmask=sh.contains(boundpoly,shsearchpoints)
    mask=np.full((len(latitudes),len(longitudes)),False)
    mask[boxlatmin:boxlatmax+1,boxlonmin:boxlonmax+1]=boxmask.reshape((boxlatmax-boxlatmin+1,boxlonmax-boxlonmin+1))
    if latdir=='dec':
        mask=np.flip(mask,axis=0)
    return mask

def center2(latitudes,longitudes,boundary,precision=1e-5):
    """
    Find a mask on a lat/lon grid from an input boundary
        
    Criteria for inclusion of a cell is whether the center of the cell
    falls within the boundary

    Functionally identical to `center()`, but utilizes matplotlib.path
    functions, which are faster (possibly due to avoidance of overhead
    in converting to shapely geometries)

    Parameters
    ----------
    latitudes : ndarray
        1D Nx1 array of latitude points (degrees)
    longitudes : ndarray
        1D Nx1 array of longitude points (degrees)
    boundary : ndarray
        2D Nx2 array of latitude/longitude points (degrees) with the
        last point equal to the first
    precision : float, default=1e-5
        value by which `boundary` is expanded, capturing indeterminate
        points that may fall on/near `boundary`

    Returns
    -------
    mask : ndarray
        2D MxN bool array where M=len(`latitudes`) and
        N=len(`longitudes`)

    See Also
    --------
    center
    nodes
    nodes2
    area

    Notes
    -----
    Regarding `precision`: Matplotlib.Path can't beat machine precision,
    and can thus give "incorrect" results for very close points or
    shapes. The `precision` input errs on being more inclusive, in
    particular to capture points falling directly on a boundary. A
    decent rule is to set the `precision` value as high as you can
    without impeding the accuracy. For instance, the default of 1e-5
    (degrees) translates to roughly 1m precision at the equator. The
    buffer can be negated by setting this input very low (to machine
    precision).
    """
    latdir=gcg.clatdir(latitudes)
    if latdir=='dec':
        latitudes=np.flip(latitudes)
    boxlatmin, boxlatmax, boxlonmin, boxlonmax = gcmu.bbox(latitudes,longitudes,boundary)
    boundpoly=mplp.Path(boundary)
    ysp, xsp = np.meshgrid(latitudes[boxlatmin:boxlatmax+1],longitudes[boxlonmin:boxlonmax+1], indexing='ij')
    searchpoints=np.hstack((ysp.reshape((-1,1)), xsp.reshape((-1,1))))
    boxmask=boundpoly.contains_points(searchpoints,radius=2*precision)
    mask=np.full((len(latitudes),len(longitudes)),False)
    mask[boxlatmin:boxlatmax+1,boxlonmin:boxlonmax+1]=boxmask.reshape((boxlatmax-boxlatmin+1,boxlonmax-boxlonmin+1))
    if latdir=='dec':
        mask=np.flip(mask,axis=0)
    return mask

def nodes(latitudes,longitudes,boundary,nodes=2,precision=1e-5):
    """
    Find a mask on a lat/lon grid from an input boundary
        
    Criteria for inclusion of a cell is whether a given number of cell
    nodes (corners) fall within the boundary 

    Parameters
    ----------
    latitudes : ndarray
        1D Nx1 array of latitude points (degrees)
    longitudes : ndarray
        1D Nx1 array of longitude points (degrees)
    boundary : ndarray
        2D Nx2 array of latitude/longitude points (degrees) with the
        last point equal to the first
    nodes : int, default=2
        number of cell nodes (corners) to use as a criteria for
        inclusion (1-4)
    precision : float, default=1e-5
        value by which `boundary` is expanded, capturing indeterminate
        points that may fall on/near `boundary`

    Returns
    -------
    mask : ndarray
        2D MxN bool array where M=len(`latitudes`) and
        N=len(`longitudes`)

    See Also
    --------
    center
    center2
    nodes2
    area

    Notes
    -----
    Regarding `precision`: Shapely can't beat machine precision, and can
    thus give "incorrect" results for very close points or shapes. The
    `precision` input errs on being more inclusive, in particular to
    capture points falling directly on a boundary. A decent rule is to
    set the `precision` value as high as you can without impeding the
    accuracy. For instance, the default of 1e-5 (degrees) translates to
    roughly 1m precision at the equator. The buffer can be negated by
    setting this input very low (to machine precision).
    """
    if nodes<1:
        warnings.warn('WARNING - valid input for nodes is 1-4, nodes<1 will result in all cells being selected')
    if nodes>4:
        warnings.warn('WARNING - valid input for nodes is 1-4, nodes>4 will result in no cells being selected')
    latdir=gcg.clatdir(latitudes)
    if latdir=='dec':
        latitudes=np.flip(latitudes)
    boxlatmin, boxlatmax, boxlonmin, boxlonmax = gcmu.bbox(latitudes,longitudes,boundary)
    latgrdspc=gcg.spacing(latitudes)
    longrdspc=gcg.spacing(longitudes)
    boxmask=np.full((boxlatmax-boxlatmin+1,boxlonmax-boxlonmin+1),False)
    boundpoly=sh.polygons(boundary).buffer(precision)
    yspMM, xspMM = np.meshgrid(latitudes[boxlatmin:boxlatmax+1]-latgrdspc/2,longitudes[boxlonmin:boxlonmax+1]-longrdspc/2, indexing='ij')
    yspPM, xspPM = np.meshgrid(latitudes[boxlatmin:boxlatmax+1]+latgrdspc/2,longitudes[boxlonmin:boxlonmax+1]-longrdspc/2, indexing='ij')
    yspMP, xspMP = np.meshgrid(latitudes[boxlatmin:boxlatmax+1]-latgrdspc/2,longitudes[boxlonmin:boxlonmax+1]+longrdspc/2, indexing='ij')
    yspPP, xspPP = np.meshgrid(latitudes[boxlatmin:boxlatmax+1]+latgrdspc/2,longitudes[boxlonmin:boxlonmax+1]+longrdspc/2, indexing='ij')
    spMM=np.hstack((yspMM.reshape((-1,1)),xspMM.reshape((-1,1))))
    spPM=np.hstack((yspPM.reshape((-1,1)),xspPM.reshape((-1,1))))
    spMP=np.hstack((yspMP.reshape((-1,1)),xspMP.reshape((-1,1))))
    spPP=np.hstack((yspPP.reshape((-1,1)),xspPP.reshape((-1,1))))
    searchpoints=np.stack((spMM,spPM,spMP,spPP)).reshape(-1,2)
    shsearchpoints=sh.points(searchpoints)
    nodesinmask=sh.contains(boundpoly,shsearchpoints).reshape(4,-1)
    boxmask=(nodesinmask.sum(axis=0)>=nodes).reshape(boxmask.shape)
    mask=np.full((len(latitudes),len(longitudes)),False)
    mask[boxlatmin:boxlatmax+1,boxlonmin:boxlonmax+1]=boxmask
    if latdir=='dec':
        mask=np.flip(mask,axis=0)
    return mask

def nodes2(latitudes,longitudes,boundary,nodes=2,precision=1e-5):
    """
    Find a mask on a lat/lon grid from an input boundary
        
    Criteria for inclusion of a cell is whether a given number of cell
    nodes (corners) fall within the boundary 

    Functionally identical to `nodes()`, but utilizes matplotlib.path
    functions, which are faster (possibly due to avoidance of overhead
    in converting to shapely geometries)

    Parameters
    ----------
    latitudes : ndarray
        1D Nx1 array of latitude points (degrees)
    longitudes : ndarray
        1D Nx1 array of longitude points (degrees)
    boundary : ndarray
        2D Nx2 array of latitude/longitude points (degrees) with the
        last point equal to the first
    nodes : int, default=2
        number of cell nodes (corners) to use as a criteria for
        inclusion (1-4)
    precision : float, default=1e-5
        value by which `boundary` is expanded, capturing indeterminate
        points that may fall on/near `boundary`

    Returns
    -------
    mask : ndarray
        2D MxN bool array where M=len(`latitudes`) and
        N=len(`longitudes`)

    See Also
    --------
    center
    center2
    nodes
    area

    Notes
    -----
    Regarding `precision`: Matplotlib.Path can't beat machine precision,
    and can thus give "incorrect" results for very close points or
    shapes. The `precision` input errs on being more inclusive, in
    particular to capture points falling directly on a boundary. A
    decent rule is to set the `precision` value as high as you can
    without impeding the accuracy. For instance, the default of 1e-5
    (degrees) translates to roughly 1m precision at the equator. The
    buffer can be negated by setting this input very low (to machine
    precision).
    """
    if nodes<1:
        warnings.warn('WARNING - valid input for nodes is 1-4, nodes<1 will result in all cells being selected')
    if nodes>4:
        warnings.warn('WARNING - valid input for nodes is 1-4, nodes>4 will result in no cells being selected')
    latdir=gcg.clatdir(latitudes)
    if latdir=='dec':
        latitudes=np.flip(latitudes)
    boxlatmin, boxlatmax, boxlonmin, boxlonmax = gcmu.bbox(latitudes,longitudes,boundary)
    latgrdspc=gcg.spacing(latitudes)
    longrdspc=gcg.spacing(longitudes)
    boxmask=np.full((boxlatmax-boxlatmin+1,boxlonmax-boxlonmin+1),False)
    boundpoly=mplp.Path(boundary)
    yspMM, xspMM = np.meshgrid(latitudes[boxlatmin:boxlatmax+1]-latgrdspc/2,longitudes[boxlonmin:boxlonmax+1]-longrdspc/2, indexing='ij')
    yspPM, xspPM = np.meshgrid(latitudes[boxlatmin:boxlatmax+1]+latgrdspc/2,longitudes[boxlonmin:boxlonmax+1]-longrdspc/2, indexing='ij')
    yspMP, xspMP = np.meshgrid(latitudes[boxlatmin:boxlatmax+1]-latgrdspc/2,longitudes[boxlonmin:boxlonmax+1]+longrdspc/2, indexing='ij')
    yspPP, xspPP = np.meshgrid(latitudes[boxlatmin:boxlatmax+1]+latgrdspc/2,longitudes[boxlonmin:boxlonmax+1]+longrdspc/2, indexing='ij')
    spMM=np.hstack((yspMM.reshape((-1,1)),xspMM.reshape((-1,1))))
    spPM=np.hstack((yspPM.reshape((-1,1)),xspPM.reshape((-1,1))))
    spMP=np.hstack((yspMP.reshape((-1,1)),xspMP.reshape((-1,1))))
    spPP=np.hstack((yspPP.reshape((-1,1)),xspPP.reshape((-1,1))))
    searchpoints=np.stack((spMM,spPM,spMP,spPP)).reshape(-1,2)
    nodesinmask=boundpoly.contains_points(searchpoints,radius=2*precision).reshape(4,-1)
    boxmask=(nodesinmask.sum(axis=0)>=nodes).reshape(boxmask.shape)
    mask=np.full((len(latitudes),len(longitudes)),False)
    mask[boxlatmin:boxlatmax+1,boxlonmin:boxlonmax+1]=boxmask
    if latdir=='dec':
        mask=np.flip(mask,axis=0)
    return mask

def area(latitudes,longitudes,boundary,area=0.5):
    """
    Find a mask on a lat/lon grid from an input boundary
        
    Criteria for inclusion of a cell is whether the area of the cell
    enclosed by the boundary is greater than some fraction between 0 and
    1

    Parameters
    ----------
    latitudes : ndarray
        1D Nx1 array of latitude points (degrees)
    longitudes : ndarray
        1D Nx1 array of longitude points (degrees)
    boundary : ndarray
        2D Nx2 array of latitude/longitude points (degrees) with the
        last point equal to the first
    area : float, default=0.5
        the fraction of the cell enclosed by the boundary to use as a
        criteria for inclusion (0-1)

    Returns
    -------
    mask : ndarray
        2D MxN bool array where M=len(`latitudes`) and
        N=len(`longitudes`)

    See Also
    --------
    center
    center2
    nodes
    nodes2
    """
    if area>1:
        warnings.warn('WARNING - valid input for area is 0-1, area>1 will result in no cells being selected')
    if area<=0:
        warnings.warn('WARNING - valid input for area is 0-1, area<=0 will result in all cells being selected')
    latdir=gcg.clatdir(latitudes)
    if latdir=='dec':
        latitudes=np.flip(latitudes)
    boxlatmin, boxlatmax, boxlonmin, boxlonmax = gcmu.bbox(latitudes,longitudes,boundary)
    latgrdspc=gcg.spacing(latitudes)
    longrdspc=gcg.spacing(longitudes)
    boxmask=np.full((boxlatmax-boxlatmin+1,boxlonmax-boxlonmin+1),False)
    boundpoly=sh.polygons(boundary)
    for la in np.arange(boxlatmin,boxlatmax+1,1):
        for lo in np.arange(boxlonmin,boxlonmax+1,1):
            LL=[latitudes[la]-latgrdspc/2,longitudes[lo]-longrdspc/2]
            HL=[latitudes[la]+latgrdspc/2,longitudes[lo]-longrdspc/2]
            LH=[latitudes[la]-latgrdspc/2,longitudes[lo]+longrdspc/2]
            HH=[latitudes[la]+latgrdspc/2,longitudes[lo]+longrdspc/2]
            cell=sh.polygons([LL,HL,HH,LH])
            if boundpoly.intersection(cell).area>=latgrdspc*longrdspc*area:
                boxmask[la-boxlatmin,lo-boxlonmin]=True
    mask=np.full((len(latitudes),len(longitudes)),False)
    mask[boxlatmin:boxlatmax+1,boxlonmin:boxlonmax+1]=boxmask
    if latdir=='dec':
        mask=np.flip(mask,axis=0)
    return mask

