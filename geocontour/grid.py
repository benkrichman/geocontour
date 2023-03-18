"""
Functions for operations involving grid values and dimensions
=============================================================
"""
import sys
import numpy as np
import geocontour.check as gcc

def spacing(dimension):
    """
    Calculate the grid spacing for a given input dimension

    Parameters
    ----------
    dimension : ndarray
        1D Nx1 array of longitude or latitude points (degrees)

    Returns
    -------
    spacing : float 
        value specifying the spacing of the input `dimension` (degrees)
    """
    gcc.cdim(dimension)
    spacing=abs(np.diff(dimension)[0])
    return spacing

def lonlens(latitudes,lonspacing=1):
    """
    Calculate the lengths of a degree (default) of longitude over a
    range of latitudes

    Parameters
    ----------
    latitudes : ndarray
        1D Nx1 array of latitude points (degrees)
    lonspacing : float, default=1
        value specifying longitude spacing
        [e.g. for lengths of a half degree of longitude at the given
        `latitudes`, `lonspacing` = 0.5]

    Returns
    -------
    lonlens : ndarray
        1D Nx1 array of longitude lengths (m)

    See Also
    --------
    latlens
    lonlen
    latlen

    References
    ----------
    Osborne, P. (2013). The Mercator Projections. Zenodo.
    <https://doi.org/10.5281/zenodo.35392>
    """
    gcc.cdim(latitudes)
    a=6378137
    b=6356752.3142
    e=np.sqrt((a**2-b**2)/a**2)
    latrads=latitudes*np.pi/180
    lonlengths=lonspacing*np.pi*a*np.cos(latrads)/(180*np.sqrt(1-(e*np.sin(latrads))**2))
    return lonlengths

def latlens(latitudes):
    """
    Calculate the lengths of a defined range of latitudes

    Parameters
    ----------
    latitudes : ndarray
        1D Nx1 array of latitude points (degrees)

    Returns
    -------
    latlens : ndarray
        1D Nx1 array of latitude lengths (m)

    See Also
    --------
    lonlens
    lonlen
    latlen

    References
    ----------
    Osborne, P. (2013). The Mercator Projections. Zenodo.
    <https://doi.org/10.5281/zenodo.35392>
    """
    gcc.cdim(latitudes)
    a=6378137
    b=6356752.3142
    e=np.sqrt((a**2-b**2)/a**2)
    latrads=latitudes*np.pi/180
    latspacing=spacing(latitudes)
    latlengths=latspacing*np.pi*a*(1-e**2)/(180*(1-(e*np.sin(latrads))**2)**(3/2))
    return latlengths

def lonlen(latitude):
    """
    Calculate the length(s) of a degree of longitude at the input
    latitude(s)

    Parameters
    ----------
    latitude : ndarray/float
        float or 1D Nx1 array of latitude point(s) (degrees)

    Returns
    -------
    lonlen : ndarray/float
        float or 1D Nx1 array of longitude length(s) (m)

    See Also
    --------
    lonlens
    latlens
    latlen

    References
    ----------
    Osborne, P. (2013). The Mercator Projections. Zenodo.
    <https://doi.org/10.5281/zenodo.35392>
    """
    a=6378137
    b=6356752.3142
    e=np.sqrt((a**2-b**2)/a**2)
    latrad=latitude*np.pi/180
    lonlength=np.pi*a*np.cos(latrad)/(180*np.sqrt(1-(e*np.sin(latrad))**2))
    return lonlength

def latlen(latitude):
    """
    Calculate the length(s) of a degree of latitude at the input
    latitude(s)

    Parameters
    ----------
    latitude : ndarray/float
        float or 1D Nx1 array of latitude point(s) (degrees)

    Returns
    -------
    lonlen : ndarray/float
        float or 1D Nx1 array of latitude length(s) (m)

    See Also
    --------
    lonlens
    latlens
    lonlen

    References
    ----------
    Osborne, P. (2013). The Mercator Projections. Zenodo.
    <https://doi.org/10.5281/zenodo.35392>
    """
    a=6378137
    b=6356752.3142
    e=np.sqrt((a**2-b**2)/a**2)
    latrad=latitude*np.pi/180
    latlength=np.pi*a*(1-e**2)/(180*(1-(e*np.sin(latrad))**2)**(3/2))
    return latlength

def areas(latitudes,longitudes,units=1):
    """
    Calculate the cell areas of a grid defined by a range of latitudes
    and longitudes

    Parameters
    ----------
    latitudes : ndarray
        1D Nx1 array of latitude points (degrees)
    longitudes : ndarray
        1D Nx1 array of longitude points (degrees)
    units : float, default=1 (m^2) 
        unit multiplier for areas
        [e.g. for km^2, mult = 1000m x 1000m = 1e6]

    Returns
    -------
    areas : ndarray
        2D MxN array of areas (default in m^2) where
        M=len(`latitudes`) and N=len(`longitudes`)

    See Also
    --------
    lonlens
    latlens

    References
    ----------
    Osborne, P. (2013). The Mercator Projections. Zenodo.
    <https://doi.org/10.5281/zenodo.35392>
    """
    gcc.cdim(latitudes)
    gcc.cdim(longitudes)
    lonspacing=spacing(longitudes)
    latlengths=latlens(latitudes)
    lonlengths=lonlens(latitudes,lonspacing=lonspacing)
    grdareas=np.repeat((latlengths*lonlengths)[:,np.newaxis],len(longitudes),axis=1)/units
    return grdareas

def clonrng(longitudes):
    """
    Find the range of a set of longitude points

    - negative (-180 to 180)
    - positive (0 to 360)
    - indeterminate (0 to 180)

    Parameters
    ----------
    longitudes : ndarray
        1D Nx1 array of longitude points (degrees)

    Returns
    -------
    longituderange : {'neg', 'pos', 'ind'} 
        descriptor for the range of the input `longitudes`
        
            ``neg``
                -180 to 180 degrees
            ``pos``
                0 to 360 degrees
            ``ind``
                0 to 180 degrees

    See Also
    --------
    switchlon
    switchind
    """
    if len(longitudes.shape)!=1:
        sys.exit('ERROR - Longitude input is not 1-dimensional')
    if (longitudes < 0).any() and (longitudes > 180).any():
        sys.exit('ERROR - Longitude input has values less than 0 and greater than 180, range should either be -180 to 180 or 0 to 360')
    elif (longitudes < 0).any():
        longituderange='neg'
    elif (longitudes > 180).any():
        longituderange='pos'
    else:
        longituderange='ind'
    return longituderange

def clatdir(latitudes):
    """
    Find the directionality of a set of latitude points
        - increasing (lowest to highest)
        - decreasing (highest to lowest)

    Parameters
    ----------
    latitudes : ndarray
        1D Nx1 array of latitude points (degrees)

    Returns
    -------
    latitudedirection : {'inc', 'dec'} 
        descriptor for the direction of the input `latitudes`

            ``inc``
                lowest to highest latitudes
            ``dec``
                highest to lowest latitudes
    """
    gcc.cdim(latitudes)
    indlatmin=latitudes.argmin()
    indlatmax=latitudes.argmax()
    if indlatmin==0 and indlatmax==len(latitudes)-1:
        latitudedirection='inc'
    elif indlatmin==len(latitudes)-1 and indlatmax==0:
        latitudedirection='dec'
    else:
        sys.exit('ERROR - Inconsistency in input latitude array, not ordered small to large or large to small')
    return latitudedirection

def switchlon(longitudes,outrange,print_output=False):
    """
    Switch a set of longitude points place between negative (-180 to
    180) and positive (0 to 360)

    Parameters
    ----------
    longitudes : ndarray
        1D Nx1 array of longitude points (degrees)
    outrange : {'pos', 'neg'} 
        select range for longitude outputs 
        [i.e. positive (0 to 360) or negative (-180 to 180)] 
    print_output : bool, default=False
        flag to print out whether input was altered

    Returns
    -------
    switchlon : ndarray
        1D Nx1 array of longitude points equal to input `longitudes`,
        but altered (if necessary) to be consistent with selected output
        range

    See Also
    --------
    clonrng
    switchind
    """
    if outrange!='pos' and outrange!='neg':
        sys.exit('ERROR - Unrecognized longitude output range. Valid inputs are \'pos\' or \'neg\'')
    inprange=clonrng(longitudes)
    if inprange==outrange or inprange=='ind':
        switchlon=longitudes
        if print_output==True:
            print('Input already has values consistent with selected output range')
    elif inprange=='neg' and outrange=='pos':
        switchlon=longitudes+360*(longitudes<0)
        if print_output==True:
            print('Switched negative input to positive output')
    elif inprange=='pos' and outrange=='neg':
        switchlon=longitudes-360*(longitudes>180)
        if print_output==True:
            print('Switched positive input to negative output')
    else:
        sys.exit('ERROR - Function should not reach this point, check source')
    return switchlon

def switchind(longitudes):
    """
    Find the index where a longitude array crosses either 0 or 180
    degrees

    Parameters
    ----------
    longitudes : ndarray
        1D Nx1 array of longitude points (degrees)

    Returns
    -------
    switchind : int 
        index for where input `longitudes` crosses 0 or 180 degrees, or
        0 if no such point exists

    See Also
    --------
    clonrng
    switchlon
    """
    gcc.cdim(longitudes)
    inprange=clonrng(longitudes)
    if inprange=='ind':
        switchind=0
    elif inprange=='neg':
        switchind=np.nonzero(longitudes>=0)[0][0]
    elif inprange=='pos':
        switchind=np.nonzero(longitudes>180)[0][0]
    return switchind
