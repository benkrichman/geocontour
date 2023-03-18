"""
Functions for checking the properties of various structures utilized by
=======================================================================
geocontour
==========
"""
import sys
import warnings
import numpy as np

def cdim(dimension,exit_on_error=True):
    """
    Check input dimension array for 1-dimensionality and regular
    spacing

    Parameters
    ----------
    dimension : ndarray
        1D Nx1 array of longitude or latitude points (degrees)
    exit_on_error : bool, default=True 
        flag to exit if an error is encountered

    Returns
    -------
    None

    See Also
    --------
    cboundary
    cmask
    ccontour
    cgeocontour
    """
    if len(dimension.shape)!=1:
        sys.exit('ERROR - Input dimension (latitude or longitude) is not 1-dimensional')
    spacing=np.diff(dimension)
    if not np.all(spacing==spacing[0]):
        if not np.allclose(spacing,spacing[0],rtol=1e-10,atol=1e-10):
            if exit_on_error:
                sys.exit('ERROR - Input dimension (latitude or longitude) has irregular spacing')
            else:
                warnings.warn('WARNING - Input dimension (latitude or longitude) has irregular spacing. Output may be inaccurate')

def cboundary(boundary):
    """
    Check array of boundary points for 2-dimensionality and proper
    ordering

    Parameters
    ----------
    boundary : ndarray 
        2D Nx2 array of latitude/longitude points (degrees) with the
        last point equal to the first

    Returns
    -------
    None

    See Also
    --------
    cdim
    cmask
    ccontour
    cgeocontour

    Notes
    -----
    Will check whether columns (lat/lon) are ordered correctly but
    **CAN'T GUARANTEE THIS**

    """
    if len(boundary.shape)!=2:
        sys.exit('ERROR - Boundary input is not 2-dimensional')
    if boundary.shape[1]!=2:
        sys.exit('ERROR - Boundary input is not oriented correctly, should be Nx2')
    if boundary.shape[0]<3:
        sys.exit('ERROR - Boundary input has 2 or less points, can not form boundary')
    if (boundary[:,0] > 90).any() or (boundary[:,0] < -90).any():
        sys.exit('ERROR - Boundary column 1 has values less than -90 or greater than 90, column 1 should be latitude')
    if (boundary[:,1] < -180).any() or (boundary[:,1] > 360).any():
        sys.exit('ERROR - Boundary column 2 has values less than -180 or greater than 360, column 2 should be longitude')
    if (boundary[:,1] < 0).any() and (boundary[:,1] > 180).any():
        sys.exit('ERROR - Boundary column 2 has values less than 0 and greater than 180, column 2 should be longitude and range should either be -180 to 180 or 0 to 360')
    if not (boundary[0,:]==boundary[-1,:]).all():
        sys.exit('ERROR - First and last boundary points are not equal, boundary should close')
    if (np.diff(boundary[:,1]) > 300).any():
        warnings.warn('WARNING - at least one boundary longitude span over 300 deg in length, boundary may cross meridian/date line\n    Suggestion: use geocontour.grid.switchlon() to swap boundary range')

def cmask(mask,latitudes=None,longitudes=None):
    """
    Check mask array for correct data type and dimensionality, and
    optionally size

    Parameters
    ----------
    mask : ndarray
        2D MxN bool array where M=len(`latitudes`) and
        N=len(`longitudes`)
    latitudes : ndarray, optional
        1D Nx1 array of latitude points (degrees)
    longitudes : ndarray, optional
        1D Nx1 array of longitude points (degrees)

    Returns
    -------
    None

    See Also
    --------
    cdim
    cboundary
    ccontour
    cgeocontour
    """
    if len(mask.shape)!=2:
        sys.exit('ERROR - Mask input is not 2-dimensional')
    if mask.dtype!=bool:
        sys.exit('ERROR - Mask data type is not bool')
    if latitudes is not None and longitudes is not None:
        cdim(latitudes)
        cdim(longitudes)
        if mask.shape[0]!=len(latitudes):
            sys.exit('ERROR - Mask dimension 0 differs from length of latitude array')
        if mask.shape[1]!=len(longitudes):
            sys.exit('ERROR - Mask dimension 1 differs from length of longitude array')

def ccontour(contour,latitudes=None,longitudes=None):
    """
    Check contour for repeating cells, closure, and connectivity, and
    optionally lat/lon range

    Parameters
    ----------
    contour : ndarray
        2D Nx2 array of ordered latitude/longitude points (degrees)
        describing the edge of a mask
    latitudes : ndarray, optional
        1D Nx1 array of latitude points (degrees)
    longitudes : ndarray, optional
        1D Nx1 array of longitude points (degrees)

    Returns
    -------
    None

    See Also
    --------
    cdim
    cboundary
    cmask
    cgeocontour
    """
    contourdiff=np.diff(contour,axis=0)
    if latitudes is not None and longitudes is not None:
        cdim(longitudes)
        cdim(latitudes)
        lonspc=abs(np.diff(longitudes)[0])
        latspc=abs(np.diff(latitudes)[0])
        spacing=np.array([latspc,lonspc])
        if (contour[:,0]<latitudes.min()).any() or (contour[:,0]>latitudes.max()).any():
            sys.exit('Input contour exceeds latitude range')
        if (contour[:,1]<longitudes.min()).any() or (contour[:,1]>longitudes.max()).any():
            sys.exit('Input contour exceeds longitude range')
        if not ((abs(contourdiff)==0)+(abs(contourdiff)==spacing)+(abs(contourdiff)==np.linalg.norm(spacing))).all():
            sys.exit('ERROR - Contour does not seem to be connected, ensure contour is continuous (8-connected or 4-connected)')
    else:
        spacing=np.array([1,1])
        if not ((abs(contourdiff)==0)+(abs(contourdiff)==spacing)+(abs(contourdiff)==np.linalg.norm(spacing))).all():
            sys.exit('ERROR - Contour does not seem to be connected, ensure contour is continuous (8-connected or 4-connected)')
    if (contourdiff==0).all(axis=1).any():
        sys.exit('ERROR - Input contour has consecutive repeating cells')
    if ((contour[0,:]-contour[-1,:])!=0).any():
        sys.exit('ERROR - Input contour does not close: last_point!=first_point')

def cgeocontour(geocontour,latitudes,longitudes):
    """
    Check geocontour for latitude/longitude range and dimension

    Parameters
    ----------
    geocontour : ndarray
        3D Nx2x5 array defining a list of N contour cells (column 1),
        their edge points (columns 2,3), segment lengths (column 4), and
        outward unit vectors (column 5)
    latitudes : ndarray
        1D Nx1 array of latitude points (degrees)
    longitudes : ndarray
        1D Nx1 array of longitude points (degrees)

    Returns
    -------
    None

    See Also
    --------
    cdim
    cboundary
    cmask
    ccontour
    """
    cdim(longitudes)
    cdim(latitudes)
    lonspc=abs(np.diff(longitudes)[0])
    latspc=abs(np.diff(latitudes)[0])
    spacing=np.array([latspc,lonspc])
    if (geocontour[:,0,0]<latitudes.min()).any() or (geocontour[:,0,0]>latitudes.max()).any():
        sys.exit('Input geocontour exceeds latitude range')
    if (geocontour[:,1,0]<longitudes.min()).any() or (geocontour[:,1,0]>longitudes.max()).any():
        sys.exit('Input geocontour exceeds longitude range')

