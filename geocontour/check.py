import sys
import warnings
import numpy as np

def cdim(dimension,exit_on_error=True):
    """
    Checks an input dimension array for 1-dimensionality and regular spacing

    Inputs (Required):
        dimension - A numpy array of longitude or latitude points (degrees)

    Inputs (Optional):
        exit_on_error - Boolean for exiting if an error is encountered, default=True

    Outputs:
        none
    """
    if len(dimension.shape)!=1:
        sys.exit('ERROR - Input dimension (latitude or longitude) is not 1-dimensional')
    spacing=np.diff(dimension)
    if not np.all(spacing==spacing[0]):
        if not np.allclose(spacing,spacing[0],rtol=1e-10,atol=1e-10):
            if exit_on_error:
                sys.exit('ERROR - Input dimension (latitude or longitude) has irregular spacing')
            else:
                warnings.warn('WARNING - Input dimension (latitude or longitude) has irregular spacing. Output may be innacurate')

def cboundary(boundary):
    """
    Checks a list of boundary points for 2-dimensionality and proper ordering
        Will check whether columns (lat/lon) are ordered correctly but CAN'T GUARANTEE THIS

    Inputs(Required):
        boundary - A 2-d Nx2 numpy array of latitude/longitude points (degrees) with the last point equal to the first

    Outputs:
        none
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
        sys.exit('ERROR - Boundary column 2 has values less than 0 and greater than 180, column 2 should be longitude')
    if not (boundary[0,:]==boundary[-1,:]).all():
        sys.exit('ERROR - First and last boundary points are not equal, boundary should close')
    if (np.diff(boundary[:,1]) > 300).any():
        warnings.warn('WARNING - at least one boundary longitude span over 300 deg in length, boundary may cross meridian/date line\n    Suggestion: use geocontour.grid.switchlon() to swap boundary range')

def cmask(mask,latitudes=None,longitudes=None):
    """
    Checks a mask for correct data type and dimensionality, and size if optional latitudes and longitudes are provided

    Inputs (Required):
        mask - A 2-d boolean numpy array of dimension MxN where M=len(latitudes) and N=len(longitudes)

    Inputs (Optional):
        latitudes - An evenly spaced numpy array of latitude points (degrees)
        longitudes - An evenly spaced numpy array of latitude points (degrees)

    Outputs:
        none
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
    Check contour for repeating cells, closure, and connectivity, and latitude/longitude range if optional latitudes and longitudes are provided

    Inputs (Required):
        contour - A 2-d Nx2 numpy array of ordered latitude/longitude points (degrees) describing the edge of a mask
        latitudes - An evenly spaced numpy array of latitude points (degrees)
        longitudes - An evenly spaced numpy array of longitude points (degrees)

    Outputs:
        none
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

    Inputs (Required):
        geocontour - A 3-d Nx2x5 numpy array defining a list of N contour cells and their edge points, lengths, and outward unit vectors
        latitudes - An evenly spaced numpy array of latitude points (degrees)
        longitudes - An evenly spaced numpy array of longitude points (degrees)

    Outputs:
        none
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

