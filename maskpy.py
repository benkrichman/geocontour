import sys
import warnings
import numpy as np
from scipy.signal import convolve
import shapely.geometry as shg
import matplotlib.pyplot as plt
import matplotlib.path as mplp
import matplotlib.colors as mplc
try:
    import cartopy as cp
    cp_exists='y'
except:
    cp_exists='n'

###Check Functions

def checkdim(dimension,exit_on_error=True):
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

def checkboundary(boundary):
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
        warnings.warn('WARNING - at least one boundary span over 300 deg in length, boundary may cross meridian/date line')

def checkmask(mask,latitudes=None,longitudes=None):
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
        if mask.shape[0]!=len(latitudes):
            sys.exit('ERROR - Mask dimension 0 differs from length of latitude array')
        if mask.shape[1]!=len(longitudes):
            sys.exit('ERROR - Mask dimension 1 differs from length of longitude array')

def checkcontour(contour,latitudes=None,longitudes=None):
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
        lonspc=gridspacing(longitudes)
        latspc=gridspacing(latitudes)
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

def checkgeocontour(geocontour,latitudes,longitudes):
    """
    Check geocontour for latitude/longitude range and dimension

    Inputs (Required):
        geocontour - A 3-d Nx2x5 numpy array defining a list of N contour cells and their edge points, lengths, and outward unit vectors
        latitudes - An evenly spaced numpy array of latitude points (degrees)
        longitudes - An evenly spaced numpy array of longitude points (degrees)

    Outputs:
        none
    """
    lonspc=gridspacing(longitudes)
    latspc=gridspacing(latitudes)
    spacing=np.array([latspc,lonspc])
    if (geocontour[:,0,0]<latitudes.min()).any() or (geocontour[:,0,0]>latitudes.max()).any():
        sys.exit('Input geocontour exceeds latitude range')
    if (geocontour[:,1,0]<longitudes.min()).any() or (geocontour[:,1,0]>longitudes.max()).any():
        sys.exit('Input geocontour exceeds longitude range')

###Grid Functions

def gridspacing(dimension):
    """
    Returns the grid spacing for a given input dimension
    
    Inputs (Required):
        dimension - An evenly spaced numpy array of longitude or latitude points (degrees)
    
    Outputs:
        gridspacing - A positive float indicating the spacing of the input dimension (degrees)
    """
    checkdim(dimension)
    spacing=abs(np.diff(dimension)[0])
    return spacing

def longitudelengths(latitudes,spacing=1):
    """
    Returns the lengths of a degree (default) of longitude over a range of latitudes 
    Source: https://doi.org/10.5281/ZENODO.35392

    Inputs (Required):
        latitudes - An evenly spaced numpy array of latitude points (degrees)

    Inputs (Optional):
        spacing - A float/int specifying longitude spacing
            e.g. for lengths of a half degree of longitude at the given latitudes, spacing=0.5

    Outputs:
        longitudelengths - A numpy array of longitude lengths (m)
    """
    checkdim(latitudes)
    a=6378137
    b=6356752.3142
    e=np.sqrt((a**2-b**2)/a**2)
    latrads=latitudes*np.pi/180
    lonlengths=spacing*np.pi*a*np.cos(latrads)/(180*np.sqrt(1-(e*np.sin(latrads))**2))
    return lonlengths

def latitudelengths(latitudes):
    """
    Returns the grid lengths of a defined range of latitudes 
    Source: https://doi.org/10.5281/ZENODO.35392

    Inputs (Required):
        latitudes - An evenly spaced numpy array of latitude points (degrees)

    Outputs:
        latitudelengths - A numpy array of latitude lengths (m)
    """
    checkdim(latitudes)
    a=6378137
    b=6356752.3142
    e=np.sqrt((a**2-b**2)/a**2)
    latrads=latitudes*np.pi/180
    spacing=gridspacing(latitudes)
    latlengths=spacing*np.pi*a*(1-e**2)/(180*(1-(e*np.sin(latrads))**2)**(3/2))
    return latlengths

def longitudelength(latitude):
    """
    Returns the length of a degree of longitude at the input latitude
    Source: https://doi.org/10.5281/ZENODO.35392

    Inputs (Required):
        latitude - A latitude point or array of points (degrees)

    Outputs:
        longitudelength - The length of a degree of longitude at the input latitude(s) (m)
    """
    a=6378137
    b=6356752.3142
    e=np.sqrt((a**2-b**2)/a**2)
    latrad=latitude*np.pi/180
    lonlength=np.pi*a*np.cos(latrad)/(180*np.sqrt(1-(e*np.sin(latrad))**2))
    return lonlength

def latitudelength(latitude):
    """
    Returns the length of a degree of latitude at the input latitude
    Source: https://doi.org/10.5281/ZENODO.35392

    Inputs (Required):
        latitude - A latitude point or array of points (degrees)

    Outputs:
        latitudelength - The length of a degree of latitude at the input latitude(s) (m)
    """
    a=6378137
    b=6356752.3142
    e=np.sqrt((a**2-b**2)/a**2)
    latrad=latitude*np.pi/180
    latlength=np.pi*a*(1-e**2)/(180*(1-(e*np.sin(latrad))**2)**(3/2))
    return latlength

def gridareas(latitudes,longitudes,units=1):
    """
    Returns the cell areas of a grid defined by a range of latitudes and longitudes

    Inputs (Required):
        latitudes - An evenly spaced numpy array of latitude points (degrees)
        longitudes - An evenly spaced numpy array of longitude points (degrees)

    Inputs (Optional):
        units - unit multiplier for outputs, default=1 (m^2) 
            e.g. for km^2, mult=1000m x 1000m=1e6

    Outputs:
        gridareas - A numpy array of areas, with dimension 0: latitude and dimension 1: longitude (default in m^2)
    """
    checkdim(latitudes)
    checkdim(longitudes)
    lonspacing=gridspacing(longitudes)
    latlengths=latitudelengths(latitudes)
    lonlengths=longitudelengths(latitudes,spacing=lonspacing)
    grdareas=np.repeat((latlengths*lonlengths)[:,np.newaxis],len(longitudes),axis=1)/units
    return grdareas

def checklongituderange(longitudes):
    """
    Returns a descriptor for the range of a set of longitude points 
        negative (-180 to 180), positive (0 to 360), or indeterminate (0 to 180) range

    Inputs (Required):
        longitudes - A numpy array of longitude points (degrees)

    Outputs:
        longituderange ('neg'/'pos'/'ind') - A descriptor for the longitude range of the input
    """
    if len(longitudes.shape)!=1:
        sys.exit('ERROR - Longitude input is not 1-dimensional')
    if (longitudes < 0).any() and (longitudes > 180).any():
        sys.exit('ERROR - Longitude input has values less than 0 and greater than 180')
    elif (longitudes < 0).any():
        longituderange='neg'
    elif (longitudes > 180).any():
        longituderange='pos'
    else:
        longituderange='ind'
    return longituderange

def checklatitudedirection(latitudes):
    """
    Returns a descriptor for the direction of a set of latitude points
        increasing or decreasing

    Inputs (Required):
        latitudes - A numpy array of latitude points (degrees)

    Outputs:
        latitudedirection ('inc'/'dec') - A descriptor for the direction of the input
    """
    checkdim(latitudes)
    indlatmin=latitudes.argmin()
    indlatmax=latitudes.argmax()
    if indlatmin==0 and indlatmax==len(latitudes)-1:
        latitudedirection='inc'
    elif indlatmin==len(latitudes)-1 and indlatmax==0:
        latitudedirection='dec'
    else:
        sys.exit('ERROR - Inconsistency in input latitude array, not ordered small to large or large to small')
    return latitudedirection

def switchlongitudes(longitudes,outrange,print_output='n'):
    """
    Returns a set of longitude points switched in place between negative (-180 to 180) and positive (0 to 360)

    Inputs (Required):
        longitudes - A numpy array of longitude points (degrees)
        outrange - ('pos'/'neg') Selection of range for longitude outputs: negative (-180 to 180) or positive (0 to 360)

    Inputs(Optional):
        print_output ('y'/'n') - flag to print out whether input was altered, default='n'

    Outputs:
        switchlongitudes - A numpy array of longitudes, altered if not consistent with selected output range
    """
    if outrange!='pos' and outrange!='neg':
        sys.exit('ERROR - Unrecognized longitude output range. Valid inputs are \'pos\' or \'neg\'')
    inprange=checklongituderange(longitudes)
    if inprange==outrange or inprange=='ind':
        switchlongitudes=longitudes
        if print_output=='y':
            print('Input already has values consistent with selected output range')
    elif inprange=='neg' and outrange=='pos':
        switchlongitudes=longitudes+360*(longitudes<0)
        if print_output=='y':
            print('Switched negative input to positive output')
    elif inprange=='pos' and outrange=='neg':
        switchlongitudes=longitudes-360*(longitudes>180)
        if print_output=='y':
            print('Switched positive input to negative output')
    else:
        sys.exit('ERROR - Function should not reach this point, check source')
    return switchlongitudes

def findswitchindex(longitudes):
    """
    Returns the index where a longitude array either crosses 0 or 180 degrees

    Inputs (Required):
        longitudes - An evenly spaced numpy array of longitude points (degrees)

    Outputs:
        switchindex - An index for where array crosses 0 or 180 degrees, or 0 if no such point exists
    """
    checkdim(longitudes)
    inprange=checklongituderange(longitudes)
    if inprange=='ind':
        switchindex=0
    elif inprange=='neg':
        switchindex=np.nonzero(longitudes>=0)[0][0]
    elif inprange=='pos':
        switchindex=np.nonzero(longitudes>180)[0][0]
    return switchindex

###Mask Fucntions

def maskboxset(latitudes,longitudes,boundary):
    """
    Checks input dimensions (lat/lon) against input boundary and returns min/max indicies of bounding box

    Inputs (Required):
        latitudes - An evenly spaced numpy array of latitude points (degrees)
        longitudes - An evenly spaced numpy array of longitude points (degrees)
        boundary - A 2-d Nx2 numpy array of latitude/longitude points (degrees)

    Outputs:
        boxlatmin - The minimum bounding box latitude index
        boxlatmax - The maximum bounding box latitude index
        boxlonmin - The minimum bounding box longitude index
        boxlonmax - The maximum bounding box longitude index
    """
    checkdim(latitudes)
    checkdim(longitudes)
    checkboundary(boundary)
    latdir=checklatitudedirection(latitudes)
    loninput=checklongituderange(longitudes)
    boundloninput=checklongituderange(boundary[:,1])
    if (loninput=='neg' and boundloninput=='pos') or (loninput=='pos' and boundloninput=='neg'):
        sys.exit('ERROR - Longitude input range is '+loninput+' and boundary longitude range is '+boundloninput)
    boundlatmin=boundary.min(axis=0)[0]
    boundlatmax=boundary.max(axis=0)[0]
    boundlonmin=boundary.min(axis=0)[1]
    boundlonmax=boundary.max(axis=0)[1]
    if latdir=='inc':
        if boundlatmin<=latitudes.min():
            warnings.warn('WARNING - boundary latitude minimum on or outside of minimum latitude')
            boxlatmin=0
        else:
            boxlatmin=((latitudes-boundlatmin)<0).nonzero()[0][-1]
        if boundlatmax>=latitudes.max():
            warnings.warn('WARNING - boundary latitude maximum on or outside of maximum latitude')
            boxlatmax=len(latitudes)-1
        else:
            boxlatmax=((latitudes-boundlatmax)>0).nonzero()[0][0]
    elif latdir=='dec':
        if boundlatmin<=latitudes.min():
            warnings.warn('WARNING - boundary latitude minimum on or outside of minimum latitude')
            boxlatmin=len(latitudes)-1
        else:
            boxlatmin=((latitudes-boundlatmin)<0).nonzero()[0][0]
        if boundlatmax>=latitudes.max():
            warnings.warn('WARNING - boundary latitude maximum on or outside of maximum latitude')
            boxlatmax=0
        else:
            boxlatmax=((latitudes-boundlatmax)>0).nonzero()[0][-1]
    if boundlonmin<=longitudes.min():
        warnings.warn('WARNING - boundary longitude minimum on or outside of minimum longitude')
        boxlonmin=0
    else:
        boxlonmin=((longitudes-boundlonmin)<0).nonzero()[0][-1]
    if boundlonmax>=longitudes.max():
        warnings.warn('WARNING - boundary longitude maximum on or outside of maximum longitude')
        boxlonmax=len(longitudes)-1
    else:
        boxlonmax=((longitudes-boundlonmax)>0).nonzero()[0][0]
    return boxlatmin, boxlatmax, boxlonmin, boxlonmax

def masksearchcenter(latitudes,longitudes,boundary,precision=1e-5):
    """
    Returns a mask over a range of input latitudes and longitudes determined by an input boundary
        Critera for inclusion of a cell is whether the center of the cell falls within the boundary

    Inputs (Required):
        latitudes - An evenly spaced numpy array of latitude points (degrees)
        longitudes - An evenly spaced numpy array of longitude points (degrees)
        boundary - A 2-d Nx2 numpy array of latitude/longitude points (degrees)

    Inputs (Optional):
        precision - Passed to shapely to increase boundary polygon area, default=1e-5
            Shapely can't beat machine precision, and can thus give "incorrect" results for very close points or shapes. This input errs on being more inclusive, in particular to capture points falling directly on a boundary. A decent rule is to set the precision value as high as you can without impeding the accuracy. For instance, the default of 1e-5 (degrees) translates to roughly 1m precision at the equator. The buffer can be negated by setting this input very low (to machine precision).

    Outputs:
        mask - A 2-d boolean numpy array of dimension MxN where M=len(latitudes) and N=len(longitudes)
    """
    latdir=checklatitudedirection(latitudes)
    if latdir=='dec':
        latitudes=np.flip(latitudes)
    boxlatmin, boxlatmax, boxlonmin, boxlonmax = maskboxset(latitudes,longitudes,boundary)
    boxmask=np.full((boxlatmax-boxlatmin+1,boxlonmax-boxlonmin+1),False)
    boundpoly=shg.Polygon(boundary).buffer(precision)
    for la in np.arange(boxlatmin,boxlatmax+1,1):
        for lo in np.arange(boxlonmin,boxlonmax+1,1):
            center=shg.Point(latitudes[la],longitudes[lo])
            boxmask[la-boxlatmin,lo-boxlonmin]=boundpoly.contains(center)
    mask=np.full((len(latitudes),len(longitudes)),False)
    mask[boxlatmin:boxlatmax+1,boxlonmin:boxlonmax+1]=boxmask
    if latdir=='dec':
        mask=np.flip(mask,axis=0)
    return mask

def masksearchcenter2(latitudes,longitudes,boundary):
    """
    Returns a mask over a range of input latitudes and longitudes determined by an input boundary
        Critera for inclusion of a cell is whether the center of the cell falls within the boundary
        Functionally matches masksearchcenter(), but utilizes matplotlib.path functions, which are probably optimized and thus is roughly 2.5*sqrt(N) faster for N points, though lacks a "precision" buffer input
    Source:
        https://stackoverflow.com/questions/50847827/how-can-i-select-the-pixels-that-fall-within-a-contour-in-an-image-represented-b
        https://stackoverflow.com/questions/16625507/checking-if-a-point-is-inside-a-polygon/23453678#23453678
        https://stackoverflow.com/questions/36399381/whats-the-fastest-way-of-checking-if-a-point-is-inside-a-polygon-in-python
        https://matplotlib.org/stable/api/path_api.html#matplotlib.path.Path.contains_point

    Inputs (Required):
        latitudes - An evenly spaced numpy array of latitude points (degrees)
        longitudes - An evenly spaced numpy array of longitude points (degrees)
        boundary - A 2-d Nx2 numpy array of latitude/longitude points (degrees)

    Outputs:
        mask - A 2-d boolean numpy array of dimension MxN where M=len(latitudes) and N=len(longitudes)
    """
    latdir=checklatitudedirection(latitudes)
    if latdir=='dec':
        latitudes=np.flip(latitudes)
    boxlatmin, boxlatmax, boxlonmin, boxlonmax = maskboxset(latitudes,longitudes,boundary)
    boundpoly=mplp.Path(boundary)
    ysp, xsp = np.meshgrid(latitudes[boxlatmin:boxlatmax+1],longitudes[boxlonmin:boxlonmax+1], indexing='ij')
    searchpoints=np.hstack((ysp.reshape((-1,1)), xsp.reshape((-1,1))))
    boxmask=boundpoly.contains_points(searchpoints)
    mask=np.full((len(latitudes),len(longitudes)),False)
    mask[boxlatmin:boxlatmax+1,boxlonmin:boxlonmax+1]=boxmask.reshape((boxlatmax-boxlatmin+1,boxlonmax-boxlonmin+1))
    if latdir=='dec':
        mask=np.flip(mask,axis=0)
    return mask

def masksearchnodes(latitudes,longitudes,boundary,nodes=2,precision=1e-5):
    """
    Returns a mask over a range of input latitudes and longitudes determined by an input boundary
        Critera for inclusion of a cell is whether a given number (default=2) of cell nodes (corners) fall within the boundary 

    Inputs (Required):
        latitudes - An evenly spaced numpy array of latitude points (degrees)
        longitudes - An evenly spaced numpy array of longitude points (degrees)
        boundary - A 2-d Nx2 numpy array of latitude/longitude points (degrees)

    Inputs (Optional):
        nodes - The number of cell nodes (corners) to use as a criteria for inclusion (1-4)
        precision - Passed to shapely to increase boundary polygon area, default=1e-5
            Shapely can't beat machine precision, and can thus give "incorrect" results for very close points or shapes. This input errs on being more inclusive, in particular to capture points falling directly on a boundary. A decent rule is to set the precision value as high as you can without impeding the accuracy. For instance, the default of 1e-5 (degrees) translates to roughly 1m precision at the equator. The buffer can be negated by setting this input very low (to machine precision).

    Outputs:
        mask - A 2-d boolean numpy array of dimension MxN where M=len(latitudes) and N=len(longitudes)
    """
    if nodes<1:
        warnings.warn('WARNING - valid input for nodes is 1-4, nodes<1 will result in all cells being selected')
    if nodes>4:
        warnings.warn('WARNING - valid input for nodes is 1-4, nodes>4 will result in no cells being selected')
    latdir=checklatitudedirection(latitudes)
    if latdir=='dec':
        latitudes=np.flip(latitudes)
    boxlatmin, boxlatmax, boxlonmin, boxlonmax = maskboxset(latitudes,longitudes,boundary)
    latgrdspc=gridspacing(latitudes)
    longrdspc=gridspacing(longitudes)
    boxmask=np.full((boxlatmax-boxlatmin+1,boxlonmax-boxlonmin+1),False)
    boundpoly=shg.Polygon(boundary).buffer(precision)
    for la in np.arange(boxlatmin,boxlatmax+1,1):
        for lo in np.arange(boxlonmin,boxlonmax+1,1):
            nodeLL=shg.Point(latitudes[la]-latgrdspc/2,longitudes[lo]-longrdspc/2)
            nodeHL=shg.Point(latitudes[la]+latgrdspc/2,longitudes[lo]-longrdspc/2)
            nodeLH=shg.Point(latitudes[la]-latgrdspc/2,longitudes[lo]+longrdspc/2)
            nodeHH=shg.Point(latitudes[la]+latgrdspc/2,longitudes[lo]+longrdspc/2)
            nodesinmask=np.array([boundpoly.contains(nodeLL),boundpoly.contains(nodeHL),boundpoly.contains(nodeLH),boundpoly.contains(nodeHH)])
            if nodesinmask.sum()>=nodes:
                boxmask[la-boxlatmin,lo-boxlonmin]=True
    mask=np.full((len(latitudes),len(longitudes)),False)
    mask[boxlatmin:boxlatmax+1,boxlonmin:boxlonmax+1]=boxmask
    if latdir=='dec':
        mask=np.flip(mask,axis=0)
    return mask

def masksearchnodes2(latitudes,longitudes,boundary,nodes=2,precision=1e-5):
    """
    Returns a mask over a range of input latitudes and longitudes determined by an input boundary
        Critera for inclusion of a cell is whether a given number (default=2) of cell nodes (corners) fall within the boundary 
        Functionally matches masksearchnodes(), but utilizes matplotlib.path functions, though speed is similar to the shapely implementation

    Inputs (Required):
        latitudes - An evenly spaced numpy array of latitude points (degrees)
        longitudes - An evenly spaced numpy array of longitude points (degrees)
        boundary - A 2-d Nx2 numpy array of latitude/longitude points (degrees)

    Inputs (Optional):
        nodes - The number of cell nodes (corners) to use as a criteria for inclusion (1-4)

    Outputs:
        mask - A 2-d boolean numpy array of dimension MxN where M=len(latitudes) and N=len(longitudes)
    """
    if nodes<1:
        warnings.warn('WARNING - valid input for nodes is 1-4, nodes<1 will result in all cells being selected')
    if nodes>4:
        warnings.warn('WARNING - valid input for nodes is 1-4, nodes>4 will result in no cells being selected')
    latdir=checklatitudedirection(latitudes)
    if latdir=='dec':
        latitudes=np.flip(latitudes)
    boxlatmin, boxlatmax, boxlonmin, boxlonmax = maskboxset(latitudes,longitudes,boundary)
    latgrdspc=gridspacing(latitudes)
    longrdspc=gridspacing(longitudes)
    boxmask=np.full((boxlatmax-boxlatmin+1,boxlonmax-boxlonmin+1),False)
    boundpoly=mplp.Path(boundary)
    for la in np.arange(boxlatmin,boxlatmax+1,1):
        for lo in np.arange(boxlonmin,boxlonmax+1,1):
            nodeLL=[latitudes[la]-latgrdspc/2,longitudes[lo]-longrdspc/2]
            nodeHL=[latitudes[la]+latgrdspc/2,longitudes[lo]-longrdspc/2]
            nodeLH=[latitudes[la]-latgrdspc/2,longitudes[lo]+longrdspc/2]
            nodeHH=[latitudes[la]+latgrdspc/2,longitudes[lo]+longrdspc/2]
            nodesinmask=boundpoly.contains_points(np.array([nodeLL,nodeHL,nodeLH,nodeHH]))
            if nodesinmask.sum()>=nodes:
                boxmask[la-boxlatmin,lo-boxlonmin]=True
    mask=np.full((len(latitudes),len(longitudes)),False)
    mask[boxlatmin:boxlatmax+1,boxlonmin:boxlonmax+1]=boxmask
    if latdir=='dec':
        mask=np.flip(mask,axis=0)
    return mask

def masksearcharea(latitudes,longitudes,boundary,area=0.5):
    """
    Returns a mask over a range of input latitudes and longitudes determined by an input boundary
        Critera for inclusion of a cell is whether the area of the cell enclosed by the boundary is greater than some fraction (default=0.5) 

    Inputs (Required):
        latitudes - An evenly spaced numpy array of latitude points (degrees)
        longitudes - An evenly spaced numpy array of longitude points (degrees)
        boundary - A 2-d Nx2 numpy array of latitude/longitude points (degrees)

    Inputs (Optional):
        area - The fraction of cell area enclosed by the boundary to use as a criteria for inclusion (0-1)

    Outputs:
        mask - A 2-d boolean numpy array of dimension MxN where M=len(latitudes) and N=len(longitudes)
    """
    if area>1:
        warnings.warn('WARNING - valid input for area is 0-1, area>1 will result in no cells being selected')
    if area<=0:
        warnings.warn('WARNING - valid input for area is 0-1, area<=0 will result in all cells being selected')
    latdir=checklatitudedirection(latitudes)
    if latdir=='dec':
        latitudes=np.flip(latitudes)
    boxlatmin, boxlatmax, boxlonmin, boxlonmax = maskboxset(latitudes,longitudes,boundary)
    latgrdspc=gridspacing(latitudes)
    longrdspc=gridspacing(longitudes)
    boxmask=np.full((boxlatmax-boxlatmin+1,boxlonmax-boxlonmin+1),False)
    boundpoly=shg.Polygon(boundary)
    for la in np.arange(boxlatmin,boxlatmax+1,1):
        for lo in np.arange(boxlonmin,boxlonmax+1,1):
            LL=[latitudes[la]-latgrdspc/2,longitudes[lo]-longrdspc/2]
            HL=[latitudes[la]+latgrdspc/2,longitudes[lo]-longrdspc/2]
            LH=[latitudes[la]-latgrdspc/2,longitudes[lo]+longrdspc/2]
            HH=[latitudes[la]+latgrdspc/2,longitudes[lo]+longrdspc/2]
            cell=shg.Polygon([LL,HL,HH,LH])
            if boundpoly.intersection(cell).area>=latgrdspc*longrdspc*area:
                boxmask[la-boxlatmin,lo-boxlonmin]=True
    mask=np.full((len(latitudes),len(longitudes)),False)
    mask[boxlatmin:boxlatmax+1,boxlonmin:boxlonmax+1]=boxmask
    if latdir=='dec':
        mask=np.flip(mask,axis=0)
    return mask

def maskedgecells(mask,latitudes=None,longitudes=None,connectivity=8):
    """
    Returns a mask of only the edge cells, and if latitudes and longitudes are provided also returns an array of the edge cells 
    
    Inputs (Required):
        mask - A 2-d boolean numpy array of dimension MxN where M=len(latitudes) and N=len(longitudes)

    Inputs (Optional):
        latitudes - An evenly spaced numpy array of latitude points (degrees)
        longitudes - An evenly spaced numpy array of longitude points (degrees)
        connectivity (4/8) - The connectivity parameter for testing edge cells, default=8
            4 tests only lateral neighbors while 8 tests lateral and diagonal neighbors

    Outputs:
        edgemask - A 2-d boolean numpy array of the same dimensions as mask input
        edgecells - A 2-d Nx2 numpy array of latitude/longitude points (degrees) of edge cells (unordered), where N is number of edge cells 
            Only output if optional latitude and longitude arguments are provided
    """
    if latitudes is not None and longitudes is not None:
        checkmask(mask,latitudes,longitudes)
    else:
        checkmask(mask)
    kernel=np.ones((3,3),dtype='int')
    kernel[1,1]=0
    if connectivity==8:
        pass
    elif connectivity==4:
        kernel[0,0]=0
        kernel[2,0]=0
        kernel[0,2]=0
        kernel[2,2]=0
    else:
        sys.exit('ERROR - connectivity='+str(connectivity)+' is not a valid selection, valid selections are 4/8')
    neighborcount=convolve(mask.astype('int'),kernel,mode='same')
    edgemask=(neighborcount<connectivity)*mask
    if latitudes is not None and longitudes is not None:
        edgecells=np.vstack((latitudes[edgemask.nonzero()[0]],longitudes[edgemask.nonzero()[1]])).T
        return edgemask, edgecells
    else:
        return edgemask

def maskvertexpoints(mask,latitudes,longitudes):
    """
    Returns the vertex points of all cells in the input mask, and the vertex points of only the mask edge

    Inputs (Required):
        mask - A 2-d boolean numpy array of dimension MxN where M=len(latitudes) and N=len(longitudes)
        latitudes - An evenly spaced numpy array of latitude points (degrees)
        longitudes - An evenly spaced numpy array of longitude points (degrees)

    Outputs:
        vertexpoints - A 2-d Nx2 numpy array of latitude/longitude points (degrees) of all vertices of mask cells
        edgevertexpoints - A 2-d Nx2 numpy array of latitude/longitude points (degrees) of all vertices of cells at the mask edge (8-connected)
    """
    checkmask(mask,latitudes,longitudes)
    vertexmask=np.full((tuple(np.array(mask.shape)+1)),0)
    maskint=mask.astype('int')
    vertexmask[:-1,:-1]+=maskint
    vertexmask[:-1,1:]+=maskint
    vertexmask[1:,:-1]+=maskint
    vertexmask[1:,1:]+=maskint
    latspc=gridspacing(latitudes)
    lonspc=gridspacing(longitudes)
    latdir=checklatitudedirection(latitudes)
    if latdir=='inc':
        vertexlatitudes=np.append(latitudes-latspc/2,latitudes[-1]+latspc/2)
    elif latdir=='dec':
        vertexlatitudes=np.append(latitudes+latspc/2,latitudes[-1]-latspc/2)
    vertexlongitudes=np.append(longitudes-lonspc/2,longitudes[-1]+lonspc/2)
    vertexpoints=np.vstack((vertexlatitudes[vertexmask.nonzero()[0]],vertexlongitudes[vertexmask.nonzero()[1]])).T
    edgevertexmask=((vertexmask<4)*(vertexmask>0))
    edgevertexpoints=np.vstack((vertexlatitudes[edgevertexmask.nonzero()[0]],vertexlongitudes[edgevertexmask.nonzero()[1]])).T
    return vertexpoints, edgevertexpoints

def findneighbors(cell,connectivity=8,direction='cw'):
    """
    Returns the neighbors of a cell, with selected connectivity and direction

    Inputs (Required):
        cell - A 1x2 numpy array describing the index of the cell

    Inputs (Optional):
        connectivity (4/8) - The connectivity parameter for testing edge cells, default=8
            4 tests only lateral neighbors while 8 tests lateral and diagonal neighbors
        direction ('cw'/'ccw') - A string selecting the direction of the returned neighbors, default='cw'

    Outputs:
        neighbors - An 8x2 or 4x2 numpy array of the neighboring cell indices for the input cell 
    """
    if connectivity==8:
        neighborrange=np.array([[-1,0],[-1,1],[0,1],[1,1],[1,0],[1,-1],[0,-1],[-1,-1]])
    elif connectivity==4:
        neighborrange=np.array([[-1,0],[0,1],[1,0],[0,-1]])
    else:
        sys.exit('ERROR - connectivity='+str(connectivity)+' is not a valid selection, valid selections are 4/8')
    if direction=='cw':
        pass
    elif direction=='ccw':
        neighborrange=neighborrange[::-1,:]
    else:
        sys.exit('ERROR - direction=\''+direction+'\' is not a valid selection, valid selections are \'cw\'/\'ccw\'')
    neighbors=cell+neighborrange
    return neighbors

def checkconnectivity(mask,checkcells='full',connectivity=8):
    """
    Returns whether a mask or its inverse are connected

    Inputs (Required):
        mask - A 2-d boolean numpy array of dimension MxN where M=len(latitudes) and N=len(longitudes)

    Inputs (Optional):
        checkcells ('full'/'empty') - A string selecting the mask cells to test (True/False), default='full'
        connectivity (4/8) - The connectivity parameter for testing edge cells, default=8
            4 tests only lateral neighbors while 8 tests lateral and diagonal neighbors
    
    Outputs:
        connected (True/False) - A boolean describing whether the input mask is connected under the input conditions
    """
    checkmask(mask)
    if checkcells=='full':
        maskcells=np.argwhere(mask==True)
    elif checkcells=='empty':
        maskcells=np.argwhere(mask==False)
    else:
        sys.exit('ERROR checkcells=\''+checkcells+'\' is not a valid selection, valid selections are \'full\'/\'empty\'')
    startcell=maskcells[0]
    checked=[]
    tocheck=[]
    tocheck.append(startcell)
    while len(tocheck)>0:
        checked.append(tocheck[0])
        cellneighbors=findneighbors(tocheck.pop(0),connectivity=connectivity)
        for k in cellneighbors:
            if (k==maskcells).all(axis=1).any() and not any(np.array_equal(k,x) for x in checked) and not any(np.array_equal(k,x) for x in tocheck):
                tocheck.append(k)
    if len(checked)==len(maskcells):
        connected=True
    else:
        connected=False
    return connected

###Contour Functions

def findstartcell(mask,searchdir='ru'):
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
    checkmask(mask)
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

def parsestartinput(start,buffermask):
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

def setstopfunction(stop,startvisits,startcell,startorientation):
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

def cleancontour(contourcells,searchcells,latitudes=None,longitudes=None,closecontour=True,remcontourrepeat=True,remsearchrepeat=False):
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
        latspc=gridspacing(latitudes)
        lonspc=gridspacing(longitudes)
        latdir=checklatitudedirection(latitudes)
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

def contourtracesquare(mask,latitudes=None,longitudes=None,direction='cw',start='auto',stop='either',startvisits=3,checkconn=False,remcontourrepeat=True,remsearchrepeat=False,closecontour=True):
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
        checkmask(mask,latitudes,longitudes)
    else:
        checkmask(mask)
    buffermask=np.full((tuple(np.array(mask.shape)+2)),False)
    buffermask[1:-1,1:-1]=mask
    if checkconn:
        fullconnectivity=checkconnectivity(buffermask,checkcells='full',connectivity=4)
        emptyconnectivity=checkconnectivity(buffermask,checkcells='empty',connectivity=4)
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
        startcell,startorientation=parsestartinput(start,buffermask)
    elif start=='auto':
        startcell,startorientation=findstartcell(buffermask,searchdir=searchdir)
    else:
        sys.exit('ERROR - start=\''+start+'\' is not a valid selection, valid selections are \'auto\' or a 2x2 array describing start cell index and start orientation')
    checkbreak=setstopfunction(stop,startvisits,startcell,startorientation)
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
    contour,contoursearch=cleancontour(contourcells,searchcells,latitudes=latitudes,longitudes=longitudes,closecontour=closecontour,remcontourrepeat=remcontourrepeat,remsearchrepeat=remsearchrepeat)
    return contour,contoursearch

def contourtracemooreneighbor(mask,latitudes=None,longitudes=None,direction='cw',start='auto',stop='either',startvisits=3,remcontourrepeat=True,remsearchrepeat=False,closecontour=True):
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
        checkmask(mask,latitudes,longitudes)
    else:
        checkmask(mask)
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
        startcell,startorientation=parsestartinput(start,buffermask)
    elif start=='auto':
        if stop=='either' or stop=='Elisoff':
            for ct,searchdir in enumerate(searchdirstrings):
                startcell,startorientation=findstartcell(buffermask,searchdir=searchdir)
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
            startcell,startorientation=findstartcell(buffermask,searchdir='ru')
    else:
        sys.exit('ERROR - start=\''+start+'\' is not a valid selection, valid selections are \'auto\' or a 2x2 array describing start cell index and start orientation')
    checkbreak=setstopfunction(stop,startvisits,startcell,startorientation)
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
            neighbors=findneighbors(cell,connectivity=8,direction=direction)
            nextneighborindex=(np.nonzero(((cell-orientation)==neighbors).all(axis=1))[0][0]+1)%8
        else:
            nextneighborindex=(nextneighborindex+1)%8
        orientation=neighbors[nextneighborindex]-cell
        cell=neighbors[nextneighborindex]
        if (cell==startcell).all():
            Nvisits+=1
        breakloop=checkbreak(cell,orientation,Nvisits)
    contour,contoursearch=cleancontour(contourcells,searchcells,latitudes=latitudes,longitudes=longitudes,closecontour=closecontour,remcontourrepeat=remcontourrepeat,remsearchrepeat=remsearchrepeat)
    return contour,contoursearch

def contourtraceimprovedmooreneighbor(mask,latitudes=None,longitudes=None,direction='cw',start='auto',stop='either',startvisits=3,remcontourrepeat=True,remsearchrepeat=False,closecontour=True):
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
        checkmask(mask,latitudes,longitudes)
    else:
        checkmask(mask)
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
        startcell,startorientation=parsestartinput(start,buffermask)
    elif start=='auto':
        if stop=='either' or stop=='Elisoff':
            for ct,searchdir in enumerate(searchdirstrings):
                startcell,startorientation=findstartcell(buffermask,searchdir=searchdir)
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
            startcell,startorientation=findstartcell(buffermask,searchdir='ru')
    else:
        sys.exit('ERROR - start=\''+start+'\' is not a valid selection, valid selections are \'auto\' or a 2x2 array describing start cell index and start orientation')
    checkbreak=setstopfunction(stop,startvisits,startcell,startorientation)
    cell=startcell
    orientation=startorientation
    searchcells=[]
    contourcells=[]
    searchcells.append(cell)
    contourcells.append(cell)
    neighbors=findneighbors(cell,connectivity=8,direction=direction)
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
            neighbors=findneighbors(cell,connectivity=8,direction=direction)
            nextneighborindex=(np.nonzero(((cell-orientation)==neighbors).all(axis=1))[0][0]+1)%8
        else:
            nextneighborindex=(nextneighborindex+1)%8
        if (cell==startcell).all():
            Nvisits+=1
        breakloop=checkbreak(cell,orientation,Nvisits)
    contour,contoursearch=cleancontour(contourcells,searchcells,latitudes=latitudes,longitudes=longitudes,closecontour=closecontour,remcontourrepeat=remcontourrepeat,remsearchrepeat=remsearchrepeat)
    return contour,contoursearch

def contourtracepavlidis(mask,latitudes=None,longitudes=None,direction='cw',start='auto',stop='Nvisits',startvisits=1,remcontourrepeat=True,remsearchrepeat=False,closecontour=True):
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
        checkmask(mask,latitudes,longitudes)
    else:
        checkmask(mask)
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
        startcell,startorientation=parsestartinput(start,buffermask)
    elif start=='auto':
        for ct,searchdir in enumerate(searchdirstrings):
            startcell,startorientation=findstartcell(buffermask,searchdir=searchdir)
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
    checkbreak=setstopfunction(stop,startvisits,startcell,startorientation)
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
    contour,contoursearch=cleancontour(contourcells,searchcells,latitudes=latitudes,longitudes=longitudes,closecontour=closecontour,remcontourrepeat=remcontourrepeat,remsearchrepeat=remsearchrepeat)
    return contour,contoursearch

def contourtraceimprovedpavlidis(mask,latitudes=None,longitudes=None,direction='cw',start='auto',stop='Nvisits',startvisits=1,remcontourrepeat=True,remsearchrepeat=False,closecontour=True):
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
        checkmask(mask,latitudes,longitudes)
    else:
        checkmask(mask)
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
        startcell,startorientation=parsestartinput(start,buffermask)
    elif start=='auto':
        for ct,searchdir in enumerate(searchdirstrings):
            startcell,startorientation=findstartcell(buffermask,searchdir=searchdir)
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
    checkbreak=setstopfunction(stop,startvisits,startcell,startorientation)
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
    contour,contoursearch=cleancontour(contourcells,searchcells,latitudes=latitudes,longitudes=longitudes,closecontour=closecontour,remcontourrepeat=remcontourrepeat,remsearchrepeat=remsearchrepeat)
    return contour,contoursearch

def contourtracefastrepresentative(mask,latitudes=None,longitudes=None,direction='cw',start='auto',stop='either',startvisits=4,remcontourrepeat=True,remsearchrepeat=False,closecontour=True):
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
        checkmask(mask,latitudes,longitudes)
    else:
        checkmask(mask)
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
        startcell,startorientation=parsestartinput(start,buffermask)
    elif start=='auto':
        startcell,startorientation=findstartcell(buffermask,searchdir=searchdir)
    else:
        sys.exit('ERROR - start=\''+start+'\' is not a valid selection, valid selections are \'auto\' or a 2x2 array describing start cell index and start orientation')
    checkbreak=setstopfunction(stop,startvisits,startcell,startorientation)
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
    contour,contoursearch=cleancontour(contourcells,searchcells,latitudes=latitudes,longitudes=longitudes,closecontour=closecontour,remcontourrepeat=remcontourrepeat,remsearchrepeat=remsearchrepeat)
    return contour,contoursearch

###Geocontour Functions

def geocontour(contour,latitudes,longitudes,connecttype='cell',simplify=False):
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
    checkcontour(contour,latitudes,longitudes)
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
    lonlengths=longitudelength(geocontour[:,0,0])
    latlengths=latitudelength(geocontour[:,0,0])
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

###Plot Functions

def plotdatasize(axobj=None,mult=1,axis='y',plottype='line'):
    """ 
    Returns a value for linewidth/markersize or size that is scaled to plotted data units 
        Note that this function must be used after any figure/axis alterations (axes limits, aspect ratio, etc.)
        Note that scaling is to data units and not axis size - if a y axis ranges from 0 to 3 and mult=1, the produced line will be 1 unit thick, or 1/3 the height of the y axis

    Inputs (Required):
        none

    Inputs (Optional):
        axobj - axes object to scale to, default=plt.gca() (current pyplot axes object)
        mult - A float multiplier for output, default=1
            e.g. to get a linewidth 1/10 the scale of the plotted units, mult=0.1
        axis ('x','y','xy') - axis to use for data scaling, default='y'
            Note that if x and y axes do not have an equal aspect ratio (e.g. axobj.set_aspect('equal')), 'xy' will attempt to average the scaling retrieved by 'x' and 'y', and will produce a warning if unequal
        plottype ('line'/'scatter') - plot type for use with output, default='line'
            Note that for lines and line markers (pyplot.plot) linewidth and markersize scale linearly, while for s (pyplot.scatter), markers scale with the square of size
            Note that using output for s will be correct in a scatter only if linewidth=0 (marker edge)

    Outputs:
        plotdatasize - linewidth/markersize/s to be used as a direct input to pyplot.plot, pyplot.scatter, etc.        
    """
    if axobj is None:
        axobj=plt.gca()
    fig=axobj.get_figure()
    pointwidth=72*fig.bbox_inches.width*axobj.get_position().width
    pointheight=72*fig.bbox_inches.height*axobj.get_position().height
    widthrange=np.diff(axobj.get_xlim())[0]
    heightrange=np.diff(axobj.get_ylim())[0]
    plotdatawidth=pointwidth/widthrange
    plotdataheight=pointheight/heightrange
    if axis=='y':
        plotdatasize=plotdataheight
    elif axis=='x':
        plotdatasize=plotdatawidth
    elif axis=='xy':
        if plotdataheight!=plotdatawidth:
            if not np.allclose(plotdataheight,plotdatawidth):
                warnings.warn('WARNING - x and y axes not scaled equally, output will not be precise')
        plotdatasize=(plotdataheight+plotdatawidth)/2
    else:
        sys.exit('ERROR - \''+axis+'\' is invalid input, select \'x\',\'y\', or \'xy\'')
    if plottype=='line':
        return mult*plotdatasize
    elif plottype=='scatter':
        return (mult*plotdatasize)**2
    else:
        sys.exit('ERROR - \''+plottype+'\' is invalid input, select \'line\' or \'scatter\'')

def plot(latitudes,longitudes,boundary=None,mask=None,contour=None,contoursearch=None,geocontour=None,vertices=None,boundingbox='all',buffer='off',grid='on',cells='default',showcontour='on',startcell='on',contourarrows='on',contoursearcharrows='on',geocontourvectors='on',emptycellcolor='lightgrey',fullcellcolor='sandybrown',boundarycolor='tab:blue',contourcolor='olivedrab',contoursearchcolor='firebrick',geocontourcolor='olivedrab',vertexcolor='tab:cyan',gridcolor='black',lw_boundary='auto',lw_contour='auto',lw_contoursearch='auto',lw_geocontour='auto',mw_contourarrows='auto',mw_contoursearcharrows='auto',mw_vertices='auto',features=None,title=None,outname='plot',outdpi='auto'):
    """
    Plots any/all maskpy-created elements: boundary, mask, contour, contoursearch, geocontour, vertices

    Inputs (Required):
        latitudes - An evenly spaced numpy array of latitude points (degrees)
        longitudes - An evenly spaced numpy array of longitude points (degrees)

    Inputs (Optional):
        boundary - A 2-d Nx2 numpy array of latitude/longitude points (degrees)
        mask - A 2-d boolean numpy array of dimension MxN where M=len(latitudes) and N=len(longitudes)
        contour - A 2-d Nx2 numpy array of ordered latitude/longitude points (degrees) describing the contour trace of a mask
        contoursearch - A 2-d Nx2 numpy array of ordered latitude/longitude points (degrees) describing the cells searched during contour tracing
        geocontour - A 3-d Nx2x5 numpy array defining a list of N contour cells and their edge points, lengths, and outward unit vectors
        vertices - A 2-d Nx2 numpy array of latitude/longitude points (degrees)
        boundingbox ('all'/'boundary'/'mask'/'contour'/'contoursearch'/'geocontour') - A string denoting the plot element(s) to be used as x and y axis limits, default='all' (uses the largest bounds from all provided elements)
        buffer ('on'/'off') - A string for creating a 1 cell buffer on the bounding box, default='off'
        grid ('on'/'off') - A string for showing/hiding the cell grid, default='on'
        cells ('default'/'none'/'mask'/'maskedge-4'/'maskedge-8'/'contour'/'geocontour') - A string denoting the plot element to be used for filling the cell grid, default='default'
            'default' uses first provided element, in this order: 'geocontour','contour','mask','none'
            'none' - no filled grid cells
            'mask' - all mask cells
            'maskedge-4' - 4-connected mask edge cells
            'maskedge-8' - 8-connected mask edge cells
            'contour' - contour cells
            'geocontour' - geocontour cells
        showcontour ('on'/'off') - A string for showing/hiding the contour, default='on'
            allows showing contour cells without line plot of contour
        startcell ('on'/'off') - A string for showing/hiding the startcell (contour or contoursearch), default='on'
        contourarrows ('on'/'off') - A string for showing/hiding directional arrows on the contour, default='on'
        contoursearcharrows ('on'/'off') - A string for showing/hiding directional arrows on the contoursearch, default='on'
        geocontourvectors ('on'/'off') - A string for showing/hiding outward normal vectors on the geocontour, default='on'
        Colors: all accept any matplotlib predefined, hex, or rgba array
            emptycellcolor - color for unmasked cells, default='lightgray'
            fullcellcolor - color for masked cells, default='sandybrown'
            boundarycolor - color for boundary, default='tab:blue'
            contourcolor - color for contour, default='olivedrab'
            contoursearchcolor - color for contoursearch, default='firebrickred'
            geocontourcolor - color for geocontour, default='olivedrab'
            vertexcolor - color for vertices, default='tab:cyan'
            gridcolor - color for grid, default='black'
        linewidths/markerwidths: a float setting fraction of cell covered by lines/markers (e.g. lw=0.5 means lines will be half as wide as grid cells, while mw=2 means markers will be as wide as 2 grid cells)
            lw_boundary - boundary linewidth, default='auto' (0.1)
            lw_contour - contour linewidth, default='auto' (0.1)
            lw_contoursearch - contoursearch linewidth, default='auto' (0.1)
            lw_geocontour - geocontour linewidth, default='auto' (0.1)
            mw_contourarrows - contour arrow markerwidth, default='auto' (0.5)
            mw_contoursearcharrows - contoursearch arrow markerwidth, default='auto' (0.5)
            mw_vertices - vertex markerwidth, default='auto' (0.4)
        features (None/'natural'/'borders') - display Earth features (if cartopy is installed), default=None
            'natural' displays coastlines, ocean, and lakes/rivers
            'borders' displays national and state/province level boundaries
        title - A string used as the plot title, default=None
        outname - A string used as the filename/path for the saved image, default='plot'
        outdpi - dpi (resolution) of the saved image, default='auto'
            'auto' scales dpi high enough to see features when zooming into a single grid cell, floor of 150
            for very large grids, 'auto' may set dpi high enough that pyplot will hang on some systems - setting dpi manually can avoid this if encountered

    Outputs:
        none

    Examples of common use cases:
    1) Plot a mask and the boundary used to create it:
        plot(latitudes,longitudes,boundary=<boundary>,mask=<mask>,lw_boundary=0.2)
    2) Plot a contour
        plot(latitudes,longitudes,countour=<countour>,boundingbox='contour',buffer='on')
    3) Plot a contoursearch overlaying contour cells
        plot(latitudes,longitudes,contour=<contour>,contoursearch=<contoursearch>,showcontour='off',cells='contour',boundingbox='contoursearch',buffer='on')
    4) Plot a geocontour overlaid onto a map projection with natural features
        plot(latitudes,longitudes,geocontour=<geocontour>,cells='geocontour',features='natural')
    """
    latspc=gridspacing(latitudes)
    lonspc=gridspacing(longitudes)
    gridlatmin=latitudes.min()-latspc/2
    gridlatmax=latitudes.max()+latspc/2
    gridlonmin=longitudes.min()-lonspc/2
    gridlonmax=longitudes.max()+lonspc/2
    latdir=checklatitudedirection(latitudes)
    boundingarray=np.array([[[gridlatmin,gridlatmax],[gridlonmin,gridlonmax]]])
    if boundary is None:
        if boundingbox=='boundary':
            sys.exit('ERROR - Can not use \'boundary\' for boundingbox if no boundary input provided')
    else:
        checkboundary(boundary)
        loninput=checklongituderange(longitudes)
        boundloninput=checklongituderange(boundary[:,1])
        if (loninput=='neg' and boundloninput=='pos') or (loninput=='pos' and boundloninput=='neg'):
            sys.exit('ERROR - Longitude input range is '+loninput+' and boundary longitude range is '+boundloninput)
        boundlatmin=np.floor((boundary.min(axis=0)[0]-gridlatmin)/latspc)*latspc+gridlatmin
        boundlatmax=np.ceil((boundary.max(axis=0)[0]-gridlatmax)/latspc)*latspc+gridlatmax
        boundlonmin=np.floor((boundary.min(axis=0)[1]-gridlonmin)/lonspc)*lonspc+gridlonmin
        boundlonmax=np.ceil((boundary.max(axis=0)[1]-gridlonmax)/lonspc)*lonspc+gridlonmax
        boundingarray=np.append(boundingarray,[[[boundlatmin,boundlatmax],[boundlonmin,boundlonmax]]],axis=0)
    if mask is None:
        if boundingbox=='mask':
            sys.exit('ERROR - Can not use \'mask\' for boundingbox if no mask input provided')
        if cells=='mask' or cells=='maskedge-8' or cells=='maskedge-4':
            sys.exit('ERROR - Can not use \''+cells+'\' for cells if no mask input provided')
    else:
        checkmask(mask,latitudes,longitudes)
        masklatitudes=latitudes[mask.sum(axis=1)>0]
        masklongitudes=longitudes[mask.sum(axis=0)>0]
        masklatmin=masklatitudes.min()-latspc/2
        masklatmax=masklatitudes.max()+latspc/2
        masklonmin=masklongitudes.min()-lonspc/2
        masklonmax=masklongitudes.max()+lonspc/2
        boundingarray=np.append(boundingarray,[[[masklatmin,masklatmax],[masklonmin,masklonmax]]],axis=0)
    if contour is None:
        if boundingbox=='contour':
            sys.exit('ERROR - Can not use \'contour\' for boundingbox if no contour input provided')
        if cells=='contour':
            sys.exit('ERROR - Can not use \'contour\' for cells if no contour input provided')
    else:
        checkcontour(contour,latitudes,longitudes)
        contlatmin=contour[:,0].min()-latspc/2
        contlatmax=contour[:,0].max()+latspc/2
        contlonmin=contour[:,1].min()-lonspc/2
        contlonmax=contour[:,1].max()+lonspc/2
        boundingarray=np.append(boundingarray,[[[contlatmin,contlatmax],[contlonmin,contlonmax]]],axis=0)
    if contoursearch is None:
        if boundingbox=='contoursearch':
            sys.exit('ERROR - Can not use \'contoursearch\' for boundingbox if no contoursearch input provided')
    else:
        contsrchlatmin=contoursearch[:,0].min()-latspc/2
        contsrchlatmax=contoursearch[:,0].max()+latspc/2
        contsrchlonmin=contoursearch[:,1].min()-lonspc/2
        contsrchlonmax=contoursearch[:,1].max()+lonspc/2
        boundingarray=np.append(boundingarray,[[[contsrchlatmin,contsrchlatmax],[contsrchlonmin,contsrchlonmax]]],axis=0)
    if geocontour is None:
        if boundingbox=='geocontour':
            sys.exit('ERROR - Can not use \'geocontour\' for boundingbox if no geocontour input provided')
        if cells=='geocontour':
            sys.exit('ERROR - Can not use \'geocontour\' for cells if no geocontour input provided')
    else:
        checkgeocontour(geocontour,latitudes,longitudes)
        geocontlatmin=geocontour[:,0,0].min()-latspc/2
        geocontlatmax=geocontour[:,0,0].max()+latspc/2
        geocontlonmin=geocontour[:,1,0].min()-lonspc/2
        geocontlonmax=geocontour[:,1,0].max()+lonspc/2
        boundingarray=np.append(boundingarray,[[[geocontlatmin,geocontlatmax],[geocontlonmin,geocontlonmax]]],axis=0)
    if boundingbox=='grid':
        ylimmin=gridlatmin
        ylimmax=gridlatmax
        xlimmin=gridlonmin
        xlimmax=gridlonmax
    elif boundingbox=='boundary':
        ylimmin=boundlatmin
        ylimmax=boundlatmax
        xlimmin=boundlonmin
        xlimmax=boundlonmax
    elif boundingbox=='mask':
        ylimmin=masklatmin
        ylimmax=masklatmax
        xlimmin=masklonmin
        xlimmax=masklonmax
    elif boundingbox=='contour':
        ylimmin=contlatmin
        ylimmax=contlatmax
        xlimmin=contlonmin
        xlimmax=contlonmax
    elif boundingbox=='contoursearch':
        ylimmin=contsrchlatmin
        ylimmax=contsrchlatmax
        xlimmin=contsrchlonmin
        xlimmax=contsrchlonmax
    elif boundingbox=='geocontour':
        ylimmin=geocontlatmin
        ylimmax=geocontlatmax
        xlimmin=geocontlonmin
        xlimmax=geocontlonmax
    elif boundingbox=='all':
        ylimmin=boundingarray.min(axis=0)[0,0]
        ylimmax=boundingarray.max(axis=0)[0,1]
        xlimmin=boundingarray.min(axis=0)[1,0]
        xlimmax=boundingarray.max(axis=0)[1,1]
    else:
        sys.exit('ERROR - boundingbox=\''+boundingbox+'\' is not a valid selection, valid selections are \'grid\'/\'boundary\'/\'mask\'/\'contour\'/\'contoursearch\'/\'geocontour\'/\'all\'')
    if buffer=='on':
        ylimmin-=latspc
        ylimmax+=latspc
        xlimmin-=lonspc
        xlimmax+=lonspc
    if cells=='default':
        if geocontour is not None:
            cells='geocontour'
        elif contour is not None:
            cells='contour'
        elif mask is not None:
            cells='mask'
        else:
            cells='none'
    if cells=='none':
        pltmask=np.full((len(latitudes),len(longitudes)),0)
    elif cells=='mask':
        pltmask=mask.astype('int')
    elif cells=='maskedge-8':
        pltmask=maskedgecells(mask,connectivity=8).astype('int')
    elif cells=='maskedge-4':
        pltmask=maskedgecells(mask,connectivity=4).astype('int')
    elif cells=='contour':
         latinds=latitudes.argsort()[np.searchsorted(latitudes,contour[:,0],sorter=latitudes.argsort())]
         loninds=np.searchsorted(longitudes,contour[:,1])
         pltmask=np.full((len(latitudes),len(longitudes)),0)
         pltmask[latinds,loninds]=1
    elif cells=='geocontour':
         latinds=latitudes.argsort()[np.searchsorted(latitudes,geocontour[:,0,0],sorter=latitudes.argsort())]
         loninds=np.searchsorted(longitudes,geocontour[:,1,0])
         pltmask=np.full((len(latitudes),len(longitudes)),0)
         pltmask[latinds,loninds]=1
    else:
        sys.exit('ERROR - cells=\''+cells+'\' is not a valid selection, valid selections are \'none\'/\'mask\'/\'maskedge-8\'/\'maskedge-4\'/\'contour\'/\'geocontour\'')
    if latdir=='inc':
        org='lower'
    elif latdir=='dec':
        org='upper'
    ext=[gridlonmin,gridlonmax,gridlatmin,gridlatmax]
    cmp=mplc.ListedColormap([emptycellcolor,fullcellcolor])
    plt.ioff()
    yminor=np.arange(ylimmin,ylimmax+latspc,latspc)
    xminor=np.arange(xlimmin,xlimmax+lonspc,lonspc)
    ymajor=(yminor[:-1]+latspc/2)[::np.ceil(len(yminor)/7).astype('int')]
    xmajor=(xminor[:-1]+lonspc/2)[::np.ceil(len(xminor)/7).astype('int')]
    fig=plt.figure()
    if features is None:
        ax=fig.add_subplot(111)
    else:
        if cp_exists=='y':
            PRO=cp.crs.PlateCarree()
            ax=fig.add_subplot(111,projection=PRO)
        else:
            sys.exit('ERROR - Could not import cartopy, features are only usable with cartopy')
    if title is not None:
        ax.set_title(title)
    ax.set_xlabel('Longitude (deg)')
    ax.set_ylabel('Latitude (deg)')
    ax.set_ylim((ylimmin,ylimmax))
    ax.set_xlim((xlimmin,xlimmax)) 
    ax.imshow(pltmask,aspect='equal',interpolation='none',extent=ext,origin=org,cmap=cmp,zorder=-1)
    pdw=plotdatasize(ax,axis='xy',mult=(latspc+lonspc)/2,plottype='line')
    pds=plotdatasize(ax,axis='xy',mult=(latspc+lonspc)/2,plottype='scatter')
    if features is None:
        ax.set_yticks(ymajor,major=True)
        ax.set_xticks(xmajor,major=True)
        ax.set_yticks(yminor,minor=True)
        ax.set_xticks(xminor,minor=True)
        if grid=='on':
            ax.grid(which='minor',linestyle=(0,(1,1)),color=gridcolor,linewidth=0.05*pdw,zorder=3)
    else:
        ax.set_yticks(ymajor)
        ax.set_xticks(xmajor)
        if grid=='on':
            gl=ax.gridlines(xlocs=xminor,ylocs=yminor,linestyle=(0,(1,1)),color=gridcolor,linewidth=0.05*pdw,zorder=5)
        if features=='natural':
            water=np.array([0.7,0.75,0.95])
            ax.add_feature(cp.feature.NaturalEarthFeature(category='physical',name='ocean',scale='10m'),color=water,edgecolor='none',linewidth=0,alpha=0.5,zorder=6)
            ax.add_feature(cp.feature.NaturalEarthFeature('physical','land','10m'),facecolor='None',edgecolor='black',alpha=0.7,linewidth=0.01*pdw,zorder=7)
            ax.add_feature(cp.feature.NaturalEarthFeature('physical','lakes','10m'),facecolor=water,edgecolor='black',linewidth=0.01*pdw,alpha=0.7,zorder=8)
            ax.add_feature(cp.feature.NaturalEarthFeature('physical','rivers_lake_centerlines','10m'),facecolor='none',edgecolor=water,linewidth=0.05*pdw,alpha=0.7,zorder=9)
        elif features=='borders':
            ax.add_feature(cp.feature.NaturalEarthFeature('cultural','admin_0_countries','10m'),facecolor='None',edgecolor='black',alpha=0.7,linewidth=0.1*pdw,zorder=6)
            ax.add_feature(cp.feature.NaturalEarthFeature('cultural','admin_1_states_provinces_lines','10m'),facecolor='None',edgecolor='black',alpha=0.7,linewidth=0.1*pdw,zorder=6)
        else:
            sys.exit('ERROR - features=\''+features+'\' is not a valid selection, valid selections are \'natural\'/\'borders\'')
    if boundary is not None:
        if lw_boundary=='auto':
            lw_boundary=0.1*pdw
        else:
            lw_boundary=lw_boundary*pdw
        ax.plot(boundary[:,1],boundary[:,0],color=boundarycolor,linewidth=lw_boundary,zorder=10)
    if vertices is not None:
        if mw_vertices=='auto':
            mw_vertices=pds*0.4**2
        else:
            mw_vertices=pds*mw_vertices**2
        ax.scatter(vertices[:,1],vertices[:,0],c=vertexcolor,s=mw_vertices,linewidth=0,zorder=10)
    if contoursearch is not None:
        if lw_contoursearch=='auto':
            lw_contoursearch=0.1*pdw
        else:
            lw_contoursearch=lw_contoursearch*pdw
        ax.plot(contoursearch[:,1],contoursearch[:,0],color=contoursearchcolor,linewidth=lw_contoursearch,zorder=11)
        if mw_contoursearcharrows=='auto':
            mw_contoursearcharrows=pds*0.5**2
        else:
            mw_contoursearcharrows=pds*mw_contoursearcharrows**2
        if startcell=='on':
            ax.scatter(contoursearch[0,1],contoursearch[0,0],marker='o',c=contoursearchcolor,s=mw_contoursearcharrows,linewidth=0,zorder=11)
        if contoursearcharrows=='on':
            contoursearchdiff=np.diff(contoursearch,axis=0)
            contoursearchpointrotation=-180/np.pi*np.arctan2(contoursearchdiff[:,1],contoursearchdiff[:,0])
            contoursearchpointlocation=contoursearch[:-1]+2*contoursearchdiff/5
            for ct,k in enumerate(contoursearchpointlocation):
                if (contoursearchdiff[ct,:]==0).all():
                    ax.scatter(k[1],k[0],marker=(4,0,0),c=contoursearchcolor,s=mw_contoursearcharrows,linewidth=0,zorder=11)
                else:
                    ax.scatter(k[1],k[0],marker=(3,0,contoursearchpointrotation[ct]),c=contoursearchcolor,s=mw_contoursearcharrows,linewidth=0,zorder=11)
    if contour is not None:
        if showcontour=='on':
            if lw_contour=='auto':
                lw_contour=0.1*pdw
            else:
                lw_contour=lw_contour*pdw
            ax.plot(contour[:,1],contour[:,0],color=contourcolor,linewidth=lw_contour,zorder=12)
            if mw_contourarrows=='auto':
                mw_contourarrows=pds*0.5**2
            else:
                mw_contourarrows=pds*mw_contourarrows**2
            if startcell=='on':
                ax.scatter(contour[0,1],contour[0,0],marker='o',c=contourcolor,s=mw_contourarrows,linewidth=0,zorder=12)
            if contourarrows=='on':
                contourdiff=np.diff(contour,axis=0)
                contourpointrotation=-180/np.pi*np.arctan2(contourdiff[:,1],contourdiff[:,0])
                contourpointlocation=contour[:-1]+2*contourdiff/5
                for ct,k in enumerate(contourpointlocation):
                    ax.scatter(k[1],k[0],marker=(3,0,contourpointrotation[ct]),c=contourcolor,s=mw_contourarrows,linewidth=0,zorder=12)
    if geocontour is not None:
        if lw_geocontour=='auto':
            lw_geocontour=0.1*pdw
        else:
            lw_geocontour=lw_geocontour*pdw
        for k in geocontour:
            ax.plot(k[1,1:3],k[0,1:3],color=geocontourcolor,linewidth=lw_geocontour,zorder=14)
        if geocontourvectors=='on':
            quivlocs=geocontour[:,:,1]+(geocontour[:,:,2]-geocontour[:,:,1])/2
            ax.quiver(quivlocs[:,1],quivlocs[:,0],geocontour[:,1,4],geocontour[:,0,4],facecolor=geocontourcolor,edgecolor='black',scale=3/(latspc+lonspc),scale_units='xy',width=0.15*(latspc+lonspc)/2,headwidth=3,headlength=3,headaxislength=2.75,units='xy',zorder=13,linewidth=pdw/150)
    if outdpi=='auto':
        outdpi=36*72/pdw
        if outdpi<150:
            outdpi=150
    fig.savefig(outname,dpi=outdpi,bbox_inches='tight')
    plt.close(fig)

###Save Functions

def save(latitudes,longitudes,boundary=None,mask=None,contour=None,contoursearch=None,geocontour=None,vertices=None,outname='save',outtype='np',maskouttxt='off',outformat='%8.3f',outdelim=' '):
    """
    Saves any/all maskpy-created elements: boundary, mask, contour, contoursearch, geocontour, vertices

    Inputs (Required):
        latitudes - An evenly spaced numpy array of latitude points (degrees)
        longitudes - An evenly spaced numpy array of longitude points (degrees)

    Inputs (Optional):
        boundary - A 2-d Nx2 numpy array of latitude/longitude points (degrees)
        mask - A 2-d boolean numpy array of dimension MxN where M=len(latitudes) and N=len(longitudes)
        contour - A 2-d Nx2 numpy array of ordered latitude/longitude points (degrees) describing the contour trace of a mask
        contoursearch - A 2-d Nx2 numpy array of ordered latitude/longitude points (degrees) describing the cells searched during contour tracing
        geocontour - A 3-d Nx2x5 numpy array defining a list of N contour cells and their edge points, lengths, and outward unit vectors
        outname - A string used as the filename/path for the saved elements, default='save'
        outtype ('np'/'xyz') - A string for selecting output filetype, default='np'
            'np': A numpy binary file that stores each element as a separate array/object
            'xyz': An xyz format text file (lat, lon, [data])
        maskouttxt ('on'/'off') - A string for choosing to also save a mask as a (1/0) text file, default='off'
        outformat - Formatting string for latitude and longitude values (in .xyz file), default='%8.3f'
        outdelim - Delimiter string for columns in output (in .xyz file), default=' '

    Outputs:
        none
    """
    if outtype=='xyz':
        np.savetxt(outname+'_latitudes.xyz',latitudes,fmt=outformat,delimiter=outdelim)
        np.savetxt(outname+'_longitudes.xyz',longitudes,fmt=outformat,delimiter=outdelim)
        if boundary is not None:
            checkboundary(boundary)
            np.savetxt(outname+'_boundary.xyz',boundary,fmt=outformat,delimiter=outdelim)
        if mask is not None:
            checkmask(mask,latitudes,longitudes)
            ysp, xsp = np.meshgrid(latitudes,longitudes,indexing='ij')
            xyzout=np.hstack((xsp.reshape((-1,1)),ysp.reshape((-1,1)),mask.reshape((-1,1))))
            np.savetxt(outname+'_mask.xyz',xyzout,fmt=[outformat, outformat, '%3d'],delimiter=outdelim)
            if maskouttxt=='on':
                np.savetxt(outname+'_mask.txt',mask.astype('int'),fmt='%1d',delimiter=outdelim)
        if contour is not None:
            checkcontour(contour,latitudes,longitudes)
            np.savetxt(outname+'_contour.xyz',contour,fmt=outformat,delimiter=outdelim)
        if contoursearch is not None:
            np.savetxt(outname+'_contoursearch.xyz',contoursearch,fmt=outformat,delimiter=outdelim)
        if geocontour is not None:
            checkgeocontour(geocontour,latitudes,longitudes)
            header1='cell'.center(16)+'entry'.center(19)+'exit'.center(18)+'length'.center(25)+'normal vector'.center(25)
            header2='lat'.center(8)+'lon'.center(9)+'lat'.center(9)+'lon'.center(9)+'lat'.center(9)+'lon'.center(9)+'degrees'.center(13)+'meters'.center(12)+'y'.center(13)+'x'.center(12)
            np.savetxt(outname+'_geocontour.xyz',geocontour.reshape((-1,10),order='F'),fmt=[outformat,outformat,outformat,outformat,outformat,outformat,'%12.7f','%12d','%12.7f','%12.7f'],delimiter=outdelim,header=header1+'\n'+header2)
        if vertices is not None:
            np.savetxt(outname+'_vertices.xyz',vertices,fmt=outformat,delimiter=outdelim)
    elif outtype=='np':
        np.savez(outname,latitudes=latitudes,longitudes=longitudes,boundary=boundary,mask=mask,contour=contour,contoursearch=contoursearch,geocontour=geocontour,vertices=vertices)
    else:
        sys.exit('ERROR - outtype=\''+outtype+'\' is not a valid selection, valid selections are \'np\'/\'xyz\'')

