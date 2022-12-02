import sys
import numpy as np
import geocontour.check as gcc

def spacing(dimension):
    """
    Returns the grid spacing for a given input dimension

    Inputs (Required):
        dimension - An evenly spaced numpy array of longitude or latitude points (degrees)

    Outputs:
        spacing - A positive float indicating the spacing of the input dimension (degrees)
    """
    gcc.cdim(dimension)
    spacing=abs(np.diff(dimension)[0])
    return spacing

def lonlens(latitudes,lonspacing=1):
    """
    Returns the lengths of a degree (default) of longitude over a range of latitudes
    Source: https://doi.org/10.5281/ZENODO.35392

    Inputs (Required):
        latitudes - An evenly spaced numpy array of latitude points (degrees)

    Inputs (Optional):
        lonspacing - A float/int specifying longitude spacing
            e.g. for lengths of a half degree of longitude at the given latitudes, lonspacing=0.5

    Outputs:
        lonlens - A numpy array of longitude lengths (m)
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
    Returns the grid lengths of a defined range of latitudes
    Source: https://doi.org/10.5281/ZENODO.35392

    Inputs (Required):
        latitudes - An evenly spaced numpy array of latitude points (degrees)

    Outputs:
        latlens - A numpy array of latitude lengths (m)
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
    Returns the length of a degree of longitude at the input latitude
    Source: https://doi.org/10.5281/ZENODO.35392

    Inputs (Required):
        latitude - A latitude point or array of points (degrees)

    Outputs:
        lonlen - The length of a degree of longitude at the input latitude(s) (m)
    """
    a=6378137
    b=6356752.3142
    e=np.sqrt((a**2-b**2)/a**2)
    latrad=latitude*np.pi/180
    lonlength=np.pi*a*np.cos(latrad)/(180*np.sqrt(1-(e*np.sin(latrad))**2))
    return lonlength

def latlen(latitude):
    """
    Returns the length of a degree of latitude at the input latitude
    Source: https://doi.org/10.5281/ZENODO.35392

    Inputs (Required):
        latitude - A latitude point or array of points (degrees)

    Outputs:
        latlen - The length of a degree of latitude at the input latitude(s) (m)
    """
    a=6378137
    b=6356752.3142
    e=np.sqrt((a**2-b**2)/a**2)
    latrad=latitude*np.pi/180
    latlength=np.pi*a*(1-e**2)/(180*(1-(e*np.sin(latrad))**2)**(3/2))
    return latlength

def areas(latitudes,longitudes,units=1):
    """
    Returns the cell areas of a grid defined by a range of latitudes and longitudes

    Inputs (Required):
        latitudes - An evenly spaced numpy array of latitude points (degrees)
        longitudes - An evenly spaced numpy array of longitude points (degrees)

    Inputs (Optional):
        units - unit multiplier for outputs, default=1 (m^2)
            e.g. for km^2, mult=1000m x 1000m=1e6

    Outputs:
        areas - A numpy array of areas, with dimension 0: latitude and dimension 1: longitude (default in m^2)
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

def clatdir(latitudes):
    """
    Returns a descriptor for the direction of a set of latitude points
        increasing or decreasing

    Inputs (Required):
        latitudes - A numpy array of latitude points (degrees)

    Outputs:
        latitudedirection ('inc'/'dec') - A descriptor for the direction of the input
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

def switchlon(longitudes,outrange,print_output='n'):
    """
    Returns a set of longitude points switched in place between negative (-180 to 180) and positive (0 to 360)

    Inputs (Required):
        longitudes - A numpy array of longitude points (degrees)
        outrange - ('pos'/'neg') Selection of range for longitude outputs: negative (-180 to 180) or positive (0 to 360)

    Inputs(Optional):
        print_output ('y'/'n') - flag to print out whether input was altered, default='n'

    Outputs:
        switchlon - A numpy array of longitudes, altered if not consistent with selected output range
    """
    if outrange!='pos' and outrange!='neg':
        sys.exit('ERROR - Unrecognized longitude output range. Valid inputs are \'pos\' or \'neg\'')
    inprange=clonrng(longitudes)
    if inprange==outrange or inprange=='ind':
        switchlon=longitudes
        if print_output=='y':
            print('Input already has values consistent with selected output range')
    elif inprange=='neg' and outrange=='pos':
        switchlon=longitudes+360*(longitudes<0)
        if print_output=='y':
            print('Switched negative input to positive output')
    elif inprange=='pos' and outrange=='neg':
        switchlon=longitudes-360*(longitudes>180)
        if print_output=='y':
            print('Switched positive input to negative output')
    else:
        sys.exit('ERROR - Function should not reach this point, check source')
    return switchlon

def switchind(longitudes):
    """
    Returns the index where a longitude array either crosses 0 or 180 degrees

    Inputs (Required):
        longitudes - An evenly spaced numpy array of longitude points (degrees)

    Outputs:
        switchind - An index for where array crosses 0 or 180 degrees, or 0 if no such point exists
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
