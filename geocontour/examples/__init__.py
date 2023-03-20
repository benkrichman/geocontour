"""
Functions serving as example use cases for geocontour
=====================================================
"""
import numpy as np
import geocontour.masksearch as gcms
import geocontour.contourtrace as gcct
from geocontour import build
import geocontour.output as gco
import importlib.resources as ilr

def small():
    """
    Run a small scale example of geocontour processing using mock data
    and save resulting plots to run directory

    1. find mask using area criteria (0.5) and plot boundary/mask
    2. trace contour using improved pavlidis algorithm and plot
       resultant contour and contour search
    3. compute geocontour from contour and plot, using simplify option

    Returns
    -------
    None

    See Also
    --------
    large

    Examples
    --------
    >>> import numpy as np
    >>> import geocontour.masksearch as gcms
    >>> import geocontour.contourtrace as gcct
    >>> from geocontour import build
    >>> import geocontour.output as gco
    >>> import importlib.resources as ilr
    >>> boundfil=ilr.files('geocontour').
    ... joinpath('data/generic_boundary.npz')
    >>> boundfil
    PosixPath('path/to/your/geocontour/geocontour/data/generic_boundary.npz')
    >>> data=np.load(boundfil)
    >>> latitudes=data['latitudes']
    >>> latitudes
    array([ 3.25,  2.75,  2.25,  1.75,  1.25,  0.75,  0.25, -0.25,
        -0.75, -1.25, -1.75, -2.25, -2.75, -3.25])
    >>> longitudes=data['longitudes']
    >>> longitudes
    array([-3.25, -2.75, -2.25, -1.75, -1.25, -0.75, -0.25,  0.25,
        0.75,  1.25,  1.75,  2.25,  2.75,  3.25])
    >>> boundary=data['boundary']
    >>> boundary
    array([[-1.3, -2.4],
       [ 1. ,  0. ],
       [ 4.9,  0. ],
       [ 1.9,  4.6],
       [-0.7,  2.5],
       [-1.1,  3.5],
       [-1.8,  3.5],
       [-1.8,  2.5],
       [-1.2,  2.5],
       [-1.2,  2.1],
       [-1.3, -2.4]])
    >>> mask=gcms.area(latitudes, longitudes, boundary)
    >>> contour,contoursearch=gcct.pavlidis_imp(mask, latitudes,
    ... longitudes)
    >>> geocontour=build(contour, latitudes, longitudes)
    >>> geocontour_simp=build(contour, latitudes, longitudes,
    ... simplify=True)
    >>> gco.plot(latitudes, longitudes, boundary=boundary, mask=mask,
    ... title='Example Mask and Boundary',
    ... outname='example_small_boundary+mask', outdpi='indep',
    ... transp=True)
    >>> gco.plot(latitudes, longitudes, mask=mask,
    ... contoursearch=contoursearch, title='Example Contour Search',
    ... outname='example_small_contoursearch', outdpi='indep',
    ... transp=True)
    >>> gco.plot(latitudes, longitudes, contour=contour,
    ... cells='contour', title='Example Contour',
    ... outname='example_small_contour', outdpi='indep', transp=True)
    >>> gco.plot(latitudes, longitudes, geocontour=geocontour,
    ... buffer=True, title='Example Geocontour',
    ... outname='example_small_geocontour', outdpi='indep', transp=True)
    >>> gco.plot(latitudes, longitudes, geocontour=geocontour_simp,
    ... buffer=True, title='Example Geocontour - Simplified',
    ... outname='example_small_geocontour_simp', outdpi='indep',
    ... transp=True)
    """
    print('\nRunning small scale geocontour example\nExample images will be output to current directory')
    boundfil=ilr.files('geocontour').joinpath('data/generic_boundary.npz')
    data=np.load(boundfil)
    latitudes=data['latitudes']
    longitudes=data['longitudes']
    boundary=data['boundary']
    mask=gcms.area(latitudes,longitudes,boundary)
    contour,contoursearch=gcct.pavlidis_imp(mask,latitudes,longitudes)
    geocontour=build(contour,latitudes,longitudes)
    geocontour_simp=build(contour,latitudes,longitudes,simplify=True)
    print('  Example figure saved as \'example_small_boundary+mask\'')
    gco.plot(latitudes,longitudes,boundary=boundary,mask=mask,title='Example Mask and Boundary',outname='example_small_boundary+mask',outdpi='indep',transp=True)
    print('  Example figure saved as \'example_small_contoursearch\'')
    gco.plot(latitudes,longitudes,mask=mask,contoursearch=contoursearch,title='Example Contour Search',outname='example_small_contoursearch',outdpi='indep',transp=True)
    print('  Example figure saved as \'example_small_contour\'')
    gco.plot(latitudes,longitudes,contour=contour,cells='contour',title='Example Contour',outname='example_small_contour',outdpi='indep',transp=True)
    print('  Example figure saved as \'example_small_geocontour\'')
    gco.plot(latitudes,longitudes,geocontour=geocontour,buffer=True,title='Example Geocontour',outname='example_small_geocontour',outdpi='indep',transp=True)
    print('  Example figure saved as \'example_small_geocontour_simp\'')
    gco.plot(latitudes,longitudes,geocontour=geocontour_simp,buffer=True,title='Example Geocontour - Simplified',outname='example_small_geocontour_simp',outdpi='indep',transp=True)

def large():
    """
    Run a large scale example of geocontour processing using the
    Mississippi River Basin boundary and save resulting plots to run
    directory

    1. find mask using area criteria (0.5) and plot boundary/mask
    2. trace contour using improved pavlidis algorithm and plot
       resultant contour and contour search
    3. compute geocontour from contour and plot, using simplify option
    4. plot geocontour with cartopy background options (borders and
       physical features) - will error and exit if cartopy not installed

    Returns
    -------
    None

    See Also
    --------
    small

    Examples
    --------
    >>> import numpy as np
    >>> import geocontour.masksearch as gcms
    >>> import geocontour.contourtrace as gcct
    >>> from geocontour import build
    >>> import geocontour.output as gco
    >>> import importlib.resources as ilr
    >>> boundfil=ilr.files('geocontour').
    ... joinpath('data/MissBasin_boundary.npz')
    >>> data=np.load(boundfil)
    >>> latitudes=data['latitudes']
    >>> longitudes=data['longitudes']
    >>> boundary=data['boundary']
    >>> mask=gcms.area(latitudes, longitudes, boundary)
    >>> contour, contoursearch=gcct.pavlidis_imp(mask, latitudes,
    ... longitudes, direction='ccw')
    >>> geocontour=build(contour, latitudes, longitudes, simplify=True)
    >>> gco.plot(latitudes, longitudes, boundary=boundary, mask=mask,
    ... title='Example Mask and Boundary\\nMississippi River Basin',
    ... outname='example_large_boundary+mask', transp=True)
    >>> gco.plot(latitudes, longitudes, mask=mask,
    ... contoursearch=contoursearch, 
    ... title='Example Contour Search\\nMississippi River Basin',
    ... outname='example_large_contoursearch', transp=True)
    >>> gco.plot(latitudes, longitudes, contour=contour,
    ... cells='contour',
    ... title='Example Contour\\nMississippi River Basin',
    ... outname='example_large_contour', transp=True)
    >>> gco.plot(latitudes, longitudes, geocontour=geocontour,
    ... title='Example Geocontour\\nMississippi River Basin',
    ... outname='example_large_geocontour', transp=True)
    >>> gco.plot(latitudes, longitudes, geocontour=geocontour,
    ... title='Example Geocontour\\nMississippi River Basin',
    ... outname='example_large_geocontour+natfeat', features='natural',
    ... transp=True)
    >>> gco.plot(latitudes, longitudes, geocontour=geocontour,
    ... title='Example Geocontour\\nMississippi River Basin',
    ... outname='example_large_geocontour+bordfeat', features='borders',
    ... transp=True)
    """
    print('\nRunning large scale geocontour example\nExample images will be output to current directory')
    boundfil=ilr.files('geocontour').joinpath('data/MissBasin_boundary.npz')
    data=np.load(boundfil)
    latitudes=data['latitudes']
    longitudes=data['longitudes']
    boundary=data['boundary']
    mask=gcms.area(latitudes,longitudes,boundary)
    contour,contoursearch=gcct.pavlidis_imp(mask,latitudes,longitudes,direction='ccw')
    geocontour=build(contour,latitudes,longitudes,simplify=True)
    print('  Example figure saved as \'example_large_boundary+mask\'')
    gco.plot(latitudes,longitudes,boundary=boundary,mask=mask,title='Example Mask and Boundary\nMississippi River Basin',outname='example_large_boundary+mask',transp=True)
    print('  Example figure saved as \'example_large_contoursearch\'')
    gco.plot(latitudes,longitudes,mask=mask,contoursearch=contoursearch,title='Example Contour Search\nMississippi River Basin',outname='example_large_contoursearch',transp=True)
    print('  Example figure saved as \'example_large_contour\'')
    gco.plot(latitudes,longitudes,contour=contour,cells='contour',title='Example Contour\nMississippi River Basin',outname='example_large_contour',transp=True)
    print('  Example figure saved as \'example_large_geocontour\'')
    gco.plot(latitudes,longitudes,geocontour=geocontour,title='Example Geocontour\nMississippi River Basin',outname='example_large_geocontour',transp=True)
    print('  Example figure saved as \'example_large_geocontour+natfeat\'')
    gco.plot(latitudes,longitudes,geocontour=geocontour,title='Example Geocontour\nMississippi River Basin',outname='example_large_geocontour+natfeat',features='natural',transp=True)
    print('  Example figure saved as \'example_large_geocontour+bordfeat\'')
    gco.plot(latitudes,longitudes,geocontour=geocontour,title='Example Geocontour\nMississippi River Basin',outname='example_large_geocontour+bordfeat',features='borders',transp=True)
