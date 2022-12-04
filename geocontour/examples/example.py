import numpy as np
import geocontour.masksearch as gcms
import geocontour.contourtrace as gcct
from geocontour import build
import geocontour.output as gco
import importlib.resources as ilr

def small():
    """
    Runs a small scale example of geocontour processing using mock data, saves resulting plots to run directory
        1. find mask using area criteria (0.5) and plot boundary/mask
        2. trace contour using improved pavlidis algorithm and plot resultant contour and contour search
        3. compute geocontour from contour and plot, using simplify option

    Inputs:
        none

    Outputs:
        none (plots saved to run directory)
    """
    print('\nRunning small scale geocontour example\nExample images will be output to current directory')
    boundfil=ilr.files('geocontour').joinpath('data/generic_boundary.npz')
    data=np.load(boundfil)
    latitudes=data['latitudes']
    longitudes=data['longitudes']
    boundary=data['boundary']
    mask=gcms.area(latitudes,longitudes,boundary)
    contour,contoursearch=gcct.pavlidis_imp(mask,latitudes,longitudes)
    geocontour=build(contour,latitudes,longitudes,simplify=True)
    print('  Example figure saved as \'example_small_boundary+mask\'')
    gco.plot(latitudes,longitudes,boundary=boundary,mask=mask,title='Example Mask and Boundary',outname='example_small_boundary+mask',outdpi='indep')
    print('  Example figure saved as \'example_small_contoursearch\'')
    gco.plot(latitudes,longitudes,mask=mask,contoursearch=contoursearch,title='Example Contour Search',outname='example_small_contoursearch',outdpi='indep')
    print('  Example figure saved as \'example_small_contour\'')
    gco.plot(latitudes,longitudes,contour=contour,cells='contour',title='Example Contour',outname='example_small_contour',outdpi='indep')
    print('  Example figure saved as \'example_small_geocontour\'')
    gco.plot(latitudes,longitudes,geocontour=geocontour,buffer='on',title='Example Geocontour',outname='example_small_geocontour',outdpi='indep')

def large():
    """
    Runs a large scale example of geocontour processing using the Mississippi River Basin boundary, saves resulting plots to run directory
        1. find mask using area criteria (0.5) and plot boundary/mask
        2. trace contour using improved pavlidis algorithm and plot resultant contour and contour search
        3. compute geocontour from contour and plot, using simplify option
        4. plot geocontour with cartopy background options (borders and physical features) - will error and exit if cartopy not installed

    Inputs:
        none

    Outputs:
        none (plots saved to run directory)
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
    gco.plot(latitudes,longitudes,boundary=boundary,mask=mask,title='Example Mask and Boundary\nMississippi River Basin',outname='example_large_boundary+mask')
    print('  Example figure saved as \'example_large_contoursearch\'')
    gco.plot(latitudes,longitudes,mask=mask,contoursearch=contoursearch,title='Example Contour Search\nMississippi River Basin',outname='example_large_contoursearch')
    print('  Example figure saved as \'example_large_contour\'')
    gco.plot(latitudes,longitudes,contour=contour,cells='contour',title='Example Contour\nMississippi River Basin',outname='example_large_contour')
    print('  Example figure saved as \'example_large_geocontour\'')
    gco.plot(latitudes,longitudes,geocontour=geocontour,title='Example Geocontour\nMississippi River Basin',outname='example_large_geocontour')
    print('  Example figure saved as \'example_large_geocontour+natfeat\'')
    gco.plot(latitudes,longitudes,geocontour=geocontour,title='Example Geocontour\nMississippi River Basin',outname='example_large_geocontour+natfeat',features='natural')
    print('  Example figure saved as \'example_large_geocontour+bordfeat\'')
    gco.plot(latitudes,longitudes,geocontour=geocontour,title='Example Geocontour\nMississippi River Basin',outname='example_large_geocontour+bordfeat',features='borders')

if __name__=='__main__':
    small()
    large()
