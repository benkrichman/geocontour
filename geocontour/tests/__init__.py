"""
Functions for testing of the geocontour package
===============================================
"""
import geocontour as gc
import numpy as np
import importlib.resources as ilr

def full():
    """
    Run all user-facing geocontour functions with test data,
    printing/saving results

    This is a useful first step to try after installation to determine
    whether everything is working properly.

    Returns
    -------
    None

    See Also
    --------
    geocontour.tests.timing
    """
    np.set_printoptions(precision=3)
    boundfil=ilr.files('geocontour').joinpath('data/generic_boundary.npz')
    data=np.load(boundfil)
    latitudes=data['latitudes']
    longitudes=data['longitudes']
    boundary=data['boundary']
    #Test Grid Funcs
    print('\n\n***TEST GRID FUNCS***\n')
    print('\n\nTest geocontour.grid.spacing()\n')
    latspc=gc.grid.spacing(latitudes)
    lonspc=gc.grid.spacing(longitudes)
    print('latspc='+str(latspc)+'\nlonspc='+str(lonspc))
    print('\n\nTest geocontour.grid.(lon/lat)lens()\n')
    lonlens=gc.grid.lonlens(latitudes)
    latlens=gc.grid.latlens(latitudes)
    print('latlens:')
    print(latlens)
    print('lonlens:')
    print(lonlens)
    print('\n\nTest geocontour.grid.(lon/lat)len()\n')
    lonlen=gc.grid.lonlen(latitudes)
    latlen=gc.grid.latlen(latitudes)
    print('latlen:')
    print(latlen)
    print('lonlen:')
    print(lonlen)
    print('\n\nTest geocontour.grid.areas()\n')
    areas=gc.grid.areas(latitudes,longitudes)
    print('areas:')
    print(areas)
    print('\n\nTest geocontour.grid.clonrng()\n')
    lonrng=gc.grid.clonrng(longitudes)
    print('lonrng='+lonrng)
    print('\n\nTest geocontour.grid.clatdir()\n')
    latdir=gc.grid.clatdir(latitudes)
    print('latdir='+latdir)
    print('\n\nTest geocontour.grid.switchlon()\n')
    swclon=gc.grid.switchlon(longitudes,'pos',print_output=True)
    print('swclon='+str(swclon))
    print('\n\nTest geocontour.grid.switchind()\n')
    swcind=gc.grid.switchind(longitudes)
    print('swcind='+str(swcind))
    #Test Mask Search Funcs
    print('\n\n***TEST MASK SEARCH FUNCS***\n')
    print('\n\nTest geocontour.masksearch.center()\n')
    mask=gc.masksearch.center(latitudes,longitudes,boundary)
    print('mask:')
    print(mask)
    print('\n\nTest geocontour.masksearch.center2()\n')
    mask=gc.masksearch.center2(latitudes,longitudes,boundary)
    print('mask:')
    print(mask)
    print('\n\nTest geocontour.masksearch.nodes()\n')
    mask=gc.masksearch.nodes(latitudes,longitudes,boundary)
    print('mask:')
    print(mask)
    print('\n\nTest geocontour.masksearch.nodes2()\n')
    mask=gc.masksearch.nodes2(latitudes,longitudes,boundary)
    print('mask:')
    print(mask)
    print('\n\nTest geocontour.masksearch.area()\n')
    mask=gc.masksearch.area(latitudes,longitudes,boundary)
    print('mask:')
    print(mask)
    #Test Mask Util Funcs
    print('\n\n***TEST MASK UTIL FUNCS***\n')
    print('\n\nTest geocontour.maskutil.bbox()\n')
    box1,box2,box3,box4=gc.maskutil.bbox(latitudes,longitudes,boundary)
    print('bbox='+str(box1)+str(box2)+str(box3)+str(box4))
    print('\n\nTest geocontour.maskutil.edge()\n')
    edgemask,edgecells=gc.maskutil.edge(mask,latitudes=latitudes,longitudes=longitudes)
    print('edgemask:')
    print(edgemask)
    print('edgecells:')
    print(edgecells)
    print('\n\nTest geocontour.maskutil.vertex()\n')
    vertexpoints,edgevertexpoints=gc.maskutil.vertex(mask,latitudes,longitudes)
    print('vertexpoints:')
    print(vertexpoints)
    print('edgevertexpoints:')
    print(edgevertexpoints)
    print('\n\nTest geocontour.maskutil.neighbors()\n')
    neighbors=gc.maskutil.neighbors(np.array([1,1]))
    print('neighbors:')
    print(neighbors)
    print('\n\nTest geocontour.maskutil.conn()\n')
    conn=gc.maskutil.conn(mask)
    print('connectivity='+str(conn))
    #Test Contour Util Funcs
    #*3/4 are internal use and run during contour tracing
    print('\n\n***TEST CONTOUR UTIL FUNCS***')
    print('   *3/4 contour util funcs are internal and run during contour tracing\n')
    print('\n\nTest geocontour.contourutil.findstart()\n')
    start=gc.contourutil.findstart(mask)
    print('start='+str(start))
    #Test Contour Trace Funcs
    print('\n\n***TEST CONTOUR TRACE FUNCS***\n')
    print('\n\nTest geocontour.contourtrace.square()\n')
    contour,contoursearch=gc.contourtrace.square(mask,latitudes=latitudes,longitudes=longitudes)
    print('contour:')
    print(contour)
    print('\n\nTest geocontour.contourtrce.moore()\n')
    contour,contoursearch=gc.contourtrace.moore(mask,latitudes=latitudes,longitudes=longitudes)
    print('contour:')
    print(contour)
    print('\n\nTest geocontour.contourtrace.moore_imp()\n')
    contour,contoursearch=gc.contourtrace.moore_imp(mask,latitudes=latitudes,longitudes=longitudes)
    print('contour:')
    print(contour)
    print('\n\nTest geocontour.contourtrace.pavlidis()\n')
    contour,contoursearch=gc.contourtrace.pavlidis(mask,latitudes=latitudes,longitudes=longitudes)
    print('contour:')
    print(contour)
    print('\n\nTest geocontour.contourtrace.pavlidis_imp()\n')
    contour,contoursearch=gc.contourtrace.pavlidis_imp(mask,latitudes=latitudes,longitudes=longitudes)
    print('contour:')
    print(contour)
    print('\n\nTest geocontour.contourtrace.MSBF()\n')
    contour,contoursearch=gc.contourtrace.MSBF(mask,latitudes=latitudes,longitudes=longitudes)
    print('contour:')
    print(contour)
    print('\n\nTest geocontour.contourtrace.ISBF()\n')
    contour,contoursearch=gc.contourtrace.ISBF(mask,latitudes=latitudes,longitudes=longitudes)
    print('contour:')
    print(contour)
    print('\n\nTest geocontour.contourtrace.TSR()\n')
    contour,contoursearch=gc.contourtrace.TSR(mask,latitudes=latitudes,longitudes=longitudes)
    print('contour:')
    print(contour)
    #Test Geocontour Funcs
    print('\n\n***TEST GEOCONTOUR FUNCS***\n')
    print('\n\nTest geocontour.build()\n')
    geocontour=gc.build(contour,latitudes,longitudes,simplify=True)
    print('geocontour:')
    print(geocontour)
    #Test Check Funcs
    print('\n\n***TEST CHECK FUNCS***\n')
    print('\n\nTest geocontour.check.cdim()\n')
    gc.check.cdim(latitudes)
    gc.check.cdim(longitudes)
    print('\n\nTest geocontour.check.cboundary()\n')
    gc.check.cboundary(boundary)
    print('\n\nTest geocontour.check.cmask()\n')
    gc.check.cmask(mask,latitudes,longitudes)
    print('\n\nTest geocontour.check.ccontour()\n')
    gc.check.ccontour(contour,latitudes,longitudes)
    print('\n\nTest geocontour.check.cgeocontour()\n')
    gc.check.cgeocontour(geocontour,latitudes,longitudes)
    #Test Output Funcs
    print('\n\n***TEST OUTPUT FUNCS***\n')
    print('\n\nTest geocontour.output.plot()\n')
    gc.output.plot(latitudes,longitudes,boundary=boundary,mask=mask,contour=contour,geocontour=geocontour,title='geocontour_test_plot',outname='geocontour_test_plot',transp=True)
    print('Plot created and saved to run directory')
    print('\n\nTest geocontour.output.save()\n')
    gc.output.save(latitudes,longitudes,boundary=boundary,mask=mask,contour=contour,contoursearch=contoursearch,geocontour=geocontour,outname='geocontour_test_output')
    print('Output saved in .npz to run directory')
