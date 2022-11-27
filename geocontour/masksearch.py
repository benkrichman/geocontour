import warnings
import numpy
import shapely.geometry
import matplotlib.path
import geocontour.grid
import geocontour.maskutil

def center(latitudes,longitudes,boundary,precision=1e-5):
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
    latdir=geocontour.grid.checklatitudedirection(latitudes)
    if latdir=='dec':
        latitudes=numpy.flip(latitudes)
    boxlatmin, boxlatmax, boxlonmin, boxlonmax = geocontour.maskutil.boxset(latitudes,longitudes,boundary)
    boxmask=numpy.full((boxlatmax-boxlatmin+1,boxlonmax-boxlonmin+1),False)
    boundpoly=shapely.geometry.Polygon(boundary).buffer(precision)
    for la in numpy.arange(boxlatmin,boxlatmax+1,1):
        for lo in numpy.arange(boxlonmin,boxlonmax+1,1):
            center=shapely.geometry.Point(latitudes[la],longitudes[lo])
            boxmask[la-boxlatmin,lo-boxlonmin]=boundpoly.contains(center)
    mask=numpy.full((len(latitudes),len(longitudes)),False)
    mask[boxlatmin:boxlatmax+1,boxlonmin:boxlonmax+1]=boxmask
    if latdir=='dec':
        mask=numpy.flip(mask,axis=0)
    return mask

def center2(latitudes,longitudes,boundary):
    """
    Returns a mask over a range of input latitudes and longitudes determined by an input boundary
        Critera for inclusion of a cell is whether the center of the cell falls within the boundary
        Functionally matches geocontour.masksearch.center(), but utilizes matplotlib.path functions, which are probably optimized and thus is roughly 2.5*sqrt(N) faster for N points, though lacks a "precision" buffer input
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
    latdir=geocontour.grid.checklatitudedirection(latitudes)
    if latdir=='dec':
        latitudes=numpy.flip(latitudes)
    boxlatmin, boxlatmax, boxlonmin, boxlonmax = geocontour.maskutil.boxset(latitudes,longitudes,boundary)
    boundpoly=matplotlib.path.Path(boundary)
    ysp, xsp = numpy.meshgrid(latitudes[boxlatmin:boxlatmax+1],longitudes[boxlonmin:boxlonmax+1], indexing='ij')
    searchpoints=numpy.hstack((ysp.reshape((-1,1)), xsp.reshape((-1,1))))
    boxmask=boundpoly.contains_points(searchpoints)
    mask=numpy.full((len(latitudes),len(longitudes)),False)
    mask[boxlatmin:boxlatmax+1,boxlonmin:boxlonmax+1]=boxmask.reshape((boxlatmax-boxlatmin+1,boxlonmax-boxlonmin+1))
    if latdir=='dec':
        mask=numpy.flip(mask,axis=0)
    return mask

def nodes(latitudes,longitudes,boundary,nodes=2,precision=1e-5):
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
    latdir=geocontour.grid.checklatitudedirection(latitudes)
    if latdir=='dec':
        latitudes=numpy.flip(latitudes)
    boxlatmin, boxlatmax, boxlonmin, boxlonmax = geocontour.maskutil.boxset(latitudes,longitudes,boundary)
    latgrdspc=geocontour.grid.gridspacing(latitudes)
    longrdspc=geocontour.grid.gridspacing(longitudes)
    boxmask=numpy.full((boxlatmax-boxlatmin+1,boxlonmax-boxlonmin+1),False)
    boundpoly=shapely.geometry.Polygon(boundary).buffer(precision)
    for la in numpy.arange(boxlatmin,boxlatmax+1,1):
        for lo in numpy.arange(boxlonmin,boxlonmax+1,1):
            nodeLL=shapely.geometry.Point(latitudes[la]-latgrdspc/2,longitudes[lo]-longrdspc/2)
            nodeHL=shapely.geometry.Point(latitudes[la]+latgrdspc/2,longitudes[lo]-longrdspc/2)
            nodeLH=shapely.geometry.Point(latitudes[la]-latgrdspc/2,longitudes[lo]+longrdspc/2)
            nodeHH=shapely.geometry.Point(latitudes[la]+latgrdspc/2,longitudes[lo]+longrdspc/2)
            nodesinmask=numpy.array([boundpoly.contains(nodeLL),boundpoly.contains(nodeHL),boundpoly.contains(nodeLH),boundpoly.contains(nodeHH)])
            if nodesinmask.sum()>=nodes:
                boxmask[la-boxlatmin,lo-boxlonmin]=True
    mask=numpy.full((len(latitudes),len(longitudes)),False)
    mask[boxlatmin:boxlatmax+1,boxlonmin:boxlonmax+1]=boxmask
    if latdir=='dec':
        mask=numpy.flip(mask,axis=0)
    return mask

def nodes2(latitudes,longitudes,boundary,nodes=2,precision=1e-5):
    """
    Returns a mask over a range of input latitudes and longitudes determined by an input boundary
        Critera for inclusion of a cell is whether a given number (default=2) of cell nodes (corners) fall within the boundary 
        Functionally matches geocontour.masksearch.nodes(), but utilizes matplotlib.path functions, though speed is similar to the shapely implementation

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
    latdir=geocontour.grid.checklatitudedirection(latitudes)
    if latdir=='dec':
        latitudes=numpy.flip(latitudes)
    boxlatmin, boxlatmax, boxlonmin, boxlonmax = geocontour.maskutil.boxset(latitudes,longitudes,boundary)
    latgrdspc=geocontour.grid.gridspacing(latitudes)
    longrdspc=geocontour.grid.gridspacing(longitudes)
    boxmask=numpy.full((boxlatmax-boxlatmin+1,boxlonmax-boxlonmin+1),False)
    boundpoly=matplotlib.path.Path(boundary)
    for la in numpy.arange(boxlatmin,boxlatmax+1,1):
        for lo in numpy.arange(boxlonmin,boxlonmax+1,1):
            nodeLL=[latitudes[la]-latgrdspc/2,longitudes[lo]-longrdspc/2]
            nodeHL=[latitudes[la]+latgrdspc/2,longitudes[lo]-longrdspc/2]
            nodeLH=[latitudes[la]-latgrdspc/2,longitudes[lo]+longrdspc/2]
            nodeHH=[latitudes[la]+latgrdspc/2,longitudes[lo]+longrdspc/2]
            nodesinmask=boundpoly.contains_points(numpy.array([nodeLL,nodeHL,nodeLH,nodeHH]))
            if nodesinmask.sum()>=nodes:
                boxmask[la-boxlatmin,lo-boxlonmin]=True
    mask=numpy.full((len(latitudes),len(longitudes)),False)
    mask[boxlatmin:boxlatmax+1,boxlonmin:boxlonmax+1]=boxmask
    if latdir=='dec':
        mask=numpy.flip(mask,axis=0)
    return mask

def area(latitudes,longitudes,boundary,area=0.5):
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
    latdir=geocontour.grid.checklatitudedirection(latitudes)
    if latdir=='dec':
        latitudes=numpy.flip(latitudes)
    boxlatmin, boxlatmax, boxlonmin, boxlonmax = geocontour.maskutil.boxset(latitudes,longitudes,boundary)
    latgrdspc=geocontour.grid.gridspacing(latitudes)
    longrdspc=geocontour.grid.gridspacing(longitudes)
    boxmask=numpy.full((boxlatmax-boxlatmin+1,boxlonmax-boxlonmin+1),False)
    boundpoly=shapely.geometry.Polygon(boundary)
    for la in numpy.arange(boxlatmin,boxlatmax+1,1):
        for lo in numpy.arange(boxlonmin,boxlonmax+1,1):
            LL=[latitudes[la]-latgrdspc/2,longitudes[lo]-longrdspc/2]
            HL=[latitudes[la]+latgrdspc/2,longitudes[lo]-longrdspc/2]
            LH=[latitudes[la]-latgrdspc/2,longitudes[lo]+longrdspc/2]
            HH=[latitudes[la]+latgrdspc/2,longitudes[lo]+longrdspc/2]
            cell=shapely.geometry.Polygon([LL,HL,HH,LH])
            if boundpoly.intersection(cell).area>=latgrdspc*longrdspc*area:
                boxmask[la-boxlatmin,lo-boxlonmin]=True
    mask=numpy.full((len(latitudes),len(longitudes)),False)
    mask[boxlatmin:boxlatmax+1,boxlonmin:boxlonmax+1]=boxmask
    if latdir=='dec':
        mask=numpy.flip(mask,axis=0)
    return mask

