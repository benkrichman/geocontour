import warnings
import numpy as np
import shapely.geometry as shg
import matplotlib.path as mplp
import geocontour.grid as gcg
import geocontour.maskutil as gcmu

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
    latdir=gcg.clatdir(latitudes)
    if latdir=='dec':
        latitudes=np.flip(latitudes)
    boxlatmin, boxlatmax, boxlonmin, boxlonmax = gcmu.bbox(latitudes,longitudes,boundary)
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

def center2(latitudes,longitudes,boundary,precision=1e-5):
    """
    Returns a mask over a range of input latitudes and longitudes determined by an input boundary
        Critera for inclusion of a cell is whether the center of the cell falls within the boundary
        Functionally matches geocontour.masksearch.center(), but utilizes matplotlib.path functions, which are probably optimized and capable of being vectorized, and thus is roughly <how much??> faster for N cells

    Inputs (Required):
        latitudes - An evenly spaced numpy array of latitude points (degrees)
        longitudes - An evenly spaced numpy array of longitude points (degrees)
        boundary - A 2-d Nx2 numpy array of latitude/longitude points (degrees)

    Inputs (Optional):
        precision - Passed to matplotlib.Path to increase boundary polygon area, default=1e-5
            Matplotlib.Path can't beat machine precision, and can thus give "incorrect" results for very close points or shapes. This input errs on being more inclusive, in particular to capture points falling directly on a boundary. A decent rule is to set the precision value as high as you can without impeding the accuracy. For instance, the default of 1e-5 (degrees) translates to roughly 1m precision at the equator. The buffer can be negated by setting this input very low (to machine precision).

    Outputs:
        mask - A 2-d boolean numpy array of dimension MxN where M=len(latitudes) and N=len(longitudes)
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
    latdir=gcg.clatdir(latitudes)
    if latdir=='dec':
        latitudes=np.flip(latitudes)
    boxlatmin, boxlatmax, boxlonmin, boxlonmax = gcmu.bbox(latitudes,longitudes,boundary)
    latgrdspc=gcg.spacing(latitudes)
    longrdspc=gcg.spacing(longitudes)
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

def nodes2(latitudes,longitudes,boundary,nodes=2,precision=1e-5):
    """
    Returns a mask over a range of input latitudes and longitudes determined by an input boundary
        Critera for inclusion of a cell is whether a given number (default=2) of cell nodes (corners) fall within the boundary 
        Functionally matches geocontour.masksearch.nodes(), but utilizes matplotlib.path functions, which are probably optimized and capable of being vectorized, and thus is roughly <how much??> faster for N cells

    Inputs (Required):
        latitudes - An evenly spaced numpy array of latitude points (degrees)
        longitudes - An evenly spaced numpy array of longitude points (degrees)
        boundary - A 2-d Nx2 numpy array of latitude/longitude points (degrees)

    Inputs (Optional):
        nodes - The number of cell nodes (corners) to use as a criteria for inclusion (1-4)
        precision - Passed to matplotlib.Path to increase boundary polygon area, default=1e-5
            Matplotlib.Path can't beat machine precision, and can thus give "incorrect" results for very close points or shapes. This input errs on being more inclusive, in particular to capture points falling directly on a boundary. A decent rule is to set the precision value as high as you can without impeding the accuracy. For instance, the default of 1e-5 (degrees) translates to roughly 1m precision at the equator. The buffer can be negated by setting this input very low (to machine precision).

    Outputs:
        mask - A 2-d boolean numpy array of dimension MxN where M=len(latitudes) and N=len(longitudes)
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
    latdir=gcg.clatdir(latitudes)
    if latdir=='dec':
        latitudes=np.flip(latitudes)
    boxlatmin, boxlatmax, boxlonmin, boxlonmax = gcmu.bbox(latitudes,longitudes,boundary)
    latgrdspc=gcg.spacing(latitudes)
    longrdspc=gcg.spacing(longitudes)
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

