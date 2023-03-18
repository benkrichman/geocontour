"""
Functions for visualizing and saving output produced by geocontour
==================================================================
"""
import sys
import warnings
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mplc
import datascale as ds
import geocontour.grid as gcg
import geocontour.check as gcc
import geocontour.maskutil as gcmu
import geocontour.contourutil as gccu
try:
    import cartopy as cp
    cp_exists='y'
except:
    cp_exists='n'

def plot(latitudes,longitudes,boundary=None,mask=None,contour=None,contoursearch=None,geocontour=None,vertices=None,boundingbox='all',buffer=False,grid=True,cells='default',showcontour=True,startcell=True,contourarrows=True,contoursearcharrows=True,fancycontoursearch=True,contoursearch_contraction=0.2,contoursearch_shift=0.25,geocontourvectors=True,emptycellcolor='lightgrey',fullcellcolor='sandybrown',boundarycolor='tab:blue',contourcolor='olivedrab',contoursearchcolor='firebrick',geocontourcolor='olivedrab',vertexcolor='tab:cyan',gridcolor='black',lw_boundary=0.1,lw_contour=0.1,lw_contoursearch='auto',lw_geocontour=0.1,mw_contourarrows=0.5,mw_contoursearcharrows='auto',mw_vertices=0.4,features=None,title=None,labs=True,outname='plot',outdpi='high',transp=False):
    """
    Plot any/all geocontour-created elements:

    - `boundary`
    - `mask`
    - `contour`
    - `contoursearch`
    - `geocontour`
    - `vertices`

    Parameters
    ----------
    latitudes : ndarray
        1D Nx1 array of latitude points (degrees)
    longitudes : ndarray
        1D Nx1 array of longitude points (degrees)
    boundary : ndarray, optional 
        2D Nx2 array of latitude/longitude points (degrees) with the
        last point equal to the first
    mask : ndarray, optional
        2D MxN bool array where M=len(`latitudes`) and
        N=len(`longitudes`)
    contour : ndarray, optional
        2D Nx2 array of ordered latitude/longitude points (degrees)
        describing the edge of a mask
    contoursearch : ndarray, optional
        2D Nx2 array of ordered latitude/longitude points (degrees)
        describing the cells searched during contour tracing
    geocontour : ndarray, optional
        3D Nx2x5 array defining a list of N contour cells (column 1),
        their edge points (columns 2,3), segment lengths (column 4), and
        outward unit vectors (column 5)
    vertices : ndarray, optional
        2D Nx2 array of latitude/longitude points (degrees)
    boundingbox : {'all', 'boundary', 'mask', 'contour',\
                  'contoursearch', 'geocontour'}, default='all' 
        select the plot element(s) to be used as x and y axis limits
        [default value 'all' uses the largest bounds from all provided
        elements]
    buffer : bool, default=True
        create a 1 cell buffer on the bounding box
    grid : bool, default=True
        show/hide the cell grid
    cells : {'default', 'none', 'mask', 'maskedge-4', 'maskedge-8',\
                'contour', 'geocontour'}, default='default'
        the plot element to be used for filling the cell grid

            ``default``
                use first provided element, in this order: 'geocontour',
                'contour', 'mask', 'none'
            ``none``
                no filled grid cells
            ``mask``
                all `mask` cells
            ``maskedge-4``
                4-connected `mask` edge cells
            ``maskedge-8``
                8-connected `mask` edge cells
            ``contour``
                `contour` cells
            ``geocontour``
                `geocontour` cells
    showcontour : bool, default=True
        show/hide the `contour` (if provided)
        [allows showing contour cells without line plot of `contour`]
    startcell : bool, default=True
        show/hide the startcell (`contour` or `contoursearch`, if
        provided)
    contourarrows : bool, default=True
        show/hide directional arrows on the `contour` (if provided)
    contoursearcharrows : bool, default=True
        show/hide directional arrows on the `contoursearch` (if
        provided)
    fancycontoursearch : bool, default=True
        plot the `contoursearch` (if provided) in a cleaner and more
        easily followed format
    contoursearch_contraction : float, default=0.2
        value determining how much `contoursearch` "shrinks" towards
        `contour` cell edges (range 0 to 0.5)
    contoursearch_shift : float, default=0.25
        value determining how much `contoursearch` shifts to avoid
        doubling back on itself (range 0 to 0.5)
    geocontourvectors : bool, default=True
        show/hide outward normal vectors on the `geocontour` (if
        provided)
    emptycellcolor : matplotlib_color, default='lightgray'
        color for unmasked cells
    fullcellcolor : matplotlib_color, default='sandybrown'
        color for masked cells
    boundarycolor : matplotlib_color, default='tab:blue'
        color for `boundary` (if provided)
    contourcolor : matplotlib_color, default='olivedrab'
        color for `contour` (if provided)
    contoursearchcolor : matplotlib_color, default='firebrickred'
        color for `contoursearch` (if provided)
    geocontourcolor : matplotlib_color, default='olivedrab'
        color for `geocontour` (if provided)
    vertexcolor : matplotlib_color, default='tab:cyan'
        color for `vertices` (if provided)
    gridcolor : matplotlib_color, default='black'
        color for grid
    lw_boundary : float, default=0.1
        `boundary` linewidth in units of grid cells (see Notes for more
        detail)
    lw_contour : float, default=0.1
        `contour` linewidth in units of grid cells (see Notes for more
        detail)
    lw_contoursearch : {'auto', float}, default='auto'
        `contoursearch` linewidth in units of grid cells (see Notes for
        more detail)
        [i.e. 'auto' specifies 0.075 if `fancycontoursearch` = True and
        0.1 if `fancycontoursearch` = False]
    lw_geocontour : float, default=0.1
        `geocontour` linewidth in units of grid cells (see Notes for
        more detail) 
    mw_contourarrows : float, default=0.5
        `contour` arrow markerwidth in units of grid cells (see Notes
        for more detail)
    mw_contoursearcharrows : {'auto', float}, default='auto'
        `contoursearch` arrow markerwidth in units of grid cells (see
        Notes for more detail)
        [i.e. 'auto' specifies 0.35 if `fancycontoursearch` = True and
        0.5 if `fancycontoursearch` = False]
    mw_vertices : float, default=0.4
        `vertex` markerwidth in units of grid cells (see Notes for more
        detail)
    features : {None, 'natural', 'borders'}, default=None
        display geographic features (if cartopy is installed)

            ``natural``
                display coastlines, ocean, and lakes/rivers
            ``borders``
                display national and state/province level boundaries
    title : str, default=None
        the plot title
    labs : bool, default=True
        x/y plot labels on or off [e.g. off if no associated lat/lon]
    outname : str, default='plot'
        filename/path for the saved image
    outdpi : {'high', 'low', 'indep', float}, default='high'
        dpi of the saved image

            ``high``
                scale dpi high enough (36 pixels per grid cell) to see
                features when zooming into a single grid cell (dpi floor
                of 100)
                [for very large grids, 'auto' may set dpi high enough
                that pyplot will hang on some systems - setting `outdpi`
                = 'low' or entering desired dpi manually can avoid this
                if encountered]
            ``low``
                scale dpi to 5 pixels per grid cell (dpi floor of 100)
            ``indep``
                set output dpi to 200 regardless of grid size/spacing
            `float`
                any other entry must be a numerical value specifying
                desired pixels per grid cell
    transp : bool, default=False
        set exterior of plot to transparent with text/labels/ticks/frame
        set to 'dimgray' for contrast against light and dark backgrounds
        (e.g. for environments in which light/dark mode may not be
        constant)

    Returns
    -------
    None

    See Also
    --------
    geocontour.contourutil.fancysearch
    datascale.plotdatasize
    datascale.plotdatadpi

    Notes
    -----
    - Some inputs will only be utilized in the presence of others. For
      example (not exhaustive), `contourcolor` will have no effect if no
      `contour` is provided, `mw_contourarrows` will have no effect if
      `contourarrows` = False, and `contoursearch_contraction` will have
      no effect if `fancycontoursearch` = False.
    - All color inputs may be provided as matplotlib pre-defined colors
      (e.g. 'red' or 'r'), hex RGB(A) strings (e.g. 'ff0000'), or RGB(A)
      tuples (e.g. (1,0,0)).
    - Inputs for linewidths and markerwidths are given in units of grid
      cells. For example, a linewidth of 0.5 means the line will be
      half as wide as grid cells for the given plot and a markerwidth of
      0.5 means the same. The `datascale` package is utilized for this
      functionality and the intent is to have plot features all scale
      and be interpretable on the plot as a default, since the primary
      function of plot is diagnostic. Example or publication figures may
      require user-altered line and marker widths, though adjustment
      should be simple as the inputs are all in units of the plotted
      grid cells.

    Examples
    --------
    1. Plot a mask and the boundary used to create it
        plot(latitudes, longitudes, boundary=<boundary>, mask=<mask>,
        lw_boundary=0.2)
    2. Plot a contour
        plot(latitudes, longitudes, countour=<countour>,
        boundingbox='contour', buffer=True)
    3. Plot a contoursearch overlaying contour cells
        plot(latitudes, longitudes, contour=<contour>,
        contoursearch=<contoursearch>, showcontour=False,
        cells='contour', boundingbox='contoursearch', buffer=True)
    4. Plot a geocontour overlaid onto a map projection with natural features
        plot(latitudes, longitudes, geocontour=<geocontour>,
        cells='geocontour', features='natural')
    """
    latspc=gcg.spacing(latitudes)
    lonspc=gcg.spacing(longitudes)
    gridlatmin=latitudes.min()-latspc/2
    gridlatmax=latitudes.max()+latspc/2
    gridlonmin=longitudes.min()-lonspc/2
    gridlonmax=longitudes.max()+lonspc/2
    latdir=gcg.clatdir(latitudes)
    boundingarray=np.array([[[gridlatmin,gridlatmax],[gridlonmin,gridlonmax]]])
    if boundary is None:
        if boundingbox=='boundary':
            sys.exit('ERROR - Can not use \'boundary\' for boundingbox if no boundary input provided')
    else:
        gcc.cboundary(boundary)
        loninput=gcg.clonrng(longitudes)
        boundloninput=gcg.clonrng(boundary[:,1])
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
        gcc.cmask(mask,latitudes,longitudes)
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
        gcc.ccontour(contour,latitudes,longitudes)
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
        gcc.cgeocontour(geocontour,latitudes,longitudes)
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
    if buffer==True:
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
        pltmask=gcmu.edge(mask,connectivity=8).astype('int')
    elif cells=='maskedge-4':
        pltmask=gcmu.edge(mask,connectivity=4).astype('int')
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
    if transp==True:
        plt.rcParams.update({
            "axes.labelcolor": "dimgray",
            "axes.edgecolor": "dimgray",
            "axes.facecolor": "white",
            "xtick.color": "dimgray",
            "ytick.color": "dimgray",
            "text.color": "dimgray",})
        fig.patch.set_visible(False)
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
    if labs==True:
        ax.set_xlabel('Longitude (deg)')
        ax.set_ylabel('Latitude (deg)')
    ax.set_ylim((ylimmin,ylimmax))
    ax.set_xlim((xlimmin,xlimmax)) 
    ax.imshow(pltmask,aspect='equal',interpolation='none',vmin=0,vmax=1,extent=ext,origin=org,cmap=cmp,zorder=-1)
    pdw=ds.plotdatasize(ax,axis='xy',mult=(latspc+lonspc)/2,plottype='line')
    pds=ds.plotdatasize(ax,axis='xy',mult=(latspc+lonspc)/2,plottype='scatter')
    if features is None:
        ax.set_yticks(ymajor)
        ax.set_xticks(xmajor)
        ax.set_yticks(yminor,minor=True)
        ax.set_xticks(xminor,minor=True)
        if grid==True:
            ax.grid(which='minor',linestyle=(0,(1,1)),color=gridcolor,linewidth=0.05*pdw,zorder=3)
    else:
        ax.set_yticks(ymajor)
        ax.set_xticks(xmajor)
        if grid==True:
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
        lw_boundary=lw_boundary*pdw
        ax.plot(boundary[:,1],boundary[:,0],color=boundarycolor,linewidth=lw_boundary,zorder=10)
    if vertices is not None:
        mw_vertices=pds*mw_vertices**2
        ax.scatter(vertices[:,1],vertices[:,0],c=vertexcolor,s=mw_vertices,linewidth=0,zorder=10)
    if contoursearch is not None:
        if fancycontoursearch==True:
            autolwcs=0.075
            automwcs=0.35
            contoursearch=gccu.fancysearch(contoursearch,contraction=contoursearch_contraction,shift=contoursearch_shift)
        else:
            autolwcs=0.1
            automwcs=0.5
        if lw_contoursearch=='auto':
            lw_contoursearch=autolwcs*pdw
        else:
            lw_contoursearch=lw_contoursearch*pdw
        if mw_contoursearcharrows=='auto':
            mw_contoursearcharrows=pds*automwcs**2
        else:
            mw_contoursearcharrows=pds*mw_contoursearcharrows**2
        ax.plot(contoursearch[:,1],contoursearch[:,0],color=contoursearchcolor,linewidth=lw_contoursearch,zorder=11)
        if startcell==True:
            ax.scatter(contoursearch[0,1],contoursearch[0,0],marker='o',c=contoursearchcolor,s=mw_contoursearcharrows,linewidth=0,zorder=11)
        if contoursearcharrows==True:
            contoursearchdiff=np.diff(contoursearch,axis=0)
            contoursearchpointrotation=-180/np.pi*np.arctan2(contoursearchdiff[:,1],contoursearchdiff[:,0])
            contoursearchpointlocation=contoursearch[:-1]+2*contoursearchdiff/5
            for ct,k in enumerate(contoursearchpointlocation):
                if (contoursearchdiff[ct,:]==0).all():
                    ax.scatter(k[1],k[0],marker=(4,0,0),c=contoursearchcolor,s=mw_contoursearcharrows,linewidth=0,zorder=11)
                else:
                    ax.scatter(k[1],k[0],marker=(3,0,contoursearchpointrotation[ct]),c=contoursearchcolor,s=mw_contoursearcharrows,linewidth=0,zorder=11)
    if contour is not None:
        if showcontour==True:
            lw_contour=lw_contour*pdw
            ax.plot(contour[:,1],contour[:,0],color=contourcolor,linewidth=lw_contour,zorder=12)
            mw_contourarrows=pds*mw_contourarrows**2
            if startcell==True:
                ax.scatter(contour[0,1],contour[0,0],marker='o',c=contourcolor,s=mw_contourarrows,linewidth=0,zorder=12)
            if contourarrows==True:
                contourdiff=np.diff(contour,axis=0)
                contourpointrotation=-180/np.pi*np.arctan2(contourdiff[:,1],contourdiff[:,0])
                contourpointlocation=contour[:-1]+2*contourdiff/5
                for ct,k in enumerate(contourpointlocation):
                    ax.scatter(k[1],k[0],marker=(3,0,contourpointrotation[ct]),c=contourcolor,s=mw_contourarrows,linewidth=0,zorder=12)
    if geocontour is not None:
        lw_geocontour=lw_geocontour*pdw
        for k in geocontour:
            ax.plot(k[1,1:3],k[0,1:3],color=geocontourcolor,linewidth=lw_geocontour,zorder=14)
        if geocontourvectors==True:
            quivlocs=geocontour[:,:,1]+(geocontour[:,:,2]-geocontour[:,:,1])/2
            ax.quiver(quivlocs[:,1],quivlocs[:,0],geocontour[:,1,4],geocontour[:,0,4],facecolor=geocontourcolor,edgecolor='black',scale=3/(latspc+lonspc),scale_units='xy',width=0.15*(latspc+lonspc)/2,headwidth=3,headlength=3,headaxislength=2.75,units='xy',zorder=13,linewidth=pdw/150)
    if outdpi=='indep':
        outdpi=200
    elif outdpi=='high':
        outdpi=ds.plotdatadpi(ax,mult=72/(latspc+lonspc),axis='xy')
        if outdpi<100:
            warnings.warn('WARNING - calculated outdpi too low, geocontour is setting outdpi to 100')
            outdpi=100
    elif outdpi=='low':
        outdpi=ds.plotdatadpi(ax,mult=10/(latspc+lonspc),axis='xy')
        if outdpi<100:
            warnings.warn('WARNING - calculated outdpi too low, geocontour is setting outdpi to 100')
            outdpi=100
    else:
        if type(outdpi) is not float and type(outdpi) is not int:
            sys.exit('ERROR - outdpi input other than \'high\'/\'low\'/\'indep\' must be a numerical value')
        else:
            outdpi=ds.plotdatadpi(ax,mult=outdpi*2/(latspc+lonspc),axis='xy')
    fig.savefig(outname,dpi=outdpi,bbox_inches='tight')
    plt.close(fig)

def save(latitudes,longitudes,boundary=None,mask=None,contour=None,contoursearch=None,geocontour=None,vertices=None,outname='save',outtype='np',maskouttxt=False,outformat='%8.3f',outdelim=' '):
    """
    Save any/all geocontour-created elements:

    - `boundary`
    - `mask`
    - `contour`
    - `contoursearch`
    - `geocontour`
    - `vertices`

    Parameters
    ----------
    latitudes : ndarray
        1D Nx1 array of latitude points (degrees)
    longitudes : ndarray
        1D Nx1 array of longitude points (degrees)
    boundary : ndarray, optional 
        2D Nx2 array of latitude/longitude points (degrees) with the
        last point equal to the first
    mask : ndarray, optional
        2D MxN bool array where M=len(`latitudes`) and
        N=len(`longitudes`)
    contour : ndarray, optional
        2D Nx2 array of ordered latitude/longitude points (degrees)
        describing the edge of a mask
    contoursearch : ndarray, optional
        2D Nx2 array of ordered latitude/longitude points (degrees)
        describing the cells searched during contour tracing
    geocontour : ndarray, optional
        3D Nx2x5 array defining a list of N contour cells (column 1),
        their edge points (columns 2,3), segment lengths (column 4), and
        outward unit vectors (column 5)
    vertices : ndarray, optional
        2D Nx2 array of latitude/longitude points (degrees)
    outname : str, default='save'
        the filename/path for the saved elements
    outtype : {'np', 'xyz'}, default='np'
        select the output filetype

            ``np``
                A numpy binary file that stores each element as a
                separate array/object
            ``xyz``
                An xyz format text file (lat, lon, [data])
    maskouttxt : bool, default=False
        select to also save a mask as a (1/0) text file
    outformat : format_str, default='%8.3f'
        formatting string for latitude and longitude values (in .xyz
        file)
    outdelim : str, default=' '
        delimiter string for columns in output (in .xyz file)

    Returns
    -------
    None
    """
    if outtype=='xyz':
        np.savetxt(outname+'_latitudes.txt',latitudes,fmt=outformat,delimiter=outdelim)
        np.savetxt(outname+'_longitudes.txt',longitudes,fmt=outformat,delimiter=outdelim)
        if boundary is not None:
            gcc.cboundary(boundary)
            np.savetxt(outname+'_boundary.xyz',boundary,fmt=outformat,delimiter=outdelim)
        if mask is not None:
            gcc.cmask(mask,latitudes,longitudes)
            ysp, xsp = np.meshgrid(latitudes,longitudes,indexing='ij')
            xyzout=np.hstack((xsp.reshape((-1,1)),ysp.reshape((-1,1)),mask.reshape((-1,1))))
            np.savetxt(outname+'_mask.xyz',xyzout,fmt=[outformat, outformat, '%3d'],delimiter=outdelim)
            if maskouttxt==True:
                np.savetxt(outname+'_mask.txt',mask.astype('int'),fmt='%1d',delimiter=outdelim)
        if contour is not None:
            gcc.ccontour(contour,latitudes,longitudes)
            np.savetxt(outname+'_contour.xyz',contour,fmt=outformat,delimiter=outdelim)
        if contoursearch is not None:
            np.savetxt(outname+'_contoursearch.xyz',contoursearch,fmt=outformat,delimiter=outdelim)
        if geocontour is not None:
            gcc.cgeocontour(geocontour,latitudes,longitudes)
            header1='cell'.center(16)+'entry'.center(19)+'exit'.center(18)+'length'.center(25)+'normal vector'.center(25)
            header2='lat'.center(8)+'lon'.center(9)+'lat'.center(9)+'lon'.center(9)+'lat'.center(9)+'lon'.center(9)+'degrees'.center(13)+'meters'.center(12)+'y'.center(13)+'x'.center(12)
            np.savetxt(outname+'_geocontour.xyz',geocontour.reshape((-1,10),order='F'),fmt=[outformat,outformat,outformat,outformat,outformat,outformat,'%12.7f','%12d','%12.7f','%12.7f'],delimiter=outdelim,header=header1+'\n'+header2)
        if vertices is not None:
            np.savetxt(outname+'_vertices.xyz',vertices,fmt=outformat,delimiter=outdelim)
    elif outtype=='np':
        np.savez(outname,latitudes=latitudes,longitudes=longitudes,boundary=boundary,mask=mask,contour=contour,contoursearch=contoursearch,geocontour=geocontour,vertices=vertices)
    else:
        sys.exit('ERROR - outtype=\''+outtype+'\' is not a valid selection, valid selections are \'np\'/\'xyz\'')

