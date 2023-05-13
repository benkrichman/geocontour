"""
Functions for comparatively testing speed of geocontour implementations
=======================================================================
"""
import sys
import numpy as np
import geocontour.masksearch as gcms
import geocontour.maskutil as gcmu
import geocontour.contourtrace as gcct
import timeit as tm
import matplotlib.pyplot as plt
import importlib.resources as ilr

def masksearch(numtests=10,runspertest=1,boundname='generic_boundary',spacing=[1.75,1.5,1.25,1,0.75,0.5,0.25,0.2],stat='min',plot=True,logax=False):
    """
    Test the timing of all mask search functions using timeit

    Parameters
    ----------
    numtests : int, default=10
        number of tests to run
    runspertest : int, default=1
        number of function calls per test
    boundname : str, default='generic_boundary'
        the name of the boundary file to use (see
        geocontour/geocontour/data/)
    spacing : array_like, default=[1.75,1.5,1.25,1,0.75,0.5,0.25,0.2]
        grid spacing value(s) to use for testing
        [i.e. if a list or numpy array is provided, testing will be run
        and results will populate for every value]
    stat : {'mean', 'median', 'min', 'max'}, default='min'
        select the statistical method for combining results from all
        tests for a given function
    plot : bool, default=True
        select to plot results
    logax : bool, default=False
        select to use log axes (if plotting)

    Returns
    -------
    results : ndarray
        2D Nx6 array of timing results in units of seconds where N is
        equal to len(`spacing`)
 
            column1: number of cells searched

            column2: center() timing

            column3: center2() timing

            column4: nodes() timing

            column5: nodes2() timing

            column6: area() timing

    See Also
    --------
    contourtrace
    geocontour.masksearch

    Notes
    -----
    The `numtests` parameter determines the number of tests that are
    later combined via the `stat` parameter. The `runspertest` parameter
    determines the number of times a function will be called in a given
    test. For example, if `numtests` = 10 and `runspertest` = 5, 50
    total function calls will be made and result in 10 test results of
    the total timing of 5 consecutive function calls.
    """
    if type(spacing) is float or type(spacing) is int:
        spacing=[spacing]
    elif type(spacing) is not np.ndarray and type(spacing) is not list:
        sys.exit('ERROR - spacing input must be single numeric (float or int) or list of numerics (numpy array or list)')
    boundfil=ilr.files('geocontour').joinpath('data/'+boundname+'.npz')
    data=np.load(boundfil)
    boundary=data['boundary']
    datalat=data['latitudes']
    datalon=data['longitudes']
    minlat=min(datalat)
    maxlat=max(datalat)
    minlon=min(datalon)
    maxlon=max(datalon)
    
    output=[]
    for sz in spacing:
        lons=np.arange(minlon,maxlon+sz,sz)
        lats=np.arange(minlat,maxlat+sz,sz)
        boxlatmin, boxlatmax, boxlonmin, boxlonmax = gcmu.bbox(lats,lons,boundary)
        numcells=((boxlatmax-boxlatmin)*(boxlonmax-boxlonmin))
        centerfunc=(lambda: gcms.center(lats,lons,boundary))
        TM=tm.Timer(centerfunc)
        centertimes=TM.repeat(numtests,runspertest)
        center2func=(lambda: gcms.center2(lats,lons,boundary))
        TM=tm.Timer(center2func)
        center2times=TM.repeat(numtests,runspertest)
        nodesfunc=(lambda: gcms.nodes(lats,lons,boundary))
        TM=tm.Timer(nodesfunc)
        nodestimes=TM.repeat(numtests,runspertest)
        nodes2func=(lambda: gcms.nodes2(lats,lons,boundary))
        TM=tm.Timer(nodes2func)
        nodes2times=TM.repeat(numtests,runspertest)
        areafunc=(lambda: gcms.area(lats,lons,boundary))
        TM=tm.Timer(areafunc)
        areatimes=TM.repeat(numtests,runspertest)
        if stat=='mean':
            centertime=np.mean(centertimes)
            center2time=np.mean(center2times)
            nodestime=np.mean(nodestimes)
            nodes2time=np.mean(nodes2times)
            areatime=np.mean(areatimes)
        elif stat=='median':
            centertime=np.median(centertimes)
            center2time=np.median(center2times)
            nodestime=np.median(nodestimes)
            nodes2time=np.median(nodes2times)
            areatime=np.median(areatimes)
        elif stat=='min':
            centertime=np.min(centertimes)
            center2time=np.min(center2times)
            nodestime=np.min(nodestimes)
            nodes2time=np.min(nodes2times)
            areatime=np.min(areatimes)
        elif stat=='max':
            centertime=np.max(centertimes)
            center2time=np.max(center2times)
            nodestime=np.max(nodestimes)
            nodes2time=np.max(nodes2times)
            areatime=np.max(areatimes)
        else:
            sys.exit('ERROR - stat=\''+stat+'\' is not a valid selection, valid selections are \'mean\'/\'median\'/\'min\'/\'max\'')
        output.append([numcells,centertime,center2time,nodestime,nodes2time,areatime])
        print('\n'+str(numcells)+' cells searched:\n   '+stat+' center time: '+str(centertime)+'\n   '+stat+' center2 time: '+str(center2time)+'\n   '+stat+' nodes time: '+str(nodestime)+'\n   '+stat+' nodes2 time: '+str(nodes2time)+'\n   '+stat+' area time: '+str(areatime))
    
    results=np.array(output)
    if plot==True:
        fig=plt.figure()
        ax=fig.add_subplot(1,1,1)
        ax.plot(results[:,0],results[:,1],label='center')
        ax.plot(results[:,0],results[:,2],label='center2')
        ax.plot(results[:,0],results[:,3],label='nodes')
        ax.plot(results[:,0],results[:,4],label='nodes2')
        ax.plot(results[:,0],results[:,5],label='area')
        ax.grid()
        if logax==True:
            ax.set_xscale('log')
            ax.set_yscale('log')
        ax.legend()
        ax.set_title('Mask Search Times\n'+stat+' of '+str(numtests)+' tests of '+str(runspertest)+' calls each')
        ax.set_ylabel('time (s)')
        ax.set_xlabel('cells searched')
        plt.savefig('test_masksearch_times')
        plt.close()
        print('\n\nTiming figure saved as \'test_masksearch_times\'')
    return results

def contourtrace(numtests=10,runspertest=1,boundname='generic_boundary',spacing=[1.75,1.5,1.25,1,0.75,0.5,0.25,0.2,0.15,0.125,0.1,0.09,0.08,0.07,0.065,0.06],stat='min',plot=True,logax=False):
    """
    Test the timing of all contour trace functions using timeit

    Parameters
    ----------
    numtests : int, default=10
        number of tests to run
    runspertest : int, default=1
        number of function calls per test
    boundname : str, default='generic_boundary'
        the name of the boundary file to use (see
        geocontour/geocontour/data/)
    spacing : array_like, default=[1.75,1.5,1.25,1,0.75,0.5,0.25,0.2]
        grid spacing value(s) to use for testing
        [i.e. if a list or numpy array is provided, testing will be run
        and results will populate for every value]
    stat : {'mean', 'median', 'min', 'max'}, default='min'
        select the statistical method for combining results from all
        tests for a given function
    plot : bool, default=True
        select to plot results
    logax : bool, default=False
        select to use log axes (if plotting)

    Returns
    -------
    results : ndarray
        2D Nx9 array of timing results in units of seconds where N is
        equal to len(`spacing`)
 
            column1: number of cells in mask

            column2: square() timing

            column3: moore() timing

            column4: moore_imp() timing

            column5: pavlidis() timing

            column6: pavlidis_imp() timing

            column7: MSBF() timing

            column8: ISBF() timing

            column9: TSR() timing

    See Also
    --------
    masksearch
    geocontour.contourtrace

    Notes
    -----
    The `numtests` parameter determines the number of tests that are
    later combined via the `stat` parameter. The `runspertest` parameter
    determines the number of times a function will be called in a given
    test. For example, if `numtests` = 10 and `runspertest` = 5, 50
    total function calls will be made and result in 10 test results of
    the total timing of 5 consecutive function calls.
    """
    if type(spacing) is float or type(spacing) is int:
        spacing=[spacing]
    elif type(spacing) is not np.ndarray and type(spacing) is not list:
        sys.exit('ERROR - spacing input must be single numeric (float or int) or list of numerics (numpy array or list)')
    boundfil=ilr.files('geocontour').joinpath('data/'+boundname+'.npz')
    data=np.load(boundfil)
    boundary=data['boundary']
    datalat=data['latitudes']
    datalon=data['longitudes']
    minlat=min(datalat)
    maxlat=max(datalat)
    minlon=min(datalon)
    maxlon=max(datalon)
    
    output=[]
    for sz in spacing:
        lons=np.arange(minlon,maxlon+sz,sz)
        lats=np.arange(minlat,maxlat+sz,sz)
        mask=gcms.center2(lats,lons,boundary)
        numcells=mask.sum()
        squarefunc=(lambda: gcct.square(mask,lats,lons))
        TM=tm.Timer(squarefunc)
        squaretimes=TM.repeat(numtests,runspertest)
        moorefunc=(lambda: gcct.moore(mask,lats,lons))
        TM=tm.Timer(moorefunc)
        mooretimes=TM.repeat(numtests,runspertest)
        moore_impfunc=(lambda: gcct.moore_imp(mask,lats,lons))
        TM=tm.Timer(moore_impfunc)
        moore_imptimes=TM.repeat(numtests,runspertest)
        pavlidisfunc=(lambda: gcct.pavlidis(mask,lats,lons))
        TM=tm.Timer(pavlidisfunc)
        pavlidistimes=TM.repeat(numtests,runspertest)
        pavlidis_impfunc=(lambda: gcct.pavlidis_imp(mask,lats,lons))
        TM=tm.Timer(pavlidis_impfunc)
        pavlidis_imptimes=TM.repeat(numtests,runspertest)
        MSBFfunc=(lambda: gcct.MSBF(mask,lats,lons))
        TM=tm.Timer(MSBFfunc)
        MSBFtimes=TM.repeat(numtests,runspertest)
        ISBFfunc=(lambda: gcct.ISBF(mask,lats,lons))
        TM=tm.Timer(ISBFfunc)
        ISBFtimes=TM.repeat(numtests,runspertest)
        TSRfunc=(lambda: gcct.TSR(mask,lats,lons))
        TM=tm.Timer(TSRfunc)
        TSRtimes=TM.repeat(numtests,runspertest)
        if stat=='mean':
            squaretime=np.mean(squaretimes)
            mooretime=np.mean(mooretimes)
            moore_imptime=np.mean(moore_imptimes)
            pavlidistime=np.mean(pavlidistimes)
            pavlidis_imptime=np.mean(pavlidis_imptimes)
            MSBFtime=np.mean(MSBFtimes)
            ISBFtime=np.mean(ISBFtimes)
            TSRtime=np.mean(TSRtimes)
        elif stat=='median':
            squaretime=np.median(squaretimes)
            mooretime=np.median(mooretimes)
            moore_imptime=np.median(moore_imptimes)
            pavlidistime=np.median(pavlidistimes)
            pavlidis_imptime=np.median(pavlidis_imptimes)
            MSBFtime=np.median(MSBFtimes)
            ISBFtime=np.median(ISBFtimes)
            TSRtime=np.median(TSRtimes)
        elif stat=='min':
            squaretime=np.min(squaretimes)
            mooretime=np.min(mooretimes)
            moore_imptime=np.min(moore_imptimes)
            pavlidistime=np.min(pavlidistimes)
            pavlidis_imptime=np.min(pavlidis_imptimes)
            MSBFtime=np.min(MSBFtimes)
            ISBFtime=np.min(ISBFtimes)
            TSRtime=np.min(TSRtimes)
        elif stat=='max':
            squaretime=np.max(squaretimes)
            mooretime=np.max(mooretimes)
            moore_imptime=np.max(moore_imptimes)
            pavlidistime=np.max(pavlidistimes)
            pavlidis_imptime=np.max(pavlidis_imptimes)
            MSBFtime=np.max(MSBFtimes)
            ISBFtime=np.max(ISBFtimes)
            TSRtime=np.max(TSRtimes)
        else:
            sys.exit('ERROR - stat=\''+stat+'\' is not a valid selection, valid selections are \'mean\'/\'median\'/\'min\'/\'max\'')
        output.append([numcells,squaretime,mooretime,moore_imptime,pavlidistime,pavlidis_imptime,MSBFtime,ISBFtime,TSRtime])
        print('\n'+str(numcells)+' cells in mask:\n   '+stat+' square time: '+str(squaretime)+'\n   '+stat+' moore time: '+str(mooretime)+'\n   '+stat+' moore_imp time: '+str(moore_imptime)+'\n   '+stat+' pavlidis time: '+str(pavlidistime)+'\n   '+stat+' pavlidis_imp time: '+str(pavlidis_imptime)+'\n   '+stat+' MSBF time: '+str(MSBFtime)+'\n   '+stat+' ISBF time: '+str(ISBFtime)+'\n   '+stat+' TSR time: '+str(TSRtime))
    
    results=np.array(output)
    if plot==True:
        fig=plt.figure()
        ax=fig.add_subplot(1,1,1)
        ax.plot(results[:,0],results[:,1],label='square')
        ax.plot(results[:,0],results[:,2],label='moore')
        ax.plot(results[:,0],results[:,3],label='moore_imp')
        ax.plot(results[:,0],results[:,4],label='pavlidis')
        ax.plot(results[:,0],results[:,5],label='pavlidis_imp')
        ax.plot(results[:,0],results[:,6],label='MSBF')
        ax.plot(results[:,0],results[:,7],label='ISBF')
        ax.plot(results[:,0],results[:,8],label='TSR')
        ax.grid()
        if logax==True:
            ax.set_xscale('log')
            ax.set_yscale('log')
        ax.legend()
        ax.set_title('Contour Trace Times\n'+stat+' of '+str(numtests)+' tests of '+str(runspertest)+' calls each')
        ax.set_ylabel('time (s)')
        ax.set_xlabel('cells in mask')
        plt.savefig('test_contourtrace_times')
        plt.close()
        print('\n\nTiming figure saved as \'test_contourtrace_times\'')
    return results

