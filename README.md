<h1><img align="left" src="https://github.com/benkrichman/geocontour/raw/main/images/icon_geocontour.png" width="130" height="130">geocontour</h1>

Utilities for masking, contour tracing, and geocontour construction for flux calculations from gridded geographic data.

\
[![DOI](https://zenodo.org/badge/550241733.svg)](https://zenodo.org/badge/latestdoi/550241733)
[![Downloads](https://pepy.tech/badge/geocontour)](https://pepy.tech/project/geocontour)
[![PyPI version](https://badge.fury.io/py/geocontour.svg)](https://badge.fury.io/py/geocontour)

## Installation

```bash
pip install geocontour
```

or

```bash
pip install git+https://github.com/benkrichman/geocontour.git@main
```

To run a full test of internal functions (minus cartopy features):

```python
geocontour.tests.full()
```

## Citation

If geocontour played a significant role in your work and you would like to cite it, the following is suggested (APA):

Krichman, B. (2023). *geocontour* (Version 1.2.1) [Computer Software]. https://doi.org/10.5281/zenodo.7707058

Bibtex:

```latex
@software{geocontour,
author={Krichman, Benjamin},
doi={10.5281/zenodo.7707058},
license={MIT},
month={3},
year={2023},
title={{geocontour}},
url={https://github.com/benkrichman/geocontour},
version={1.2.1},
note = {Computer Software}
}
```

## Features

### Masks

Selectable/tunable criteria for masks created from input boundary coordinates
- cell center (multiple methods with variable precision)
- node ratio (multiple methods with variable precision)
- area ratio

Useful mask operators
- return mask connectivity (and null connectivity)
- return mask edge cells
- return mask vertex points

### Contours

Implements 4 existing algorithms for contour tracing, and two improvements on known algorithms
- square tracing [^IPP][^Toussaint]
- moore neighbor tracing [^IPP][^Toussaint]
- improved moore neighbor tracing (capturing inside corners)
- pavlidis tracing [^IPP][^Pavlidis]
- improved pavlidis tracing (capturing inside corners)
- fast representative tracing [^FRT]

Tuning of contours created from tracing input masks
- trace direction
- selectable and adjustable stopping conditions
- automatic or manual selection of starting cell
- selectable connection type (cell to cell or cell edge to cell center)
- simplification of output contour (removal of repeating cells)
- selectable contour closure
- usable for an associated lat/lon grid or on a non-specified grid

Useful contour operators
- return full search path for a contour trace
- return cell neighbors with connectivity and directional input
- return starting cell for contour tracing and check that starting cells work for a given algorithm
- return visually improved contour search

### Geocontours

From an input contour, create a closed geospatial contour with calculated segment lengths and outward unit vectors (for example: useful in calculating flux across a bounding surface from a geospatial data set)

Options for tuning criteria of geocontours created from input contours
- selectable connection type (cell to cell or cell edge to cell center)
- optionally simplify geocontours at the cell level to shorten and improve compute times in practical applications

### Timing

Timing modules for easy comparison between mask search methods or contour tracing algorithms using timeit. 

Note that in mask search and contour tracing care has been taken to implement algorithms in a fast and efficient manner through utilization of shapely and matplotlib builtins and through numpy vectorization where possible. However, not everything is speed optimized where optimization would necessitate significantly more complexity or utilization of external low level libraries or custom functions. The timing modules exist for intercomparison amongst methods, but also for giving users a reasonable expectation of performance.

### Visualization

Easy and semi-automated plotting function for visualization of boundaries/masks/contours/contour searches/geocontours
- buffers
- grid overlay
- mask/contour cell visibility
- directional indicators for contours and contour searches
- outward unit vector indicators for geocontours
- automatic calculation of feature size and output resolution
- display of natural features or political boundaries (optional with cartopy installed)
- selectable marker/line/arrow/cell size/color/style
- optional transparency mode for presentation/publication use

## Example Use Case

\*to reconstruct these examples use (or view)
```python
geocontour.examples.small()
geocontour.examples.large()
```

### mask search

Given a series of lat/lon points constituting a geographical boundary, and a set of gridded data on a lat/lon grid, find an appropriate mask to select gridded data within the boundary:

Use the 'area' approach to mask calculation, defaulting to selection of all cells for which 50% or greater falls withing the boundary. Note that boundary falls outside gridded data bounds at some points and those cells inside the boundary but outside the gridded data bounds are not included in the mask.
```python
mask=geocontour.masksearch.area(latitudes,longitudes,boundary)
geocontour.output.plot(latitudes,longitudes,boundary=boundary,mask=mask,title='Example Mask and Boundary',outname='example_small_boundary+mask',outdpi='indep')
```

<img src=https://github.com/benkrichman/geocontour/raw/main/images/example_small_boundary%2Bmask.png width="450" height="450">

### contour trace

Given the previously calculated mask, find the outer edge using a contour tracing algorithm:

Use the improved Pavlidis algorithm to trace the contour. Note that the contoursearch plot shows the start cell as a circle, directional arrows for each segment, and diamonds where cells are consecutively and repeatedly searched. In the case of the pavlidis algorithm these diamonds show where the orientation turned 90 degrees. Similarly the contour plot uses a circle to mark the start cell and arrows to signify direction.
```python
contour,contoursearch=geocontour.contourtrace.pavlidis_imp(mask,latitudes,longitudes)
geocontour.output.plot(latitudes,longitudes,mask=mask,contoursearch=contoursearch,title='Example Contour Search',outname='example_small_contoursearch',outdpi='indep')
geocontour.output.plot(latitudes,longitudes,contour=contour,cells='contour',title='Example Contour',outname='example_small_contour',outdpi='indep')
```

<p float="middle">
<img src=https://github.com/benkrichman/geocontour/raw/main/images/example_small_contoursearch.png width="400" height="400"/>
<img src=https://github.com/benkrichman/geocontour/raw/main/images/example_small_contour.png width="400" height="400"/>
</p>

### construct geocontour

Given the previously calculated contour, construct the geocontour to determine contour segment lengths and outward normal vectors:

Use the build function of geocontour to construct the geocontour. Note that in the second plot the 'simplify' option is used, combining cells with multiple visits into single segments exactly equal to the vector combination of segments in the cell. The directional information contained in the contour has been discarded, and in the case of simplification may not be extractable from the geocontour.
```python
geocontour=geocontour.build(contour,latitudes,longitudes)
geocontour_simp=geocontour.build(contour,latitudes,longitudes,simplify=True)
geocontour.output.plot(latitudes,longitudes,geocontour=geocontour,buffer='on',title='Example Geocontour',outname='example_small_geocontour',outdpi='indep')
geocontour.output.plot(latitudes,longitudes,geocontour=geocontour_simp,buffer='on',title='Example Geocontour - Simplified',outname='example_small_geocontour_simp',outdpi='indep')
```

<p float="middle">
<img src=https://github.com/benkrichman/geocontour/raw/main/images/example_small_geocontour.png width="400" height="400"/>
<img src=https://github.com/benkrichman/geocontour/raw/main/images/example_small_geocontour_simp.png width="400" height="400"/>
</p>

### project geocontour against map features

Given a large geocontour (in this case, the Mississippi River Basin) project against natural features and political borders (requires cartopy):

```python
geocontour.output.plot(latitudes,longitudes,geocontour=geocontour,title='Example Geocontour\nMississippi River Basin',outname='example_large_geocontour+natfeat',features='natural')
```
<img src=https://github.com/benkrichman/geocontour/raw/main/images/example_large_geocontour%2Bnatfeat.png width="800">

```python
geocontour.output.plot(latitudes,longitudes,geocontour=geocontour,title='Example Geocontour\nMississippi River Basin',outname='example_large_geocontour+bordfeat',features='borders')
```

<img src=https://github.com/benkrichman/geocontour/raw/main/images/example_large_geocontour%2Bbordfeat.png width="800">

## Function Overview

\*to see full function documentation use
```python
help(geocontour.module.function)
```

### check

#### geocontour.check.cdim()
Checks an input dimension array for 1-dimensionality and regular spacing

#### geocontour.check.cboundary()
Checks a list of boundary points for 2-dimensionality and proper ordering

#### geocontour.check.cmask()
Checks a mask for correct data type and dimensionality, and size if optional latitudes and longitudes are provided

#### geocontour.check.ccontour()
Check contour for repeating cells, closure, and connectivity, and latitude/longitude range if optional latitudes and longitudes are provided

#### geocontour.check.cgeocontour()
Check geocontour for latitude/longitude range and dimension

### grid

#### geocontour.grid.spacing()
Returns the grid spacing for a given input dimension

#### geocontour.grid.lonlens()
Returns the lengths of a degree (default) of longitude over a range of latitudes [^Osborne]

#### geocontour.grid.latlens()
Returns the grid lengths of a defined range of latitudes [^Osborne]

#### geocontour.grid.lonlen()
Returns the length of a degree of longitude at the input latitude [^Osborne]

#### geocontour.grid.latlen()
Returns the length of a degree of latitude at the input latitude [^Osborne]

#### geocontour.grid.areas()
Returns the cell areas of a grid defined by a range of latitudes and longitudes

#### geocontour.grid.clonrng()
Returns a descriptor for the range of a set of longitude points
  * negative (-180 to 180), positive (0 to 360), or indeterminate (0 to 180) range

#### geocontour.grid.clatdir()
Returns a descriptor for the direction of a set of latitude points (increasing or decreasing)

#### geocontour.grid.switchlon()
Returns a set of longitude points switched in place between negative (-180 to 180) and positive (0 to 360)

#### geocontour.grid.switchind()
Returns the index where a longitude array either crosses 0 or 180 degrees

### masksearch

#### geocontour.masksearch.center()
Returns a mask over a range of input latitudes and longitudes determined by an input boundary
  - Critera for inclusion of a cell is whether the center of the cell falls within the boundary

#### geocontour.masksearch.center2()
Returns a mask over a range of input latitudes and longitudes determined by an input boundary
  - Critera for inclusion of a cell is whether the center of the cell falls within the boundary
  - Functionally matches geocontour.masksearch.center(), but utilizes matplotlib.path functions, which are faster (possibly due to avoidance of overhead in converting to shapely geometries)

#### geocontour.masksearch.nodes()
Returns a mask over a range of input latitudes and longitudes determined by an input boundary
  - Critera for inclusion of a cell is whether a given number (default=2) of cell nodes (corners) fall within the boundary

#### geocontour.masksearch.nodes2()
Returns a mask over a range of input latitudes and longitudes determined by an input boundary
  - Critera for inclusion of a cell is whether a given number (default=2) of cell nodes (corners) fall within the boundary 
  - Functionally matches geocontour.masksearch.nodes(), but utilizes matplotlib.path functions, which are faster (possibly due to avoidance of overhead in converting to shapely geometries)

#### geocontour.masksearch.area()
Returns a mask over a range of input latitudes and longitudes determined by an input boundary
  - Critera for inclusion of a cell is whether the area of the cell enclosed by the boundary is greater than some fraction (default=0.5) 

### maskutil

#### geocontour.maskutil.bbox()
Checks input dimensions (lat/lon) against input boundary and returns min/max indicies of bounding box

#### geocontour.maskutil.edge()
Returns a mask of only the edge cells, and if latitudes and longitudes are provided also returns an array of the edge cells

#### geocontour.maskutil.vertex()
Returns the vertex points of all cells in the input mask, and the vertex points of only the mask edge

#### geocontour.maskutil.neighbors()
Returns the neighbors of a cell, with selected connectivity and direction

#### geocontour.maskutil.conn()
Returns whether a mask or its inverse are connected

### contourtrace

#### geocontour.contourtrace.square()
Returns the contour trace of a mask input using the square tracing algorithm [^IPP][^Toussaint]

#### geocontour.contourtrace.moore()
Returns the contour trace of a mask input using the Moore neighbor tracing algorithm [^IPP][^Toussaint]

#### geocontour.contourtrace.moore_imp()
Returns the contour trace of a mask input using an improved Moore neighbor tracing algorithm 
  - Captures inside corners missed by Moore neighbor tracing

#### geocontour.contourtrace.pavlidis()
Returns the contour trace of a mask input using the Pavlidis tracing algorithm [^IPP][^Pavlidis]

#### geocontour.contourtrace.pavlidis_imp()
Returns the contour trace of a mask input using an improved Pavlidis tracing algorithm
  - Captures inside corners missed by Pavlidis tracing

#### geocontour.contourtrace.TSR()
Returns the contour trace of a mask input using two-step representative tracing [^FRT]

### contourutil

#### geocontour.contourutil.findstart()
Returns a starting cell for a contour, given a mask and a search criteria

#### geocontour.contourutil.parsestart()
Checks start input for contour tracing
  - Mainly used internally for contour trace functions

#### geocontour.contourutil.setstop()
Returns a stopping function for use in contour tracing while loop
  - Mainly used internally for contour trace functions

#### geocontour.contourutil.clean()
Returns a cleaned contour that will pass checks
  - Mainly used internally for contour trace functions

#### geocontour.contourutil.fancysearch()
Returns a contoursearch that is visually more easy to follow

### geocontour

#### geocontour.build()
Returns a geocontour from a contour input

### output

#### geocontour.output.plot()
Plots any/all geocontour-created elements: boundary, mask, contour, contoursearch, geocontour, vertices

#### geocontour.output.save()
Saves any/all geocontour-created elements: boundary, mask, contour, contoursearch, geocontour, vertices

### tests

#### geocontour.tests.test.full()
Runs all user-facing geocontour functions with test data, printing/saving results

#### geocontour.tests.timing.masksearch()
Tests the timing of all mask search functions (geocontour.masksearch) using timeit

#### geocontour.tests.timing.contourtrace()
Tests the timing of all contour trace functions (geocontour.contourtrace) using timeit

### examples

#### geocontour.examples.small()
Runs a small scale example of geocontour processing using mock data, saves resulting plots to run directory
  - find mask using area criteria (0.5) and plot boundary/mask
  - trace contour using improved pavlidis algorithm and plot resultant contour and contour search
  - compute geocontour from contour and plot, using simplify option

#### geocontour.examples.large()
Runs a large scale example of geocontour processing using the Mississippi River Basin boundary, saves resulting plots to run directory
  - find mask using area criteria (0.5) and plot boundary/mask
  - trace contour using improved pavlidis algorithm and plot resultant contour and contour search
  - compute geocontour from contour and plot, using simplify option
  - plot geocontour with cartopy background options (borders and physical features) - will error and exit if cartopy not installed
  - plots will be large (dpi is auto-calculated to create enough resolution to zoom in for diagnostic use - this setting can be changed for quick plotting)


[^IPP]:Ghuneim, A.G. (2000). *Contour Tracing*. McGill University. <https://www.imageprocessingplace.com/downloads_V3/root_downloads/tutorials/contour_tracing_Abeer_George_Ghuneim/alg.html>

[^Toussaint]:Toussaint, G.T. (2010). *Grids Connectivity and Contour Tracing* [Lesson Notes]. McGill University. <http://www-cgrl.cs.mcgill.ca/~godfried/teaching/mir-reading-assignments/Chapter-2-Grids-Connectivity-Contour-Tracing.pdf>

[^Pavlidis]:Pavlidis, T. (1982) Algorithms for Graphics and Image Processing. Computer Science Press, New York, NY. <https://doi.org/10.1007/978-3-642-93208-3>

[^FRT]:Seo, J., Chae, S., Shim, J., Kim, D., Cheong, C., & Han, T.-D. (2016). Fast Contour-Tracing Algorithm Based on a Pixel-Following Method for Image Sensors. Sensors, 16(3), 353. <https://doi.org/10.3390/s16030353>

[^Osborne]:Osborne, P. (2013). The Mercator Projections. Zenodo. <https://doi.org/10.5281/zenodo.35392>
