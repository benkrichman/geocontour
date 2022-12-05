<h1><img align="left" src="https://github.com/benkrichman/geocontour/raw/main/images/icon_geocontour.png" width="130" height="130">geocontour</h1>

Utilities for masking, contour tracing, and geocontour construction for flux calculations from gridded geographic data.

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

## Features

### Masks

Selectable criteria for masks created from input boundary coordinates
- cell center
- area ratio
- node ratio

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

### Geocontours

From an input contour, create a closed geospatial contour with calculated segment lengths and outward unit vectors (for example: useful in calculating flux across a bounding surface from a geospatial data set)

Options for tuning criteria of geocontours created from input contours
- selectable connection type (cell to cell or cell edge to cell center)
- optionally simplify geocontours at the cell level to shorten and improve compute times in practical applications

### Visualization

Easy and semi-automated plotting function for visualization of boundaries/masks/contours/contour searches/geocontours
- buffers
- grid overlay
- mask cell visibility
- directional indicators for contours and contour searches
- outward unit vector indicators for geocontours
- automatic calculation of feature size and output resolution
- display of natural features or political boundaries (optional with cartopy installed)
- selectable marker/line/arrow/cell size/color/style

## Example Use Case

\*to reconstruct these examples use (or view)
```python
geocontour.examples.small()
geocontour.examples.large()
```

### mask search

Given a series of lat/lon points constituting a geographical boundary, and a set of gridded data on a lat/lon grid, find an appropriate mask to select gridded data within the boundary:

Use the 'area' approach to mask calculation, defaulting to selection of all cells for which 50% or greater falls withing the boundary. Note that boundary falls outside gridded data bounds at some points.
```python
mask=geocontour.masksearch.area(latitudes,longitudes,boundary)
```

![geocontour.masksearch.area() example](https://github.com/benkrichman/geocontour/raw/main/images/example_small_boundary%2Bmask.png?raw=true)

<img src=https://github.com/benkrichman/geocontour/raw/main/images/example_small_boundary%2Bmask.png width="800" height="800">

### contour trace

Given the previously calculated mask, find the outer edge using a contour tracing algorithm:

Use the improved Pavlidis algorithm to trace the contour.
```python
contour,contoursearch=geocontour.contourtrace.pavlidis_imp(mask,latitudes,longitudes)
```
![geocontour.contourtrace.pavlidis_imp() example contour](https://github.com/benkrichman/geocontour/raw/main/images/example_small_contoursearch.png?raw=true)

![geocontour.contourtrace.pavlidis_imp() example search](https://github.com/benkrichman/geocontour/raw/main/images/example_small_contour.png?raw=true)

### construct geocontour

Given the previously calculated contour, construct the geocontour to determine contour segment lengths and outward normal vectors:

Use the build function of geocontour to construct the geocontour. Note that the 'simplify' option is used, combining cells with multiple visits into single segments.
```python
geocontour=geocontour.build(contour,latitudes,longitudes,simplify=True)
```

![geocontour.build() example](https://github.com/benkrichman/geocontour/raw/main/images/example_small_geocontour.png?raw=true)

### project geocontour against map features

Given a large geocontour (in this case, the Mississippi River Basin) project against natural features and political borders:

```python
geocontour.output.plot(latitudes,longitudes,geocontour=geocontour,features='natural')
```
![geocontour.output.plot() nat example](https://github.com/benkrichman/geocontour/raw/main/images/example_large_geocontour%2Bnatfeat.png?raw=true)

```python
geocontour.output.plot(latitudes,longitudes,geocontour=geocontour,features='borders')
```
![geocontour.output.plot() bord example](https://github.com/benkrichman/geocontour/raw/main/images/example_large_geocontour%2Bbordfeat.png?raw=true)

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
  - Functionally matches geocontour.masksearch.center(), but utilizes matplotlib.path functions, which are probably optimized and thus is roughly 2.5\*sqrt(N) faster for N points, though lacks a "precision" buffer input

#### geocontour.masksearch.nodes()
Returns a mask over a range of input latitudes and longitudes determined by an input boundary
  - Critera for inclusion of a cell is whether a given number (default=2) of cell nodes (corners) fall within the boundary

#### geocontour.masksearch.nodes2()
Returns a mask over a range of input latitudes and longitudes determined by an input boundary
  - Critera for inclusion of a cell is whether a given number (default=2) of cell nodes (corners) fall within the boundary 
  - Functionally matches geocontour.masksearch.nodes(), but utilizes matplotlib.path functions, though speed is similar to the shapely implementation

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

### geocontour

#### geocontour.build()
Returns a geocontour from a contour input

### output

#### geocontour.output.plot()
Plots any/all maskpy-created elements: boundary, mask, contour, contoursearch, geocontour, vertices

#### geocontour.output.save()
Saves any/all maskpy-created elements: boundary, mask, contour, contoursearch, geocontour, vertices

### tests

#### geocontour.tests.full()
Runs all user-facing geocontour functions with test data, printing/saving results

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


[^IPP]:Contour Tracing Algorithms, [Pattern Recognition Project](https://www.imageprocessingplace.com/downloads_V3/root_downloads/tutorials/contour_tracing_Abeer_George_Ghuneim/alg.html)
[^Toussaint]:Grids Connectivity and Contour Tracing, [Lesson Notes](http://www-cgrl.cs.mcgill.ca/~godfried/teaching/mir-reading-assignments/Chapter-2-Grids-Connectivity-Contour-Tracing.pdf)
[^Pavlidis]:Algorithms for Graphics and Image Processing, [doi:10.1007/978-3-642-93208-3](https://link.springer.com/book/10.1007/978-3-642-93208-3)
[^FRT]:Fast Contour-Tracing Algorithm Based on a Pixel-Following Method for Image Sensors, [doi:10.3390/s16030353](https://www.mdpi.com/1424-8220/16/3/353)
[^Kovalevsky]:Other Source to Note: [Vladimir Kovalevsky](http://www.kovalevsky.de/index.htm)
[^Osborne]:The Mercator Projections [doi:10.5281/ZENODO.35392](https://zenodo.org/record/35392)
