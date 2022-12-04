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

### mask search

Given a series of lat/lon points constituting a geographical boundary, and a set of gridded data on a lat/lon grid, find an appropriate mask to select gridded data within the boundary:

Use the 'area' approach to mask calculation, defaulting to selection of all cells for which 50% or greater falls withing the boundary. Note that boundary falls outside gridded data bounds at some points.
```python
mask=geocontour.masksearch.area(latitudes,longitudes,boundary)
```

![geocontour.masksearch.area() example](https://github.com/benkrichman/geocontour/raw/main/images/example_small_boundary%2Bmask.png?raw=true)

### contour trace

Given the previously calculated mask, find the outer edge using a contour tracing algorithm:

Use the improved Pavlidis algorithm to trace the contour.
```python
contour,contoursearch=geocontour.contourtrace.pavlidis_imp(mask,latitudes,longitudes)
```
![geocontour.contourtrace.pavlidis_imp() example contour](https://github.com/benkrichman/geocontour/raw/main/images/example_small_contoursearch.png?raw=true)

![geocontour.contourtrace.pavlidis_imp() example search](https://github.com/benkrichman/geocontour/raw/main/images/example_small_contour.png?raw=true)


##Function Overview

[^IPP]:Contour Tracing Algorithms, [Pattern Recognition Project](https://www.imageprocessingplace.com/downloads_V3/root_downloads/tutorials/contour_tracing_Abeer_George_Ghuneim/alg.html)
[^Toussaint]:Grids Connectivity and Contour Tracing, [Lesson Notes](http://www-cgrl.cs.mcgill.ca/~godfried/teaching/mir-reading-assignments/Chapter-2-Grids-Connectivity-Contour-Tracing.pdf)
[^Pavlidis]:Algorithms for Graphics and Image Processing, [doi:10.1007/978-3-642-93208-3](https://link.springer.com/book/10.1007/978-3-642-93208-3)
[^FRT]:Fast Contour-Tracing Algorithm Based on a Pixel-Following Method for Image Sensors, [doi:10.3390/s16030353](https://www.mdpi.com/1424-8220/16/3/353)
[^Kovalevsky]:Other Source to Note: [Vladimir Kovalevsky](http://www.kovalevsky.de/index.htm)
