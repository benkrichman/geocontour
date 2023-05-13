<h1><img align="left" src="https://github.com/benkrichman/geocontour/raw/main/images/icon_geocontour.png" width="130" height="130">geocontour</h1>

Utilities for masking, contour tracing, and geocontour construction with gridded geographic data.

\
[![DOI](https://zenodo.org/badge/550241733.svg)](https://zenodo.org/badge/latestdoi/550241733)
[![Downloads](https://pepy.tech/badge/geocontour)](https://pepy.tech/project/geocontour)
[![PyPI version](https://badge.fury.io/py/geocontour.svg)](https://badge.fury.io/py/geocontour)
[![Documentation Status](https://readthedocs.org/projects/geocontour/badge/?version=latest)](https://geocontour.readthedocs.io/en/latest/?badge=latest)


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
## Documentation

Find the full documentation hosted [here](https://geocontour.readthedocs.io/en/latest/index.html).

## Citation

If geocontour played a significant role in your work and you would like to cite it, the following is suggested (APA):

Krichman, B. (2023). *geocontour* (Version 1.3.0) [Computer Software]. https://doi.org/10.5281/zenodo.7402714

Bibtex:

```latex
@software{geocontour,
author={Krichman, Benjamin},
doi={10.5281/zenodo.7402714},
license={MIT},
month={5},
year={2023},
title={{geocontour}},
url={https://github.com/benkrichman/geocontour},
version={1.3.0},
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

Implements 6 existing algorithms for contour tracing, and two improvements on known algorithms
- square tracing (a.k.a. simple boundary follower/Papert's turtle algorithm) [^3][^4][^5][^7][^8]
- moore neighbor tracing [^3][^8]
- improved moore neighbor tracing (capturing inside corners)
- pavlidis tracing [^3][^6]
- improved pavlidis tracing (capturing inside corners)
- modified simple boundary follower [^1][^2][^4]
- improved simple boundary follower [^1][^2][^7]
- fast representative tracing [^7]

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
<img src=https://github.com/benkrichman/geocontour/raw/main/images/example_small_contoursearch.png width="340" height="340"/>
<img src=https://github.com/benkrichman/geocontour/raw/main/images/example_small_contour.png width="340" height="340"/>
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
<img src=https://github.com/benkrichman/geocontour/raw/main/images/example_small_geocontour.png width="340" height="340"/>
<img src=https://github.com/benkrichman/geocontour/raw/main/images/example_small_geocontour_simp.png width="340" height="340"/>
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

[^1]: Cheong, C.-H., & Han, T.-D. (2006). Improved Simple Boundary Following Algorithm. Journal of KIISE: Software and Applications, 33(4), 427–439. <https://koreascience.kr/article/JAKO200622219415761.pdf>

[^2]: Cheong, C.-H., Seo, J., & Han, T.-D. (2006). Advanced Contour Tracing Algorithms based on Analysis of Tracing Conditions. Proceedings of the 33rd KISS Fall Conference, 33, 431–436. <https://koreascience.kr/article/CFKO200614539217302.pdf>

[^3]: Ghuneim, A.G. (2000). *Contour Tracing*. McGill University. <https://www.imageprocessingplace.com/downloads_V3/root_downloads/tutorials/contour_tracing_Abeer_George_Ghuneim/alg.html>

[^4]: Gose, E., Johnsonbaugh, R., & Jost, S. (1996). Pattern recognition and image analysis. Prentice Hall PTR.

[^5]: Papert, S. (1973). Uses of Technology to Enhance Education (No. 298; Artificial Intelligence). Massachusetts Institute of Technology. <https://dspace.mit.edu/handle/1721.1/6213>

[^6]: Pavlidis, T. (1982) Algorithms for Graphics and Image Processing. Computer Science Press, New York, NY. <https://doi.org/10.1007/978-3-642-93208-3>

[^7]: Seo, J., Chae, S., Shim, J., Kim, D., Cheong, C., & Han, T.-D. (2016). Fast Contour-Tracing Algorithm Based on a Pixel-Following Method for Image Sensors. Sensors, 16(3), 353. <https://doi.org/10.3390/s16030353>

[^8]: Toussaint, G.T. (2010). *Grids Connectivity and Contour Tracing* [Lesson Notes]. McGill University. <http://www-cgrl.cs.mcgill.ca/~godfried/teaching/mir-reading-assignments/Chapter-2-Grids-Connectivity-Contour-Tracing.pdf>

