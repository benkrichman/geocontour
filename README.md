<h1><img align="left" src="https://github.com/benkrichman/geocontour/raw/main/images/icon_geocontour.png" width="200" height="200">geocontour</h1>

Utilities for masking, contour tracing, and geocontour construction for flux calculations from gridded geographic data.


## Masks

Options for tuning criteria of masks created from input boundary coordinates
- cell center
- area ratio
- node ratio

Useful mask operators
- return mask connectivity (and null connectivity)
- return mask edge cells
- return mask vertex points

## Contours

Implements 4 existing algorithms for contour tracing, and two improvements on known algorithms
- square tracing
- moore neighbor tracing
- improved moore neighbor tracing (capturing inside corners)
- pavlidis tracing
- improved pavlidis tracing (capturing inside corners)
- fast representative tracing (see doi:10.3390/s16030353)

Options for tuning critera of contours created from tracing input masks
- trace direction
- selectable and tunable stopping conditions
- automatic or manual selection of starting cell
- selectable connection type (cell to cell or cell edge to cell center)
- simplification of output contour (removal of repeating cells)
- selectable contour closure
- usable for an associated lat/lon grid or on a non- specified grid

Useful contour operators
- return full search path for a contour trace
- return cell neighbors with connectivity and directional input
- return starting cell for contour tracing and check that starting cells work for a given algorithm

## Geocontours

From an input contour, create a closed geospatial contour with calculated segment lengths and outward unit vectors (for example: useful in calculating flux across a bounding surface from a geospatial data set)

Options for tuning criteria of geocontours created from input contours
- selectable connection type (cell to cell or cell edge to cell center)
- simplify geocontours at the cell level to shorten and improve compute times in practical applications

## Visualization

Easy and semi-automated plotting function for visualization of boundaries/masks/contours/contour searches/geocontours
- buffers
- grid overlay
- mask cell visibility
- directional indicators for contours and contour searches
- outward unit vector indicators for geocontours
- automatic calculation of output resolution
- display of natural features or political boundaries (optional with cartopy installed)
- selectable marker/line/arrow/cell size/color/style



