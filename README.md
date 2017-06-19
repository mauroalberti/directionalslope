# DirectionalSlope

DirectionalSlope is a Python plugin for QGIS, for calculating directional slopes from DEM. 


## Aim

The slope of a topographic surface is usually calculated as the maximum gradient of a surface at a given point. In other cases, we could be interested in the slope along a constant orientation in space, i.e., directional slope. It is typically used to highlight structures perpendicular to the direction of analysis. A generalization of this concept is the slope along directions that may vary in space, based on external parameters, such as the flows of wind or ice. 

This algorithm calculates the slope of a surface, both directional and steepest. Directions may be uniform or variable in space. 

![alt text](/ims/directionalslope_gui.png "DirectionaSlope plugin interface")

*Fig. 1. DirectionaSlope plugin interface.*


## Methodology

The estimation of the directional slope for surfaces sampled along grids relies on the adaptation of the analytical techniques for continuous and differentiable surfaces. Usually we derive gradient by applying a 3x3 kernel to each cell (except those at the border). The figure below illustrates the adopted convention for cell indexing. 

![alt text](/ims/matrix-grid_small.jpg "Kernel for slope calculation")

*Fig. 2. 3x3 kernel for slope calculation.*


To determine the directional slope, we must first calculate the directional gradients along the two principal axes, x and y. From the values of dz / dx and dz / dy, one can determine the slope of the plane that approximates the surface at the local cell (i, j) using the formula (Neteler & Mitasova, 2008, eq. A.27; see also Geospatial Analysis - a comprehensive guide):


*directional slope (alpha) = arctan [(dz/dx) * sin(alpha) + (dz/dy) * cos(alpha)]*

where alpha represents the orientation along which we calculate the slope. This value is the angle that the direction makes with the top of the map and varies from 0° to 360° clockwise.

Various methods are available for the estimation of the gradients along the two axes. A detailed analysis can be found in Burrough & McDonnell (1998), in the "First and higher order derivatives of a continuous surface" section of the "Spatial analysis using continuous fields" chapter.

Two methods are available in the module: the Zevenbergen & Thorne (1987) and the Horn (1981) methods. The former considers only two cells for each gradient calculation, while the Horn (1981), suitable for rough surfaces, considers six cell values. 


### Zevenbergen & Thorne (1987) method

Gradients along the two axes are:

*dz/dx = [ Elev(i, j + 1) - Elev(i, j - 1) ] / 2 cell_size*

*dz/dy = [ Elev(i - 1, j) - Elev(i + 1, j) ] / 2 cell_size*

In the case of cells at the grid edge, this method can be replaced by the first difference method, where the central cell itself provides the missing data. 

![alt text](/ims/simple_1off_small.jpg "Zevenbergen & Thorne (1987) method")

*Fig. 3. Zevenbergen & Thorne (1987) method.*

### Horn (1981) method

The gradient estimation along an axis derives from the values of six cells in a 3x3 kernel. The weight applied to each cell depends on its position relative to the central cell. This method is the most suitable for rough surfaces.

Formulas of the gradients are (Burrough & McDonnell, 1998, eqs. 8.9 and 8.10):

*dz/dx = [ (Elev(i-1, j + 1) + 2 Elev(i, j + 1) + Elev(i+1, j + 1)) - (Elev(i-1, j - 1) + 2 Elev(i, j - 1) + Elev(i+1, j - 1)] / 8 cell_size*


*dz/dy = [ (Elev(i - 1, j-1) + 2 Elev(i - 1, j) + Elev(i - 1, j+1)) - (Elev(i + 1, j-1) + 2 Elev(i + 1, j) + Elev(i + 1, j+1)) ] / 8 cell_size*
           
   
![alt text](/ims/horn_small.jpg "Horn (1981) method")

*Fig. 4. Horn (1981) method.*


## Data processing

### Input

Required input data is the DEM, the source for the slope calculation. It must be projected (i.e., not in Long-Lat coordinates). The unit for vertical distances has to be the same as for the horizontal ones.

#### Uniform direction

For uniform direction analysis, the user can provide direction values in three ways:
 
 a) a single value (in degrees), e.g. 300;
 
 b) two values (start and stop), separated by a forward slash, e.g. 80/90: the slopes will be calculated for the range comprised between start and stop values, with a step of 1° (if start is less than stop) or -1° (if stop is less than start);
 
 c) three values, separated by a forward slash, e.g., 0/350/10: they correspond to start, stop and step values. A set of slopes will be automatically calculated with the given range and with the provided step. Each result file will have a suffix in its name describing the specific direction value.

#### Spatially-changing directions

For variable orientation analysis, a grid where each cell stores a specific direction (in decimal degrees) is required. Its geographical extent and cell size of the direction grid may be different from those of the DEM. 


### Output

The output consists of one or more grids, in ESRI ascii format, with the same geographical parameters (i.e., cell resolution, lower-left corner) as the input DEM, and slope values expressed in decimal degrees, positive for downwards orientations, and negative for upwards orientations (the opposite with respect to the previous versions). Maximum slope is obviously always positive.


The output files are saved in the folder defined by the user, with the name of each file will be made up by the basename and a suffix related to the specific calculation parameters (both used slope method and analysis type parameters):

#### Slope methods

  'ZT': Zevenbergen & Thorne (1987)
  
  'H': Horn (1981)
  
  
#### Analysis types

  'm': maximum slope
  
  'ud': uniform directions, followed by the specific angle value
  
  'vd': variable directions


### Examples

  'slope01_H_ud_180.0.asc': 'slope01' is the basename, 'H' means Horn (1981) method, 'ud' is the uniform directions option, with a value of 180.0 degrees.
  
  'result_ZT_vd.asc': 'result' is the basename, 'ZT' means Zevenbergen & Thorne (1987) method, 'vd' is the variable directions option. 
  
  
## Case study

Note: the sign convention for slope followed in this example refers to the previous versions of the module, inverted with respect to the current one.

We derive the directional slope of a glacial surface along the ice flow directions for the Reeves Glacier (Victoria Land, East Antarctica). The ice flow data derives from a PhD thesis by D. Biscaro. 

![alt text](/ims/reeves_fluxes.png "Glacial fluxes in the grounding zone of the Reeves Glacier (Antarctica)")

*Fig. 5. Glacial fluxes in the grounding zone of the Reeves Glacier (Antarctica). Data D. Biscaro.*

DEM data derives from the Radarsat Ramp DEM, that covers the whole of Antarctica with a nominal cell resolution of 200 m. Data were processed in Grass and then exported as ESRI Grid ASCII files. We use the Zevenbergen & Thorne (1987), since we are dealing with smooth surfaces. The results of the module were then imported in Grass and then exported into vtk format, to be visualized in Paraview. A profile of the directional slope is also presented. 

![alt text](/ims/mappa_slope_ink_small.jpg "Directional slopes in the grounding zone of the Reeves Glacier (Antarctica)")

*Fig. 6. Directional slopes in the grounding zone of the Reeves Glacier (Antarctica).*

Results are between 5.43° upward and 19.01° downward. There is a clear distinction between two main zones. The bulk of the slope variability occurs the ice grounded area (at the left). The slopes are about 0° in the floating ice zone (at the right). The separation between the two zones corresponds to a band of slopes reaching the steepest downward values, of about 19°, and that nearly corresponds to the grounding line, were the continental ice detaches from the bedrock and begins floating in the sea. We note also in the grounded area (bottom-right) some zones were the ice flows upwards: these situations are likely linked to subglacial obstacles that constrain the ice movements upwards. 


![alt text](/ims/profilo_small.jpg "Directional slope profile along A-B")

*Fig. 7. Directional slope profile along A-B (A: left; B: right).*






