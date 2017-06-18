

import os
from math import *
import numpy as np
import gdal
from gdalconst import *


# 2D Array coordinates 
class ArrCoord:

    def __init__(self, ival, jval):
        self.i = ival
        self.j = jval

    def get_i(self):
        """get i value"""

        return self.i

    def get_j(self):
        """get j value"""

        return self.j


class Point:
    """3D point class"""

    def __init__(self, x, y, z=0.0):
        self.x = x
        self.y = y
        self.z = z

    def distance(self, pt):
        """calculate euclidean distance between two points"""

        return sqrt((self.x - pt.x)**2 + (self.y - pt.y)**2 + (self.z - pt.z)**2)

    def movedby(self, sx, sy):
        """create a new point shifted by given amount"""

        return Point( self.x + sx , self.y + sy)


class GDALParameters:
    """GDAL raster parameters"""

    def __init__(self): 
        self.driver_shortname = None
        self.datatype = None
        self.nodatavalue = None
        self.projection = None
        self.topleftX = None
        self.topleftY = None
        self.pixsizeEW = None
        self.pixsizeNS = None
        self.rows = None
        self.cols = None
        self.rotationA = None
        self.rotationB = None

    def set_driverShortName(self, raster_driver):
        """set driver format"""

        self.driver_shortname = raster_driver

    def get_driverShortName(self):
        """get driver format"""

        return self.driver_shortname

    def set_dataType(self, data_type):
        """set data type"""

        self.datatype = data_type

    def get_dataType(self):
        """get data type"""

        return self.datatype

    def set_noDataValue(self, nodataval):
        """set nodata value"""

        self.nodatavalue = nodataval

    def get_noDataValue(self):
        """get nodata value"""

        return self.nodatavalue

    def set_projection(self, proj):
        """set projection"""

        self.projection = proj

    def get_projection(self):
        """get projection"""

        return self.projection

    def set_topLeftX(self, topleftX):
        """set topleftX value"""

        self.topleftX = topleftX

    def get_topLeftX(self):
        """get topleftX value"""

        return self.topleftX

    def set_topLeftY(self, topleftY):
        """set topleftY value"""

        self.topleftY = topleftY

    def get_topLeftY(self):
        """get topleftY value"""

        return self.topleftY

    def set_pixSizeEW(self, pixsizeEW):
        """set pixsizeEW value"""

        self.pixsizeEW = pixsizeEW

    def get_pixSizeEW(self):
        """get pixsizeEW value"""

        return self.pixsizeEW

    def set_pixSizeNS(self, pixsizeNS):
        """set pixsizeNS value"""

        self.pixsizeNS = pixsizeNS

    def get_pixSizeNS(self):
        """get pixsizeNS value"""

        return self.pixsizeNS

    def set_rows(self, rows):
        """set rows number"""

        self.rows = rows

    def get_rows(self):
        """get rows number"""
        return self.rows

    def set_cols(self, cols):
        """set cols number"""

        self.cols = cols

    def get_cols(self):
        """get cols number"""

        return self.cols

    def set_rotationA(self, rotationA):
        """set rotationA value"""

        self.rotationA = rotationA

    def get_rotationA(self):
        """get rotationA value"""

        return self.rotationA

    def set_rotationB(self, rotationB):
        """set rotationB value"""

        self.rotationB = rotationB

    def get_rotationB(self):
        """get rotationB value"""

        return self.rotationB

    def check_params(self):
        """check di absence of axis rotations or pixel size differences"""
        
        # set tolerance value

        diff_tolerance = 1e-06 
        
        # check if pixel size can be considered the same in the two axis directions

        if abs(abs(self.pixsizeEW) - abs(self.pixsizeNS))/abs(self.pixsizeNS) > diff_tolerance :
            raise RasterParametersErrors('Pixel sizes in x and y directions are different in raster')
            
        # check for the absence of axis rotations

        if abs(self.rotationA) > diff_tolerance or abs(self.rotationB) > diff_tolerance:
            raise RasterParametersErrors('There should be no axis rotation in raster')
        
        return


class RasterParametersErrors(Exception):
    """exception for raster parameters"""

    pass  


class OutputErrors(Exception):
    """exception for output errors"""

    pass


class SlopeErrors(Exception):
    """exception for slope calculation errors"""

    pass
    

class SpatialDomain:
    """Rectangular spatial domain class"""

    def __init__(self, pt_init, pt_end):

        self.pt_init = pt_init # lower-left corner of the domain, class: Point
        self.pt_end = pt_end # top-right corner of the domain, class: Point

    def get_start_point(self):
        """get start point"""

        return self.pt_init

    def get_end_point(self):
        """get end point"""

        return self.pt_end

    def get_xrange(self):
        """get x range of spatial domain"""

        return self.pt_end.x-self.pt_init.x

    def get_yrange(self):
        """get y range of spatial domain"""

        return self.pt_end.y-self.pt_init.y

    def get_zrange(self):
        """get z range of spatial domain"""

        return self.pt_end.z-self.pt_init.z

    def get_horiz_area(self):
        """get horizontal area of spatial domain"""

        return self.get_xrange()*self.get_yrange()


class Grid:

    def __init__(self):
        
        self.grid_domain = None	
        self.grid_data = None			
        self.nodatavalue = None

    def set_grid_domain(self, pt_init, pt_end):
        """set grid domain"""
        
        self.grid_domain = SpatialDomain(pt_init, pt_end)

    def get_grid_domain(self):
        """get grid domain"""
        
        return self.grid_domain

    def set_grid_data(self, data_array):
        """set grid data"""
        
        self.grid_data = data_array

    def get_grid_data(self):
        """get grid data"""
        
        return self.grid_data

    def set_nodatavalue(self, nodata):
        """set nodata value"""
        
        self.nodatavalue = nodata

    def get_nodatavalue(self):
        """get nodata value"""
        
        return self.nodatavalue

    def get_ylines_num(self):
        """get row number of the grid domain"""
        
        return np.shape(self.grid_data)[0]		

    def get_xlines_num(self):
        """column number of the grid domain"""
        
        return np.shape(self.grid_data)[1]		

    def get_zlines_num(self):
        """z line number of the grid domain"""
        
        return np.shape(self.grid_data)[2]	

    def get_cellsize_x(self):
        """returns the cell size of the gridded dataset in the x direction"""
        
        return self.grid_domain.get_xrange()/float(self.get_xlines_num())

    def get_cellsize_y(self):
        """returns the cell size of the gridded dataset in the y direction"""
        
        return self.grid_domain.get_yrange()/float(self.get_ylines_num())

    def get_cellsize_z(self):
        """returns the cell size of the gridded dataset in the z direction"""
        
        return self.grid_domain.get_zrange()/float(self.get_zlines_num())

    def get_cellsize_horiz_mean(self):
        """returns the mean horizontal cell size"""
        
        return (self.get_cellsize_x()+self.get_cellsize_y())/2.0

    def geog2gridcoord(self, curr_Pt):
        """converts from geographic to grid coordinates"""
        
        currArrCoord_grid_i = (self.get_grid_domain().get_end_point().y - curr_Pt.y)/self.get_cellsize_y()
        currArrCoord_grid_j = (curr_Pt.x - self.get_grid_domain().get_start_point().x)/self.get_cellsize_x()
        
        return ArrCoord(currArrCoord_grid_i, currArrCoord_grid_j)

    def grid2geogcoord(self, currArrCoord):
        """converts from grid to geographic coordinates"""
        
        currPt_geogr_y = self.get_grid_domain().get_end_point().y - currArrCoord.i*self.get_cellsize_y()
        currPt_geogr_x = self.get_grid_domain().get_start_point().x + currArrCoord.j*self.get_cellsize_x()
        
        return Point(currPt_geogr_x, currPt_geogr_y)

    def write_esrigrid(self, outgrid_fn):
        """writes ESRI ascii grid"""

        outgrid_fn = str(outgrid_fn)
        
        # null value in output slope grid (ESRI format)

        esri_nullvalue = -99999 

        # checking existance of output slope grid

        if os.path.exists(outgrid_fn):
          os.remove(outgrid_fn)

        try:
            outputgrid = open(outgrid_fn, 'w')  # create the output ascii file
        except:
            raise OutputErrors('Unable to open output grid: ' + outgrid_fn)
       
        if outputgrid is None:
            raise OutputErrors('Unable to create output grid: ' + outgrid_fn)

        # writes header of grid ascii file

        outputgrid.write('NCOLS %d\n' % self.get_xlines_num())
        outputgrid.write('NROWS %d\n' % self.get_ylines_num())
        outputgrid.write('XLLCORNER %.8f\n' % self.grid_domain.get_start_point().x)
        outputgrid.write('YLLCORNER %.8f\n' % self.grid_domain.get_start_point().y)
        outputgrid.write('CELLSIZE %.8f\n' % self.get_cellsize_horiz_mean())
        outputgrid.write('NODATA_VALUE %d\n' % esri_nullvalue)

        esrigrid_outvalues = np.where(np.isnan(self.grid_data), esri_nullvalue, self.grid_data)
        
        # output of results

        for i in range(0, self.get_ylines_num()):
                for j in range(0, self.get_xlines_num()):
                        outputgrid.write('%.8f ' % (esrigrid_outvalues[i,j]))
                outputgrid.write('\n')
                
        outputgrid.close()

        return True


def horn_gradient(z):
    """calculate x and y gradients according to Horn (1981) method"""

    nrows, ncols = z.shape

    z_00 = z[0:nrows-2, 0:ncols-2]
    z_10 = z[1:nrows-1, 0:ncols-2]
    z_20 = z[2:nrows, 0:ncols-2]

    z_02 = z[0:nrows-2, 2:ncols]
    z_12 = z[1:nrows-1, 2:ncols]
    z_22 = z[2:nrows, 2:ncols]
    
    dz_dx = (z_02+2.0*z_12+z_22-z_00-2.0*z_10-z_20)/8.0
    
    z_00 = z[0:nrows-2, 0:ncols-2]
    z_01 = z[0:nrows-2, 1:ncols-1]
    z_02 = z[0:nrows-2, 2:ncols]

    z_20 = z[2:nrows, 0:ncols-2]
    z_21 = z[2:nrows, 1:ncols-1]
    z_22 = z[2:nrows, 2:ncols]
    
    dz_dy = (z_20+2.0*z_21+z_22-z_00-2.0*z_01-z_02)/8.0	
    
    return (dz_dy, dz_dx)


def maxslope(dz_dy, dz_dx):
    """calculate maximum slope values (array, values in degrees)"""

    p = dz_dy*dz_dy + dz_dx*dz_dx

    return np.arctan(np.sqrt(p))*180.0/pi 


def directionalslope(dz_dx, dz_dy, direction):
    """calculate directional slope values (array, values in degrees)"""
    
    direction_rad = direction*pi / 180.0
    tan_alpha = ((dz_dx*sin(direction_rad)) + (dz_dy*cos(direction_rad)))

    return np.arctan(tan_alpha)*(-180.0)/pi  # output in degrees - downward slopes: positive, upward slopes, negative


def dirslope_from_directions(dz_dx, dz_dy, direction_array):
    """calculate directional slope values from direction array (array, values in degrees)"""
    
    directions_rad = direction_array*pi / 180.0
    tan_alpha = ((dz_dx*np.sin(directions_rad)) + (dz_dy*np.cos(directions_rad)))

    return np.arctan(tan_alpha)*(-180.0)/pi  # output in degrees - downward slopes: positive, upward slopes, negative


def get_unifdir_params(text):
    """reads uniform direction parameters and converts them into a list of values"""

    # error - no input

    if text == '':
        raise ValueError('No input for uniform direction orientations')
        
    # initialize list of list of orientations

    values = []
    
    # tokenize by ';'

    tokens = text.split(';')
     
    # populate list of list

    for token in tokens:
        
        # tokenize by '/'

        curr_range = token.split('/') 
        
        # convert to float and populate sublist

        value_range = []
        for value in curr_range:
            try:
                value_range.append(float(value))                
            except ValueError:
                raise ValueError('Format of input values is not correct')
        
        # populate list

        values.append(value_range)
 
    # initialize result list

    result = []
    
    # populates list of direction values

    for value_range in values:
        
        if len(value_range) == 1:
            value_range.append(value_range[0])
            value_range.append(1.0)            
        elif len(value_range) <= 3:
            if len(value_range) == 2:
                if value_range[0] > value_range[1]:
                    value_range.append(-1.0)
                else:
                    value_range.append(1.0)
                    
            if value_range[0] > value_range[1] and value_range[2] >= 0.0:
                raise ValueError('Increment for uniform direction analysis cannot be positive or zero when start value is larger than stop value')

            if value_range[0] < value_range[1] and value_range[2] <= 0.0:
                raise ValueError('Increment for uniform direction analysis cannot be negative or zero when start value is smaller than stop value')
        else:
            raise ValueError('Too many numbers in input')

        values_array = np.arange(value_range[0], value_range[1]+value_range[2]/2, value_range[2])
        
        for value in values_array:
            result.append(value)
                
    return result


def read_raster_band(source_file):
    """read input raster band based on GDAL"""

    # GDAL register

    gdal.AllRegister
    
    # open raster file and check operation success

    raster_data = gdal.Open(str(source_file), GA_ReadOnly)    
    if raster_data is None:
        raise IOError('Unable to open raster') 

    # initialize DEM parameters

    raster_params = GDALParameters()
    
    # get driver type for current raster

    raster_params.set_driverShortName(raster_data.GetDriver().ShortName)

    # get current raster projection

    raster_params.set_projection(raster_data.GetProjection())   

    # get row and column numbers

    raster_params.set_rows(raster_data.RasterYSize)
    raster_params.set_cols(raster_data.RasterXSize)
    
    # get and check number of raster bands - it must be one

    raster_bands = raster_data.RasterCount
    if raster_bands > 1:
        raise TypeError('More than one raster band in raster')
    
    # set critical grid values from geotransform array

    raster_params.set_topLeftX(raster_data.GetGeoTransform()[0])
    raster_params.set_pixSizeEW(raster_data.GetGeoTransform()[1])
    raster_params.set_rotationA(raster_data.GetGeoTransform()[2])
    raster_params.set_topLeftY(raster_data.GetGeoTransform()[3])
    raster_params.set_rotationB(raster_data.GetGeoTransform()[4])
    raster_params.set_pixSizeNS(raster_data.GetGeoTransform()[5])
 
    # get single band

    band = raster_data.GetRasterBand(1)
    
    # get no data value for current band

    raster_params.set_noDataValue(band.GetNoDataValue())

    # read data from band

    grid_values = band.ReadAsArray(0,0,raster_params.get_cols(),raster_params.get_rows())
    if grid_values is None:
        raise IOError('Unable to read data from raster')
     
    # transform data into numpy array

    data = np.asarray(grid_values) 

    # if nodatavalue exists, set null values to NaN in numpy array

    if raster_params.get_noDataValue() is not None:
        data = np.where(abs(data-raster_params.get_noDataValue())> 1e-05, data, np.NaN) 

    return raster_params, data


def calculate_paired_values(dem_grid, vardir_grid):
    """find nearest values from vardir_grid for each cell center of dem_grid (return Grid)"""

    # inizialite grid

    paired_vardir_grid = Grid()
    
    # set grid values

    paired_vardir_grid.set_nodatavalue( vardir_grid.get_nodatavalue())
    paired_vardir_grid.set_grid_domain( dem_grid.get_grid_domain().get_start_point().movedby(dem_grid.get_cellsize_x(), dem_grid.get_cellsize_y()), \
                                        dem_grid.get_grid_domain().get_end_point().movedby(dem_grid.get_cellsize_x()*(-1), dem_grid.get_cellsize_y()*(-1)))
 
    # create array that will store the vardir values

    paired_values = np.zeros((dem_grid.get_ylines_num()-2, dem_grid.get_xlines_num()-2))*np.NaN
    
    for i in xrange(paired_values.shape[0]):
        for j in xrange(paired_values.shape[1]): 
            
            # offset between grid dem cell corner and paired_values cell center, in grid coordinates
            offset_dem_array = 1.5
            
            # geographic coordinates of dem grid cell center
            currGeogCoordPt = dem_grid.grid2geogcoord(ArrCoord(i+offset_dem_array,j+offset_dem_array))
            
            # grid coordinates in vardir grid
            currVarDirArrCoor = vardir_grid.geog2gridcoord(currGeogCoordPt)
            
            if currVarDirArrCoor.get_i() >= 0.0 and currVarDirArrCoor.get_i() < vardir_grid.get_ylines_num() and \
              currVarDirArrCoor.get_j() >= 0.0 and currVarDirArrCoor.get_j() < vardir_grid.get_xlines_num():
                paired_i = int(currVarDirArrCoor.get_i())
                paired_j = int(currVarDirArrCoor.get_j())
                paired_values[i,j] = vardir_grid.get_grid_data()[paired_i, paired_j]
                 
    paired_vardir_grid.set_grid_data( paired_values)

    return paired_vardir_grid


def create_grid(grid_params, grid_data):
    
    # inizialite grid

    grid = Grid()
    
    # define lower-left corner as initial point

    pt_init = Point( grid_params.get_topLeftX(), grid_params.get_topLeftY() - abs(grid_params.get_pixSizeNS())*grid_params.get_rows()) 
    
    # define top-right corner as end point

    pt_end = Point( grid_params.get_topLeftX() + abs(grid_params.get_pixSizeEW())*grid_params.get_cols(), grid_params.get_topLeftY())  
    
    # set grid values

    grid.set_nodatavalue( grid_params.nodatavalue)
    grid.set_grid_domain( pt_init, pt_end)
    grid.set_grid_data( grid_data)

    return grid	                    


def methods_calculation(gradients, directional_methods, analysis_parameters, output_parameters, method_string):
    """calculates results based on the chosen methods and parameters"""

    # extracts parameters

    dz_dy, dz_dx = gradients
    maxslope_analysis, unifdir_analysis, vardir_analysis = directional_methods
    dem_grid, unif_direction_values, vardir_grid = analysis_parameters
    output_grid_params, output_path, output_suffix = output_parameters
 
    # maximum slope case

    if maxslope_analysis:                    
        
        outfilepath = output_path + method_string + "_m" + output_suffix 
        slope_array = maxslope(dz_dy, dz_dx)                                   
        if slope_array is not None:
            try:
                create_grid(output_grid_params, slope_array).write_esrigrid(outfilepath) 
            except (OutputErrors):
                raise OutputErrors('Unable to create output file: ' + outfilepath)
        else:
            raise SlopeErrors('Unable to calculate maximum slope')

    # uniform direction case

    if unifdir_analysis:
        
        for direction in unif_direction_values:
            
            outfilepath = output_path + method_string + "_ud_" + str(direction) + output_suffix  
            slope_array = directionalslope(dz_dx, dz_dy, direction)
            if slope_array is not None: 
                try:
                    create_grid(output_grid_params, slope_array).write_esrigrid(outfilepath) 
                except OutputErrors:
                    raise OutputErrors('Unable to create output file: ' + outfilepath)
            else:
                raise SlopeErrors('Unable to calculate uniform direction slope')
 
    # variable direction case

    if vardir_analysis:             

        outfilepath = output_path + method_string + "_vd" + output_suffix                
        slope_array = dirslope_from_directions(dz_dx, dz_dy, vardir_grid.grid_data)                                   
        if slope_array is not None: 
            try:
                create_grid(output_grid_params, slope_array).write_esrigrid(outfilepath) 
            except OutputErrors:
                raise OutputErrors('Unable to create output file: ' + outfilepath)
        else:
            raise SlopeErrors('Unable to calculate variable direction slope')
        
    return
     
