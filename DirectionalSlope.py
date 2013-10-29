"""
/***************************************************************************
 DirectionalSlope
                                 A QGIS plugin
 Calculates directional slope
                              -------------------
        begin                : 2011-10-25
        version              : 2013-10-24 (1.2.0), QGIS 2.0 compatible
        copyright            : (C) 2011-2013 by Mauro Alberti - www.malg.eu
        email                : alberti.m65@gmail.com
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 3 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
"""


# Import the PyQt and QGIS libraries
from PyQt4.QtCore import *
from PyQt4.QtGui import *
from qgis.core import *

# Initialize Qt resources from file resources.py
import resources

# Import the code for the dialog
from DirectionalSlope_dialog import DirectionalSlopeDialog

# Import of module analysis classes
from DirectionalSlope_classes import *


class DirectionalSlope:

    
    def __init__(self, iface):
        # constructor
        
        # Save reference to the QGIS interface
        self.iface = iface

    
    def initGui(self):
        # GUI inizialization
        
        # Create action that will start plugin configuration
        self.action = QAction(QIcon(":/plugins/DirectionalSlope/icon.png"), \
            "DirectionalSlope", self.iface.mainWindow())
        
        # connect the action to the run method
        self.action.triggered.connect( self.run )

        # Add toolbar button and menu item
        self.iface.addToolBarIcon(self.action)
        self.iface.addPluginToMenu("DirectionalSlope", self.action)
    
    
    def unload(self):
        # Remove the plugin menu item and icon
        
        self.iface.removePluginMenu("DirectionalSlope",self.action)
        self.iface.removeToolBarIcon(self.action)


    
    def run(self):
        # run method 

        # create the dialog
        dlg = DirectionalSlopeDialog()

        # append loaded raster layers to combo boxes	
        self.layermap = QgsMapLayerRegistry.instance().mapLayers() 

        for (name,layer) in self.layermap.iteritems():
            if layer.type() == QgsMapLayer.RasterLayer: 
                dlg.InputDem_combobox.addItem(layer.name())		
                dlg.VariableOrientations_raster_cb.addItem(layer.name())               
  
        # show the dialog
        dlg.show()    

        # user choice
        result = dlg.exec_() 
 
        # process data
        if result == 1:

            # get analysis parameters from GUI
            
            dem_name = dlg.InputDem_combobox.currentText()

            zt_method = dlg.ZevenThornMeth_chbox.isChecked()
            h_method = dlg.HornMeth_chbox.isChecked()              
            
            maxslope_analysis = dlg.MaximumSlope_chbox.isChecked()
            unifdir_analysis = dlg.UniformDirection_chbox.isChecked()
            vardir_analysis = dlg.VariableOrientations_chbox.isChecked()

            unifdir_inputtext = dlg.UniformDirection_value.text()            
            vardirectgrid_name = dlg.VariableOrientations_raster_cb.currentText()
            
            output_folder = dlg.ResultFolder_linedit.text()
            output_basename = dlg.ResultBasename_linedit.text()
            
            
            # verify input parameters
            
            if dem_name is None or dem_name == '':
                QMessageBox.critical(self.iface.mainWindow(), "DEM input", "No DEM defined")   
                return  

            if (not (zt_method or h_method)):
                QMessageBox.critical(self.iface.mainWindow(), "Slope methods", "No method defined")   
                return 
                
            if (not (maxslope_analysis or unifdir_analysis or vardir_analysis)):
                QMessageBox.critical(self.iface.mainWindow(), "Slope analysis types", "No option chosen")   
                return  

            if unifdir_analysis and unifdir_inputtext == '':
                QMessageBox.critical(self.iface.mainWindow(), "Slope analysis types", "No input for uniform direction analysis")   
                return 
                
            if output_folder is None or output_folder == '':
                QMessageBox.critical(self.iface.mainWindow(), "Output folder", "No output folder defined")   
                return                 
        
            if output_basename is None or output_basename == '':
                QMessageBox.critical(self.iface.mainWindow(), "Output basename", "No output name defined")   
                return                
                
                
            # define DEM input file 
            
            for (name,layer) in self.layermap.iteritems():
                if layer.name() == dem_name:       
                    dem_layer = layer
                    break                    
            dem_source = dem_layer.source()
            
            
            # read DEM parameters and array data
            
            try:
                input_grid_params, dem_array = read_raster_band(dem_source)
            except (IOError, TypeError), e:            
                QMessageBox.critical(self.iface.mainWindow(), "DEM input error", str(e))   
                return


            # check geometric parameters of DEM
            
            try:
                input_grid_params.check_params()
            except (Raster_Parameters_Errors), e:
                QMessageBox.critical(self.iface.mainWindow(), "DEM parameters error", str(e))   
                return                
        
        
            # create DEM grid
            
            dem_grid = create_grid(input_grid_params, dem_array)        


            # preprocessing for uniform direction analysis
            
            unif_direction_values = None
            if unifdir_analysis:
                head_unifdirinput_error = 'Error in uniform direction input'
                try:
                    unif_direction_values = get_unifdir_params(unifdir_inputtext)
                except ValueError, e:
                    QMessageBox.critical(self.iface.mainWindow(), head_unifdirinput_error, str(e))   
                    return                    
                if len(unif_direction_values) == 0:
                    QMessageBox.critical(self.iface.mainWindow(), head_unifdirinput_error, "Unable to read uniform direction data")   
                    return   
 
 
            # preprocessing for variable direction analysis
            
            paired_vardir_grid = None
            if vardir_analysis:
                # read orientation data
                if vardirectgrid_name is None:
                    QMessageBox.critical(self.iface.mainWindow(), "Variable directional analysis", "No raster provided")   
                    return                   
                for (name,layer) in self.layermap.iteritems():
                    if layer.name() == vardirectgrid_name:       
                        vardirectgrid_layer = layer
                        break
                vardirectgrid_source = vardirectgrid_layer.source()
                
                try:
                    vardir_params, vardir_array = read_raster_band(vardirectgrid_source)
                except (IOError, TypeError), e:            
                    QMessageBox.critical(self.iface.mainWindow(), "Direction grid input", str(e))   
                    return 
 
                try:
                    vardir_params.check_params()
                except (Raster_Parameters_Errors), e:
                    QMessageBox.critical(self.iface.mainWindow(), "Direction grid parameters", str(e))   
                    return
           
                vardir_grid = create_grid(vardir_params, vardir_array) 
                 
                try:
                    paired_vardir_grid = calculate_paired_values(dem_grid, vardir_grid)
                except:
                    QMessageBox.critical(
					self.iface.mainWindow(), 
					"Variable direction grid", 
					"Unable to calculate paired_vardir_grid. Please report error to alberti.m65@gmail.com" 
					)
   
                    return 

                
            
            ## START OF SLOPE CALCULATIONS
            
            
            # output grid parameters
            
            output_grid_params = GDALParameters()
            output_grid_params.set_driverShortName(input_grid_params.get_driverShortName())
            output_grid_params.set_dataType(input_grid_params.get_dataType())
            output_grid_params.set_noDataValue(input_grid_params.get_noDataValue())
            output_grid_params.set_projection(input_grid_params.get_projection())
            output_grid_params.set_topLeftX(input_grid_params.get_topLeftX() + abs(input_grid_params.get_pixSizeEW()))
            output_grid_params.set_topLeftY(input_grid_params.get_topLeftY() - abs(input_grid_params.get_pixSizeNS()))
            output_grid_params.set_pixSizeEW(input_grid_params.get_pixSizeEW())
            output_grid_params.set_pixSizeNS(input_grid_params.get_pixSizeNS())
            output_grid_params.set_rows(input_grid_params.get_rows() - 2)
            output_grid_params.set_cols(input_grid_params.get_cols() - 2)
            output_grid_params.set_rotationA(input_grid_params.get_rotationA())
            output_grid_params.set_rotationB(input_grid_params.get_rotationB())
           
           
            # output file path and suffix
            
            output_path = os.path.join(str(output_folder), str(output_basename))
            output_suffix = ".asc"


            # general parameters
            
            directional_methods = ( maxslope_analysis, unifdir_analysis, vardir_analysis )
            analysis_values = ( dem_grid, unif_direction_values, paired_vardir_grid )
            output_parameters = ( output_grid_params, output_path, output_suffix )
   
            
            # Zevenbergen and Thorne (1987) method analysis 
            
            if zt_method:
                
                method_string = "_ZT"
                
                dz_dy, dz_dx = np.gradient(dem_grid.grid_data)               
                dz_dx = dz_dx[1:-1,1:-1]/dem_grid.get_cellsize_horiz_mean()
                dz_dy = (dz_dy[1:-1,1:-1]/dem_grid.get_cellsize_horiz_mean())*(-1)

                try:
                    methods_calculation((dz_dy, dz_dx), directional_methods, analysis_values, output_parameters, method_string)
                except (Slope_Errors, Output_Errors), e:
                    QMessageBox.critical(self.iface.mainWindow(), "Calculation error", str(e))
                    return


            # Horn (1981) method analysis 
            if h_method:
                
                method_string = "_H"
                
                dz_dy, dz_dx = horn_gradient(dem_grid.grid_data)               
                dz_dx = dz_dx/dem_grid.get_cellsize_horiz_mean()
                dz_dy = (dz_dy/dem_grid.get_cellsize_horiz_mean())*(-1)
                
                try:
                    methods_calculation((dz_dy, dz_dx), directional_methods, analysis_values, output_parameters, method_string)
                except (Slope_Errors, Output_Errors), e:
                    QMessageBox.critical(self.iface.mainWindow(), "Calculation error", str(e))
                    return            

            # all done
            QMessageBox.information(
				self.iface.mainWindow(), 
				"Module output", 
				"Processings completed.\nResults saved in folder: "+str(output_folder)
				)






