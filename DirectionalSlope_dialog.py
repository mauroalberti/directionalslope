# -*- coding: utf-8 -*-


# Import PyQt libraries
from PyQt4.QtCore import *
from PyQt4.QtGui import *


# create the dialog 
class DirectionalSlopeDialog( QDialog ):


    def __init__( self ):

        super( DirectionalSlopeDialog, self ).__init__()
        
        self.initialize() 
                
        self.setupUI()
        
        # self.linkQGIS()
        
        # self.setConnections()   
        

    def initialize( self ):

        self.setWindowTitle( "DirectionalSlope v. 1.1.1" )
               
    
    def setupUI( self ):

        self.tabs_layout = QVBoxLayout()
        self.setLayout( self.tabs_layout )
        
        self.tabWidget = QTabWidget( )
        self.tabs_layout.addWidget( self.tabWidget )        
        
        self.setupCalculationTab()
        self.tabWidget.addTab( self.calc_tab, 'Directional slope' )
        
        self.setupHelpTab()
        self.tabWidget.addTab( self.help_tab, 'Help' )
        
        self.setupAboutTab()
        self.tabWidget.addTab( self.about_tab, 'About' )       


        self.tabWidget.setCurrentIndex(0)

        QObject.connect(self.btn_ok_cancel, SIGNAL( "accepted()" ), self.accept )
        QObject.connect(self.btn_ok_cancel, SIGNAL( "rejected()" ), self.reject )

        QMetaObject.connectSlotsByName(self)

        self.resize(565, 534)
        
        

    def setupCalculationTab(self):
        
        self.calc_tab = QWidget()        
        self.calc_layout = QVBoxLayout()
        self.calc_tab.setLayout( self.calc_layout )        
        
        self.setupInputDEM()
        self.setupMethods()        
        self.setupAnalysisTypes()
        self.setupResultFiles()     
        self.setupChoice()
         

            
        
    def setupInputDEM(self):

        # input DEM group box
        self.InputDEM_groupbox = QGroupBox( "Input DEM")
        self.InputDEM_layout = QGridLayout()        
        self.InputDEM_groupbox.setLayout( self.InputDEM_layout ) 

        self.InputDEM_layout.addWidget( QLabel(self.tr("DEM name")) , 0, 0, 1, 1  ) 
        self.InputDem_combobox = QComboBox()
        self.InputDEM_layout.addWidget( self.InputDem_combobox, 0, 1, 1, 3 )
        
        self.calc_layout.addWidget( self.InputDEM_groupbox )
        

    def setupMethods(self):
                
        # method group box
        self.Methods_groupbox = QGroupBox( "Methods")
        self.Methods_layout = QHBoxLayout()        
        self.Methods_groupbox.setLayout( self.Methods_layout ) 
        
        self.ZevenThornMeth_chbox = QCheckBox("Zevenbergen and Thorne (1987)")

        self.HornMeth_chbox = QCheckBox("Horn (1981)")

        self.Methods_layout.addWidget( self.ZevenThornMeth_chbox )
        self.Methods_layout.addWidget( self.HornMeth_chbox )

        self.calc_layout.addWidget( self.Methods_groupbox )
        

    def setupAnalysisTypes(self):
        
        # method group box
        self.AnalysisTypes_groupbox = QGroupBox( "Slope analysis types" )
        self.AnalysisTypes_layout = QGridLayout() 
        self.AnalysisTypes_groupbox.setLayout( self.AnalysisTypes_layout )

        self.UniformDirection_chbox = QCheckBox( "Sl. along uniform directions" )
        self.AnalysisTypes_layout.addWidget( self.UniformDirection_chbox , 0, 0)       
       
        self.UniformDirection_value = QLineEdit(self.AnalysisTypes_groupbox)
        self.UniformDirection_value.setPlaceholderText("e.g., 90; 170/180; 250/290/5")        
        self.AnalysisTypes_layout.addWidget( self.UniformDirection_value , 0, 1)


        self.VariableOrientations_chbox = QCheckBox( "Sl. along directions defined in grid " )
        self.AnalysisTypes_layout.addWidget( self.VariableOrientations_chbox , 1, 0)
        
        self.VariableOrientations_raster_cb = QComboBox()
        self.AnalysisTypes_layout.addWidget( self.VariableOrientations_raster_cb , 1, 1)        

        self.MaximumSlope_chbox = QCheckBox( "Maximum sl." )
        self.AnalysisTypes_layout.addWidget( self.MaximumSlope_chbox , 2, 0)        

        self.calc_layout.addWidget( self.AnalysisTypes_groupbox )
 

    def setsaveDirectory(self):
        
        dirName = QFileDialog.getExistingDirectory( self, 
                                                      self.tr("Select save directory"),
                                                      QDir.currentPath(), 
                                                      QFileDialog.ShowDirsOnly|
                                                      QFileDialog.ReadOnly )
        
        if dirName.isEmpty(): return
                
        self.ResultFolder_linedit.setText( dirName )
              
         
    def setupResultFiles(self):

        self.Results_groupbox = QGroupBox( "Output rasters" )
        self.Results_layout = QGridLayout() 
        self.Results_groupbox.setLayout( self.Results_layout )        
 
        self.Results_layout.addWidget( QLabel(self.tr("Save results in folder")) , 0, 0)
         
        self.ResultFolder_linedit = QLineEdit()
        self.Results_layout.addWidget( self.ResultFolder_linedit , 0, 1)
        
        self.ResultFolder_QPushButton = QPushButton(self.tr("Choose directory"))
        self.Results_layout.addWidget( self.ResultFolder_QPushButton , 0, 2)
        
        QObject.connect(self.ResultFolder_QPushButton, SIGNAL( "clicked()" ), self.setsaveDirectory )
 
        self.Results_layout.addWidget( QLabel(self.tr("Result file basename")), 1, 0)
                        
        self.ResultBasename_linedit = QLineEdit()
        self.ResultBasename_linedit.setPlaceholderText("e.g., slope01")
        self.Results_layout.addWidget( self.ResultBasename_linedit, 1, 1, 1, 2)
                
        self.calc_layout.addWidget( self.Results_groupbox )
 

    def setupChoice(self):        

        self.btn_ok_cancel = QDialogButtonBox(self.calc_tab)
        self.btn_ok_cancel.setOrientation(Qt.Horizontal)
        self.btn_ok_cancel.setStandardButtons(QDialogButtonBox.Cancel|QDialogButtonBox.Ok)
        
        self.calc_layout.addWidget( self.btn_ok_cancel )
                
        
 
    def setupHelpTab(self):

        help_html = """
        <html><head><title>Directional slope</title></head>
        <body>
        <h2>Directional slope</h2>
        <p>
                This module calculates the topographic slope from a DEM, both directional and steepest (maximum). 
                Directional slope can be calculated along orientations that are uniform or variable in space. 
                In the latter case, the user must provide them in an additional grid. 
                Data is assumed to be in a planar spatial frame (i.e., already projected, not as long-lat values), 
                with distance unit being the same for both horizontal and vertical coordinates (e.g., meters).
        </p>
        <h3>Methods</h3>
        <p>
                The methods used to estimate the slope value are two: Zevenbergen &amp; Thorne (1987), and Horn (1981). 
                The former is best suited for rough surfaces, while the latter performs better for smooth surfaces (Jones, 1997, cited in Burrough &amp; McDonnell, 1998).
        </p>
        <p>
        Three options for slope calculation are available: 
        <br />- maximum slope: steepest angle slope;
        <br />- slope along uniform directions: slope along directions that have a constant value in space;
        <br />- slope along variable directions: directions along which slopes will be calculated can vary in space.
        </p>
        <h3>Input</h3>
        <p>
                The axis frame for analysis is x - eastwards and y - northwards. 
                Direction angles are measured clockwise starting form the North, from 0&deg; to 360&deg;. 
                Therefore, 90&deg; is eastwards, 180&deg; southwards, 270&deg; westwards and 36&deg; again northwards.
        </p>
        <h4>Parameters for slope along uniform directions</h4>
        <p>
        The user can provide direction values in three ways:
        <br />a) entering a single value (in degrees), e.g. 300;
        <br />b) two values (start and stop), separated by a forward slash, e.g. 80/90: 
                the slopes will be calculated for the range comprised between start and stop values, with a step of 1&deg; (if start is less than stop) 
                or -1&deg; (if stop is less than start);
        <br />c) three values, separated by a forward slash, e.g., 0/350/10: they correspond to start, stop and step values. 
        A set of slopes will be automatically calculated with the given range and with the provided step. 
        Each result file will have a suffix in its name describing the specific direction value.
        </p>
        <h4>Parameters for slope along variable directions</h4>
        <p>
        The orientations along which to calculate slope are provided in an additional direction grid, storing values in decimal degrees. 
        </p>
        <h3>Output</h3>
        <p>
        Calculated slope values are in decimal degrees, positive for downwards directional slopes and negative for upwards ones. Maximum slope is obviously always positive.
        </p>
        <p>
        The output is stored as a set of ESRI ASCII grid files in the folder defined by the user, with the same geographical parameters (i.e., cell resolution, lower-left corner) as the input DEM. The name of each file will be made up by the basename and a suffix related to the specific calculation parameters (both used slope method and analysis type parameters):
        <br /><br />Slope methods
        <i>
        <br /> 'ZT': Zevenbergen &amp; Thorne (1987)
        <br /> 'H': Horn (1981)
        </i>
        <br /><br />Analysis types
        <i>
        <br /> 'm': maximum slope
        <br /> 'ud': uniform directions, plus the specific angle used
        <br /> 'vd': variable directions
        </i>
        </p>
        <h4>Examples</h4>
        <p>
        'slope01_H_ud_180.0.asc': 'slope01' is the basename, 'H' means Horn (1981) method, 'ud' is the uniform directions option, with a value of 180.0 degrees.
        </p>
        <p>
        'result_ZT_vd.asc': 'result' is the basename, 'ZT' means Zevenbergen & Thorne (1987) method, 'vd' is the variable directions option.
        </p>
        </body></hmtl>
        """

        self.help_tab = QWidget( )
        self.help_layout = QVBoxLayout()
        self.help_tab.setLayout( self.help_layout ) 

        self.help_textBrwsr = QTextBrowser(self.help_tab)

        self.help_textBrwsr.setHtml( help_html )


        self.help_layout.addWidget( self.help_textBrwsr )
        
        
 
    def setupAboutTab(self):       

        self.about_tab = QWidget( )
        self.about_layout = QVBoxLayout()
        self.about_tab.setLayout( self.about_layout ) 
        
        
        self.about_textBrwsr = QTextBrowser(self.about_tab)
        

        self.about_textBrwsr.setHtml("<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-weight:400; font-style:normal;\">\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-weight:600;\">Directional Slope </span></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span> Vers. 1.1.1 - February, 16, 2013.</span></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span> Created by Mauro Alberti - </span><a href=\"mailto:alberti.m65@gmail.com\"><span style=\" text-decoration: underline; color:#0000ff;\">alberti.m65@gmail.com</span></a></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span> License: GPL vers. 3.</span></p>\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span>See </span><a href=\"http://www.malg.eu/\"><span style=\" text-decoration: underline; color:#0000ff;\">www.malg.eu</span></a><span> for updates.</span></p></body></html>" )



        self.about_textBrwsr.setOpenExternalLinks(True)
        
        self.about_layout.addWidget( self.about_textBrwsr )
 



 


