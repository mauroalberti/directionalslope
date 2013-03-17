
from DirectionalSlope import DirectionalSlope

def name():
    return "Directional Slope"

def description():
    return "Calculates directional slope along uniform or variable orientations"

def version():
    return "1.1.1"

def icon():
    return "icon.png"

def qgisMinimumVersion():
    return "1.8"

def classFactory(iface):
    return DirectionalSlope(iface)
