from __future__ import absolute_import

from .DirectionalSlope import DirectionalSlope


def classFactory(iface):
    return DirectionalSlope(iface)
