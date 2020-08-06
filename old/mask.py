# myles stokowski
# 07.03.2020
# mask.py
# quick script for masking the OID rasters by the landcover road mask

import glob
import numpy as np
import pandas as pd
try: 
    from osgeo import gdal
    from osgeo.gdalnumeric import * 
    from osgeo.gdalconst import *
except ImportError:
    sys.exit('ERROR: cannot find GDAL/OGR modules')

def mask(road_raster, landcover_mask, dsttif):
    road = gdal.Open(road_raster, gdal.GA_ReadOnly)
    road_array = BandReadAsArray(road.GetRasterBand(1))

    landcover = gdal.Open(landcover_mask, gdal.GA_ReadOnly)
    landcover_array = BandReadAsArray(landcover.GetRasterBand(1))

    masked = -9999 * (landcover_array == 0) + road_array * (landcover_array == 1)

    # create new output raster
    datatype = gdal.GDT_Float32
    drv = gdal.GetDriverByName('GTiff')
    
    dst = drv.Create(dsttif, road.RasterXSize, road.RasterYSize, 1, datatype, options=['COMPRESS=DEFLATE'])
    dst.SetProjection(road.GetProjectionRef())
    dst.SetGeoTransform(road.GetGeoTransform())

    band = dst.GetRasterBand(1)
    band.SetNoDataValue(-9999)
    band.WriteArray(masked)

    band = None
    dst = None
    src = None
    ref = None