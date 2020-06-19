"""
Abigail Stone
6/18/2020
vector2raster.py

Converts a .shp file to a .tif raster using the extent/resolution of a reference .tif
"""

import sys
import argparse
try:
    from osgeo import ogr, osr, gdal
except ImportError:
    sys.exit('ERROR: cannot find GDAL/OGR modules')



def vector2raster(inputshp, outtif, reftif, ops):
    datatype = gdal.GDT_Float32
    burnVal = 0

    # Read files
    ref = gdal.Open(reftif, gdal.GA_ReadOnly)
    shapefile = ogr.Open(inputshp)
    shp_layer = shapefile.GetLayer()

    # create output layer
    out = gdal.GetDriverByName('GTiff').Create(
    outtif, ref.RasterXSize, ref.RasterYSize, 1, datatype, options=['COMPRESS=DEFLATE']
    )
    out.SetProjection(ref.GetProjectionRef())
    out.SetGeoTransform(ref.GetGeoTransform())

    # Rasterize!
    band = out.GetRasterBand(1)
    band.SetNoDataValue(1)
    gdal.RasterizeLayer(out, [1], shp_layer, burn_values=[burnVal], options=ops)

    # close datasets
    band = None
    out = None
    ref = None
    shapefile = None



if __name__ == "__main__":
    DESCRIPTION = "Converts a .shp file to a .tif raster file"

    parser = argparse.ArgumentParser(description=DESCRIPTION)

    parser.add_argument("shp_in",
                        help="Path to shapefile containing the layer to be \
                        rasterized (required)")

    parser.add_argument("out_path",
                        help="Destination path for output .tif file (required)")

    parser.add_argument("ref",
                        help="Reference .tif file of footprint (required)")

    parser.add_argument("-a", "--attribute",
                        help="attribute from .shp file to rasterize")

    args = parser.parse_args()

    ops = []
    if args.attribute:
        ops.append('ATTRIBUTE='+args.attribute)

    vector2raster(args.shp_in, args.out_path, args.ref, ops)
