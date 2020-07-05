"""
Myles Stokowski
6.19.2020

outputs new .tif with data from input .tif and projection/alignment from
reference .tif

conceived use case is converting 1m VA landcover layer to 3m
resolution and same extent/alignment as SAR data
"""
import sys
import argparse
try:
    from osgeo import gdal
except ImportError:
    sys.exit('ERROR: cannot find GDAL/OGR modules')

def convert_resolution(srctif, dsttif, reftif):
    """
    output dsttif raster with data from srctif, resolution/extent from reftif
    """
    datatype = gdal.GDT_Byte
    drv = gdal.GetDriverByName('GTiff')

    # read files
    ref = gdal.Open(reftif, gdal.GA_ReadOnly)
    src = gdal.Open(srctif, gdal.GA_ReadOnly)

    # create output layer
    dst = drv.Create(dsttif, ref.RasterXSize, ref.RasterYSize, 1, datatype, \
        options=['COMPRESS=DEFLATE'])
    dst.SetProjection(ref.GetProjectionRef())
    dst.SetGeoTransform(ref.GetGeoTransform())

    src_proj = src.GetProjection()
    ref_proj = ref.GetProjection()

    # fill with data from source
    band = dst.GetRasterBand(1)
    band.SetNoDataValue(-9999)
    gdal.ReprojectImage(src, dst, src_proj, ref_proj, gdal.GRA_NearestNeighbour)

    # close datasets
    band = None
    dst = None
    src = None
    ref = None

if __name__ == "__main__":
    DESCRIPTION = "Reprojects data from a .tif to the dimensions of a reference .tif"

    parser = argparse.ArgumentParser(description=DESCRIPTION)

    parser.add_argument("srctif",
                        help="path to source .tif containing data (required)")

    parser.add_argument("dsttif",
                        help="destination path to output .tif (required)")

    parser.add_argument("reftif",
                        help="path to reference .tif containing target \
                                resolution (required)")

    args = parser.parse_args()

    convert_resolution(args.srctif, args.dsttif, args.reftif)
