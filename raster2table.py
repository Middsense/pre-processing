"""
Myles Stokowski
6.22.20

Functionality: Takes in a 'mask' and 'data' raster and outputs a  DataFrame
with the values from both rasters for each pixel where the neither raster is
'nodata'. Then, provides options for how to summarize the 'data' raster
for each value/group in the 'mask' raster.

Middsense use case: calculate average SAR amplitude for the ~100 pixels under
each segment of road, uniquely identified by the road OBJECTID
data raster = a SAR amplitude raster
mask raster = an OBJECTID road raster

Documentation/variable names are largely generalized, but sometimes refer directly
to the Middsense data
"""

import argparse
import numpy as np
import pandas as pd
import time
try:
    from osgeo import gdal
    from osgeo.gdalnumeric import *
    from osgeo.gdalconst import *
except ImportError:
    sys.exit('ERROR: cannot find GDAL/OGR modules')

def raster2table(datatif, masktif):
    """
    opens the rasters using GDAL, converts to numpy arrays, extracts the
    unmasked data pixels, and returns a DataFrame with these values

    datatif: path to the raster whose values we want (SAR amplitude raster)
    masktif: path to the raster used as a spatial mask (OID raster)
    """

    # read files
    data = gdal.Open(datatif, gdal.GA_ReadOnly)
    mask = gdal.Open(masktif, gdal.GA_ReadOnly)

    # convert .tif > single band > numpy array
    data_array = BandReadAsArray(data.GetRasterBand(1))
    mask_array = BandReadAsArray(mask.GetRasterBand(1))

    # create an OID and amplitude stack, filtering out nodata values
    # since we are filtering pixels, one dimension is lost, meaning that the
    # notion of location is only given by OID, not index in the array
    not_nodata = np.logical_and(mask_array != -9999, data_array != -9999)
    stack_array = np.stack((mask_array[not_nodata], data_array[not_nodata]), axis=1)

    # convert to a DataFrame, subsetting the data filename to extract only the
    # image acquisition date, which is unique for each
    df = pd.DataFrame(stack_array, columns=['OID', 'amp_' + str(datatif[-59:-51])])

    # close datasets
    data_array = None
    mask_array = None
    data = None
    mask = None

    return df

def summarize(df, method):
    """
    groups dataframe 'df' by 'OID' column and summarizes by method 'method'
    """
    grouped = df.groupby(by='OID')

    if method == 'mean':
        summarized = grouped.mean()
    elif method == 'count':
        summarized = grouped.count()
    elif method == 'median':
        summarized = grouped.median()
    else:
        print('invalid or no method chosen')

    return summarized

def merge(mask_tifs, data_tifs, out_path, method):
    """
    perform the calculation for all the (data, mask) pairs
    and merge them together into one csv
    """
    start_time = time.time()
    temp0 = []
    n = len(mask_tifs) * len(data_tifs)
    progress = 0

    print("calculating summaries on all " + str(n) + " (data, mask) pairs")

    for data in data_tifs:
        temp1 = []
        for mask in mask_tifs:
            temp1.append(summarize(raster2table(data, mask), method))
            progress += 1
            print("calculation completed for " + str(progress) + "/" + str(n) + " file pairs")
        temp0.append(pd.concat(temp1, axis=0))
    out = pd.concat(temp0, axis=1)

    out.to_csv(path_or_buf=out_path)

    print("execution time (seconds): {} \nfor {} data rasters, {} mask rasters".format( \
        str(time.time() - start_time), len(data_tifs), len(mask_tifs)))


if __name__ == "__main__":

    DESCRIPTION = "summarize data raster for each mask value"

    parser = argparse.ArgumentParser(description=DESCRIPTION)

    parser.add_argument('-m', '--mask_tifs', nargs='+',
                        help="path to mask .tif (required)")
    parser.add_argument('-d', '--data_tifs', nargs='+',
                        help="path to data .tif (required)")
    parser.add_argument('-o', '--out_path',
                        help="path to output csv (required)")
    parser.add_argument('-s', '--summarize_method',
                        help="summarize method: from {'mean', 'median', 'count', 'none'}")

    args = parser.parse_args()

    merge(args.mask_tifs, args.data_tifs, args.out_path, args.summarize_method)
