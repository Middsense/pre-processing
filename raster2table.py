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
    
    # create an OID and amplitude stack with only the pixels where OID!=1 (the nodata value)
    # since we are filtering pixels, one dimension is lost, meaning that the
    # notion of location is only given by OID, not index in the array
    # TODO OID 1 is getting masked out, need to change nodata value in the OID rasters
    not_nodata = np.logical_and(mask_array != -9999, data_array != -9999) 
    # stack = np.stack((mask_array[mask_array != 1], data_array[mask_array != 1]), axis = 1)
    stack = np.stack((mask_array[not_nodata], data_array[not_nodata]), axis = 1)
    
    # convert to a DataFrame, subsetting the data filename to extract only the
    # date the image was taken, which is unique for each
    df = pd.DataFrame(stack, columns = ['OID', 'amp_' + str(datatif[-59:-51])]) 
    
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
    grouped = df.groupby(by = 'OID')
    #return grouped
    if method == 'mean':
        summarized = grouped.mean()
    elif method == 'count':
        summarized = grouped.count()
    elif method == 'median':
        summarized = grouped.median()
        
    return summarized


if __name__ == "__main__":
    
    start_time = time.time()

    DESCRIPTION = "summarize data raster for each mask value"
    
    parser = argparse.ArgumentParser(description=DESCRIPTION)

    parser.add_argument('-m', '--mask_tifs', nargs = '+',
                        help = "path to mask .tif (required)")
    parser.add_argument('-d', '--data_tifs', nargs = '+',
                        help = "path to data .tif (required)")
    parser.add_argument('-o', '--out_path',
                        help = "path to output csv (required)")
    parser.add_argument('-s', '--summarize',
                        help = "summarize method: from {'mean', 'median', 'count', 'none'}")
    parser.add_argument('-i', '--interactive', action = 'store_true',
                        help = "flag to pause to confirm calculation")
    args = parser.parse_args()

    n = len(args.mask_tifs) * len(args.data_tifs)

    if args.interactive:
        proceed = input("run on all " + str(n) + " (data, mask) pairs? (y/n): ")
    else:
        print("calculating summaries on all " + str(n) + " (data, mask) pairs")

    if (not args.interactive) or proceed == 'y':
        temp0 = []
        progress = 0
        for data in args.data_tifs:
            temp1 = []
            for mask in args.mask_tifs:
                temp1.append(summarize(raster2table(data, mask), args.summarize))
                progress += 1
                print("calculation completed for " + str(progress) + "/" + str(n) + " file pairs") 
            temp0.append(pd.concat(temp1, axis = 0))

        out = pd.concat(temp0, axis = 1)    
        out.to_csv(path_or_buf = args.out_path)

        print("execution time (seconds): " + str(time.time() - start_time))
        print("for {} data rasters, {} mask rasters".format( \
            len(args.data_tifs), len(args.mask_tifs)))
    else: 
        print("exiting")
