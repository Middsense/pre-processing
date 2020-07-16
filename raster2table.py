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
    from osgeo.gdalnumeric import * # pylint: disable=unused-wildcard-import
    from osgeo.gdalconst import * # pylint: disable=unused-wildcard-import
except ImportError:
    sys.exit('ERROR: cannot find GDAL/OGR modules')

from glob import glob
from scipy import sparse

datatif = '../data/SAR_Amplitude_Raw/20110829_110031_20110829_110037_CSKS3-H4-0B_HH_SLC_data.tif'
masktif = '../data/road-rasters/road-oid-2011.tif'

# def raster2table(datatif, masktif):
#     """
#     opens the rasters using GDAL, converts to numpy arrays, extracts the
#     unmasked data pixels, and returns a DataFrame with these values

#     datatif: path to the raster whose values we want (SAR amplitude raster)
#     masktif: path to the raster used as a spatial mask (OID raster)
#     """

#     # read files
#     data = gdal.Open(datatif, gdal.GA_ReadOnly)
#     mask = gdal.Open(masktif, gdal.GA_ReadOnly)

#     # convert .tif > single band > numpy array
#     data_array = BandReadAsArray(data.GetRasterBand(1))
#     mask_array = BandReadAsArray(mask.GetRasterBand(1))

#     # create an OID and amplitude stack, filtering out nodata values
#     # since we are filtering pixels, one dimension is lost, meaning that the
#     # notion of location is only given by OID, not index in the array
#     not_nodata = np.logical_and(mask_array != -9999, data_array != -9999)
#     stack_array = np.stack((mask_array[not_nodata], data_array[not_nodata]), axis=1)

#     # convert to a DataFrame, subsetting the data filename to extract only the
#     # image acquisition date, which is unique for each
#     df = pd.DataFrame(stack_array, columns=['OID', 'amp_' + str(datatif[-59:-51])])

#     # close datasets
#     data_array = None
#     mask_array = None
#     data = None
#     mask = None

#     return df


def gen_all_roads_array(all_roads_raster_path):
    '''
    Parameters
    ----------
    all_roads_raster_path : String
        Path to the rasterized road layer with all roads 2011-2015

    Returns
    -------
    all_roads_dense : np.array
        dense boolean (used as a mask) matrix of pixels where there is a road
        measurement in any year 2011-2015
    '''
    all_roads = gdal.Open(all_roads_raster_path, gdal.GA_ReadOnly)
    all_roads_dense = BandReadAsArray(all_roads.GetRasterBand(1))
    all_roads_dense[all_roads_dense == -9999] = 0 #TODO could instead set nodata value to 0 in shp -> tif conversion
    all_roads_dense = all_roads_dense.astype(dtype='bool')
    # can we mask the sar amplitude rasters using a sparse matrix version of the any_rd raster?
    # or do we need to keep it as a numpy array?
    return all_roads_dense

# TODO becuase of the road segments that go partly outside of the image footprint,
# there are some pixels in the rasterized any_rd.tif that have underlying -9999
# values in the SAR tifs that result in the sparse images having ~3000 fewer pixels
# than the sparse roads. Below is a quick fix, the long term solution is probably
# to just mask out the road segments that are not within the SAR footprint , either
# when creating the road shapefiles or rasters
temp_img = gdal.Open(sar_image_list[0], gdal.GA_ReadOnly) # arbitrary image choice
temp_img_array = BandReadAsArray(temp_img.GetRasterBand(1))
any_road_dense[temp_img_array == -9999] = False
any_road_dense.sum() # check that this value matches the size of each sparse image

def gen_sparse_image(image, any_rd_array):
    '''
    generate a single sparse SAR image
    '''
    img = gdal.Open(image, gdal.GA_ReadOnly)
    img_dense = BandReadAsArray(img.GetRasterBand(1))
    img_dense[img_dense == 0] = -1
    img_dense[img_dense == -9999] = 0
    img_dense[any_road_dense == False] = 0 # TODO this should be 0 not -9999... need to first change any_rd tif
    img_sparse = sparse.coo_matrix(img_dense)
    return img_sparse

def gen_all_sparse_images():
    """
    loop through all images, generate sparse versions, add them to a list
    TODO or dictionary with key = (date)_(raw/despeck)
    """
    sparse_images = []
    for image in sar_image_list:
        sparse_images.append(gen_sparse_image(image, any_road_dense))
    return sparse_images

def gen_sparse_roads(road_raster):
    roads = gdal.Open(road_raster, gdal.GA_ReadOnly)
    roads_dense = BandReadAsArray(roads.GetRasterBand(1))
    roads_dense[roads_dense == -9999] = 0
    roads_dense[any_road_dense == False] = 0 # temporary, remove once rd segments outside of image are handled better
    roads_sparse = sparse.coo_matrix(roads_dense)
    return roads_sparse

def gen_all_sparse_roads():
    sparse_roads = []
    for roads in road_raster_list:
        sparse_roads.append(gen_sparse_roads(roads))
    return sparse_roads

img_dict = dict(zip(sar_image_list, sparse_images))
road_dict = dict(zip(road_raster_list, sparse_roads))

# for now, just for a single road, image pair...
# mask the sparse image matrix by the sparse rd matrix
# output the data of both into a pixels DF
img0 = sparse_images[0]
rd0 = sparse_roads[0]

img0_masked = img0.multiply(rd0 > 0).tocoo()

pixels = pd.DataFrame({
    'oid': rd0.data,
    'amp': img0_masked.data})

summarize(pixels, mean)
































def summarize(df, method):
    """
    groups dataframe 'df' by 'OID' column and summarizes by method 'method'
    """
    grouped = df.groupby(by='oid')

    if method == 'mean':
        summarized = grouped.mean()
    elif method == 'count':
        summarized = grouped.count()
    elif method == 'median':
        summarized = grouped.median()
    else:
        print('invalid or no method chosen')

    return summarized



def summarize2(df):
    # handle zeros by calculating percent zero (can later decide to ignore segments based on % zero)
    # handle high outliers by calculating boxplot, removing values about


    # Calculate the number of zero amplitude values for each OID
    zero_count = df.groupby('OID')['amp_20110829'].apply(lambda x: x[x==0].count()).rename('zero_count')

    # Calculate the total number of amplitude values for each OID
    total_count = df.groupby('OID')['amp_20110829'].count().rename('total_count')

    # Calculate the percentage of zero amp values for each OID
    percent_zero = (zero_count/total_count).rename('percent_zero')

    for i in [zero_count, total_count, percent_zero]:
        print(len(i))

    out = pd.concat([zero_count, total_count, percent_zero], axis=1)

    # lets see the first 20 rows with the highest % zeros
    out.sort_values(by='percent_zero', ascending=False).head(20)

    # now, remove the zero values and calculate the summary stats

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

def merge2(mask_tifs, data_tifs):
    for data in data_tifs:
        data_array = BandReadAsArray(gdal.Open(datatif, gdal.GA_ReadOnly))


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
