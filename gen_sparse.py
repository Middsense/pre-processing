"""
gen_sparse.py
"""
import time
import numpy as np
import pandas as pd
from scipy import sparse
from glob import glob

try:
    from osgeo import gdal
    from osgeo.gdalnumeric import * # pylint: disable=unused-wildcard-import
    from osgeo.gdalconst import * # pylint: disable=unused-wildcard-import
except ImportError:
    sys.exit('ERROR: cannot find GDAL/OGR modules')


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
    all_roads_dense = all_roads_dense.astype(dtype='bool') #TODO this should be done during shp -> tif
    return all_roads_dense

# def clean_all_roads_array(all_roads_dense, REFTIF):
#     # TODO becuase of the road segments that go partly outside of the image footprint,
#     # there are some pixels in the rasterized any_rd.tif that have underlying -9999
#     # values in the SAR tifs that result in the sparse images having ~3000 fewer pixels
#     # than the sparse roads. Below is a quick fix, the long term solution is probably
#     # to just mask out the road segments that are not within the SAR footprint, either
#     # when creating the road shapefiles or rasters
#     temp_img = gdal.Open(REFTIF, gdal.GA_ReadOnly) # arbitrary image choice
#     temp_img_array = BandReadAsArray(temp_img.GetRasterBand(1))
#     all_roads_dense[temp_img_array == -9999] = False
#     return all_roads_dense
#     #all_roads_dense.sum() # check that this value matches the size of each sparse image

def gen_sparse_image(image, all_roads_dense):
    '''
    generate a single sparse SAR image
    '''
    img = gdal.Open(image, gdal.GA_ReadOnly)
    img_dense = BandReadAsArray(img.GetRasterBand(1))
    img_dense[img_dense == 0] = -1 # so we don't lose the zero values when we convert to sparse matrix
    img_dense[img_dense == -9999] = 0 # so converting to sparse matrix removes nodata (-9999) values
    img_dense[all_roads_dense == False] = 0 # masks out non-road pixels from SAR image
    img_sparse = sparse.coo_matrix(img_dense) # dense np array --> sparse matrix
    return img_sparse

def gen_all_sparse_images(sar_image_list, all_roads_dense):
    """
    loop through all images, generate sparse versions, add them to a list
    TODO or dictionary with key = (date)_(raw/despeck)
    """
    sparse_images = []
    image_labels = []
    for image in sar_image_list:
        image_labels.append(image.split('/')[-1][0:8]) # TODO verify we're getting just the date
        sparse_images.append(gen_sparse_image(image, all_roads_dense))
    image_dict = dict(zip(image_labels, sparse_images))
    return image_dict

def gen_sparse_roads(road_raster, all_roads_dense):
    roads = gdal.Open(road_raster, gdal.GA_ReadOnly)
    roads_dense = BandReadAsArray(roads.GetRasterBand(1))
    roads_dense[roads_dense == -9999] = 0
    roads_sparse = sparse.coo_matrix(roads_dense)
    return roads_sparse

def gen_all_sparse_roads(road_raster_list, all_roads_dense):
    sparse_roads = []
    road_labels = []
    for road_raster in road_raster_list:
        road_labels.append(road_raster.split('/')[-1][9:13]) # TODO different slicing for masked roads
        sparse_roads.append(gen_sparse_roads(road_raster, all_roads_dense))
    road_dict = dict(zip(sparse_roads, road_labels))
    return road_dict



"""
below, some quick testing, mostly to figure out what to do with each (img, rd) pair
"""

# img_dict = dict(zip(sar_image_list, sparse_images))
# road_dict = dict(zip(road_raster_list, sparse_roads))

# for now, just for a single road, image pair...
# mask the sparse image matrix by the sparse rd matrix
# output the data of both into a pixels DF
# img0 = sparse_images[0]
# rd0 = sparse_roads[0]

# img0_masked = img0.multiply(rd0 > 0).tocoo()

# pixels = pd.DataFrame({
#     'oid': rd0.data,
#     'amp': img0_masked.data})

# summarize(pixels, mean)
