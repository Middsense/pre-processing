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
    """
    Parameters
    ----------
    all_roads_raster_path : String
        Path to the rasterized road layer with all roads 2011-2015

    Returns
    -------
    all_roads_dense : np.array
        dense boolean (used as a mask) matrix of pixels where there is a road
        measurement in any year 2011-2015
    """
    all_roads = gdal.Open(all_roads_raster_path, gdal.GA_ReadOnly)
    all_roads_dense = BandReadAsArray(all_roads.GetRasterBand(1))
    return all_roads_dense

def gen_sparse_image(image, all_roads_dense):
    """
    generate a single sparse SAR image
    """
    img = gdal.Open(image, gdal.GA_ReadOnly)
    img_dense = BandReadAsArray(img.GetRasterBand(1))

    img_dense[img_dense == 0] = -1 # so we don't lose the zero values when we convert to sparse matrix
    img_dense[img_dense == -9999] = 0 # so converting to sparse matrix removes nodata (-9999) values
    img_dense[all_roads_dense == False] = 0 # masks out non-road pixels from SAR image
    img_sparse = sparse.coo_matrix(img_dense) # dense np array --> sparse matrix

    return img_sparse

def gen_all_sparse_images(sar_image_list, all_roads_dense):
    """
    loop through all images, generate sparse versions, add them to a dictionary
    """
    sparse_images = []
    image_labels = []

    nimgs = len(sar_image_list)
    progress = 0
    for image in sar_image_list:
        image_labels.append(image.split('/')[-1][0:8])
        sparse_images.append(gen_sparse_image(image, all_roads_dense))

        progress += 1
        print('Generated '+str(progress)+' of '+str(nimgs)+' sparse SAR matrices')

    image_dict = dict(zip(image_labels, sparse_images))
    return image_dict

def gen_sparse_roads(road_raster):
    """
    Generate a single sparse road (OID) image
    """
    roads = gdal.Open(road_raster, gdal.GA_ReadOnly)
    roads_dense = BandReadAsArray(roads.GetRasterBand(1))

    roads_dense[roads_dense == -9999] = 0
    roads_sparse = sparse.coo_matrix(roads_dense)

    return roads_sparse

def gen_all_sparse_roads(road_raster_list):
    """
    Loops through all road images, generates sparse versions, and adds them to a dictionary
    """
    sparse_roads = []
    road_labels = []

    nroads = len(road_raster_list)
    progress = 0
    for road_raster in road_raster_list:
        road_labels.append(road_raster.split('/')[-1][:-4])
        sparse_roads.append(gen_sparse_roads(road_raster))

        progress+=1
        print('Generated '+str(progress)+' of '+str(nroads)+' sparse road matrices')

    road_dict = dict(zip(road_labels, sparse_roads))
    return road_dict
