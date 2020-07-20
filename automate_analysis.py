"""
automate_analysis.py
Myles Stokowski
07.03.2020

all data available in shared drives
all scripts available from middsense github repo

required inputs:
    Staunton_maint.gdb (road maintenance data)
    RITA2_Site1_Staunton_Amplitude_Image_Footprint.shp (SAR footprint)
    SAR_Amplitude/*.tif (SAR raw and/or despeck images)
    Raster_Tiles_AOI/*.tif (landcover raster tiles)

usage:
- move data files or edit input path constants to ensure script can find inputs
- call generate_csv() or generate_all_csvs()
"""

import os
import time
from glob import glob
import pandas as pd
from scipy import sparse

import vector2raster as v2r
import change_raster_resolution as crr
#import raster2table as r2t
import gen_sparse
import mask

# TODO using these libraries should probably happen just in other files
# called by automate_analysis
import geopandas as gpd
try:
    from osgeo import gdal
    from osgeo.gdalnumeric import * # pylint: disable=unused-wildcard-import
    from osgeo.gdalconst import * # pylint: disable=unused-wildcard-import
except ImportError:
    sys.exit('ERROR: cannot find GDAL/OGR modules')

start = time.time()

# define constants

# define data directory and set as working directory
# (all subsequent paths are relative paths to DATA_DIR)
DATA_DIR = '../data2/'
os.chdir(DATA_DIR)

# define input paths

# road quality data (a .gdb and the layers we want in that .gdb)
ROAD_GDB = './Staunton_Maint.gdb'
BUFFERED_SINGLE_YEAR_ROAD_LAYERS = [
    'Staunton_Maint_2011_All_Buf12_10_10', 'Staunton_Maint_2012_All_Buf12_10_10',
    'Staunton_Maint_2013_All_Buf12_10_10', 'Staunton_Maint_2014_All_Buf12_10_10',
    'Staunton_Maint_2015_All_Buf12_10_10']
CENTERLINE_SINGLE_YEAR_ROAD_LAYERS = ['Staunton_Maint_2011_All', 'Staunton_Maint_2012_All',
    'Staunton_Maint_2013_All', 'Staunton_Maint_2014_All', 'Staunton_Maint_2015_All']


# landcover data
LANDCOVER_TILE_DIR = './Raster_Tiles_AOI/'

# SAR data
FOOTPRINT = './SAR_Footprint/RITA2_Site1_Staunton_Amplitude_Image_Footprint.shp'
RAW_SAR_DIR = './SAR_Amplitude_Raw/'
DESPECK_SAR_DIR = './SAR_Amplitude_Despeckled/'
REFTIF = os.listdir(RAW_SAR_DIR)[0] # an arbitrary SAR image used only for alignment

# define output paths

# directories
ROAD_SHP_DIR = './road-shapefiles/'
ROAD_RASTER_DIR = './road-rasters/'
LANDCOVER_RASTER_DIR = './landcover-rasters/'
CSV_DIR = './csv/'
SPARSE_DIR = './sparse/'

# files
LANDCOVER_MERGED = LANDCOVER_RASTER_DIR + 'landcover-merged.tif'
LANDCOVER_REPROJECTED = LANDCOVER_RASTER_DIR + 'landcover-merged-reprojected.tif'
LANDCOVER_MASKED = LANDCOVER_RASTER_DIR + 'landcover-merged-reprojected-masked.tif'
ROAD_CSV = CSV_DIR + 'all_road_features.csv'
ALL_ROADS_SHP = ROAD_SHP_DIR + 'all_roads_within_footprint.shp'


# check if inputs are in place
def check_inputs():
    missing_inputs = []
    for inp in [ROAD_GDB, LANDCOVER_TILE_DIR, FOOTPRINT, RAW_SAR_DIR, DESPECK_SAR_DIR]:
        if not os.path.exists(inp):
            missing_inputs.append(inp)
    if len(missing_inputs) > 0:
        print('exiting, missing input(s): ')
        print(*missing_inputs, sep="\n")
    else:
        print('all inputs in place')

check_inputs()

# 00
# TODO this should probably get its own file or at least method
# TODO need to handle creating an all_roads.shp (and raster also)
# read in road maintenance .gdb
# select only features completely within footprint
# write a .shp file for each year (2011-2015)
# write a single .shp with all the road features (all years)
# write a single .csv with all the road features (all years)
if os.path.exists(ROAD_SHP_DIR):
    print('skipping .gdb to .shp conversion, {} directory already exists'.format(ROAD_SHP_DIR))
else:
    print('converting ' + ROAD_GDB + ' to shapefiles and csv')
    os.mkdir(ROAD_SHP_DIR)
    os.mkdir(CSV_DIR)
    footprint = gpd.read_file(filename=FOOTPRINT)

    out_gdfs = []

    for layer in BUFFERED_SINGLE_YEAR_ROAD_LAYERS + CENTERLINE_SINGLE_YEAR_ROAD_LAYERS:

        # read directly from .gdb with geopandas
        in_layer = gpd.read_file(filename=ROAD_GDB, layer=layer)

        # keep only road segments completely within border
        join = gpd.sjoin(in_layer, footprint, how='left', op='within')
        join = join.loc[join['index_right'] == 0]

        outshp = layer + '_within_footprint.shp'
        join.to_file(ROAD_SHP_DIR + outshp)

        # only output a merged all_roads csv and shapefile for the buffered roads
        if layer[-12:] == '_Buf12_10_10':
            out_gdfs.append(join)

    # merge all features from all years
    all_road_features = pd.concat(out_gdfs)
    all_road_features.to_file(ROAD_SHP_DIR + 'all_roads_within_footprint.shp')
    all_road_features.drop(columns='geometry').to_csv(ROAD_CSV) # IRI .csv


# 01 rasterize .shp files, burn in OID value
if os.path.exists(ROAD_RASTER_DIR):
    print('skipping .shp to .tif conversion, {} directory already exists'.format(\
        ROAD_RASTER_DIR))
else:
    print('rasterizing .shp road layers')
    os.mkdir(ROAD_RASTER_DIR)

    for shp in glob(ROAD_SHP_DIR + '*.shp'):
        OUTTIF = shp.replace(ROAD_SHP_DIR, ROAD_RASTER_DIR).replace('.shp', '.tif')
        v2r.vector2raster(shp, OUTTIF, RAW_SAR_DIR + REFTIF, ['ATTRIBUTE=OBJECTID'])

# 02 create all roads boolean raster (a mask for creating sparse SAR matrices)
if os.path.exists(ROAD_RASTER_DIR + 'all_roads_bool.tif'):
    print('skipping creating road mask "all_roads_bool.tif"')
else:
    print('creating road mask "all_roads_bool.tif"')
    all_roads_raster = ROAD_RASTER_DIR+'all_roads_within_footprint.tif'
    out_mask = ROAD_RASTER_DIR+'all_roads_bool.tif'
    MASK_COMMAND = 'gdal_calc.py -A {} --calc="A != -9999" --outfile {}'.format(all_roads_raster, out_mask) +\
    ' --co "COMPRESS=DEFLATE" --type=Byte --co="NBITS=1" --NoDataValue=0 --format Gtiff'
    os.system(MASK_COMMAND)

# 10 setup a landcover raster output folder
if not os.path.exists(LANDCOVER_RASTER_DIR):
    os.mkdir(LANDCOVER_RASTER_DIR)

# 11 merge landcover raster tiles into one .tif
if os.path.exists(LANDCOVER_MERGED):
    print('skipping raster merge, {} already exists'.format(LANDCOVER_MERGED))
else:
    print('merging landcover raster tiles')
    merge_tiles = glob(LANDCOVER_TILE_DIR + '*.tif')
    WARP_COMMAND = 'gdalwarp -r near -co COMPRESS=DEFLATE --config '\
        + 'GDAL_CACHEMAX 5000 -wm 5000 -cutline ' + FOOTPRINT + ' ' \
        + ' '.join(merge_tiles) + ' ' + LANDCOVER_MERGED
    os.system(WARP_COMMAND)

# test of gdalwarp as a Python library function (rather than calling from cli with os)
# unfortunately, runs slower (at least as configured) than os.system implementation
# options = gdal.WarpOptions(cutlineDSName=FOOTPRINT, resampleAlg="near",\
#     creationOptions=["COMPRESS=DEFLATE"])
# gdal.Warp(LANDCOVER_RASTER_DIR + 'landcover-merged.tif', merge_tiles, options=options)

# 12 reproject and change resolution
if os.path.exists(LANDCOVER_REPROJECTED):
    print('skipping reprojection, {} already exists'.format(LANDCOVER_REPROJECTED))
else:
    print("reprojecting landcover raster and matching resolution")
    crr.convert_resolution(LANDCOVER_MERGED, LANDCOVER_REPROJECTED, RAW_SAR_DIR + REFTIF)

# 13 reclassify to create road mask
if os.path.exists(LANDCOVER_MASKED):
    print('skipping reclassification, {} already exists'.format(LANDCOVER_MASKED))
else:
    print('reclassifiyng landcover')
    reclass_command = 'gdal_calc.py -A ' + LANDCOVER_REPROJECTED + ' --outfile='\
        + LANDCOVER_MASKED + ' --calc="logical_or(A==21, A==22)" ' \
        + '--NoDataValue=-9999 --co "COMPRESS=DEFLATE" --type=Float32 --format Gtiff'
    os.system(reclass_command)

# 20 mask OBJECT_ID (OID) rasters with landcover raster
for road_raster in glob(ROAD_RASTER_DIR + '*.tif'):
    out = road_raster[:-4] + '-masked.tif'
    if os.path.exists(out):
        print('skipping mask, {} already exists'.format(out))
    else:
        print('creating masked road raster {}'.format(out))
        mask.mask(road_raster, LANDCOVER_MASKED, out)

##############################################################################
# at this point in pipeline, all intermediate files have been created, now we
# proceed by summarizing the SAR data by road segment into tables/csvs
##############################################################################

# 30 load the SAR images and road rasters into memory as sparse matrices

def save_sparse(sparse_dict, path):
    """
    helper function to save each element in dict as .npz files
    """
    for key in sparse_dict:
        sparse.save_npz(path+key, sparse_dict[key])

def gen_all_sparse():
    """
    Generate all sparse images as dictionaries of sparse matrices and save as .npz
    """
    if os.path.exists(SPARSE_DIR):
        print('sparse matrix directory already exists, skipping rater --> sparse matrix conversion')
    else:
        os.mkdir(SPARSE_DIR)
        os.mkdir(SPARSE_DIR+'RAW/')
        os.mkdir(SPARSE_DIR+'DESPECKLED/')
        os.mkdir(SPARSE_DIR+'roads/')

        # dense roads to mask each SAR image during conversion to sparse matrix (output is numpy array)
        all_roads_dense = gen_sparse.gen_all_roads_array(ROAD_RASTER_DIR+'all_roads_bool.tif')

        # Raw images
        raw_imgs = glob(RAW_SAR_DIR+'*.tif')
        sparse_raw = gen_sparse.gen_all_sparse_images(raw_imgs, all_roads_dense) # dict of 67 sparse matrices
        save_sparse(sparse_raw, SPARSE_DIR+'RAW/')

        # Despeckled images
        despeck_imgs = glob(DESPECK_SAR_DIR+'*.tif')
        sparse_despeck = gen_sparse.gen_all_sparse_images(despeck_imgs, all_roads_dense) #dict of 67 sparse matrices
        save_sparse(sparse_despeck, SPARSE_DIR+'DESPECKLED/')

        # sparse matrices for single year road rasters
        road_rasters = glob(ROAD_RASTER_DIR+'*.tif')
        sparse_roads = gen_sparse.gen_all_sparse_roads(road_rasters) #out: list??
        save_sparse(sparse_roads, SPARSE_DIR+'roads/')


# now that we have all the data loaded into memory, we'll do summarization calculations
# for each (road, SAR image) pair
# TODO logic to loop through the pairs/combinations

# load the sparse matrices back from files into memory
#for sar_sparse in glob(SPARSE_IMAGE_DIR)

# choose SAR {raw, despeck}
# choose road {centerlines, buffered, buffered and masked by landcover}
SPARSE_IMAGE_DIR = './Sparse/DESPECKLED/'
SPARSE_ROAD_DIR = './Sparse/Roads/'

concat_list = []

for sar_sparse in glob(SPARSE_IMAGE_DIR + '*'):
    img = sparse.load_npz(sar_sparse)
    for road_sparse in glob(SPARSE_ROAD_DIR + '*'):
        rd = sparse.load_npz(road_sparse)

        img_masked = img.multiply(rd > 0).tocoo()
        # check same size
        if img_masked.shape != rd.shape:
            print('error, sparse matrix sizes dont match!!!')

        # convert to DataFrame
        sar_col_name = sar_sparse[0:8] # slice date from SAR filename #TODO need to handle raw/despeck
        pixels = pd.DataFrame({
            'oid': rd.data,
            sar_col_name: img_masked.data
        })

        concat_list.append(pixels)



# for each of sparse matrices (img, rd):
    # img_masked = img.multiply(rd > 0).tocoo()
    # pixels = pd.DataFrame({
    #     'oid': rd.data,
    #     'amp': img_masked.data
    # })
    # summarize(pixels)

# TODO logic to get all these summary stats into a nice DataFrame with
# columns labels in format "SARdate_statistic" --> "20110829_mean"
# rows labels OID
# filenames: stats_raw/despeck_landcover/nomask_centerline/buffered.csv

# 30 summarize into .csv
def generate_csv(masked_by_landcover, despeckled):
    '''
    masked_by_landcover: {True, False}
    despeckled: {True, False}
    '''
    if masked_by_landcover:
        mask_tifs = glob(ROAD_RASTER_DIR + '*masked.tif')
        mask_tag = 'landcover_mask'
    else:
        mask_tifs = glob(ROAD_RASTER_DIR + '*[0-9].tif')
        mask_tag = 'no_mask'

    if despeckled:
        data_tifs = glob(DESPECK_SAR_DIR + '*.tif')
        despeck_tag = 'despeck'
    else:
        data_tifs = glob(RAW_SAR_DIR + '*.tif')
        despeck_tag = 'raw'

    out_path = str('{}_{}_{}_sar_amplitude_all_images_all_roads.csv'.format(\
        CSV_DIR, despeck_tag, mask_tag))
    print('starting to create ' + out_path)

    # TODO no longer using merge....
    #r2t.merge(mask_tifs, data_tifs, out_path)

# def generate_all_csvs():
#     '''
#     generate output csvs for all combinations of options
#     '''
#     for i in [False]: #[True, False]:
#         for j in [True, False]:
#             for method in ['mean', 'median', 'count']:
#                 generate_csv(i, j, method)
#
# generate_all_csvs()

print('total runtime (s): ' + str(time.time() - start))
