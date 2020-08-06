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
import pandas as pd
import numpy as np
from glob import glob
from scipy import sparse

# other Middsense files
import vector2raster as v2r
import change_raster_resolution as crr
import gen_sparse
# import mask
import summarize
import merge_clean

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

"""
OPTIONS
"""
# compute statistics with log: exp(average(log(data)))
LOG_FILTER = False

"""
WORKING DIRECTORY
"""
# define data directory and set as working directory
# (all subsequent paths are relative paths to DATA_DIR)
DATA_DIR = '../data/'
os.chdir(DATA_DIR)

"""
INPUTS
(these should exist before running)
"""
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

"""
OUTPUTS
(these will all be created as the script runs)
"""
# directories
ROAD_SHP_DIR = './road-shapefiles/'
ROAD_RASTER_DIR = './road-rasters/'
LANDCOVER_RASTER_DIR = './landcover-rasters/'
CSV_DIR = './csv/'
SPARSE_DIR = './sparse/'
PKL_DIR = './pickles/'

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

"""
00 Formatting road data
TODO this should probably get its own file or at least method


read in road maintenance .gdb
select only features completely within footprint
write a .shp file for each year (2011-2015), both for buffered and centerline roads
write a .shp with all road features (all years), both for buffered and centerline roads
write a single .csv with all the road features (all years)
"""
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
    print('skipping creating road mask "all_roads_bool.tif", file already exists')
else:
    print('creating road mask "all_roads_bool.tif"')
    all_roads_raster = ROAD_RASTER_DIR+'all_roads_within_footprint.tif'
    out_mask = ROAD_RASTER_DIR+'all_roads_bool.tif'
    MASK_COMMAND = 'gdal_calc.py -A {} --calc="A != -9999" --outfile {}'.format(all_roads_raster, out_mask) +\
    ' --co "COMPRESS=DEFLATE" --type=Byte --co="NBITS=1" --NoDataValue=0 --format Gtiff --quiet'
    os.system(MASK_COMMAND)

"""
Landcover masking steps
"""

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
for road_raster in glob(ROAD_RASTER_DIR + '*within_footprint.tif'):
    out = road_raster[:-4] + '_landcovermasked.tif'
    if os.path.exists(out):
        print('skipping mask, {} already exists'.format(out))
    else:
        print('creating masked road raster {}'.format(out))
        #mask.mask(road_raster, LANDCOVER_MASKED, out)
        MASK_COMMAND = 'gdal_calc.py -A {} -B {} --calc="-9999*(B==0) + A*(B==1)" --outfile {}'.format(road_raster, LANDCOVER_MASKED, out) +\
        ' --co "COMPRESS=DEFLATE" --type=Float32 --NoDataValue=-9999 --format Gtiff --quiet'
        os.system(MASK_COMMAND)


"""
Sparse matrix processing of images and road features

at this point in pipeline, all intermediate files have been created, now we
proceed by summarizing the SAR data by road segment into tables/csvs
"""

# 30 generate sparse matrix versions of the road rasters and SAR images
# clipped to the road pixels, and save these sparse matrices as files (.npz)

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
        print('sparse matrix directory already exists, skipping raster --> sparse matrix conversion')
    else:
        print('generating sparse matrices and saving to .npz')
        os.mkdir(SPARSE_DIR)
        os.mkdir(SPARSE_DIR+'raw/')
        os.mkdir(SPARSE_DIR+'despeckled/')
        os.mkdir(SPARSE_DIR+'roads/')

        # dense roads to mask each SAR image during conversion to sparse matrix (output is numpy array)
        all_roads_dense = gen_sparse.gen_all_roads_array(ROAD_RASTER_DIR+'all_roads_bool.tif')

        # Raw images
        print('generating raw SAR sparse matrices')
        raw_imgs = glob(RAW_SAR_DIR+'*.tif')
        sparse_raw = gen_sparse.gen_all_sparse_images(raw_imgs, all_roads_dense) # dict of 67 sparse matrices
        save_sparse(sparse_raw, SPARSE_DIR+'raw/')

        # Despeckled images
        print('generating despeckled SAR sparse matrices (hang in there)')
        despeck_imgs = glob(DESPECK_SAR_DIR+'*.tif')
        sparse_despeck = gen_sparse.gen_all_sparse_images(despeck_imgs, all_roads_dense) #dict of 67 sparse matrices
        save_sparse(sparse_despeck, SPARSE_DIR+'despeckled/')

        # sparse matrices for single year road rasters
        # this should generate 5 years * 2 (buffered/centerlines) * 2 (landcovermasked/notmasked) = 20 sparse matrices
        print('generating road OID sparse matrices')
        road_rasters = glob(ROAD_RASTER_DIR+'Staunton_Maint_[0-9]*.tif')
        sparse_roads = gen_sparse.gen_all_sparse_roads(road_rasters) #dict of sparse matrices for each road raster
        save_sparse(sparse_roads, SPARSE_DIR+'roads/')

# Generate and save all sparse matrices: raw, despeck, and roads
# gen_all_sparse()

# TODO could load into memory before looping, low priority

"""
Summarizing!
"""

sar_type = ['raw/', 'despeckled/']
buffer_type = ['All_Buf12_10_10_within_footprint', 'All_within_footprint']
mask_type = ['_landcovermasked.npz', '.npz']

# list of possible combinations of data sources (raw/despeck, buffered/not buffered, masked/unmasked)
datagroups = np.array(np.meshgrid(sar_type, buffer_type, mask_type)).T.reshape(-1, 3)
names = np.array(np.meshgrid(['raw', 'despeck'], ['buffered', 'centerline'], ['masked', 'unmasked'])).T.reshape(-1, 3)


print('Summarizing...')
# 31 Generates a .csv file for each combination of data inputs
for i, dg in enumerate(datagroups):

    out_path = str('{}{}_{}_{}.csv'.format(CSV_DIR, names[i][0], names[i][1], names[i][2]))

    # check if this summary .csv already exists
    if os.path.exists(out_path):
        print('skipping summary of {}, file already exists'.format(out_path))

    else:
        print('converting ' + ROAD_GDB + ' to shapefiles and csv')

        # For each .csv, we summarize each image over the road files for each year:
        summarylist = []

        # Each image in selected stack (raw/despeck)
        for img_path in glob(SPARSE_DIR + dg[0] + '*'):

            img = sparse.load_npz(img_path)

            concat_list = []

            # Each year of road segments matching the road type we want
            for rd_path in glob(SPARSE_DIR + 'roads/*' + dg[1] + dg[2]):

                rd = sparse.load_npz(rd_path)

                # mask SAR image to only pixels where there was an OID road segment that year
                img_masked = img.multiply(rd > 0).tocoo()
                rd_masked = rd.multiply((img_masked > 0) + (img_masked == -1)).tocoo()

                # Check that the masked files are the same size
                if img_masked.size != rd_masked.size:
                    print('error, sparse matrix sizes dont match!!!')

                # Convert to DataFrame
                sar_col_name = img_path.split('/')[-1][0:8] + '_'


                pixels = pd.DataFrame({
                    'oid':rd_masked.data,
                    'amp': img_masked.data
                })

                pixels = pixels.replace(-1, 0)

                concat_list.append(pixels)


            # concatenate year lists for this image
            all_oids = pd.concat(concat_list, axis=0)

            # Summarize !!!!
            summarized = summarize.all_metrics(all_oids, LOG_FILTER)
            summarized.columns = sar_col_name + summarized.columns

            # list of dataframes by image
            summarylist.append(summarized)


        fullsummary = pd.concat(summarylist, axis=1)
        fullsummary.index.name = 'oid'

        fullsummary.to_csv(out_path)

        print('summarized ' + out_path)

"""
Merging and cleaning the datasets
"""

if os.path.exists(PKL_DIR):
    print('skipping merge/clean, {} directory already exists'.format(\
        PKL_DIR))
else:
    print('merging and cleaning datasets')
    os.mkdir(PKL_DIR)

    # dictionary to hold datasets
    sar_datasets = {}

    for path in glob(CSV_DIR + '*masked*'):
        key = path.split('/')[-1][:-4]
        df = pd.read_csv(path, index_col='oid')
        sar_datasets[key] = df

    roads = pd.read_csv(CSV_DIR+'all_road_features.csv')

    print("merging...")

    merged_datasets = {}

    for key in sar_datasets:
        # join and merge
        joined = merge_clean.join_roads(sar_datasets[key], roads)
        merged_datasets[key] = merge_clean.clean(joined)

        # save output as .pkl
        out_path = PKL_DIR + key + '.pkl'
        merged_datasets[key].to_pickle(path=out_path)



print('total runtime (s): ' + str(time.time() - start))
