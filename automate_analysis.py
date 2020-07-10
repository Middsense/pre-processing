"""
automate_analysis.py
Myles Stokowski
07.03.2020

all data available in shared drives
all scripts available from middsense github repo

required inputs:
    Staunton_maint.gdb (road maintenance data)
    RITA2_Site1_Staunton_Amplitude_Image_Footprint.shp (SAR footprint)
    SAR_Amplitude/*.tif (SAR raw or despeck images)
    Raster_Tiles_AOI/*.tif (landcover raster tiles)

usage:
- move data files or edit input path constants to ensure script can find inputs
- call generate_csv() or generate_all_csvs()
"""

import os
import time
import glob
import vector2raster as v2r
import change_raster_resolution as crr
import raster2table as r2t
import dbf2csv as d2c
import mask

start = time.time()

# define data directory and set as working directory
# (all subsequent paths are relative paths to DATA_DIR)
DATA_DIR = '../data/'
os.chdir(DATA_DIR)

# define input paths
ROAD_GDB = './Staunton_Maint.gdb'
LANDCOVER_TILE_DIR = './Raster_Tiles_AOI/'
CUTLINE = './SAR_Footprint/RITA2_Site1_Staunton_Amplitude_Image_Footprint.shp'
RAW_SAR_DIR = './SAR_Amplitude_Raw/'
DESPECK_SAR_DIR = './SAR_Amplitude_Despeckled/'
REFTIF = os.listdir(RAW_SAR_DIR)[0] # an arbitrary SAR image used only for alignment

# define output paths
ROAD_SHP_DIR = './road-shapefiles/'
ROAD_RASTER_DIR = './road-rasters/'
LANDCOVER_RASTER_DIR = './landcover-rasters/'
LANDCOVER_MERGED = LANDCOVER_RASTER_DIR + 'landcover-merged.tif'
LANDCOVER_REPROJECTED = LANDCOVER_RASTER_DIR + 'landcover-merged-reprojected.tif'
LANDCOVER_MASKED = LANDCOVER_RASTER_DIR + 'landcover-merged-reprojected-masked.tif'
ROAD_DBF = ROAD_SHP_DIR + 'Staunton_Maint_All_Buf12_10_10.dbf'
ROAD_CSV = 'all_road_features.csv'

# check that inputs are all in right places
missing_inputs = []
for inp in [ROAD_GDB, LANDCOVER_TILE_DIR, CUTLINE, RAW_SAR_DIR]:
    if not os.path.exists(inp):
        missing_inputs.append(inp)
if len(missing_inputs) > 0:
    print('exiting, missing input(s): ')
    print(*missing_inputs, sep="\n")
else:
    print('all inputs in place')

# 00 convert road maintenance .gdb to multiple shapefiles
if os.path.exists(ROAD_SHP_DIR):
    print('skipping .gdb to .shp conversion, {} directory already exists'.format(\
        ROAD_SHP_DIR))
else:
    print('converting ' + ROAD_GDB + ' to shapefiles in ' + ROAD_SHP_DIR)
    os.mkdir(ROAD_SHP_DIR)
    GDB2SHP_COMMAND = 'ogr2ogr -f "ESRI Shapefile" ' + ROAD_SHP_DIR + ' ' + ROAD_GDB
    os.system(GDB2SHP_COMMAND)

# 01 rasterize .shp files, burn in OID value
if os.path.exists(ROAD_RASTER_DIR):
    print('skipping .shp to .tif conversion, {} directory already exists'.format(\
        ROAD_RASTER_DIR))
else:
    print('rasterizing .shp road layers')
    os.mkdir(ROAD_RASTER_DIR)
    for shp in glob.glob(ROAD_SHP_DIR + '*[0-9]_All_Buf12_10_10.shp'):
        OUTTIF = ROAD_RASTER_DIR + 'road-oid-' + str(shp[-24:-20]) + '.tif'
        v2r.vector2raster(shp, OUTTIF, REFTIF, ['ATTRIBUTE=OBJECTID'])

# 02 convert the .dbf with all (buffered) features to csv
if os.path.exists(ROAD_CSV):
    print('skipping .dbf to .csv conversion, {} already exists'.format(ROAD_CSV))
else:
    print('converting {} to {}'.format(ROAD_DBF, ROAD_CSV))
    d2c.dbf2csv(ROAD_DBF, ROAD_CSV)

# 10 setup a landcover raster output folder
if not os.path.exists(LANDCOVER_RASTER_DIR):
    os.mkdir(LANDCOVER_RASTER_DIR)

# 11 merge landcover raster tiles into one .tif
if os.path.exists(LANDCOVER_MERGED):
    print('skipping raster merge, {} already exists'.format(LANDCOVER_MERGED))
else:
    print('merging landcover raster tiles')
    merge_tiles = glob.glob(LANDCOVER_TILE_DIR + '*.tif')
    WARP_COMMAND = 'gdalwarp -r near -co COMPRESS=DEFLATE --config '\
        + 'GDAL_CACHEMAX 5000 -wm 5000 -cutline ' + CUTLINE + ' ' \
        + ' '.join(merge_tiles) + ' ' + LANDCOVER_MERGED
    os.system(WARP_COMMAND)

# test of gdalwarp as a Python library function (rather than calling from cli with os)
# unfortunately, runs slower (at least as configured) than os.system implementation
# options = gdal.WarpOptions(cutlineDSName=CUTLINE, resampleAlg="near",\
#     creationOptions=["COMPRESS=DEFLATE"])
# gdal.Warp(LANDCOVER_RASTER_DIR + 'landcover-merged.tif', merge_tiles, options=options)

# 12 reproject and change resolution
if os.path.exists(LANDCOVER_REPROJECTED):
    print('skipping reprojection, {} already exists'.format(LANDCOVER_REPROJECTED))
else:
    print("reprojecting landcover raster and matching resolution")
    crr.convert_resolution(LANDCOVER_MERGED, LANDCOVER_REPROJECTED, REFTIF)

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

for road_raster in glob.glob(ROAD_RASTER_DIR + '*[0-9].tif'):
    out = road_raster[:-4] + '-masked.tif'
    if os.path.exists(out):
        print('skipping mask, {} already exists'.format(out))
    else:
        print('creating masked road raster {}'.format(out))
        mask.mask(road_raster, LANDCOVER_MASKED, out)

# 30 summarize into .csv
def generate_csv(masked_by_landcover, despeckled, method):
    '''
    masked_by_landcover: {True, False}
    despeckled: {True, False}
    method: {'mean', 'median', 'count'}
    '''
    if masked_by_landcover:
        mask_tifs = glob.glob(ROAD_RASTER_DIR + '*masked.tif')
        mask_tag = 'landcover_mask'
    else:
        mask_tifs = glob.glob(ROAD_RASTER_DIR + '*[0-9].tif')
        mask_tag = 'no_mask'

    if despeckled:
        data_tifs = glob.glob(DESPECK_SAR_DIR + '*.tif')
        despeck_tag = 'despeck'
    else:
        data_tifs = glob.glob(RAW_SAR_DIR + '*.tif')
        despeck_tag = 'raw'

    out_path = str('{}{}_{}_{}_sar_amplitude_all_images_all_roads.csv'.format(\
        DATA_DIR, method, despeck_tag, mask_tag))
    print('starting to create ' + out_path)

    r2t.merge(mask_tifs, data_tifs, out_path, method)

def generate_all_csvs():
    '''
    generate output csvs for all combinations of options
    '''
    for i in [False]: #[True, False]:
        for j in [True, False]:
            for method in ['mean', 'median', 'count']:
                generate_csv(i, j, method)

generate_all_csvs()

print('total runtime (s): ' + str(time.time() - start))
