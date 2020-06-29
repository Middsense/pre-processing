# automate_analysis.py
# Myles Stokowski
# 6.29.2020

# Helper script that calls a mix of scripts written by Middsense team and
# command line gdal commands. This script is intended to implement the
# functionality graphed in the OneNote collab note 'Model V2'

import os
import vector2raster as v2r
import change_raster_resolution as crr
import raster2table as r2t
import time

start = time.time()

# 00 convert road maintenance .gdb to multiple .shp
path_to_shp = '../data/road-shapefiles/'
if not os.path.exists(path_to_shp):
    os.mkdir(path_to_shp)
    path_to_gdb = '../data/Staunton_Maint.gdb'
    gdb2shp_command = 'ogr2ogr -f "ESRI Shapefile" ' + path_to_shp + ' ' + path_to_gdb
    print('converting ' + path_to_gdb + ' to .shp in ' + path_to_shp)
    os.system(gdb2shp_command)
else: 
    print(path_to_shp + ' directory already exists, skipping .gdb to .shp conversion')


# 01 rasterize .shp files, burn in OID value
reftif = '../data/SAR_Amplitude/20110829_110031_20110829_110037_CSKS3-H4-0B_HH_SLC_data.tif'
shpdir = '../data/road-shapefiles/'
outdir = '../data/road-rasters/'
if not os.path.exists(outdir):
    print('rasterizing .shp road layers')
    os.mkdir(outdir)
    for shp in os.listdir(shpdir):
        #print(shp[-19:])
        if (shp[-20:] == '_All_Buf12_10_10.shp') and (shp[-21] in ['1','2','3','4','5']):
            #outtif = '../data/road-rasters/road-oid-' + str(shp[-24:-20]) + '.tif'
            outtif = outdir + 'road-oid-' + str(shp[-24:-20]) + '.tif'
            inshp = shpdir + shp
            v2r.vector2raster(inshp, outtif, reftif, ['ATTRIBUTE=OBJECTID'])
else: 
    print(outdir + ' directory already exists, skipping .shp to .tif conversion')

# 10 merge landcover rasters
merged = '../data/landcover-merged.tif'
if not os.path.exists(merged):
    print('merging landcover raster tiles')
    cutline = '../data/SAR_Footprint/RITA2_Site1_Staunton_Amplitude_Image_Footprint.shp'
    tile_dir = '../data/Raster_Tiles_AOI/'
    tiles = ''
    for tile in os.listdir(tile_dir):
        if tile[-4:] == '.tif':
            tiles += tile_dir + tile + ' '

    warp_command = 'gdalwarp -r near -co COMPRESS=DEFLATE --config GDAL_CACHEMAX 5000 -wm 5000 ' \
    + '-cutline ' + cutline + ' ' + tiles + merged
    os.system(warp_command)
else:
    print(merged + ' already exists, skipping raster merge')

# 11 reproject and change resolution
# we use the same reftif from 01
reprojected = '../data/landcover-merged-reprojected.tif'
if not os.path.exists(reprojected):
    print("reprojecting landcover raster and matching resolution")
    crr.convert_resolution(merged, reprojected, reftif)
else: 
    print(reprojected + ' already exists, skipping reprojection')

# 12 reclassify to create road mask
mask = '../data/landcover-merged-reprojected-mask.tif'
if not os.path.exists(mask):
    print('reclassifiyng landcover')
    reclass_command = 'gdal_calc.py -A ' + reprojected + ' --outfile=' + mask + \
    ' --calc="logical_or(A==21, A==22)" --NoDataValue=-9999 --co "COMPRESS=DEFLATE" ' \
        + '--type=Float32 --format Gtiff'
    os.system(reclass_command)
else:
    print(mask + ' already exists, skipping reclassification')

# 20 TODO logical and to mask OBJECT_ID (OID) raster with landcover raster

# 30 summarize amplitude within each OID road segment
# just using command line interface becuase __main__ method 
# has significant functionality built in
summarize_command = 'python raster2table.py -d ../data/SAR_Amplitude/20110829*.tif ' + \
    '-m ../data/road-rasters/*.tif -o ../data/summarized.csv -s "mean"'
os.system(summarize_command)
print('summary table created succesfully, exiting')

# 40 join average SAR with other road data by OID

print('runtime (s): ' + str(time.time() - start))
