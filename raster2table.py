# Myles Stokowski
# 6.22.20
#
# script to generate a table of average SAR amplitude within each 
# unique road segment identified by OBJECTID

import sys
import argparse
import numpy as np
import pandas as pd
try: 
    from osgeo import gdal
    from osgeo.gdalnumeric import *
    from osgeo.gdalconst import *
except ImportError:
    sys.exit('ERROR: cannot find GDAL/OGR modules')

# this chunk (before defining raster_group_means) is just for
# testing to pass the files without running from command line
# so I can run interactively
maskpath = '/home/myles/middsense/OID Road Rasters/'
datapath = '/home/myles/middsense/data_from_drive/SAR_Amplitude/'
oid = pd.DataFrame()
oid_list = [oid]
for data in os.listdir('/home/myles/middsense/data_from_drive/SAR_Amplitude/'):
    for mask in os.listdir('/home/myles/middsense/OID Road Rasters/'):
        if mask[-7:] == 'All.tif':
            #print(datapath + data, maskpath + mask)
            #pd.concat(oid_list, raster_group_means(datapath + data, maskpath + mask, ''))
            oid_list.append(raster_group_means(datapath + data, maskpath + mask, ''))
oid_concat = pd.concat(oid_list)

# which OID values are missing?
unique_list = oid_concat.sort_values('OID').index.tolist()
unique_list = list( map(int, unique_list) )
missing = []
for i in range(1,max(unique_list)):
    if i not in unique_list:
        missing.append(i)

# TODO strange, OID 1 is missing, many other values between 3926 - 14711
# no apparent pattern to whats missing, total 347 missing elts
# did some poking around at the missing features in QGIS, seems like they are values
# where the road was monitored multiple times in one year (ex. OID 14711 and 14712)
# were both monitored in 2013, and only the value 14712 comes through in the rasterized
# OID roads created with vector2raster.py
# also checked with OID 3926, same deal (overlapped with OID 6868, both in year 2011)
# TODO so... does this mean that the better approach is not ever converting to raster? 
# or some workaround while still using rasters?

datatif = '/home/myles/middsense/data_from_drive/SAR_Amplitude/20110829_110031_20110829_110037_CSKS3-H4-0B_HH_SLC_data.tif'
masktif = '/home/myles/middsense/OID Road Rasters/Staunton_OID_2011_IS.tif'

def raster_group_means(datatif, masktif, dsttif):
    """
    datatif: the raster whose values we want (SAR amplitude raster)
    masktif: the raster used as a spatial mask (OID raster)
    dsttif: output
    """
    #datatype = gdal.GDT_Float32
    #drv = gdal.GetDriverByName('GTiff')

    # read files
    data = gdal.Open(datatif, gdal.GA_ReadOnly)
    mask = gdal.Open(masktif, gdal.GA_ReadOnly)

    # convert .tif > single band > numpy array
    data_array = BandReadAsArray(data.GetRasterBand(1))
    mask_array = BandReadAsArray(mask.GetRasterBand(1))
    
    # create an OID and amplitude stack with only the pixels where OID!=1 (the nodata value)
    # since we are filtering pixels, one dimension is lost, meaning that the
    # notion of location is only given by OID, not index in the array
    # TODO also need to mask out the SAR nodata values which are skewing averages (-9999) 
    # TODO OID 1 is getting masked out, need to change nodata value in the OID rasters
    stack = np.stack((mask_array[mask_array != 1], data_array[mask_array != 1]), axis = 1)
    df = pd.DataFrame(stack, columns = ['OID', 'amp_' + str(datatif)])
    means = df.groupby(by = 'OID').mean()

    return(means)
    
    # create an OID and amplitude stack with all the pixels in the extent
    # stack_all = np.stack((mask_array, data_array), axis = 2)


    # df[0].unique() # list of OIDs
    # np.unique(mask_array[mask_array != 1]) # alternatively, via numpy
    #dst_array = (mask_array!=1) * data_array
    #mask_array_01 = mask_array != 1
    #dst_array = data_array[mask_array_01]
    #mask_array_clip = mask_array[mask_array_01]
    
    #stack = np.stack((dst_array, mask_array_clip), axis =1)
    #stack = np.stack((mask_array, data_array), axis = 2)
    #stack2 = stack[(stack!=1).nonzero()[:2]]
    #np.savetxt("bar2.csv", stack, delimiter = ",")
    #print((dst_array[0,0]))
    #numpy.savetxt("foo.csv", dst_array, delimiter = ",")

    # TODO data and mask should have the same resolution, alignment, etc.
    # could add a check before proceeding

    # create output layer
    # dst = drv.Create(dsttif, data.RasterXSize, data.RasterYSize, 1, datatype, options=['COMPRESS=DEFLATE'])
    # dst.SetProjection(data.GetProjectionRef())
    # dst.SetGeoTransform(data.GetGeoTransform())
        
    # write the out file
    #dst_band = dst.GetRasterBand(1)
    #dst_band.SetNoDataValue(-9999)
    #BandWriteArray(dst_band, dst_array)

    # close datasets
    #data_band = None
    #mask_band = None
    data = None
    mask = None
    #dst_band = None
    #dst = None
    

if __name__ == "__main__":
    DESCRIPTION = "placeholder"

    parser = argparse.ArgumentParser(description=DESCRIPTION)

    #parser.add_argument("datatif")
    #parser.add_argument("masktif")
    #parser.add_argument("dsttif")

    parser.add_argument('-m', '--mask_tifs', nargs = '*')
    #parser.add_argument("mask_tif", nargs=n_)
    parser.add_argument('-d', '--data_tifs', nargs='*')
    

    args = parser.parse_args()

    print("masks: " + str(args.mask_tifs))
    print("data: " + str(args.data_tifs))

    print("masks: ")
    for mask in args.mask_tifs:
        print(mask)
    
    print("data: ")
    for data in args.data_tifs:
        print(data)

    #OIDs = []
    #for i in range(5):
        #print(len(raster_group_means(args.data_tifs[0], args.mask_tifs[i], '')))
        #OIDs.append((raster_group_means(args.data_tifs[0], args.mask_tifs[i], '').index))
    
    #OIDs_array = np.array(OIDs)
    #np.savetxt("bar3.csv", OIDs_array, delimiter = ",")

    # TODO all the means tables together have 15930 rows, less than the 
    # 16277 range(1,16277) features in Staunton_Maint_All_Buf12_10_10.dbf
    # not sure what is being left out

    
    #raster_group_means(args.datatif, args.masktif, args.dsttif)
