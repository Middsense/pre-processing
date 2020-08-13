# pre-processing

**automate_analysis.py**

This is the main driver script that performs most of the spatial analysis. Input paths should lead to the raw data; outputs will be a series of .pkl files with additional statistics.

Two flags may be set at the beginning of the script:
* LOG_FILTER uses log-correction when computing statistics
* MEAN_NORMALIZE subtracts a pixel-wise temporal mean reference image from each image in the stack before computing statistics.

Required inputs:
* Staunton_maint.gdb (road maintenance data)
* RITA2_Site1_Staunton_Amplitude_Image_Footprint.shp (SAR footprint)
* SAR_Amplitude/\*.tif (SAR raw and/or despeck images)
* Raster_Tiles_AOI/\*.tif (landcover raster tiles)

**change_raster_resolution.py**

Helper script for changing the projection/alignment of a .tif file. Used to change 1m Virginia landcover dataset to match the 3m resolution and extent of the SAR data.

Input: Raster to realign (.tif), reference .tif from SAR dataset

Output: Realigned .tif file


**gen_sparse.py**

Functions for creating sparse matrices for image and road OID datasets.
* gen_all_roads_array returns a dense boolean mask of pixels where there is a road measurement at any point in the dataset
* gen_sparse_image creates a single sparse SAR image
* gen_all_sparse_images creates a sparse SAR image stack as a dictionary
* gen_sparse_roads creates a single sparse road OID layer
* gen_all_sparse_roads creates the whole sparse OID stack

**pixel2road.py**

Handles summary statistics at the pixel level.

Input: DataFrame of pixel level data for all images with one set of OIDs

Output: DataFrame of pixel level data with mean, median, count, standard deviation, minimum, maximum, and "boxplot" min/max


**pixelwise.py**

Creates the SAR image stack from the sparse matrices to prepare for pixelwise statistics (see pixel2road.py) The output DataFrame has the row and column identifying the pixel location in the raster, any corresponding OIDs from all the OID datasets, and the pixel's amplitude value for each image in the stack.

Input: Sparse matrices (.npz) for SAR images and road OIDs

Output: DataFrame representing the image stack


**road_stats.py**

Handles data cleaning and aggregation at the road level. Removes outliers and null values, zero amplitude values, and invalid date ranges; adds road-level statistics and finds the closest SAR acquisition date for each IRI measurements; joins road quality data and amplitude data

**vector2raster.py**

Converts a .shp file to a .tif raster using the extent/resolution of a reference .tif

Input: .shp file to be rasterized; reference .tif file for the region

Output: .tif rasterized version
