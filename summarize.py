"""
7/16/2020
summarize.py

functions to compute summary statistics on amplitude data 
"""
import numpy as np
import pandas as pd
import scipy.sparse as sp

from osgeo import ogr, osr, gdal
from osgeo.gdalnumeric import *
from osgeo.gdalconst import *

from glob import glob
import argparse

"""
Summary Functions
"""

def zero_count(df):
  """
  count number of amplitude 0s per object
  """
  return df.loc[df['amp']==0, 'oid'].value_counts()

def total_count(df):
  """
  counts number of pixels in region
  """
  return df['oid'].value_counts()

def quantiles(grouped):
  """
  computes quantiles, iqr, and min and max of iqr
  """
  # quantiles
  q1 = grouped.quantile(0.25)
  q3 = grouped.quantile(0.75)

  # IQR and 'box plot' extents
  iqr = q3 - q1
  iqr_prod = 1.5 * iqr
  min = q1 - iqr_prod
  max = q3 + iqr_prod

  # concatenate
  out = pd.concat([q1, q3, iqr, min, max], axis=1)
  out.columns = ['q1', 'q3', 'iqr', 'min', 'max']

  return out

def all_metrics(df):
  """
  compute all summaries
  """
  count = total_count(df).rename('count')
  zero = zero_count(df).rename('zero_count')

  # group by OID
  grouped = df.groupby('OID')

  # quantiles
  quant = quantiles(grouped)

  # mean and median
  mean = grouped.mean()
  median = grouped.median()

  merge = pd.concat([mean, median, count, zero, quant], axis=1, join='outer')

  return merge

"""
nothing to see here :)
"""


# if __name__ == "__main__":
#
#     DESCRIPTION = "Compute summary statistics of SAR images over road segments"
#
#     parser = argparse.ArgumentParser(description=DESCRIPTION)
#
#     parser.add_argument('-r', '--road_tifs', nargs='+',
#                         help="path to directory of road OID .tif files (required)")
#
#     parser.add_argument('-d', '--data_tifs', nargs='+',
#                         help="path to directory of SAR .tif files (required)")
#
#     parser.add_argument('-o', '--out_path',
#                         help="path to output csv (required)")
#
#
#     args = parser.parse_args()
#
#
#     data_tifs = glob(args.data_tifs)
#     road_tifs = glob(args.road_tifs)

    # TODO: CALL SPARSE MATRIX THING!


    ## there was other stuff down here if we want to call this as its on script
