"""
pixel2road.py
Created: 7/16/2020
Last Updated: 8/10/2020

functions to compute summary statistics on amplitude data
"""
import numpy as np
import pandas as pd

"""
Summary Functions
"""


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
  
  # rename columns
  q1, q3, iqr = q1.add_suffix('_q1'), q3.add_suffix('_q3'), iqr.add_suffix('_iqr')
  min, max = min.add_suffix('_bmin'), max.add_suffix('_bmax')

  return pd.concat([q1, q3, iqr, min, max], axis=1)

def all_metrics(df):
  """
  compute all summaries
  """
  
  # group by OID
  grouped = df.groupby('oid', as_index=True)

  # count pixels per road segment
  count = grouped.count().add_suffix('_count')
  zerocount = df.eq(0).groupby('oid').sum().add_suffix('_zerocount')

  # mean, median, stddev, min, max (not boxplot min/max)
  mean = grouped.mean().add_suffix('_mean')
  median = grouped.median().add_suffix('_median')
  std = grouped.std().add_suffix('_std')
  min = grouped.min().add_suffix('_min')
  max = grouped.max().add_suffix('_max')

  # combine the statistics series --> dataframe 
  merge = pd.concat([mean, median, count, zerocount, std, min, max], axis=1, join='outer')
  merge = merge.fillna(0) # fill 0s for zero_count column

  # quantiles
  quant = quantiles(grouped)

  # pixel filtering based on boxplot min/max
  maxs = quant.filter(like='bmax')
  mins = quant.filter(like='bmin')
  maxs.columns = (maxs.columns.str.rstrip('_bmax'))
  mins.columns = (mins.columns.str.rstrip('_bmin'))

  mask = df.lt(maxs) & df.gt(mins)
  
  # quick check if axes are same for mask and data DataFrames
  if not np.unique(df.index == mask.index) and np.unique(df.columns == mask.columns):
    print('error! pixel dataframe and boxplot filtering mask dataframe do not align')

  # apply mask, calculate stats on filtered pixels
  filtered = df[mask]  
  filtered_mean = filtered.groupby('oid').mean().add_suffix('_filtered_mean')
  filtered_std = filtered.groupby('oid').std().add_suffix('_filtered_std')
  
  # add the filtered stats and 0.25 and 0.75 quantiles
  merge = pd.concat([merge, quant.filter(like='q1'), quant.filter(like='q3'), \
                     filtered_mean, filtered_std], axis=1, join='outer')

  # sort columns
  merge = merge.sort_index(axis=1)

  return merge
