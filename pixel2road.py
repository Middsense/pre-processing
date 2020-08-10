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

  q1, q3, iqr, min, max = q1.add_suffix('_q1'), q3.add_suffix('_q3'), iqr.add_suffix('_iqr'), min.add_suffix('_min'), max.add_suffix('_max')

  # concatenate
  out = pd.concat([q1, q3, iqr, min, max], axis=1)

  return out

def all_metrics(df):
  """
  compute all summaries
  """
  zero = zero_count(df).rename('zero_count')

  # group by OID
  grouped = df.groupby('oid', as_index=True)

  # quantiles
  # quant = quantiles(grouped)

  # built in pandas functions
  # count
  count = grouped.count().add_suffix('_count')
  zerocount = df1.set_index('oid').eq(0).groupby('oid').sum().add_suffix('_zerocount')

  # mean, median, stddev, min, max (not boxplot min/max)
  mean = grouped.mean().add_suffix('_mean')
  median = grouped.median().add_suffix('_median')
  std = grouped.std().add_suffix('_std')
  min = grouped.min().add_suffix('_min')
  max = grouped.max().add_suffix('_max')



  # merge = pd.concat([mean, median, count, zero, quant, std, min, max], axis=1, join='outer')
  merge = pd.concat([mean, median, count, std, min, max], axis=1, join='outer')
  merge = merge.fillna(0) # fill 0s for zero_count column

  # quantiles
  quant = quantiles(grouped)
  merge = merge.concat(quant, axis=1, join='outer')

  # # join the pixels with the road-level stats
  # joined = df.join(merge, on='oid', how='left')
  #
  # # filter the pixels within the min_boxplot and max_boxplot values for that rd segment
  # filtered = joined.loc[joined['amp'].ge(joined['min_boxplot'])\
  #                     & joined['amp'].le(joined['max_boxplot'])]
  #
  # # calculate std and mean on the filtered pixels
  # filtered_mean = filtered.groupby('oid')['amp'].mean().rename('filtered_mean')
  # filtered_std = filtered.groupby('oid')['amp'].std().rename('filtered_std')
  #
  # # add the filtered stats, drop the unecessary boxplot stats (which are products
  # # of q1 and q3 and can be regenerated later)
  # merge = merge.drop(columns={'iqr', 'min_boxplot', 'max_boxplot'})
  # merge = pd.concat([merge, filtered_mean, filtered_std], axis=1, join='outer')


  return merge
