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

  q1, q3, iqr = q1.add_suffix('_q1'), q3.add_suffix('_q3'), iqr.add_suffix('_iqr')
  min, max = min.add_suffix('_bmin'), max.add_suffix('_bmax')

  # concatenate
  out = pd.concat([q1, q3, iqr, min, max], axis=1)

  return out

def all_metrics(df):
  """
  compute all summaries
  """

  # group by OID
  grouped = df.groupby('oid', as_index=True)

  # built in pandas functions
  # count
  count = grouped.count().add_suffix('_count')
  zerocount = df.set_index('oid').eq(0).groupby('oid').sum().add_suffix('_zerocount')

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

  # Mean filtering
  maxs = quant.filter(like='bmax')
  mins = quant.filter(like='bmins')

  merge = pd.concat([merge, quant, filtered_mean], axis=1, join='outer')

  # sort columns
  merge = merge.reindex(sorted(merge.columns), axis=1)

  return merge
