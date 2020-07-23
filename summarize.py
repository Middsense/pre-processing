"""
7/16/2020
summarize.py

functions to compute summary statistics on amplitude data
"""
import numpy as np
import pandas as pd

"""
Summary Functions
"""

def zero_count(df):
  """
  count number of amplitude 0s per object
  only calculates count for objects with any zeros, zero-filling occurs later
  """
  return df.loc[df['amp']==0, 'oid'].value_counts()

def total_count(df):
  """
  counts number of pixels in region
  # TODO: data type????/ - myles looked at this, don't think its worth fiddling with

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
  out.columns = ['q1', 'q3', 'iqr', 'min_boxplot', 'max_boxplot']

  return out

def all_metrics(df):
  """
  compute all summaries
  """
  count = total_count(df).rename('count')
  zero = zero_count(df).rename('zero_count')

  # group by OID
  grouped = df.groupby('oid', as_index=True)

  # quantiles
  quant = quantiles(grouped)

  # built in pandas functions
  # mean, median, stddev, min, max (not boxplot min/max)
  mean = grouped.mean().rename(columns={'amp':'mean'})
  median = grouped.median().rename(columns={'amp':'median'})
  std = grouped.std().rename(columns={'amp': 'std'})
  min = grouped.min().rename(columns={'amp': 'min'})
  max = grouped.max().rename(columns={'amp': 'max'})

  merge = pd.concat([mean, median, count, zero, quant, std, min, max], axis=1, join='outer')
  merge = merge.fillna(0) # fill 0s for zero_count column

  # mean and stddev on the data filtered by the boxplot min and max
  # amp = df.set_index('oid')['amp']
  # filtered = amp.gt(merge['min_boxplot']) & amp.lt(merge['max_boxplot'])
  # filtered_mean = filtered.groupby(filtered.index).mean().rename('filtered_mean')
  # filtered_std = filtered.groupby(filtered.index).std().rename('filtered_std')
  # merge2 = pd.concat([merge, filtered_mean, filtered_std], axis=1, join='outer')
  # merge2 = merge2.drop(columns={'iqr', 'min_boxplot', 'max_boxplot'})

  # join the pixels with the road-level stats
  joined = df.join(merge, on='oid', how='left')

  # filter the pixels within the min_boxplot and max_boxplot values for that rd segment
  filtered = joined.loc[joined['amp'].ge(joined['min_boxplot'])\
                      & joined['amp'].le(joined['max_boxplot'])]

  # calculate std and mean on the filtered pixels
  filtered_mean = filtered.groupby('oid')['amp'].mean().rename('filtered_mean')
  filtered_std = filtered.groupby('oid')['amp'].std().rename('filtered_std')

  # add the filtered stats, drop the unecessary boxplot stats (which are products
  # of q1 and q3 and can be regenerated later)
  merge = merge.drop(columns={'iqr', 'min_boxplot', 'max_boxplot'})
  merge = pd.concat([merge, filtered_mean, filtered_std], axis=1, join='outer')


  return merge
