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
  out.columns = ['q1', 'q3', 'iqr', 'min', 'max']

  return out

def all_metrics(df):
  """
  compute all summaries
  """
  count = total_count(df).rename('count')
  zero = zero_count(df).rename('zero_count')

  # group by OID
  grouped = df.groupby('oid')

  # quantiles
  quant = quantiles(grouped)

  # mean and median
  mean = grouped.mean()
  median = grouped.median()

  merge = pd.concat([mean, median, count, zero, quant], axis=1, join='outer')
  merge.rename(columns='amp': 'mean', 'amp.1': 'median')
  merge.fillna(0) # fill 0s for zero_count column

  return merge
