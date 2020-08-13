"""
pixel2road.py
Created: 7/16/2020
Last Updated: 8/10/2020

functions to compute summary statistics on amplitude data
"""
import numpy as np
import pandas as pd
import pixelwise

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

def all_metrics(df, LOG_FILTER, MEAN_NORMALIZE):
    """
    compute all summaries
    """
    # count zeros before applying log or mean normalizing
    zerocount = df.eq(0).groupby('oid').sum().add_suffix('_zerocount')
    count = df.groupby('oid').count().add_suffix('_count')

    # zeros cause problems for logs (-inf) and potentially other steps
    df = df.replace(0, np.nan)

    if LOG_FILTER:
        df = np.log(df)

    if MEAN_NORMALIZE:
        df = pixelwise.mean_normalize(df)

    # group by OID
    grouped = df.groupby('oid', as_index=True)

    # mean, median, stddev, min, max (not boxplot min/max)
    mean = grouped.mean().add_suffix('_mean')
    median = grouped.median().add_suffix('_median')
    std = grouped.std().add_suffix('_std')
    min = grouped.min().add_suffix('_min')
    max = grouped.max().add_suffix('_max')

    # combine the multiple statistics series --> a dataframe
    merge = pd.concat([mean, median, std, min, max], axis=1, join='outer')

    # quantiles
    quant = quantiles(grouped)

    # pixel filtering based on boxplot min/max
    maxs = quant.filter(like='bmax')
    mins = quant.filter(like='bmin')
    maxs.columns = (maxs.columns.str.rstrip('_bmax'))
    mins.columns = (mins.columns.str.rstrip('_bmin'))

    # create mask for values within boxplot min and max
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

    if LOG_FILTER:
        merge = np.exp(merge)

    merge = pd.concat([merge, count, zerocount], axis=1, join='outer').fillna(0)

    # sort columns
    merge = merge.sort_index(axis=1)

    return merge
