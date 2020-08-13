"""
merge_clean.py

Created: 7/6/2020
Last updated: 8/7/2020

Some data cleaning to remove outliers and null values; preparation for classification
Inputs are speciific outputs from other middsense scripts

Inputs:
    amplitude_csv: .csv with aggregated amplitude values for each road segment and image
    road_csv: .csv of all road features (from original road segment .gdb)

Mostly from the Pavement_Quality_V2 notebook
"""

import argparse
import numpy as np
import pandas as pd

from scipy import stats

from datetime import timedelta
import warnings

# dictionaries for assigning quality labels
QUALITY = {'Interstate' : {range(0, 60) : 'Excellent',
                           range(60, 100) : 'Good',
                           range(100, 140) : 'Fair',
                           range(140, 200) : 'Poor',
                           range(200, 600) : 'Very poor'},
           'Primary' : {range(0, 60) : 'Excellent',
                        range(60, 100) : 'Good',
                        range(100, 140) : 'Fair',
                        range(140, 200) : 'Poor',
                        range(200, 600) : 'Very poor'},
           'Secondary' : {range(0, 95) : 'Excellent',
                          range(95, 170) :  'Good',
                          range(170, 220) : 'Fair',
                          range(220, 280) : 'Poor',
                          range(280, 600) : 'Very poor'}
           }
QUALITY_NUM = {'Very poor' : 0, 'Poor': 1, 'Fair' : 2, 'Good' : 3, 'Excellent' : 4}


def n_closest(row, sar_dates, max_diff):
    """
    called in df.apply, returns Gaussian-weighted mean of averages within max_diff days
    """

    iri_date = row['Date_Teste']

    # convert dates to datetime format
    iri_datetime = pd.to_datetime(iri_date, format='%Y%m%d')
    sar_datetimes = pd.to_datetime(sar_dates, format='%Y%m%d')

    # finds SAR dates closest to IRI test date
    max_diff = timedelta(days=max_diff)
    closest_dates = pd.DataFrame({'diff': iri_datetime - sar_datetimes, 'sar_dates':sar_dates, 'sar_datetime':sar_datetimes})
    closest_dates = closest_dates[abs(closest_dates['diff']) < max_diff]

    if closest_dates.empty:

        return pd.Series([np.NaN, np.NaN, np.NaN])

    else:
        # closest mean and closest standard deviation
        closest_dates['closest_mean'] = closest_dates.apply(lambda row2: row[str(row2['sar_dates'])+'_mean'], axis=1)
        closest_dates['closest_std'] = closest_dates.apply(lambda row2: row[str(row2['sar_dates'])+'_std'], axis=1)

        # parameters for Gaussian
        mu, sig = 0, 30
        closest_dates['weight'] = np.exp(-1*np.power(closest_dates['diff'].dt.days - mu, 2.) / (2 * np.power(sig, 2.)))

        avg = np.average(closest_dates['closest_mean'], weights=closest_dates['weight'])
        std = np.average(closest_dates['closest_std'], weights=closest_dates['weight'])

        # Theil-Sen slopes
        ms, mi, ls, us = stats.mstats.theilslopes(closest_dates['closest_mean'], closest_dates['diff'].dt.days)
        ts_slope = ms.astype(np.float32)

        # polyfit slope
        # ignore warnings from the polyfit
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', np.RankWarning)

            m, b = np.polyfit(closest_dates['diff'].dt.days, closest_dates['closest_mean'], 1)
            slope = m.astype(np.float32)

        return pd.Series([slope, avg, std, ts_slope])


def join_roads(roads, data):
    """
    Join aggregated amplitude values to each road segment
    """
    merged = roads.join(data, how='right')

    return merged

def clean(df):
    """
    Data cleaning and road-level statistics:
    * removes invalid IRI values
    * removes road segments with pavement quality dates outside our range of interest
    * finds the SAR acquisition date closest to the IRI test date and and adds columns
      of statistics for that date
    * removes rows containing zero amplitude values
    * adds road quality labels based on IRI categories (for categorical classification labels)
    * computes the polyfit slope and Gaussian weighted average of average SAR values acquired within
      a number of days (default = 45) of IRI testing
    """
    # convert Date_Teste column to ints
    df['Date_Teste'] = df['Date_Teste'].astype(int)

    # remove invalid IRI values (0 and -1)
    df = df.loc[df['NIRI_Avg'] > 0]

    # remove road segments with pavement quality from before 2011 or after 2014
    df = df.loc[(df['Date_Teste'] < 20150000) & (df['Date_Teste'] > 20110000)]

    # add a closest date column
    SAR_DATES = np.array([c[0:8] for c in df.columns if '_filtered_mean' in c]).astype(int)
    SAR_DATES = np.sort(SAR_DATES)
    df['closest_date'] = df.apply(lambda row: np.min(np.where(SAR_DATES >= row['Date_Teste'], SAR_DATES, 99999999)), axis=1)

    # add CLOSEST mean/median
    df['closest_mean'] = df.apply(lambda row: row[str(row['closest_date'])+'_filtered_mean'], axis=1)
    df['closest_std'] = df.apply(lambda row: row[str(row['closest_date'])+'_filtered_std'], axis=1)
    df['closest_median'] = df.apply(lambda row: row[str(row['closest_date'])+'_median'], axis=1)
    df['closest_iqr'] = df.apply(lambda row: row[str(row['closest_date'])+'_q3'] - row[str(row['closest_date'])+'_q1'], axis=1)

    # mark rows containing zero amplitude values (True/False contains zeros in at least one image)
    df['zeroamp'] = df.apply(lambda row: not np.any([row[c] for c in df.columns if 'zero_count' in c]), axis=1)
    df = df.loc[df['zeroamp'] == True] # removes columns containing zero values

    # add road quality label based on IRI categories (and int version to make classification easier)
    df['quality'] = df.apply(lambda row: next((v for k, v in QUALITY[row['VDOT_Sys_I']].items() if row['NIRI_Avg'] in k), 0), axis=1)
    df['qualityINT'] = df['quality'].map(QUALITY_NUM)

    # add polyfit slope of closest_mean values within 45 days
    # add gaussian weighted average and stddev of closest_mean values within 45 days
    df[['pf_slope', 'gauss_closest_mean', 'gauss_closest_std', 'ts_slope']] = df.apply(n_closest, args=(SAR_DATES, 45), axis=1)

    return df
