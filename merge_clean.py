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
from datetime import timedelta


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

    iri_datetime = pd.to_datetime(iri_date, format='%Y%m%d')
    sar_datetimes = pd.to_datetime(sar_dates, format='%Y%m%d')

    max_diff = timedelta(days=max_diff)
    closest_dates = pd.DataFrame({'diff': iri_datetime - sar_datetimes, 'sar_dates':sar_dates, 'sar_datetime':sar_datetimes})
    closest_dates = closest_dates[abs(closest_dates['diff']) < max_diff]

    if closest_dates.empty:

        return pd.Series([np.NaN, np.NaN, np.NaN])

    else:
        closest_dates['closest_mean'] = closest_dates.apply(lambda row2: row[str(row2['sar_dates'])+'_mean'], axis=1)
        closest_dates['closest_std'] = closest_dates.apply(lambda row2: row[str(row2['sar_dates'])+'_std'], axis=1)

        # parameters for Gaussian
        mu, sig = 0, 10
        closest_dates['weight'] = np.exp(-1*np.power(closest_dates['diff'].dt.days - mu, 2.) / (2 * np.power(sig, 2.)))

        avg = np.average(closest_dates['closest_mean'], weights=closest_dates['weight'])
        std = np.average(closest_dates['closest_std'], weights=closest_dates['weight'])

        # polyfit slope
        m, b = np.polyfit(closest_dates['diff'].dt.days, closest_dates['closest_mean'], 1)
        slope = m.astype(np.float32)

        return pd.Series([slope, avg, std])


def join_roads(roads, data):
    """
    Join aggregated amplitude values to each road segment
    """
    merged = roads.join(data, how='right')

    return merged

def clean(df):
    """
    Data cleaning
    """

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
    df[['pf_slope', 'gauss_closest_mean', 'gauss_closest_std']] = df.apply(n_closest, args=(SAR_DATES, 45), axis=1)

    return df


if __name__ == "__main__":

    DESCRIPTION = "Merge road data to amplitude data and do some cleaning. Default saves to .pkl"

    parser = argparse.ArgumentParser(description=DESCRIPTION)

    parser.add_argument('road_csv',
                        help="Path to road .csv (required)")
    parser.add_argument('amplitude_csv',
                        help="Path to amplitude .csv (required)")
    parser.add_argument('out_path',
                        help="path to output .pkl or csv (required)")
    parser.add_argument('--csv', action='store_true',
                        help='Save as csv instead of .pkl')

    args = parser.parse_args()

    # Read input and index
    roads = pd.read_csv(args.road_csv, index_col='OBJECTID')
    data = pd.read_csv(args.amplitude_csv, index_col='oid')

    # join
    joined = join_roads(roads, data)

    # filtering and cleaning
    cleaned = clean(joined)

    # Output
    if args.csv:
        cleaned.to_csv(path_or_buf=args.out_path)
    else:
        cleaned.to_pickle(path=args.out_path)
