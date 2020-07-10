"""
merge_clean.py
Abigail Stone
7/6/2020

Some data cleaning to remove outliers and null values; preparation for classification
Inputs are speciific outputs from other middsense scripts

Inputs:
    amplitude_csv: .csv with aggregated amplitude values for each road segment and image
    road_csv: .csv of all road features (from original road segment .gdb)

** VERY PRELIMINARY **
This doesn't do much at the moment but eventually we'll want to remove outliers other
than 0 based on IQR calculations
"""

import argparse
import numpy as np
import pandas as pd
try:
    from osgeo import gdal
except ImportError:
    sys.exit('ERROR: cannot find GDAL/OGR modules')


def join_roads(roads, data):
    """
    Join aggregated amplitude values to each road segment
    """

    merged = roads.set_index(roads.OBJECTID).join(data.set_index(data.OID), how='right')

    return merged

def clean(df):
    """
    Data cleaning
    """

    # remove invalid IRI values (0 and -1)
    df = df.loc[df['NIRI_Avg'] > 0]

    # TODO:  eventually we'll want to remove outliers here (based on IQR stats, etc.)

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

    roads = pd.read_csv(args.road_csv)
    data = pd.read_csv(args.amplitude_csv)

    joined = join_roads(roads, data)
    cleaned = clean(joined)

    if args.csv:
        cleaned.to_csv(path_or_buf=args.out_path)
    else:
        cleaned.to_pickle(path=args.out_path)
