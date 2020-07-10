#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Gretchen Doyle
7/9/20

PURPOSE: Given a csv file with the road OBJECTID and corresponding amplitude values,
calculates the number of zeros within each unqiue OID, as well as the percentage of zeros.
Next, summarizes the sum, mean and standard deviation for each OID before and after cleaning
the data so that all zero amplitude values are removed.

NOTE: you can use a test csv file called amp_test.csv which is uploaded

"""

import pandas as pd
import numpy as np

def summarize(df):
    """

    Parameters
    ----------
    df : dataframe containing columns 'OID' and 'amp'

    Returns
    -------
    grouped: dataframe created by grouping 'df' by 'OID' and calculating the average, 
    mean, and standard deviation of the corresponding amplitude values

    """
    grouped = df.groupby('OID')['amp'].agg([np.sum, np.mean, np.std])
    return (grouped)


# Create dataframe from csv file
df = pd.read_csv("~/Desktop/amp_test.csv")

# Calculate the number of zero amplitude values for each OID
zero_count = df.groupby('OID')['amp'].apply(lambda x: x[x==0].count())

# Calculate the total amplitude values for each OID
amp_total = df.groupby('OID')['amp'].count()

# Calculate the percentage of zero amp values for each OID
for unique_value in 'OID':
    percent_zeros = zero_count/amp_total
    
#Summarize original dataset
print("Summary before data cleaning:\n")
print(summarize(df))


#Create a new dataframe with the following columns: unique OID, number of zeros, percentage of zeros
new_df = pd.concat([amp_total, zero_count, percent_zeros], axis=1)
new_df.columns = ['amp_total', 'zero_count', 'percent_zeros']


#Remove zero amp values from df
print("Summary after data cleaning:\n")
df = df.loc[df.amp!=0]


#Summarize dataset with no zero amp values
summarize(df)

