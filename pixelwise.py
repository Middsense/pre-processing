"""
pixelwise.py
08.06.2020

goal: get sparse matrices into an easily operable format so that we can perform
pixelwise operations across the 67-image stack

in: SAR and road OID sparse matrices

out: DataFrames with each row corresponding to a pixel. columns are: 
    - the row and col identifying the pixel location in the raster
    - the OID if a pixel is within a road segment in each of the road datasets
    - the pixel's amplitude value for each image
"""

import pandas as pd
import numpy as np
from scipy import sparse
from glob import glob

def gen_img_stack(img_dir):
    
    ref_sparse = sparse.load_npz(img_dir + '/20110829.npz') #arbitrary sp matrix for row and col
    
    all_images_df = pd.DataFrame({'row': ref_sparse.row,                                  
                                  'col': ref_sparse.col})

    for path in glob(img_dir + '*'):
        name = path.split('/')[-1][:-4]
        sp = sparse.load_npz(path)
        if np.all(sp.col == all_images_df['col'].values) and np.all(sp.row == all_images_df['row'].values):    
            all_images_df[name] = sp.data
        else:
            print('error: problem with sparse matrix pixel alignment')
            
    # we have leftover -1 values from the sparse matrices
    # this was a way to hold onto zero value pixels, which sparse matrices 
    # treat as nodata. here, we cast them back to 0
    all_images_df = all_images_df.replace(-1, 0)
            
    return all_images_df


def merge_img_rd_stack(img_stack, rd_dir):
    
    for path in glob(rd_dir + '*'):
        
        # handle column naming
        name = path.split('/')[-1][:-4].split('_')
        tag = 'oid_'
        tag += name[2]
        if 'Buf12' in name:
            tag += '_buffered'
        else: 
            tag += '_centerline'
        if 'landcovermasked' in name:
            tag += '_masked'
        else:
            tag += '_unmasked'
        
        # load the rd sparse matrices and convert to df
        rd_sp = sparse.load_npz(path)
        
        rd_df = pd.DataFrame({
            'row': rd_sp.row,
            'col': rd_sp.col,
            tag: rd_sp.data 
            })    
        
        # add each road df's OID column to img_stack 
        # TODO probably slow/inefficient merging
        img_stack = pd.merge(left=rd_df, right=img_stack, how='right', on=['row', 'col'])
        
    return img_stack