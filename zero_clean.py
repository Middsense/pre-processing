def load_sar_dict():
    """
    out: dictionary with all the SAR datasets as dataframes
    """
    DATA_DIR = '/content/drive/Shared drives/Remote sensing/Summer 2020/Output Attribute Tables/21Jul_filtered/'

    sar_datasets = {}

    for path in glob(DATA_DIR + '*masked*'): # just the SAR csvs
        key = path.split('/')[-1][:-4]
        df = pd.read_csv(path, index_col='oid')
        sar_datasets[key] = df
    
    return sar_datasets

def get_zero_coords(df):
  """
  helper called by clean_zeros
  in: SAR dataframe
  out: list of (column, index) tuples where zero_count > 0
  """
  zeros = df.filter(regex='_zero_count') > 0
  zero_list = []

  for i in zeros.columns:
    for j in zeros.index:
      if zeros[i][j] == True:
        zero_list.append((i[:8],j))

  return zero_list

def clean_zeros(df):
  """
  in: SAR dataframe
  out: SAR dataframe with NaN for all cells corresponding to (oid, sar date) pairs 
       with any zero value pixels. Each of these pairs is ~8 cells, one for each 
       statistic on that oid road segment for that sar image date
  """
  df2 = df # we need to make a copy bc (I think) this function changes df in place otherwise
  zero_coords = get_zero_coords(df2)

  for i in zero_coords:
    date = i[0]
    oid = i[1]
    columns = df2.columns.str.startswith(date)
    df2.loc[oid, columns] = np.NaN

  return df2

def gen_zero_filtered_datasets(sar_datasets):
  """
  in: sar_datasets dictionary of dataframes
  out: new dictionary with NaNs in place of cells affected by zero-value pixels
  """
  sar_datasets_nozeros = {}

  for key, df in sar_datasets.items():
    sar_datasets_nozeros[key] = clean_zeros(df)

  return sar_datasets_nozeros