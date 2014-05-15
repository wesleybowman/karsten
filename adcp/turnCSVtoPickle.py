import pandas as pd
from glob import glob

files = glob('*.csv')

for i in files:
    data = pd.read_csv(i, index_col=0)
    data.index = pd.to_datetime(data.index)

    name = '{0}.pkl'.format(i.strip('.csv'))
    data.to_pickle(name)
