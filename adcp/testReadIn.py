import numpy as np
import pandas as pd


# if pandas is available, best to use it.
adcp = pd.read_csv('dngrid_adcp_2012.txt')

# or use numpy cause pandas isn't on placentia
adcp = np.load('dngrid_adcp_2012.pkl')

# these are equivalent.
