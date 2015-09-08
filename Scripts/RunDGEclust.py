import sys, os
import numpy as np
import pandas as pd
import matplotlib.pylab as pl
from dgeclust import CountData, SimulationManager
from dgeclust.models import NBinomModel
 
os.chdir('/Users/HeathOBrien/Bioinformatics/Selaginella_DGE/DGEclust')

table = sys.argv[1]
species = os.path.split(table)[-1].split('count')[0]
if species == 'MOEL':
    group_names = ["%s%i" % t for t in zip([species] * 3, range(1,4))]
else:
    group_names = ["%s%i" % t for t in zip([species] * 4, range(1,5))]
    
counts = pd.read_table(table)
row_sums = counts.sum(1)
idxs = row_sums > np.percentile(row_sums, 40)   # identifies the upper 60% of the data
counts_filt = counts[idxs]  

mgr = SimulationManager()
data = CountData(counts_filt, groups=group_names)
mdl = NBinomModel(data, outdir=species + '_DGEclust')

mgr.new(data, mdl)

