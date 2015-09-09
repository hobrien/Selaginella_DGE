import sys, os
import numpy as np
import pandas as pd
from dgeclust import CountData, SimulationManager

from dgeclust.models import NBinomModel
from dgeclust import compare_groups
 
os.chdir('/Users/HeathOBrien/Bioinformatics/Selaginella_DGE/DGEclust')

table = sys.argv[1]
species = os.path.split(table)[-1].split('count')[0]
if species == 'MOEL':
    group_names = ["%s%i" % t for t in zip([species] * 3, range(1,4))]
else:
    group_names = ["%s%i" % t for t in zip([species] * 4, range(1,5))]
    
counts = pd.read_table(table, index_col=0)

mgr = SimulationManager()
data = CountData(counts, groups=group_names)
mdl = NBinomModel.load(species + '_DGEclust')

comp1 = sys.argv[2]
comp2 = sys.argv[3]
res, nsamples = compare_groups(data, mdl, group1=comp1, group2=comp2)

outfile = comp1 + comp2[-1] + '.txt'
res.sort('FDR').to_csv(outfile, sep='\t', index_label="Cluster")
