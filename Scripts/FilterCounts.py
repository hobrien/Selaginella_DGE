#!/usr/local/bin/python

"""
Filter counts by selected min and max
"""


import sys, getopt, warnings
from os import path
import pandas as pd
import numpy as np

def main(argv):
  usage = 'FilterCounts.py --min --max --perc --in inputfile --out >outfile'
  minimum = 0
  maximum = 0
  infilename = ''
  outfilename = ''
  rownames = 0
  percent = 0
  try:
     opts, args = getopt.getopt(argv,"hi:m:x:p:o:r",["in=", "out=", "min=", "max=", "perc", "rownames="])
  except getopt.GetoptError:
     print usage
     sys.exit(2)
  for opt, arg in opts:
     if opt == '-h':
        print usage
        sys.exit()
     elif opt in ("-m", "--min"):
        minimum = int(arg)
     elif opt in ("-x", "--max"):
        maximum = int(arg)
     elif opt in ("-i", "--in"):
        infilename = arg
     elif opt in ("-o", "--out"):
        outfilename = arg
     elif opt in ("-p", "--prec"):
        percent = arg
     elif opt in ("-r", "--rownames"):
        rownames = 1
  if not outfilename:
      outfilename = path.splitext(path.basename(infilename))[0] + "_filtered.txt"
  counts = pd.read_table(infilename)
  row_sums = counts.sum(1)
  if percent > 0:
      idxs = row_sums > np.percentile(row_sums, percent)
      counts = counts[idxs]
  if minimum > 0:
      idxs = row_sums > minimum
      counts = counts[idxs]
  if maximum > 0:
      idxs = row_sums < maximum
      counts = counts[idxs]
  counts.to_csv(outfilename, sep='\t')      
  
          
if __name__ == "__main__":
   main(sys.argv[1:])


