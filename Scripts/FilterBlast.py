#!/usr/local/bin/python

"""
A number of functions to add data to Selaginella database from csv files or from sequence
files"""


import sys, getopt, csv, warnings, re
import MySQLdb as mdb
from Bio import SeqIO
from os import path


def main(argv):
  usage = 'FilterBlast.py -i <inputfile> >outfile'
  infilename = ''
  database = 'SelaginellaGenomics'
  try:
     opts, args = getopt.getopt(argv,"hi:",["infile="])
  except getopt.GetoptError:
     print usage
     sys.exit(2)
  for opt, arg in opts:
     if opt == '-h':
        print usage
        sys.exit()
     elif opt in ("-i", "--infile"):
        infilename = arg
  con = mdb.connect('localhost', 'root', '', database);
  with con:
    cur = con.cursor()
    with open(infilename, 'rU') as f:
      reader=csv.reader(f,delimiter='\t')
      for row in reader:
        try:
          (qseqid, sseqid, pident, hitlen, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore) = row
        except ValueError:
          continue
        skip_row = 0
        for seqid in (qseqid, sseqid):
          seqid = seqid.replace('KRUS', 'KRAUS')
          if seqid.split('|')[0] in ('KRAUS', 'MOEL', 'UNC', 'WILD'):
            cur.execute("SELECT CorsetGroups.non_redundant FROM CorsetGroups, CodingSequences WHERE CorsetGroups.seqID = CodingSequences.seqID AND CodingSequences.geneID = %s", seqid.replace('|', '_'))
            non_redundant = cur.fetchone()
            if not non_redundant[0] == '1':
              skip_row = 1
        if not skip_row:
          print '\t'.join(row)
          
if __name__ == "__main__":
   main(sys.argv[1:])


