#!/Users/HeathOBrien/anaconda/bin/python

"""
A number of functions to add data to Selaginella database from csv files or from sequence
files"""


import sys, getopt, csv, warnings, re
import MySQLdb as mdb
from Bio import SeqIO
from os import path


def main(argv):
  usage = 'AddData -f <function> -i <inputfile> -s <sequencefile> -d database_name'
  function = ''
  infilename = ''
  seqfilename = ''
  database = 'SelaginellaGenomics'
  name = ''
  try:
     opts, args = getopt.getopt(argv,"hvf:i:s:d:n:",["function", "infile=", "seqfile=", "database=", "name="])
  except getopt.GetoptError:
     print usage
     sys.exit(2)
  for opt, arg in opts:
     if opt == '-h':
        print usage
        sys.exit()
     elif opt in ("-f", "--function"):
        function = arg
     elif opt in ("-i", "--infile"):
        infilename = arg
     elif opt in ("-s", "--seqfile"):
        seqfilename = arg
     elif opt in ("-d", "--database"):
        database = arg
     elif opt in ("-n", "--name"):
        name = arg
     elif opt in ("-v", "--verbose"):
        global verbose
        verbose = 1
        
  if function == 'ref_coding':
    infilename = 1  #don't need a infile for this function. This prevents an error in the following lines
  con = mdb.connect('localhost', 'root', '', database)
  with con:
    cur = con.cursor()
    if function == 'TAIR':
      add_TAIR(cur, infilename)
    if function == 'sequences':
      add_seqs(cur, infilename)
    if function == 'orthologs':
    #  add_ortholog_info(cur, infilename) 
      add_orthologs(cur, infilename)
    if function == 'coding':
      add_coding(cur, infilename)
    if function == 'ref_coding':
      ref_coding(cur)
    if function == 'nr':
      non_redundant(cur, infilename)
    if function == 'counts':
      add_counts(cur, infilename, name)
    if function == 'blast':
      add_blast(cur, infilename, 'Blast')
    if function == 'blastx':
      add_blast(cur, infilename, 'BlastX')
    if function == 'corset_clusters':
      add_corset_clusters(cur, infilename, name)
    if function == 'corset_counts':
      add_corset_counts(cur, infilename, name)
    if function == 'corset_nr':
      corset_nr(cur, name)
    if function == 'de':
      differential_expression(cur, infilename, name)
    if function == 'dep':
      de_pvalues(cur, infilename, name)
    if function == 'ath':
      add_Ath(cur, infilename)

def PrintCommand(command, options=()):
  if type(options) is str:
    if command.count("%s") != 1:
      sys.exit("Command requires %s options. %s supplied" % (command.count("%s"), 1))
    options = (options, "")
  elif command.count("%s") != len(options):
    sys.exit("Command requires %s options. %s supplied" % (command.count("%s"), len(options)))
  for param in options:
    command =command.replace("%s", "'" + param + "'", 1)
  return command + "\n"
  
def ExecuteCommand(cur, command, options=()):
    if verbose:
        sys.stderr.write(PrintCommand(command, options))
    try:
        cur.execute(command, options)
    except mdb.IntegrityError, e:
        warnings.warn("%s" % e)
        pass
  

def add_Ath(cur, infilename):
  with open(infilename, 'rU') as f:
    reader=csv.reader(f,delimiter='\t')
    headers = reader.next()
    for row in reader:
      try:
        (EnsemblGeneID, EnsemblTranscriptID, PercIdentity, HomologyType, AthEnsemblGeneID) = row
      except ValueError:
        warnings.warn("wrong number of columns")
        continue
      if 'EFJ' in EnsemblTranscriptID:
         EnsemblTranscriptID = EnsemblTranscriptID.replace('EFJ', 'EFJ_')
      else:
         EnsemblTranscriptID = 'EFJ_' + EnsemblTranscriptID   
      if not AthEnsemblGeneID:
         continue
      command = "INSERT INTO AthHomologs(seqID, athID, percen_identity, homology_type) VALUES(%s, %s, %s, %s)"
      options = (EnsemblTranscriptID, AthEnsemblGeneID, PercIdentity, HomologyType)
      ExecuteCommand(cur, command, options)

def de_pvalues(cur, infilename, name):
  command = "INSERT INTO Expression(clusterID, SpeciesID, sample1, sample2, Posteriors, FDR) VALUES(%s, %s, %s, %s, %s, %s)"
  with open(infilename, 'rU') as f:
    reader=csv.reader(f,delimiter='\t')
    headers = reader.next()
    basename = path.splitext(path.basename(infilename))[0]
    sample1 = basename[:-1]
    sample2 = sample1[:-1] + basename[-1]
    for row in reader:
      try:
        (clusterID, Posteriors, FDR) = row
      except ValueError:
        continue
      options = (clusterID, sample1[:-1], sample1, sample2, Posteriors, FDR)
      ExecuteCommand(cur, command, options)

def differential_expression(cur, infilename, name):
  command = "UPDATE Expression SET sample1Expr = %s, sample2Expr = %s, avgLogExpr = %s, rLogFC = %s WHERE clusterID = %s AND sample1 = %s AND sample2 = %s"
  with open(infilename, 'rU') as f:
    reader=csv.reader(f,delimiter='\t')
    headers = reader.next()
    basename = path.splitext(path.basename(infilename))[0]
    sample1 = basename[:-1]
    sample2 = sample1[:-1] + basename[-1]
    for row in reader:
      try:
        (clusterID, sample1_expression, sample2_expression, avgLogExpr, rLogFC) = row
      except ValueError:
        continue
      options = (sample1_expression, sample2_expression, avgLogExpr, rLogFC, clusterID, sample1, sample2)
      ExecuteCommand(cur, command, options)

def corset_nr(cur, speciesID):
    print "SELECT CorsetGroups.seqID, Max(Blast.hitlen) FROM Blast, CorsetGroups, CodingSequences WHERE CorsetGroups.seqID = CodingSequences.seqID AND CodingSequences.geneID = Blast.qseqid AND Blast.qseqid != Blast.sseqid AND CorsetGroups.speciesID = %s GROUP BY CorsetGroups.clusterID" % speciesID
    cur.execute("SELECT CorsetGroups.seqID, Max(Blast.hitlen) FROM Blast, CorsetGroups, CodingSequences WHERE CorsetGroups.seqID = CodingSequences.seqID AND CodingSequences.geneID = Blast.qseqid AND Blast.qseqid != Blast.sseqid AND CorsetGroups.speciesID = %s GROUP BY CorsetGroups.clusterID", speciesID)
    for (seqID, cluster_size) in cur.fetchall():
      try:
        print seqID, cluster_size
        cur.execute('UPDATE CorsetGroups SET non_redundant = 1 WHERE seqID =  %s' , (seqID))
      except mdb.IntegrityError, e:
        #warnings.warn("%s" % e)
        pass


def add_corset_clusters(cur, infilename, name):
  with open(infilename, 'rU') as f:
    reader=csv.reader(f,delimiter='\t')
    for row in reader:
      try:
        (seqID, clusterID) = row
      except ValueError:
        continue
      cur.execute("INSERT INTO CorsetGroups(speciesID, seqID, clusterID) VALUES(%s, %s, %s)", (name, seqID, clusterID))

def add_corset_counts(cur, infilename, name):
  with open(infilename, 'rU') as f:
    reader=csv.reader(f,delimiter='\t')
    for row in reader:
      try:
        (clusterID, leaf1, leaf1b, leaf2, leaf2b, leaf3, leaf3b, leaf4, leaf4b) = row
        cur.execute("INSERT INTO CorsetCounts(speciesID, clusterID, leaf1, leaf1b, leaf2, leaf2b, leaf3, leaf3b, leaf4, leaf4b) VALUES(%s, %s, %s, %s, %s, %s, %s, %s, %s, %s)", (name, clusterID, leaf1, leaf1b, leaf2, leaf2b, leaf3, leaf3b, leaf4, leaf4b))
      except ValueError:  #In the case of MOEL, where there are only 7 columns
        try:
          (clusterID, leaf1, leaf1b, leaf2, leaf2b, leaf3, leaf3b) = row
          cur.execute("INSERT INTO CorsetCounts(speciesID, clusterID, leaf1, leaf1b, leaf2, leaf2b, leaf3, leaf3b) VALUES(%s, %s, %s, %s, %s, %s, %s, %s)", (name, clusterID, leaf1, leaf1b, leaf2, leaf2b, leaf3, leaf3b))
        except ValueError:
          continue

def add_blast(cur, infilename, table):
  if table == 'BlastX':
      command = "INSERT INTO BlastX(qseqid, sseqid, pident, hitlen, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore) VALUES(%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)"
  else:
      command = "INSERT INTO Blast(qseqid, sseqid, pident, hitlen, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore) VALUES(%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)"
  with open(infilename, 'rU') as f:
    reader=csv.reader(f,delimiter='\t')
    for row in reader:
      try:
        (qseqid, sseqid, pident, hitlen, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore) = row
      except ValueError:
        continue
      qseqid = qseqid.replace('|', '_')
      sseqid = sseqid.replace('|', '_')
      qseqid = qseqid.replace('KRUS', 'KRAUS')
      sseqid = sseqid.replace('KRUS', 'KRAUS')
      options = (qseqid, sseqid, pident, hitlen, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore)
      if verbose:
         sys.stderr.write(PrintCommand(command, options))
      try:
        cur.execute(command, options)
      except mdb.IntegrityError, e:
        warnings.warn("%s" % e)
        pass
        

def add_counts(cur, infilename, name):
  infile = open(infilename, 'r')
  for line in infile.readlines():
    line = line.strip()
    (geneID, count) = line.split('\t')
    if geneID not in ('no_feature', 'ambiguous', 'too_low_aQual', 'not_aligned', 'alignment_not_unique'):
      cur.execute("SELECT * FROM Counts WHERE geneID = %s", (geneID))
      results = cur.fetchall()
      if len(results) == 1:
        #print "UPDATE Counts SET %s = '%s' WHERE geneID = '%s'" % (name, count, geneID)
        cur.execute("UPDATE counts SET %s = '%s' WHERE geneID = '%s'" % (name, count, geneID))
      elif len(results) == 0:
        #print "INSERT INTO Counts(geneID, %s) VALUES('%s', '%s')" % (name, geneID, count)
        cur.execute("INSERT INTO Counts(geneID, %s) VALUES('%s', '%s')" % (name, geneID, count))
      else:
        sys.exit("%s entries in DB for %s" % (len(results), geneID))
        
def add_aligned(cur, seqfilename):
  for sequence in SeqIO.parse(seqfilename, "fasta"):
    try:
      #print 'INSERT INTO Sequences(seqID, locus, accessionID, sequence)  VALUES(%s, %s, %s, %s)' , (seqID, locus, accessionID, sequence)
      cur.execute('UPDATE Sequences SET aligned = %s WHERE seqID = %s' , (sequence.seq.upper(), sequence.id))
    except mdb.IntegrityError, e:
      #warnings.warn("%s" % e)
      pass


"""add genbank sequences obtained using something like:
"blastdbcmd -db nt -entry_batch genbankIDs.txt >sequences.fa"

infile should have columns as follows:

seqID, locus, accessionID
"""
def genbank(cur, infilename, seqfilename):
  seqdict = SeqIO.to_dict(SeqIO.parse(seqfilename, "fasta"), key_function=get_accession)
  with open(infilename, 'rU') as f:
    reader=csv.reader(f,delimiter='\t')
    next(reader, None) 
    for row in reader:
      try:
        (accessionID, locus, seqID) = row
      except ValueError:
        warnings.warn("row length %s does not match the expected number of columns (3). Double-check delimiter" % len(row))
        continue
      sequence = seqdict[seqID].seq
      
      try:
        #print 'INSERT INTO Sequences(seqID, locus, accessionID, sequence)  VALUES(%s, %s, %s, %s)' , (seqID, locus, accessionID, sequence)
        cur.execute('INSERT INTO Sequences(seqID, locus, accessionID, sequence)  VALUES(%s, %s, %s, %s)' , (seqID, locus, accessionID, sequence))
      except mdb.IntegrityError, e:
        #warnings.warn("%s" % e)
        pass

def non_redundant(cur, infilename):
  for seq_record in SeqIO.parse(infilename, "fasta"):
    id = seq_record.id
    id = id.replace('|', '_')
    try:
      #print 'INSERT INTO Sequences(seqID, sequence, species, length)  VALUES(%s, %s, %s, %s)' % (id,seq,species,length)
      cur.execute('UPDATE OrthoGroups SET non_redundant = 1 WHERE geneID =  %s' , (id))
    except mdb.IntegrityError, e:
      #warnings.warn("%s" % e)
      pass


def add_orthologs(cur, infilename):
  infile = open(infilename, 'r')
  regex = re.compile('\d+')
  for line in infile.readlines():
    genes = line.split()
    group = genes.pop(0)
    group = regex.findall(group)[-1]
    for gene in genes:
      gene = gene.replace('|','_')
      gene = gene.replace('KRUS', 'KRAUS')
      try:
        #print "INSERT INTO OrthoGroups (geneID, orthoID) VALUES ('%s', '%s')" % (gene, group)
        cur.execute("INSERT INTO OrthoGroups (geneID, orthoID) VALUES (%s, %s)" , (gene, group))
      except mdb.IntegrityError, e:
        warnings.warn("%s" % e)
        sys.exit()
          
def add_TAIR(cur, infilename):
  infile = open(infilename, 'r')
  for line in infile.readlines():
    line = line.strip()
    fields = line.split('\t')
    locus = fields[0]
    if locus == 'Locus Identifier': #header row
      continue
    subfields = fields[2].split(';')
    description = subfields[0]
    try:
      symbol = fields[4]
    except IndexError:
      symbol = ''
    print "INSERT INTO tair (locus, description, symbol) VALUES ('%s', '%s', '%s')" % (locus, description, symbol)
    try:
      cur.execute("INSERT INTO tair (locus, description, symbol) VALUES (%s, %s, %s)", (locus, description, symbol))
    except mdb.IntegrityError, e:
      continue

def add_seqs(cur, infilename):
  for seq_record in SeqIO.parse(infilename, "fasta"):
    seq = seq_record.seq
    length = len(seq)
    #id = seq_record.id
    if 'comp' in seq_record.id:
      species = seq_record.id.split("comp")[0]
      #print species
      id = seq_record.id
    elif 'scaffold-' in seq_record.id:
      id = "_".join(seq_record.id.split("-")[1:3])
      species = seq_record.id.split("-")[1]
    elif 'EFJ' in seq_record.id:
      species = 'EFJ'
      id = seq_record.id.split(' ')[0].replace('EFJ', 'EFJ_')
    elif 'ADH' in seq_record.id:
      species = 'EFJ'
      id = 'EFJ_' + seq_record.id.split(' ')[0]
    try:
      #print 'INSERT INTO Sequences(seqID, sequence, species, length)  VALUES(%s, %s, %s, %s)' % (id,seq,species,length)
      cur.execute('INSERT INTO Sequences(seqID, sequence, species, length)  VALUES(%s, %s, %s, %s)' , (id,seq,species,length))
    except mdb.IntegrityError, e:
      #warnings.warn("%s" % e)
      pass

def add_coding(cur, infilename):
  """This will loop through a bed file of coding coordinates from transdecoder and add 
  info to the database. It is pretty straightforward EXCEPT that there's a bunch of info
  encoded in the name field that needs to be extracted"""
  with open(infilename, 'rU') as f:
    reader=csv.reader(f,delimiter='\t')
    for row in reader:
      try:
        (seqID, start, end, name, score, strand, thickStart, thickEnd, itemRgb, blockCount, blockSizes, blockStarts) = row
      except ValueError:
        continue
      species = name[3:name.find('comp')]
      geneID = species + '_' + name.split('|')[-1].split(':')[0].split('_')[0].replace('m.','')
      thickStart = int(thickStart) + 1
      if 'internal_len' in name.split('|')[-1].split(':')[1]:
        start_codon = 0
        stop_codon = 0
      elif '5prime_partial_len' in name.split('|')[-1].split(':')[1]:
        start_codon = 0
        stop_codon = 1
      elif '3prime_partial_len' in name.split('|')[-1].split(':')[1]:
        start_codon = 1
        stop_codon = 0
      elif 'complete_len' in name.split('|')[-1].split(':')[1]:
        start_codon = 1
        stop_codon = 1
      else:
        warnings.warn("could not parse %s" % name.split('|')[-1].split(':')[1])
      try:  
        #print 'INSERT INTO CodingSequences(geneID, seqID, species, start, end, strand, start_codon, stop_codon)  VALUES(%s, %s, %s, %s, %s, %s, %s, %s)' % (geneID, seqID, species, thickStart, thickEnd, strand, start_codon, stop_codon)
        cur.execute('INSERT INTO CodingSequences(geneID, seqID, species, start, end, strand, start_codon, stop_codon)  VALUES(%s, %s, %s, %s, %s, %s, %s, %s)', (geneID, seqID, species, thickStart, thickEnd, strand, start_codon, stop_codon))
      except mdb.IntegrityError, e:
        warnings.warn("%s" % e)
        pass
        
def ref_coding(cur):
  """This will loop through all non-BLUELEAF sequences (which are all CDSs) and add their
  info to the CodingSequences table (the info is the same as the Sequences table, but 
  keeping everything consistent will make life easier downstream):
  CodingSequences.geneID = Sequences.seqID
  CodingSequences.seqID = Sequences.seqID
  CodingSequences.start = 1
  CodingSequences.end = Sequences.length
  CodingSequences.strand = +
  CodingSequences.start_codon = NULL
  CodingSequences.stop_codon = NULL"""
  cur.execute("SELECT Sequences.seqid, Sequences.length, Sequences.species FROM Sequences, Species WHERE Sequences.species = Species.speciesID AND Species.source NOT LIKE 'BLUELEAF'")
  for (seqid, len, species) in cur.fetchall():
    try:
      cur.execute('INSERT INTO CodingSequences(geneID, seqID, species, start, end, strand)  VALUES(%s, %s, %s, %s, %s, %s)', (seqid, seqid, species, 1, len, '+'))
    except mdb.IntegrityError, e:
      warnings.warn("%s" % e)
      pass

"""This is an older function that I'm replacing with the orthomcl version"""        
def add_ortholog_info(cur, infilename):
  infile = open(infilename, 'r')
  for line in infile.readlines():
    line = line.strip()
    EFJ = line.split('\t')[1]
    if EFJ == 'Ensembl Transcript ID': #header row
      continue
    try:
      ATH = line.split('\t')[4]
    except IndexError:
      continue
    cur.execute("SELECT clusternum FROM sequences WHERE seqid = %s", (EFJ))
    try:
      clusternum = cur.fetchone()[0]
    except:
      continue
    cur.execute("UPDATE sequences SET clusternum = %s WHERE seqid = %s", (clusternum, ATH))
    print "UPDATE sequences SET clusternum = %s WHERE seqid = '%s'" % (clusternum, ATH)

def get_accession(record):
    """"Given a SeqRecord, return the accession number as a string.
  
    e.g. "gi|2765613|gb|Z78488.1|PTZ78488" -> "Z78488.1"
    """
    parts = record.id.split("|")
    assert len(parts) == 5 and parts[0] == "gi" and parts[2] == "gb"
    return parts[3].split('.')[0]
    
if __name__ == "__main__":
   verbose = 0
   main(sys.argv[1:])


