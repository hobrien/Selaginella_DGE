#!/Users/HeathOBrien/anaconda/bin/python

"""
Use blastX results to merge multiple coding sequences that are colinear with the same
reference gene into a single 'multi-exon' coding sequence.

Effectively, this means using the same geneID for multiple intervals and adjusting the 
coordinates to exclude non-homologous regions after frameshift mutations.

Minimum overlap is set to .1 by default, meaning that slight overlaps are not merged.
This is particularly important for organelles, where there can be overlapping ORFs
"""


import sys, getopt, warnings
import MySQLdb as mdb


def main(argv):
  usage = 'FixChimeras -d database_name -o overlap_cutoff [-s <species>|-c <contigID>] [-v]'
  species = ''
  contig = ''
  database = 'SelaginellaGenomics'
  overlap = .1
  try:
     opts, args = getopt.getopt(argv,"hvs:o:c:",["help", "verbose", "species=", "database=", "overlap=", "contig="])
  except getopt.GetoptError:
     print usage
     sys.exit(2)
  for opt, arg in opts:
     if opt in ("-h", "--help"):
        print usage
        sys.exit()
     elif opt in ("-s", "--species"):
        species = arg
     elif opt in ("-c", "--contig"):
        contig = arg
     elif opt in ("-d", "--database"):
        database = arg
     elif opt in ("-overlap", "--overlap"):
        overlap = float(arg)
     elif opt in ("-v", "--verbose"):
        global verbose
        verbose = 1

  con = mdb.connect(host='localhost', user='root', db=database)
  with con:
    cur = con.cursor()
    if contig:
        command = """SELECT CodingSequences.seqID from CodingSequences, CorsetGroups 
                           WHERE CorsetGroups.seqID = CodingSequences.seqID 
                           AND CorsetGroups.non_redundant = 1 
                           AND CodingSequences.seqID = %s 
                           GROUP BY CodingSequences.SeqID 
                           HAVING COUNT(CodingSequences.geneID) > 1"""
        options = (contig)
    elif species:
        command = """SELECT CodingSequences.seqID from CodingSequences, CorsetGroups 
                           WHERE CorsetGroups.seqID = CodingSequences.seqID 
                           AND CorsetGroups.non_redundant = 1 
                           AND CodingSequences.species = %s 
                           GROUP BY CodingSequences.SeqID 
                           HAVING COUNT(CodingSequences.geneID) > 1"""
        options = (species)
    else:
        command = """SELECT CodingSequences.seqID from CodingSequences, CorsetGroups 
                           WHERE CorsetGroups.seqID = CodingSequences.seqID 
                           AND CorsetGroups.non_redundant = 1 
                           GROUP BY CodingSequences.SeqID 
                           HAVING COUNT(CodingSequences.geneID) > 1"""
        options = ()
    if verbose:
       sys.stderr.write(PrintCommand(command, options))
    cur.execute(command, options)
    for row in cur.fetchall():
        seqID = row[0]
        Intervals = GetIntervals(seqID, overlap, cur)
        OldCoding = GetCoding(seqID, cur)
        
        for new_coding in UpdateCoding(Intervals, OldCoding, overlap):
            command = """INSERT INTO CodingSequences(geneID, seqID, species, start_pos, end_pos, strand, start_codon, 
                         stop_codon, notes) VALUES(%s, %s, %s, %s, %s, %s, %s, %s, %s);"""
            options = (new_coding['geneID'], new_coding['seqID'], new_coding['species'], new_coding['start_pos'],
                       new_coding['end_pos'], new_coding['strand'], new_coding['start_codon'], new_coding['stop_codon'],
                       new_coding['notes'])
            if verbose:
                sys.stderr.write(PrintCommand(command, options))
        

def UpdateCoding(intervals, old_coding, max_overlap):
    merged_coding_list = []
    deleted_coding = []
    #add start and stop codon annotations to blast intervals
    for x in range(len(intervals)):
        sorted_hits = intervals[x]
        for y in range(len(sorted_hits)):
            sorted_hits[y]['start_codon'] = 0
            sorted_hits[y]['stop_codon'] = 0
            sorted_hits[y]['species'] = old_coding[0]['species']
        merged_coding_list.append([])    

        #Determine which coding sequences overlap each blast interval
        for coding_sequence in old_coding:
            overlap = compare_ranges((coding_sequence['start_pos'], coding_sequence['end_pos']), 
                             (sorted_hits[0]['qstart'], sorted_hits[-1]['qend']))
            if overlap/float(sorted_hits[-1]['qend'] - sorted_hits[0]['qstart'] + 1) > max_overlap:
                deleted_coding.append(coding_sequence['geneID']) #coding sequence overlaps blast hit.
                if coding_sequence['strand'] == sorted_hits[0]['strand']:
                    merged_coding_list[x].append(coding_sequence)
                else:
                    pass #coding sequence on opposite strand of blast hit. spurious

        #check for IN FRAME start/stop codons at the ends of the first/last coding sequence
        if len(merged_coding_list[x]) > 0:
            sorted_hits = intervals[x]
            merged_coding = merged_coding_list[x]
            #compare frames
            if (sorted_hits[0]['qstart'] - merged_coding[0]['start_pos']) % 3 == 0:
                start = min(sorted_hits[0]['qstart'], merged_coding[0]['start_pos'])
            else:
                start = sorted_hits[0]['qstart']
            if (sorted_hits[0]['qend'] - merged_coding[0]['end_pos']) % 3 == 0:
                end = max(sorted_hits[-1]['qend'], merged_coding[-1]['end_pos'])
            else:
                end = sorted_hits[-1]['qend']    
            # examining start/stop codons'
            if merged_coding[0]['strand'] == 1:
                if merged_coding[0]['start_codon'] == 1:
                    # print 'in frame start found'
                    sorted_hits[0]['qstart'] = start
                    if merged_coding[0]['start_pos'] == start:
                        sorted_hits[0]['start_codon'] = 1
            if merged_coding[-1]['stop_codon'] == 1:
            # 'in frame end found'
                sorted_hits[-1]['qend'] = end    
            if merged_coding[-1]['end_pos'] == end:
                sorted_hits[-1]['stop_codon'] = 1
            else:
                if merged_coding[0]['stop_codon'] == 1:
                    sorted_hits[0]['qstart'] = start
                    if merged_coding[0]['start_pos'] == start:
                        sorted_hits[0]['stop_codon'] = 1
                if merged_coding[-1]['start_codon'] == 1:
                    sorted_hits[-1]['qend'] = end    
                    if merged_coding[-1]['end_pos'] == end:
                        sorted_hits[-1]['start_codon'] = 1
            if len(sorted_hits) == 1:
                sorted_hits[0]['notes'] = ''
                for coding in merged_coding:
                    if coding['stop_codon'] == 1:
                        if (coding['strand'] == 1 and coding['end_pos'] < sorted_hits[0]['qend'] and 
                            (sorted_hits[0]['qend'] -  coding['end_pos']) % 3 == 0
                            ):
                            sorted_hits[0]['notes'] = 'Non-sense'
                        elif (coding['strand'] == -1 and coding['start_pos'] > sorted_hits[0]['qstart'] and
                              (sorted_hits[0]['qstart'] -  coding['start_pos']) % 3 == 0
                              ):
                            sorted_hits[0]['notes'] = 'Non-sense'
            else:
                for y in range(len(sorted_hits)):
                    sorted_hits[y]['notes'] = 'frame-shift'
            for y in range(len(sorted_hits)):
                sorted_hits[y]['geneID'] = merged_coding[0]['geneID']
            intervals[x] = sorted_hits
    new_coding = []
    for interval in intervals:
        for sub_interval in interval:
            new_coding.append({'geneID': sub_interval['geneID'], 
                              'seqID': sub_interval['qseqid'], 
                              'species': sub_interval['species'],
                              'start_pos': sub_interval['qstart'],
                              'end_pos': sub_interval['qend'],
                              'strand': sub_interval['strand'],
                              'start_codon': sub_interval['start_codon'],
                              'stop_codon': sub_interval['stop_codon'],
                              'notes': sub_interval['notes']})
    for item in old_coding:
        if item['geneID'] not in deleted_coding:
            new_coding.append(item)
    return (new_coding)
    
def GetCoding(qseqid, cur):
    command = "SELECT * FROM CodingSequences WHERE seqID = %s"
    options = (qseqid)
    cur.execute(command, options)
    coding_sequences = []
    for (id, geneID, seqID, species, start, end, strand, start_codon, stop_codon, notes) in cur.fetchall():
        if strand == '+':
            strand = 1
        else:
            strand = -1
        coding_sequence = {'geneID': geneID,
                           'seqID': seqID,
                           'species': species,
                           'start_pos': int(start),
                           'end_pos': int(end),
                           'strand': int(strand),
                           'start_codon': int(start_codon),
                           'stop_codon': int(stop_codon),
                           'notes': notes}
        coding_sequences.append(coding_sequence)
    coding_sequences.sort(key=lambda k: int(k['start_pos']))
    return coding_sequences
    
def GetIntervals(qseqid, max_overlap, cur):
    col_names = 'qseqid sseqid qstart qend sstart send'
    intervals = []
    command = "SELECT DISTINCT sseqid FROM BlastX WHERE qseqid = %s ORDER BY evalue"
    options = (qseqid)
    if verbose:
       sys.stderr.write(PrintCommand(command, options))
    cur.execute(command, options)
    for hit in cur.fetchall():
        results = []
        command = """SELECT qseqid, sseqid, qstart, qend, sstart, send FROM BlastX WHERE qseqid = %s 
                    AND sseqid = %s ORDER BY evalue;"""
        options = (qseqid, hit[0])
        if verbose:
            sys.stderr.write(PrintCommand(command, options))
        cur.execute(command, options)
        for row in cur.fetchall():
            results.append(parse_blast_stats(col_names, row))
        sorted_hits = combine_hits(results)
        include = 1
        for interval in intervals:
            overlap = compare_ranges((interval[0]['qstart'], interval[-1]['qend']), 
                               (sorted_hits[0]['qstart'], sorted_hits[-1]['qend']))
            if overlap/float(sorted_hits[-1]['qend'] - sorted_hits[0]['qstart']) > max_overlap:
                include = 0
        if include:
            intervals.append(sorted_hits)  
    return intervals
    
def compare_ranges(range1, range2):
    if range1[0] >= range2[1] or range2[0] >= range1[1]:
        return 0
    elif range1[0] >= range2[0] and range1[1] <= range2[1]:
        return range1[1] - range1[0] + 1
    elif range2[0] >= range1[0] and range2[1] <= range1[1]:
        return range2[1] - range2[0] + 1
    elif range1[0] >= range2[0]:
        return range2[1] - range1[0] + 1
    elif range2[0] >= range1[0]:
        return range1[1] - range2[0] + 1
    else:
        print "how can it be that %s < %s and %s < %s AND %s < %s!" % (str(range2[0]), str(range1[0]), str(range1[0]),
                                                                       str(range2[1]), str(range2[0]), str(range1[1])) 
                                                                       
def combine_hits(results):
    """This will return a sorted list of query coordinate pairs with parts that overlap removed"""
    non_overlapping = []
    for result in results:
        include = 1
        for included_result in non_overlapping:
            if result['strand'] != included_result['strand']:
                #print "Skipping opposite strand hits"
                include = 0
                break
            if result['qstart'] >= included_result['qstart'] and result['qend'] <= included_result['qend']: 
                #print "result nested in included_result"
                include = 0
                break
            elif result['qstart'] <= included_result['qend'] and result['qend'] >= included_result['qend']:
                #print "start of result overlaps included_result"
                result['qstart'] = included_result['qend']+1
                result['sstart'] = included_result['send']+1
            elif result['qstart'] <= included_result['qstart'] and result['qend'] >= included_result['qstart']:
                #print "end of result overlaps included_result"
                result['qend'] = included_result['qstart']-1
                result['send'] = included_result['sstart']-1
            elif result['qstart'] <= included_result['qstart'] and result['qend'] >= included_result['qend']:
                #print "included_result nested in result (unlikely but formally possible)"
                include = 0
                break
            #Start with the most conservative possibility: any overlap of subjects is excluded
            if result['sstart'] > included_result['sstart'] and result['sstart'] < included_result['send']:
                #print "start of subject (%s) overlaps included result (%s)" % (str(result['sstart']),
                #                                                               str(included_result['sstart']) + ' - ' +
                #                                                               str(included_result['send']))
                include = 0
                break
            if result['send'] > included_result['sstart'] and result['send'] < included_result['send']:
                #print "end of subject (%s) overlaps included result (%s)" % (str(result['send']),
                #                                                             str(included_result['sstart']) + ' - ' +
                #                                                             str(included_result['send']))
                include = 0
                break       
        if include == 1:
            non_overlapping.append(result)
        non_overlapping.sort(key=lambda k: int(k['qstart']))
    return non_overlapping
    
def parse_blast_stats(column_names, row):
  """This will convert stats from blast hits to the correct numeric format and determine strand
  start and end coordinates are not reversed for negative strand results"""
  result = {}
  fields = column_names.split()
  if len(row) != len(fields):
    sys.exit("number of columns does not match specified file format. Please recheck column headers")
  for (x, field) in enumerate(fields):
    if field in ('qlen', 'slen', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'qframe', 'sstart', 'send', 'sframe', 'sstrand'):
        try:
            result[field] = int(row[x])
        except ValueError:
            sys.exit("'%s' not an integer" % str(row[x]))
    elif field in ('pident', 'evalue', 'bitscore'):
        try:
            result[field] = float(row[x])
        except ValueError:
            sys.exit("'%s' not numeric" % str(row[x]))
    else:  
        result[field] = row[x]
  result['strand'] = 1
  if result['sstart'] > result['send']:
      (result['sstart'], result['send']) = (result['send'], result['sstart'])
      result['strand'] = -1
  if result['qstart'] > result['qend']:
      (result['qstart'], result['qend']) = (result['qend'], result['qstart'])
      result['strand'] = -1
  
  return result 

def PrintCommand(command, options=()):
  if type(options) is str:
    if command.count("%s") != 1:
      sys.exit("Command requires %s options. %s supplied" % (command.count("%s"), 1))
    options = (options, "")
  elif command.count("%s") != len(options):
    sys.exit("Command requires %s options. %s supplied" % (command.count("%s"), len(options)))
  for param in options:
    command =command.replace("%s", "'" + str(param) + "'", 1)
  return command + "\n"

 
if __name__ == "__main__":
   verbose = 0
   main(sys.argv[1:])
