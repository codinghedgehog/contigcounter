#!/usr/bin/python
#
# Contig Counter
#
# Usage contigcounter.py <BLAST result file> [-debug]
#
# This script takes a BLAST result file and tabulates the number of top hits for each match
# across all the QUERYs in the file.
#

import sys
import re

VERSION = '1.1'

print "\nContig Counter v" + VERSION 

if (len(sys.argv) < 2):
    print "Missing filename.\n"
    print "Usage: " + sys.argv[0] + " <BLAST result file> [-debug]"
    sys.exit(1)

try:
    blastFile = open(sys.argv[1],"r")
except IOError as e:
    print "Unable to open file " + sys.argv[1]
    print "I/O error({0}): {1}".format(e.errno, e.strerror)
except:
    print "Unexpected error:", sys.exc_info()[0]
    raise

# TODO: Should use argparse here.
debugMode = False

if len(sys.argv) == 3 and sys.argv[2] == "-debug":
    debugMode=True

foundHeader = False
newQuery = False
getNextHit = False
results=dict()

# Processing the file.
# Verify that it is a BLAST result file (has a BLASTN header line)
# Then grab the first match line a QUERY section. Build up a dictionary of counts for each hit.
# TODO: Really should refactor this to use the State pattern.
for line in blastFile:
    # Skip empty lines.
    line = line.rstrip()
    if re.match("^\s*$",line):
        continue
    elif not foundHeader and re.match("^\s*BLASTN",line):
        if debugMode: print "Found file header: " + line + " -- OK"
        foundHeader = True
    elif foundHeader and not newQuery:
        queryMatch = re.match("^\s*Query= (?P<queryName>.+)$",line)
        if queryMatch:
            newQuery = True
            if debugMode: print "At query " + queryMatch.group('queryName') 
    elif foundHeader and newQuery and not getNextHit:
        if line.find("Sequences producing significant alignments:") == 0:
            getNextHit = True
            if debugMode: print "Tabulating top hit for query " + queryMatch.group('queryName') 
    elif foundHeader and newQuery and getNextHit:
        #print "Processing " + line
        hitMatch = re.match("^\s*(?P<seqstring>.+\|.+?)\s{2,}(?P<score>.+)\s{2,}(?P<evalue>.+)$",line)
        if (hitMatch):
            newQuery=False
            getNextHit=False
            seqString = hitMatch.group('seqstring').strip()
            if debugMode: print "Tally added for " + seqString 
            if seqString in results:
                results[seqString] = results[seqString] + 1
            else:
                results[seqString] = 1

if not foundHeader:
    print "*** WARNING: File does not appear to be a BLAST result file.  Nothing processed.\n"
else:
    print "\n===== FINAL REPORT =====\n"
    print "Match                                                                              #Hits"
    print "---------------------------------------------------------------------------------  -----\n"
    for result in sorted(results.items(),key=lambda x: x[1],reverse=True):
        print "{0:81}  {1}".format(result[0],str(result[1]))

print ""

# Cleanup
blastFile.close()

