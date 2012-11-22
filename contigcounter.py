#!/usr/bin/python
#
# Contig Counter
#
# Usage contigcounter.py <BLAST result file> [-debug] [-kf #[,#]]
#
# This script takes a BLAST result file and tabulates the number of top hits for each match
# across all the QUERYs in the file.
#
# Parameters:
# BLAST result file = The path and name of the input file to process.
# -debug = Debug flag for verbose reporting during processing.
# -kf #[,#] = Key field indicator.  Default is the entire BLAST description (excluding score and evalue columns).
# Each "field" is separated by spaces for purposes of counting fields.
#

import sys
import re
import argparse

VERSION = '1.4.0'

# USER EXCEPTIONS #
class GetBlastResultKeyError(Exception):
    """Custom user exception raised when there is a failure to generate a composite (or default) key from a blast result line.
       value = Error message
       innerException = Original exception caught in the function.
    """
    def __init__(self,value,innerException = None):
        self.value =value
        self.innerException = innerException

    def __str__(self):
        return repr(self.value)

# FUNCTIONS #

def get_aggregation_key(blast_desc):
    """Takes a blast result line and returns the key value used to store in the results hash"""

    # Don't want to recreate the keyArray variable each time -- do it once and store it as an attribute of this function.
    try:
        if not getattr(get_aggregation_key,"keyArray",False):
            if args.key_fields is None:
                # If no key field parameter was specified, default is the entire BLAST entry description.
                get_aggregation_key.keyArray = None
            else:
                # Otherwise certain fields are meant to be used as the composite key for tallying.
                get_aggregation_key.keyArray = args.key_fields.split(',')

        blastKey=""
        blastDescArray = blast_desc.split()
        if get_aggregation_key.keyArray is None:
            return blast_desc
        else:
            for keyField in get_aggregation_key.keyArray:
                blastKey += " " + blastDescArray[int(keyField) - 1]
            return blastKey.strip()
    except Exception as e:
        #print "*** WARNING: Failed to get key field for entry (using full line): " + blast_desc
        raise GetBlastResultKeyError(blast_desc,e)

# MAIN #

print "\nContig Counter v" + VERSION  + "\n"

# Define our parameters and parase with argparse module.
argParser = argparse.ArgumentParser()

# First argument is required and must be the file path/name to the BLAST file.
argParser.add_argument("blast_file",help="The path and name of the input BLAST file to process") 

# Optional parmeters include debug flag and key field indicator (for stats aggregation purposes)
argParser.add_argument("-debug","--debug",action="store_true",help="Verbose output for debugging")
argParser.add_argument("-kf","--key-fields",help="Specify which space-delimited fields from BLAST results to use for stats aggregation.  Default is the entire description.  Format is #[,#[,# .. ]] e.g. 1,3 for first and third \"words\" of description")

args = argParser.parse_args()

if args.debug:
    print args

blastFilename = args.blast_file

try:
    blastFile = open(blastFilename,"r")
except IOError as e:
    print "Unable to open file " + blastFilename
    print "I/O error({0}): {1}".format(e.errno, e.strerror)
    sys.exit(1)
except:
    print "Unexpected error:", sys.exc_info()[0]
    raise


# Validate key fields parameter, if provided.
if args.key_fields and not re.match("\d+(,\d+)*",args.key_fields):
    print "Invalid key field parameter.  Format is #[,#[,# ... ]], e.g. 1 or 1,2,3"
    sys.exit(1)

# Set debug mode.
debugMode = args.debug

foundHeader = False
newQuery = False
getNextHit = False
results=dict()
lineCount = 0
entryCount = 0
warningCount = 0

print "Processing " + blastFilename + "..."

# Processing the file.
# Verify that it is a BLAST result file (has a BLASTN header line)
# Then grab the first match line a QUERY section. Build up a dictionary of counts for each hit.
# TODO: Should consider refactoring this to use the State pattern.
for line in blastFile:
    lineCount += 1
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
        # Basically want to grab all but the last two space-delimited fields (i.e. the score and e-value) in the line.
        hitMatch = re.match("^\s*(?P<seqstring>.+)\s+(?P<score>\S+)\s+(?P<evalue>\S+)$",line)
        if (hitMatch):
            if debugMode: print "Match hit for seqstring: " + hitMatch.group('seqstring') 
            entryCount += 1
            newQuery=False
            getNextHit=False
            seqString = hitMatch.group('seqstring').strip()
            try:
                seqKey = get_aggregation_key(seqString)
            except GetBlastResultKeyError as e:
                seqKey = seqString
                warningCount += 1
                print "*** WARNING: Failed to get key field for entry on line {0}: {1}".format(lineCount, line)
                print "(Will use full description as key)"
                if debugMode:
                    print "Matched seqString is " + seqString
                    print "Exception details:"
                    print "Error value: " + str(e.value)
                    print "From original exception: ", e .innerException
                    print str(e.innerException)

            if debugMode: print "Tally added for " + seqKey
            if seqKey in results:
                results[seqKey] = results[seqKey] + 1
            else:
                results[seqKey] = 1

if not foundHeader:
    print "*** WARNING: File does not appear to be a BLAST result file.  Nothing processed.\n"
else:
    print "\n===== FINAL REPORT =====\n"
    print "BLAST file: " + blastFilename + "\n"
    print "Match                                                                              #Hits"
    print "---------------------------------------------------------------------------------  -----\n"
    for result in sorted(results.items(),key=lambda x: x[1],reverse=True):
        print "{0:81}  {1}".format(result[0],str(result[1]))

print ""

print "Total results tallied: " + str(entryCount)
print "Total warnings: " + str(warningCount)

print ""

# Cleanup
blastFile.close()

