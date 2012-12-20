#!/usr/bin/python
#
# Contig Counter
#
# Usage contigcounter.py <BLAST result file> [-debug] [-kf #[,#]] [-exclude <exclude file>]
#
# This script takes a BLAST result file and tabulates the number of top hits for each match
# across all the QUERYs in the file.
#
# Parameters:
# BLAST result file = The path and name of the input file to process.
# -debug = Debug flag for verbose reporting during processing.
# -kf #[,#] = Key field indicator.  Default is the entire BLAST description (excluding score and evalue columns).
# Each "field" is separated by spaces for purposes of counting fields.
# -exclude <exclude file> = File containing terms to exclude, if found in the RAW hit description.  Case insensitive.
#

import sys
import os
import re
import argparse
import string
import cStringIO
import math

VERSION = '1.5.1'

# CLASSES #
class HitResultObject(object):
    """This class represents all the tallied hits found for a given sequence string or key field."""
    def __init__(self,seqName,score,evalue):
        self.Name = seqName
        self.hits = 1
        self.scores = [score]
        self.evalues = [evalue]

    def addTally(self,score,evalue):
        self.hits += 1
        self.scores.append(score)
        self.evalues.append(evalue)

    def getTally(self):
        return self.hits

    def getScores(self):
        return self.scores

    def getEvalues(self):
        return self.evalues

    def getScoreAverage(self):
        return float(sum(self.scores))/len(self.scores)

    def getScoreSum(self):
        return sum(self.scores)

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

def numeric(numstr):
    """Returns an int or float value represented by the input string, or ValueError if no conversion is possible."""
    value=None
    try:
        value = int(numstr)
    except ValueError:
        value = float(numstr)

    return value

# MAIN #

print "\nContig Counter v" + VERSION  + "\n"

# Define our parameters and parase with argparse module.
argParser = argparse.ArgumentParser()

# First argument is required and must be the file path/name to the BLAST file.
argParser.add_argument("blast_file",help="The path and name of the input BLAST file to process") 

# Optional parmeters include debug flag, key field indicator (for stats aggregation purposes), and exclusion file.
argParser.add_argument("-debug","--debug",action="store_true",help="Verbose output for debugging")
argParser.add_argument("-kf","--key-fields",help="Specify which space-delimited fields from BLAST results to use for stats aggregation.  Default is the entire description.  Format is #[,#[,# .. ]] e.g. 1,3 for first and third \"words\" of description")
argParser.add_argument("-exclude",help="File containing terms (one expression per line) to exclude, if found in the RAW hit description (not just the key fields).  Case insensitive.")

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
    print "Unexpected error while opening BLAST file:", sys.exc_info()[0]
    raise


# Validate key fields parameter, if provided.
if args.key_fields and not re.match("\d+(,\d+)*",args.key_fields):
    print "Invalid key field parameter.  Format is #[,#[,# ... ]], e.g. 1 or 1,2,3"
    sys.exit(1)

# Set debug mode.
debugMode = args.debug

# Working variables.
foundHeader = False
newQuery = False
getNextHit = False
results={}
lineCount = 0
entryCount = 0
warningCount = 0
exclusionRegex = None
excludedResults={}

# Load exclusion expressions.
excludeFilename = args.exclude
if excludeFilename:
    print "Loading exclusion file " + excludeFilename + "..."

    try:
        excludeFile = open(excludeFilename,"r")
    except IOError as e:
        print "Unable to open exclusion file " + excludeFilename
        print "I/O error({0}): {1}".format(e.errno, e.strerror)
        sys.exit(1)
    except:
        print "Unexpected error while opening exclude file:", sys.exc_info()[0]
        raise


    # Build up the exclusion regular expression by writing it to a pseudo-file, for efficiency's sake.
    exclusionFileStr = cStringIO.StringIO()
    for excludeLine in excludeFile:
        exclusionFileStr.write("({})|".format(re.escape(excludeLine.rstrip(os.linesep))))

    exclusionRegex = exclusionFileStr.getvalue().rstrip('|')
    excludeFile.close()

if debugMode:
    print "Exclusion regex is {}".format(exclusionRegex)

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
            newQuery=False
            getNextHit=False
            seqString = hitMatch.group('seqstring').strip()
            hitScore = numeric(hitMatch.group('score').strip())
            hitEvalue = numeric(hitMatch.group('evalue').strip())
            # Apply exclusion filter.
            if exclusionRegex:
                excludeMatch = re.search(exclusionRegex,seqString,re.IGNORECASE)
                if excludeMatch:
                    if excludeMatch.group(0) in excludedResults:
                        excludedResults[excludeMatch.group(0)].addTally(hitScore,hitEvalue)
                    else:
                        excludedResults[excludeMatch.group(0)] = HitResultObject(excludeMatch.group(0),hitScore,hitEvalue)
                    if debugMode: print "Excluding hit " + seqString + " (matched on exclude filter " + excludeMatch.group(0) + ")"
                    continue
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
                results[seqKey].addTally(hitScore,hitEvalue)
            else:
                results[seqKey] = HitResultObject(seqKey,hitScore,hitEvalue)

# Cleanup
blastFile.close()


# Report
if not foundHeader:
    print "*** WARNING: File does not appear to be a BLAST result file.  Nothing processed.\n"
    sys.exit(1)
else:
    print "\n===== FINAL REPORT =====\n"
    print "BLAST file: " + blastFilename + "\n"
    print "Match                                                                #Hits  Score Sum  Score Avg"
    print "-------------------------------------------------------------------  -----  ---------  ---------\n"
    for result in sorted(results.items(),key=lambda x: x[1].getTally(),reverse=True):
        print "{:<67}  {:>5}  {:9.2f}  {:9.2f}".format(result[0],str(result[1].getTally()),result[1].getScoreSum(),result[1].getScoreAverage())

# Now find the total, average, and standard deviation of the number of reported hits.
totalHits = sum([ x[1].getTally() for x in results.items()])
avgHits = float(totalHits) / len(results)

if len(results) > 1:
    stddevHits = math.sqrt(float(sum([(x[1].getTally() - avgHits)**2 for x in results.items()])/(len(results) - 1)))
else:
    stddevHits = math.sqrt(float(sum([(x[1].getTally() - avgHits)**2 for x in results.items()])/len(results)))

print ""
print "Total hits: " + str(totalHits)
print "Average hits: {:0.2f}".format(avgHits)
print "Hit Standard Deviation: {:0.2f}".format(stddevHits)
print ""
print "Total reported results: " + str(len(results))
print "Total excluded results: " + str(len(excludedResults))
print "Total warnings: " + str(warningCount)

print ""

if excludedResults:
    print "\n===== EXCLUSION REPORT =====\n"
    print "Exclusion Filter Expression                                                        #Hits"
    print "---------------------------------------------------------------------------------  -----\n"
    for excludedResult in sorted(excludedResults.items(),key=lambda x: x[1].getTally(),reverse=True):
        print "{0:81}  {1}".format(excludedResult[0],str(excludedResult[1].getTally()))

print ""



