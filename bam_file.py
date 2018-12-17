
import numpy
import gc
import time
import random
import pysam
import copy
import os
import sys
import json
import pysam
#import sets
import coverage_file
import subprocess
import array
import numpy as np
try:
    import numpy
except ImportError:
    print('WARNING: module numpy was not loaded.')

try:
    import progressbar
except ImportError:
    print('WARNING: module progressbar was not loaded.')

import xlwt

#from matplotlib import pyplot

sys.path.append('/home/agarcia/PycharmProjects/ngscat3')

import bed_file


# TMP ='/home/fjavier/tmp/'
TMP = '/home/agarcia/tmp/'
# CHR_LENGTHS = '/data/reference_genomes/human/human_g1k_v37.genome'
CHR_LENGTHS = '/home/javi/MGP/data/reference_genomes/human/human_g1k_v37.genome'


class bam_file(pysam.Samfile):

    def __init__(self, _filename=None):
        self._nreads = None

        super().__init__()
        self._nreads = None

    def run(self, command):
        """*****************************************************************************************************************
        Task: launches a system call
        Inputs:
            command: string containing the system call.
        *****************************************************************************************************************"""

        print('CMD: ' + command)
        # Checks whether an error occurred during the execution of the system call
        try:
            subprocess.check_call(command)
        except:
            print('Some error ocurred while executing: ')
            print('', command)

    def nreads(self):
        return self.mapped + self.unmapped


    def issorted(self):
        print('Checking sorting of ' + str(self.filename))

        visitedcontigs = set()
        previouscontig = None
        previousstart = None
        previousend = None

        try:
            read = next(self)
            readsavailable = True
            sorted = True
            currend = read.pos + sum([nbases[1] for nbases in read.cigar if nbases[0] != 1]) - 1
        except StopIteration:
            readsavailable = False

        rc = 1
        #		while(readsavailable and ((previouscontig!=read.rname and (read.rname not in visitedcontigs)) or
        #								 (previousstart<read.pos or (previousstart==read.pos and previousend<=currend)))):
        while (readsavailable and ((previouscontig != read.rname and (read.rname not in visitedcontigs)) or
                                   (previousstart <= read.pos))):

            if (not rc % 1000000):
                print(str(rc) + ' reads checked')

            visitedcontigs.add(read.rname)
            previouscontig = read.rname
            previousstart = read.pos
            previousend = currend

            try:
                read = next(self)
                readsavailable = True
                currend = read.pos + sum([nbases[1] for nbases in read.cigar if nbases[0] != 1]) - 1
            except StopIteration:
                readsavailable = False
            except TypeError:
                if (read.is_unmapped):
                    print('ERROR: unmapped reads found at ' + str(self.filename))

                else:
                    print('ERROR: incorrect bam format')
                print(' Read position: ' + str(rc))
                print('	Alignment entry: ' + str(read))
                print(' Exiting.')
                sys.exit(1)

            rc += 1

        print('	Done.')

        if readsavailable:
            print('WARNING: bam not sorted.')
            print('	Check read ' + str(read))
            print('	and the read before.')

            return False

        else:
            return True

    def mappingsize(self):
        totalsize = 0

        for contig in self.header['SQ']:
            totalsize += contig['LN']

        return totalsize

    def sort_bam(self):
        """*******************************************************************************************************************************************
        JPFLORIDO.
        Task:  Sort BAM file by position and creat the corresponding index
        Inputs: BAM file to be sorted
        Outputs: Sorted file (object) with the corresponding index created
        *******************************************************************************************************************************************"""
        pid = str(os.getpid())
        print('Sorting BAM according to position...')
        sortedBAMfilename = TMP + pid + os.path.basename(self.filename.decode('utf-8')) + ".sorted"+ '.bam'
        print(sortedBAMfilename)
        #		self.run('samtools sort '+self.filename+' '+sortedBAMfilename)
        pysam.sort("-o", sortedBAMfilename, self.filename)

        pysam.index(sortedBAMfilename)

        return bam_file(sortedBAMfilename)

    # def myCoverageBed(self, target, numberReads=None, writeToFile=None, executiongranted=None, tmpdir=None,
    #                   bedGraphFile=None):
    #     """*************************************************************************************************************
    #             JPFLORIDO
    #             Task:  Custom method equivalent to bedtool's coverage bed
    #             Inputs: BAM file, target file (BED file), number of desired reads to be taken into account from the BAM
    #             file (default all) and option to write results to file the same way as bedtool does
    #             Outputs:
    #                 positionArray_intersect: a vector of bp coordinates in which coverage changes according to the
    #                 intersection of BAM file and target file
    #                 coverageArray_intersect: coverage (number of reads) for a given bp position
    #                 chromosomeCoordinates: a dictionary in which each entry is related to a chromosome and the content
    #                                        indicates the starting and ending position in positionArray_intersect and
    #                                         coverageArray_intersect arrays
    #                            finalBed: target file with overlapping exons removed and sorted (object)
    #             Other issues: it is intented to return also positions off target. However, coverage around exons changes
    #                          a lot, so there is an important increasing
    #             # in the memory used
    #             *****************************************************************************************************"""
    #
    #     global TMP
    #
    #     if (executiongranted != None):
    #         executiongranted.acquire()
    #
    #     if (tmpdir != None):
    #         TMP = tmpdir
    #
    #     pid = str(os.getpid())
    #
    #     # # Check whether BAM file is sorted
    #     # command = 'samtools view -H ' + self.filename + ' | grep SO:'
    #     # fd = os.popen(command)
    #     # outputCommand = fd.read()
    #     # # BAM not sorted -> Sort and indexing BAM
    #     # if ((fd.close() != None) or ('coordinate' not in outputCommand.split('SO:')[-1])):
    #     #     #			print 'ERROR: bam file must be sorted.'
    #     #     #			print '	Exiting.'
    #     #     #			sys.exit(1)
    #     #     sortedBam = self.sort_bam()
    #     # else:
    #     #     sortedBam = self
    #     #
    #     if (writeToFile != None):  # Results written to output file
    #         fdw = open(writeToFile, 'w')
    #     sortedBam = self
    #     positionArray_intersect = []
    #     coverageArray_intersect = []
    #     chromosomeCoordinates = {}  # Dictionary that controls, for a given chromosome,
    #     # its starting position and end position in the previous intersection arrays
    #
    #
    #     # Load target file, remove overlapping regions, sort it and load it
    #     bed = bed_file.bed_file(target)
    #     sortedBed = bed.my_sort_bed(tmpdir=TMP)
    #     # Base 1!!! # This generates a BED file in base 1 (Non-standard BED)
    #     nonOverlappingBed = sortedBed.non_overlapping_exons(1,tmpdir=TMP)
    #     finalBed = nonOverlappingBed.my_sort_bed(tmpdir=TMP)  # BED file in base 1 (Non-standard BED)
    #     finalBed.load_custom(-1)  # Load chromosome and positions in base 1....(finalBed is in base 1->Non-standard BED)
    #
    #     # Move along all chromosomes
    #     for currentChromosome in finalBed.chrs.keys():
    #         # Get all reads of current chromosome it there are reads in such chromosome
    #
    #         positionArray = []  # Stores bp positions along the current chromosome
    #         coverageArray = []  # Stores coverage values for a related position in chromosome (positionArray)
    #
    #         if (currentChromosome in sortedBam.references):
    #             # There might me information about the chromosome in the BAM header but no reads
    #             if (sortedBam.count(str(currentChromosome)) > 0):
    #                 allReads = sortedBam.fetch(str(currentChromosome))
    #
    #                 initPositions = []  # Structure that stores initial positions of each read
    #                 endPositions = []  # Structure that stores end positions of each read
    #                 for currentRead in allReads:
    #                     initPositions.append(int(currentRead.pos) + 1)#Fetch is 0-base indexing,working on 1-base indexing
    #                     endPositions.append(int(currentRead.aend))
    #
    #                 # Select a given number of reads according to numReadsDesired
    #                 if (numberReads is None):
    #                     numReadsDesired = sortedBam.nreads()
    #                 else:
    #                     numReadsDesired = numberReads
    #
    #                 numpy.random.seed(1)
    #                 selected = numpy.random.uniform(
    #                     size=len(initPositions)) <= (numReadsDesired * 1.0 / sortedBam.nreads())
    #
    #                 # Convert to numpy arrays
    #                 initPositions = numpy.array(initPositions, dtype=numpy.uint32)
    #                 endPositions = numpy.array(endPositions, dtype=numpy.uint32)
    #
    #                 selectedPositions_init = initPositions[selected]
    #                 selectedPositions_end = endPositions[selected]
    #                 selectedPositions_end += 1  # At the end of the position, the coverage counts. Sum 1 to say that at this position the coverage decreases
    #
    #                 # Sort each vector independtly
    #                 selectedPositions_init.sort()
    #                 selectedPositions_end.sort()
    #
    #                 totalPositions_init = len(selectedPositions_init)
    #                 totalPositions_end = len(selectedPositions_end)
    #
    #
    #                 # First iteration is done manually
    #                 positionArray.append(selectedPositions_init[0])  # Init position of the first read
    #                 coverageArray.append(1)  # A single read (coverage=1)
    #
    #                 indexInit = 1  # Controls index along selectedPositions_init (fist position in this array has been read)
    #                 indexEnd = 0  # Controls index along selectedPositions_end
    #
    #                 while indexEnd < totalPositions_end:  # While there are reads to be visited
    #                     if (indexInit < totalPositions_init):  # There are still reads that have to be visited
    #                         if (selectedPositions_init[indexInit] < selectedPositions_end[indexEnd]):  # If current position in init is smaller than current position in end
    #                             position = selectedPositions_init[indexInit]
    #                             indexInit += 1
    #                             partialSum = 1
    #                         elif (selectedPositions_init[indexInit] > selectedPositions_end[
    #                             indexEnd]):  # If current position in end is greater than current position in init
    #                             position = selectedPositions_end[indexEnd]
    #                             indexEnd += 1
    #                             partialSum = -1
    #                         else:  # If current position in init is equal to the position in end
    #                             position = selectedPositions_end[indexEnd]
    #                             indexInit += 1
    #                             indexEnd += 1
    #                             partialSum = 0
    #                     else:  # All starting positions for all reads have been visited
    #                         position = selectedPositions_end[indexEnd]
    #                         indexEnd += 1
    #                         partialSum = -1
    #
    #                     # Check whether position is already in the vector of positions
    #                     if (position == positionArray[-1]):  # More than a read starts or ends at the same time
    #                         coverageArray[-1] += partialSum
    #                     elif (partialSum != 0):  # If partialSum==0, then a read ends and a new read start -> do not update information
    #                         positionArray.append(position)
    #                         coverageArray.append(coverageArray[-1] + partialSum)
    #
    #                     pbarIter = indexInit + indexEnd
    #                 # pbar.update(pbarIter)
    #                 # pbar.finish()
    #
    #                 # Transform positionArray and coverageArray to numpy arrays to save memory
    #                 positionArray = numpy.array(positionArray, dtype=numpy.uint32)
    #                 coverageArray = numpy.array(coverageArray, dtype=numpy.uint16)
    #
    #                 numPositions = len(positionArray)
    #
    #                 # Create a bedgraph with coverage per position for all reads contained in the current chromosome
    #                 if (bedGraphFile != None):
    #                     extension = '.'
    #                     onlyName = os.path.basename(bedGraphFile)
    #                     components = onlyName.split(extension)
    #                     prefixFile = extension.join(components[:-1])
    #                     newFileName = prefixFile + '.' + str(currentChromosome) + '.' + components[-1]
    #                     if (len(os.path.dirname(bedGraphFile)) > 0):
    #                         newFileNameFULL = os.path.dirname(bedGraphFile) + '/' + newFileName
    #                     else:
    #                         newFileNameFULL = newFileName
    #
    #                     fdw_bedGraph = file(newFileNameFULL, 'w')
    #                     fdw_bedGraph.write('track type=bedGraph name=coverage_chr_' + str(currentChromosome) + str(self.filename)
    #                                        + '" description="coverage per position for ' + str(self.filename) + ' chromosome ' +
    #                                        str(currentChromosome) + '"\n')
    #
    #                     positionArray_bedGraph = positionArray - 1  # BED GRAPH displays in base 1, although data must use 0-base indexing (http://genome.ucsc.edu/FAQ/FAQtracks#tracks1). See format of bedgraph http://genome.ucsc.edu/goldenPath/help/bedgraph.html
    #                     coverageArray_bedGraph = coverageArray
    #                     for index in range(0, len(positionArray_bedGraph) - 1):
    #                         if (coverageArray_bedGraph[index] != 0):
    #                             fdw_bedGraph.write(str(currentChromosome) + ' ' + str(positionArray_bedGraph[index]) +' ' +
    #                                 str(positionArray_bedGraph[index + 1]) + ' ' + str(coverageArray_bedGraph[index]) +'\n')
    #
    #                     fdw_bedGraph.close()
    #             else:
    #                 numPositions = 0
    #
    #         else:  # There are no reads for the current chromosome
    #             numPositions = 0
    #
    #
    #         indexCoverage = 0
    #         nextPositionChange = 0
    #         previousIndexCoverage = indexCoverage
    #         firstPosition = len(positionArray_intersect)  # First position in the intersection vectors (positionArray_intersect and coverageArray_intersect)
    #
    #         for currentExon in finalBed.chrs[str(currentChromosome)]:
    #             exon_init = int(currentExon[0])
    #             exon_end = int(currentExon[1])
    #
    #             # Search for the first position that matches with the beginning of currentExon
    #             thereAreReads = False
    #             # firstIndexCoverage=indexCoverage
    #             readsAfterExonBeg = False  # Whether there are reads after the starting position of the exon ##NEW
    #
    #             # Check if change in coverage from the previous iteration is before exon_init. If not, move one position back so that previous coverage value is recovered
    #             if (indexCoverage > 0 and indexCoverage < numPositions and positionArray[nextPositionChange] > exon_init):
    #                 indexCoverage = indexCoverage - 1
    #
    #             if (indexCoverage < numPositions and positionArray[indexCoverage] > exon_init and positionArray[
    #                 indexCoverage] <= exon_end):  # There are reads in the exon, but the first read starts after the beginning of the exon
    #                 startPosition = indexCoverage
    #                 readsAfterExonBeg = True
    #             else:  # Reads start before the beginning of the exon
    #                 while (indexCoverage < numPositions and positionArray[indexCoverage] <= exon_init):
    #                     indexCoverage += 1
    #                 startPosition = indexCoverage - 1
    #
    #
    #
    #             if (indexCoverage != previousIndexCoverage):
    #                 thereAreReads = True
    #             # Search for the first position that matches with the end of currentExon
    #             while (indexCoverage < numPositions and positionArray[indexCoverage] <= exon_end):
    #                 indexCoverage += 1
    #             endPosition = indexCoverage - 1
    #
    #             nextPositionChange = indexCoverage  # positionArray[indexCoverage] contains next position after exon_end in which coverage changes
    #
    #             if (endPosition != startPosition or thereAreReads):  # Puede darse el caso de que un exon este en una sola coordenada de positionArray -> tenerlo en cuenta!!!!
    #                 numElements = endPosition - startPosition + 1  # Num of coverage/bases values to be interseted
    #                 if (writeToFile is None):
    #                     if (readsAfterExonBeg):
    #                         positionArray_intersect.extend([exon_init])
    #                         coverageArray_intersect.extend([0])
    #
    #                     positionArray_intersect.extend(positionArray[startPosition:(endPosition + 1)])
    #                     coverageArray_intersect.extend(coverageArray[startPosition:(endPosition + 1)])
    #
    #                     # It may happen that positionArray[startPosition] and positionArray[endPosition+1] are not equal to the beginning and end of the exon respectively
    #                     # Thus, it is forced to modify positionArray_intersect[indexIntersect] and positionArray_intersect[indexIntersect_end-1] to the beggining and end of the exon respectively
    #                     if (not readsAfterExonBeg):
    #                         positionArray_intersect[-numElements] = exon_init
    #                 # positionArray_intersect[-1]=exon_end # No es valido si la posicion final del exon tiene el mismo coverage que otra posicion anterior
    #                 else:  # Write to file base per base current exon
    #                     keys = range(exon_init, exon_end + 1)
    #                     dicExon = dict(zip(keys, [-1] * len(keys)))
    #
    #                     currentPositionArray = positionArray[startPosition:(endPosition + 1)]
    #                     if (not readsAfterExonBeg):
    #                         currentPositionArray[0] = exon_init
    #                     else:  # If reads after after exon_init, place a zero coverage to the beginning of the exon
    #                         dicExon[exon_init] = 0
    #
    #                     currentCoverageArray = coverageArray[startPosition:(endPosition + 1)]
    #                     dicExon.update(
    #                         zip(currentPositionArray, currentCoverageArray))  # Fill positions where coverage changes
    #
    #                     previousKey = sorted(dicExon.keys())[0]
    #
    #                     for currentKey in sorted(dicExon.keys()):
    #                         if (dicExon[currentKey] == -1):  # Take into account the coverage of the last position (key) that has coverage
    #                             # dicExon[currentKey]=dicExon[previousKey]
    #                             fdw.write(
    #                                 str(currentChromosome) + '\t' + str(exon_init) + '\t' + str(exon_end) + '\t' + str(
    #                                     dicExon[previousKey]) + '\n')
    #                         else:
    #                             fdw.write(
    #                                 str(currentChromosome) + '\t' + str(exon_init) + '\t' + str(exon_end) + '\t' + str(
    #                                     dicExon[currentKey]) + '\n')
    #                             previousKey = currentKey
    #
    #                 previousIndexCoverage = indexCoverage
    #             else:  # A exon is not covered -> startPosition moves until the end but endPosition could not move more. However, we want to point that, although there are no reads for
    #                 # the current exon, the initial position of the exon and its coverage (zero) are stored
    #                 if (writeToFile is None):
    #                     positionArray_intersect.append(exon_init)  # Initial position of exon
    #                     coverageArray_intersect.append(0)  # Zero coverage
    #                 else:
    #                     keys = range(exon_init, exon_end + 1)
    #                     dicExon = dict(zip(keys, [0] * len(keys)))
    #
    #                     for currentKey in sorted(dicExon.keys()):
    #                         fdw.write(str(currentChromosome) + '\t' + str(exon_init) + '\t' + str(exon_end) + '\t' + str(
    #                             dicExon[currentKey]) + '\n')
    #
    #                 indexCoverage = previousIndexCoverage  # Move backwards in positionArray to search reads for the next exon (all exons have not visited yet)
    #
    #         # Get the positions in the intersection vector for the current chromosome
    #         chromosomeCoordinates[currentChromosome] = (firstPosition, len(positionArray_intersect) - 1)
    #
    #
    #         del positionArray
    #         gc.collect()
    #
    #         del coverageArray
    #         gc.collect()
    #
    #     positionArray_intersect = numpy.array(positionArray_intersect, dtype=numpy.uint32)
    #     coverageArray_intersect = numpy.array(coverageArray_intersect, dtype=numpy.uint16)
    #
    #     if (writeToFile != None):
    #         fdw.close()
    #
    #     if (executiongranted != None):
    #         executiongranted.release()
    #
    #     return [positionArray_intersect, coverageArray_intersect, chromosomeCoordinates, finalBed]


    def myCoverageBed(self, target, numberReads=None, tmpdir=None,
                      bedGraphFile=None):
        """*************************************************************************************************************
                JPFLORIDO
                Task:  Custom method equivalent to bedtool's coverage bed
                Inputs: BAM file, target file (BED file), number of desired reads to be taken into account from the BAM
                file (default all) and option to write results to file the same way as bedtool does
                Outputs:
                    positionArray_intersect: a vector of bp coordinates in which coverage changes according to the
                    intersection of BAM file and target file
                    coverageArray_intersect: coverage (number of reads) for a given bp position
                    chromosomeCoordinates: a dictionary in which each entry is related to a chromosome and the content
                                           indicates the starting and ending position in positionArray_intersect and
                                            coverageArray_intersect arrays
                               finalBed: target file with overlapping exons removed and sorted (object)
                Other issues: it is intented to return also positions off target. However, coverage around exons changes
                             a lot, so there is an important increasing
                # in the memory used
                *****************************************************************************************************"""

        global TMP


        if (tmpdir != None):
            TMP = tmpdir

        pid = str(os.getpid())

        # # Check whether BAM file is sorted
        # command = 'samtools view -H ' + self.filename + ' | grep SO:'
        # fd = os.popen(command)
        # outputCommand = fd.read()
        # # BAM not sorted -> Sort and indexing BAM
        # if ((fd.close() != None) or ('coordinate' not in outputCommand.split('SO:')[-1])):
        #     #			print 'ERROR: bam file must be sorted.'
        #     #			print '	Exiting.'
        #     #			sys.exit(1)
        #     sortedBam = self.sort_bam()
        # else:
        #     sortedBam = self
        #
        # if (writeToFile != None):  # Results written to output file
        #     fdw = open(writeToFile, 'w')

        sortedBam = self
        positionArray_intersect = []
        coverageArray_intersect = []
        chromosomeCoordinates = {}  # Dictionary that controls, for a given chromosome,
        # its starting position and end position in the previous intersection arrays


        # Load target file, remove overlapping regions, sort it and load it
        bed = bed_file.bed_file(target)
        sortedBed = bed.my_sort_bed(tmpdir=TMP)
        # Base 1!!! # This generates a BED file in base 1 (Non-standard BED)
        nonOverlappingBed = sortedBed.non_overlapping_exons(1,tmpdir=TMP)
        finalBed = nonOverlappingBed.my_sort_bed(tmpdir=TMP)  # BED file in base 1 (Non-standard BED)
        finalBed.load_custom(-1)  # Load chromosome and positions in base 1....(finalBed is in base 1->Non-standard BED)

        # Instanciate class
        coveragefile = coverage_file.Coveragefile(self.filename)
        coveragefile.bedfilename  = target
        covertotal = array.array('I')

        idxregion = 0
        jdxregion = 0

        # Move along all chromosomes
        for currentChromosome in finalBed.chrs.keys():

            # Get all reads of current chromosome it there are reads in such chromosome
            positionArray = []  # Stores bp positions along the current chromosome
            coverageArray = []  # Stores coverage values for a related position in chromosome (positionArray)

            if (currentChromosome in sortedBam.references):
                # There might me information about the chromosome in the BAM header but no reads
                if (sortedBam.count(str(currentChromosome)) > 0):
                    allReads = sortedBam.fetch(str(currentChromosome))

                    initPositions = []  # Structure that stores initial positions of each read
                    endPositions = []  # Structure that stores end positions of each read
                    for currentRead in allReads:
                        initPositions.append(int(currentRead.pos) + 1)#Fetch is 0-base indexing,working on 1-base indexing
                        endPositions.append(int(currentRead.aend))

                    # Select a given number of reads according to numReadsDesired
                    if (numberReads is None):
                        numReadsDesired = sortedBam.nreads()
                    else:
                        numReadsDesired = numberReads

                    numpy.random.seed(1)
                    selected = numpy.random.uniform(
                        size=len(initPositions)) <= (numReadsDesired * 1.0 / sortedBam.nreads())

                    # Convert to numpy arrays
                    initPositions = numpy.array(initPositions, dtype=numpy.uint32)
                    endPositions = numpy.array(endPositions, dtype=numpy.uint32)

                    selectedPositions_init = initPositions[selected]
                    selectedPositions_end = endPositions[selected]
                    selectedPositions_end += 1  # At the end of the position, the coverage counts. Sum 1 to say that at this position the coverage decreases

                    # Sort each vector independtly
                    selectedPositions_init.sort()
                    selectedPositions_end.sort()

                    totalPositions_init = len(selectedPositions_init)
                    totalPositions_end = len(selectedPositions_end)




                    # First iteration is done manually
                    positionArray.append(selectedPositions_init[0])  # Init position of the first read
                    coverageArray.append(1)  # A single read (coverage=1)

                    indexInit = 1  # Controls index along selectedPositions_init (fist position in this array has been read)
                    indexEnd = 0  # Controls index along selectedPositions_end

                    while indexEnd < totalPositions_end:  # While there are reads to be visited
                        if (indexInit < totalPositions_init):  # There are still reads that have to be visited
                            if (selectedPositions_init[indexInit] < selectedPositions_end[indexEnd]):  # If current position in init is smaller than current position in end
                                position = selectedPositions_init[indexInit]
                                indexInit += 1
                                partialSum = 1
                            elif (selectedPositions_init[indexInit] > selectedPositions_end[
                                indexEnd]):  # If current position in end is greater than current position in init
                                position = selectedPositions_end[indexEnd]
                                indexEnd += 1
                                partialSum = -1
                            else:  # If current position in init is equal to the position in end
                                position = selectedPositions_end[indexEnd]
                                indexInit += 1
                                indexEnd += 1
                                partialSum = 0
                        else:  # All starting positions for all reads have been visited
                            position = selectedPositions_end[indexEnd]
                            indexEnd += 1
                            partialSum = -1

                        # Check whether position is already in the vector of positions
                        if (position == positionArray[-1]):  # More than a read starts or ends at the same time
                            coverageArray[-1] += partialSum
                        elif (partialSum != 0):  # If partialSum==0, then a read ends and a new read start -> do not update information
                            positionArray.append(position)
                            coverageArray.append(coverageArray[-1] + partialSum)

                        pbarIter = indexInit + indexEnd
                    # pbar.update(pbarIter)
                    # pbar.finish()

                    # Transform positionArray and coverageArray to numpy arrays to save memory
                    positionArray = numpy.array(positionArray, dtype=numpy.uint32)
                    coverageArray = numpy.array(coverageArray, dtype=numpy.uint16)

                    numPositions = len(positionArray)

                    # Create a bedgraph with coverage per position for all reads contained in the current chromosome
                    if (bedGraphFile != None):
                        extension = '.'
                        onlyName = os.path.basename(bedGraphFile)
                        components = onlyName.split(extension)
                        prefixFile = extension.join(components[:-1])
                        newFileName = prefixFile + '.' + str(currentChromosome) + '.' + components[-1]
                        if (len(os.path.dirname(bedGraphFile)) > 0):
                            newFileNameFULL = os.path.dirname(bedGraphFile) + '/' + newFileName
                        else:
                            newFileNameFULL = newFileName

                        fdw_bedGraph = open(newFileNameFULL, 'w')
                        fdw_bedGraph.write('track type=bedGraph name=coverage_chr_' + str(currentChromosome) + str(self.filename)
                                           + '" description="coverage per position for ' + str(self.filename) + ' chromosome ' +
                                           str(currentChromosome) + '"\n')

                        positionArray_bedGraph = positionArray - 1  # BED GRAPH displays in base 1, although data must use 0-base indexing (http://genome.ucsc.edu/FAQ/FAQtracks#tracks1). See format of bedgraph http://genome.ucsc.edu/goldenPath/help/bedgraph.html
                        coverageArray_bedGraph = coverageArray
                        for index in range(0, len(positionArray_bedGraph) - 1):
                            if (coverageArray_bedGraph[index] != 0):
                                fdw_bedGraph.write(str(currentChromosome) + ' ' + str(positionArray_bedGraph[index]) +' ' +
                                                   str(positionArray_bedGraph[index + 1]) + ' ' + str(coverageArray_bedGraph[index]) +'\n')

                        fdw_bedGraph.close()
                else:
                    numPositions = 0

            else:  # There are no reads for the current chromosome
                numPositions = 0

            # COVERAGE FILE GENERATION, takes 2 list below.

            # Generate chromosome object that will contain regions objects
            chromosome = coverage_file.Chromosome(currentChromosome)

            indexCoverage = 0
            nextPositionChange = 0
            previousIndexCoverage = indexCoverage
            firstPosition = len(positionArray_intersect)  # First position in the intersection vectors (positionArray_intersect and coverageArray_intersect)

            for currentExon in finalBed.chrs[str(currentChromosome)]:

                region = coverage_file.Region()
                region.start = int(currentExon[0])
                region.end = int(currentExon[1])
                covlist = []
                exon_init = int(currentExon[0])
                exon_end = int(currentExon[1])
                # lexon.append([exon_init, exon_end])


                # Search for the first position that matches with the beginning of currentExon
                thereAreReads = False
                # firstIndexCoverage=indexCoverage
                readsAfterExonBeg = False  # Whether there are reads after the starting position of the exon ##NEW

                # Check if change in coverage from the previous iteration is before exon_init. If not, move one position back so that previous coverage value is recovered
                if (indexCoverage > 0 and indexCoverage < numPositions and positionArray[nextPositionChange] > exon_init):
                    indexCoverage = indexCoverage - 1

                if (indexCoverage < numPositions and positionArray[indexCoverage] > exon_init and positionArray[
                    indexCoverage] <= exon_end):  # There are reads in the exon, but the first read starts after the beginning of the exon
                    startPosition = indexCoverage
                    readsAfterExonBeg = True
                else:  # Reads start before the beginning of the exon
                    while (indexCoverage < numPositions and positionArray[indexCoverage] <= exon_init):
                        indexCoverage += 1
                    startPosition = indexCoverage - 1



                if (indexCoverage != previousIndexCoverage):
                    thereAreReads = True
                # Search for the first position that matches with the end of currentExon
                while (indexCoverage < numPositions and positionArray[indexCoverage] <= exon_end):
                    indexCoverage += 1
                endPosition = indexCoverage - 1

                nextPositionChange = indexCoverage  # positionArray[indexCoverage] contains next position after exon_end in which coverage changes

                if (endPosition != startPosition or thereAreReads):  # Puede darse el caso de que un exon este en una sola coordenada de positionArray -> tenerlo en cuenta!!!!

                    keys = range(exon_init, exon_end + 1)
                    dicExon = dict(zip(keys, [-1] * len(keys)))

                    currentPositionArray = positionArray[startPosition:(endPosition + 1)]
                    if (not readsAfterExonBeg):
                        currentPositionArray[0] = exon_init
                    else:  # If reads after after exon_init, place a zero coverage to the beginning of the exon
                        dicExon[exon_init] = 0

                    currentCoverageArray = coverageArray[startPosition:(endPosition + 1)]
                    # Fill positions where coverage changes
                    dicExon.update(zip(currentPositionArray, currentCoverageArray))

                    previousKey = sorted(dicExon.keys())[0]

                    for currentKey in sorted(dicExon.keys()):
                        jdxregion += 1

                        if (dicExon[currentKey] == -1):  # Take into account the coverage of the last position (key) that has coverage

                            covlist.append(dicExon[previousKey])

                        else:
                            covlist.append(dicExon[currentKey])
                            previousKey = currentKey

                    covertotal.extend(covlist)
                    region.covStartIndex = idxregion
                    region.covEndIndex = jdxregion
                    region.mean = np.mean(covlist)
                    region.std = np.std(covlist)

                    #Go for the next region
                    idxregion = jdxregion

                    previousIndexCoverage = indexCoverage

                else:  # A exon is not covered -> startPosition moves until the end but endPosition could not move more. However, we want to point that, although there are no reads for
                    # the current exon, the initial position of the exon and its coverage (zero) are stored

                    keys = range(exon_init, exon_end + 1)
                    dicExon = dict(zip(keys, [0] * len(keys)))

                    for currentKey in sorted(dicExon.keys()):
                        jdxregion += 1

                        covlist.append(dicExon[currentKey])


                    covertotal.extend(covlist)
                    region.covStartIndex = idxregion
                    region.covEndIndex = jdxregion
                    region.mean = np.mean(covlist)
                    region.std = np.std(covlist)
                    region.zeropos = [i for i,x in enumerate(covlist) if x == 0]
                    idxregion = jdxregion
                    # Move backwards in positionArray to search reads for the next exon (all exons have not visited yet)
                    indexCoverage = previousIndexCoverage

                chromosome.regions.append(region)



            coveragefile.chromosomes.append(chromosome)
        coveragefile.coverages = np.array(covertotal ,dtype= np.uint32)

        del positionArray
        gc.collect()

        del coverageArray
        gc.collect()



        return [coveragefile, finalBed]

    def myReadsOnTarget(self, target):
        """*******************************************************************************************************************************************
        JPFLORIDO
        Task:  Custom method equivalent to bedtool's intersect bed. It is made to count the number of reads on/off target and the number of reads on/off target that start and end at the same position
        Inputs: BAM file (self) and target file
        Outputs:
            dicOnTarget: dictionary containing the number of reads on target per chromosome
            dicOffTarget: dictionary containing the number of reads off target per chromosome
            dicTotal: dictionary containing the number of total reads (on + off target) per chromosome
            duplicatesOnTarget: array containing number of reads on target that start and end at the same position at different levels of reads per start/end position; duplicatesOnTarget[0] = number of reads that start/end at a unique position;  duplicatesOnTarget[1] = number of reads that start and end at the same position (at most two reads per start/end positions); uplicatesOnTarget[n] = no.of reads that start and end at the same position (at most n+1 reads per start/end positions)
            duplicatesOffTarget: array containing number of reads off target that start and end at the same position at different levels of reads per start/end position; duplicatesOffTarget[0] = number of reads that start/end at a unique position;  duplicatesOffTarget[1] = number of reads that start and end at the same position (at most two reads per start/end positions); uplicatesOffTarget[n] = no.of reads that start and end at the same position (at most n+1 reads per start/end positions)


        Other issues: if a read targets two regions in the BED file, the read is count twice. Pysam works in 0-base indexing, so, the overlapping regions
        # will be removed from the target file as they are -> no conversion to real zero or one-base indexing
        *******************************************************************************************************************************************"""

        global TMP

        # Load target file, remove overlapping regions, sort it and load it
        bed = bed_file.bed_file(target)
        sortedBed = bed.my_sort_bed(tmpdir=TMP)
        nonOverlappingBed = sortedBed.non_overlapping_exons(1, tmpdir=TMP)  # Base-1 indexing
        finalBed = nonOverlappingBed.my_sort_bed(tmpdir=TMP)  # BED file in base 1 (Non-standard BED)
        finalBed.load_custom(-1)  # Load chromosome and positions as they are

        dicOnTarget = {}
        dicOffTarget = {}
        dicTotalReads = {}
        countDuplicatesOnTarget = {}
        countDuplicatesOffTarget = {}
        previousChromosome = '0'  # False chromosome
        dicOnTargetChr = {}
        dicOffTargetChr = {}
        readsOnTarget = 0  # fjavier: Overall number of reads on target
        readsOnTargetChr = 0
        readsOffTargetChr = 0
        thereisInfo = False
        duplicatesOnTarget = []  # Contains reads that appear once, two times, three times...on target
        duplicatesOffTarget = []  # Contains reads that appear once, two times, three times...off target

        for currentRead in self:
            # Get chromosome info from currentRead
            try:
                currentChromosome = self.getrname(currentRead.tid)
            except ValueError:
                if (currentRead.tid < 0):
                    print("\nPlease, check that BAM file has only mapped reads")
                    exit()
            if (currentChromosome != previousChromosome):
                #				print "Examining Chromosome "+str(currentChromosome)+'... \n'
                if (thereisInfo):  # A new chromosome is started. Store information about previous chromosome
                    #					print "Saving information of chromosome "+ str(previousChromosome)+'\n'
                    readsOnTarget += readsOnTargetChr
                    dicOnTarget[str(previousChromosome)] = readsOnTargetChr  # Count reads on target for the current chromosome
                    dicOffTarget[str(previousChromosome)] = readsOffTargetChr  # Count reads off target for the current chromosome
                    dicTotalReads[str(previousChromosome)] = readsOnTargetChr + readsOffTargetChr  # Count all reads for the current chromosome

                    # Get counts of each read for the current chromosome and accumulate them in two dictionaries
                    countReadsOnTarget = dicOnTargetChr.values()
                    if (len(countReadsOnTarget) > 0):
                        for count in range(min(countReadsOnTarget), max(countReadsOnTarget) + 1):
                            if (count in countDuplicatesOnTarget):
                                countDuplicatesOnTarget[count] += sum([1 for x in countReadsOnTarget if x == count])
                            else:
                                countDuplicatesOnTarget[count] = sum([1 for x in countReadsOnTarget if x == count])

                    countReadsOffTarget = dicOffTargetChr.values()
                    if (len(countReadsOffTarget) > 0):
                        for count in range(min(countReadsOffTarget), max(countReadsOffTarget) + 1):
                            if (count in countDuplicatesOffTarget):
                                countDuplicatesOffTarget[count] += sum([1 for x in countReadsOffTarget if x == count])
                            else:
                                countDuplicatesOffTarget[count] = sum([1 for x in countReadsOffTarget if x == count])

                    del dicOnTargetChr
                    gc.collect()

                    del dicOffTargetChr
                    gc.collect()

                if (currentChromosome in finalBed.chrs.keys()):  # Check whether chromosome related to the read is in the BED file
                    targets_chr = numpy.array(finalBed.chrs[str(currentChromosome)])  # all rows, 1st column -> exon_init positions; all rows, 2nd column -> exon_end positions
                    numRegions = len(targets_chr[:, 0])

                dicOnTargetChr = {}
                dicOffTargetChr = {}
                readsOnTargetChr = 0
                readsOffTargetChr = 0
                startPos = 0
                previousChromosome = currentChromosome
                thereisInfo = False

            if (currentChromosome in finalBed.chrs.keys()):  # Go through the process only if currentchromosome is present in the region
                thereisInfo = True
                currentReadInit = int(currentRead.pos) + 1  # BED has been transformed to base 1
                currentReadEnd = int(currentRead.aend)
                currentPos = startPos
                stop = False

                while ((not stop) and currentPos < numRegions):
                    if ((currentReadInit >= targets_chr[currentPos, 0] and currentReadInit <= targets_chr[currentPos, 1])
                            or(currentReadEnd >= targets_chr[currentPos, 0] and currentReadEnd <= targets_chr[currentPos, 1])):

                        # Read on target: read init inside the region or read end inside the region

                        if ((currentReadInit, currentReadEnd) in dicOnTargetChr):  # Read on target
                            dicOnTargetChr[(currentReadInit, currentReadEnd)] += 1
                        else:
                            dicOnTargetChr[(currentReadInit, currentReadEnd)] = 1

                        readsOnTargetChr += 1

                        stop = True
                        startPos = currentPos
                    else:  # Check to the next region, unless a read ends before the current region
                        if (currentReadEnd < targets_chr[currentPos, 0]):
                            # A read ends before the current region -> stop searching. Read off target
                            if ((currentReadInit, currentReadEnd) in dicOffTargetChr):
                                dicOffTargetChr[(currentReadInit, currentReadEnd)] += 1
                            else:
                                dicOffTargetChr[(currentReadInit, currentReadEnd)] = 1

                            readsOffTargetChr += 1
                            stop = True
                        else:
                            currentPos += 1

                if (currentPos == numRegions):  # Read off target -> we moved along all regions and none of them mathes with the read
                    if ((currentReadInit, currentReadEnd) in dicOffTargetChr):
                        dicOffTargetChr[(currentReadInit, currentReadEnd)] += 1
                    else:
                        dicOffTargetChr[(currentReadInit, currentReadEnd)] = 1

                    readsOffTargetChr += 1

        if (thereisInfo):
            # Store information for the last chromosome
            #			print "Saving information of chromosome "+ str(currentChromosome)+'\n'
            readsOnTarget += readsOnTargetChr
            dicOnTarget[
                str(currentChromosome)] = readsOnTargetChr  # Count reads on target for the current chromosome
            dicOffTarget[
                str(currentChromosome)] = readsOffTargetChr  # Count reads off target for the current chromosome
            dicTotalReads[str(
                currentChromosome)] = readsOnTargetChr + readsOffTargetChr  # Count all reads for the current chromosome

            # Get counts of each read for the current chromosome and accumulate them in two dictionaries
            countReadsOnTarget = list(dicOnTargetChr.values())
            if (len(countReadsOnTarget) > 0):
                for count in range(min(countReadsOnTarget), max(countReadsOnTarget) + 1):
                    if (count in countDuplicatesOnTarget):
                        countDuplicatesOnTarget[count] += sum([1 for x in countReadsOnTarget if x == count])
                    else:
                        countDuplicatesOnTarget[count] = sum([1 for x in countReadsOnTarget if x == count])

            countReadsOffTarget = list(dicOffTargetChr.values())
            if (len(countReadsOffTarget) > 0):
                for count in range(min(countReadsOffTarget), max(countReadsOffTarget) + 1):
                    if (count in countDuplicatesOffTarget):
                        countDuplicatesOffTarget[count] += sum([1 for x in countReadsOffTarget if x == count])
                    else:
                        countDuplicatesOffTarget[count] = sum([1 for x in countReadsOffTarget if x == count])

        # Store information about duplicates in an array

        if (len(countDuplicatesOnTarget.keys()) > 0):
            duplicatesOnTarget = numpy.zeros(max(countDuplicatesOnTarget.keys()))
            for currentKey in sorted(countDuplicatesOnTarget.keys()):
                duplicatesOnTarget[currentKey - 1] = countDuplicatesOnTarget[currentKey] * currentKey
            # It is multiplied by currentKey, since countDuplicatesOnTarget[currentKey] has the number of start/end positions that have currentKey reads. We want the number of total reads.

        if (len(countDuplicatesOffTarget.keys()) > 0):
            duplicatesOffTarget = numpy.zeros(max(countDuplicatesOffTarget.keys()))
            for currentKey in sorted(countDuplicatesOffTarget.keys()):
                duplicatesOffTarget[currentKey - 1] = countDuplicatesOffTarget[currentKey] * currentKey

        # Check whether there are no reads for a given chromosome contained in BED file

        for currentChromosome in finalBed.chrs.keys():
            if (not currentChromosome in dicOnTarget.keys()):  # If chromosome is not in the dictionary is because there are no reads for that chromosome
                dicOnTarget[str(currentChromosome)] = 0
                dicTotalReads[str(currentChromosome)] = 0
        # FJAVIER: dicOffTarget is not actually needed yet, that is why it is not being returned
        return [readsOnTarget, dicOnTarget, dicTotalReads, duplicatesOnTarget, duplicatesOffTarget]


    def reads_on_target(self, bed, outdir, bamlist=[], legend=None, maxduplicates = 5, executiongranted=None, onoff_status=None,
                        duplicates_status=None,
                        retonduplicates=None, retoffduplicates=None, enrichment=None, percontarget=None, tmpdir=None,
                        warnthreshold=80):
        """************************************************************************************************************************************************************
        Task: Print reads on traget and off target
        Inputs:
            bed: bed file with capture coordinates
            outdir: Output folder
            bamlist: list of bam_file objects representing bam files which are also wanted to be analyzed.
            legend: list of strings containing the label to tag each bam file in the xls and png files that will be created.

            onoff_status: multiprocessing.Value object to return whether the number of on-target reads is extremely low (False) or not (True)
            duplicates_status: multiprocessing.Value object to return whether the number of duplicated reads on-target is greater than the number of duplicated
                off-target (False) or not (True).
            enrichment: multiprocessing.Array objecto to return the enrichment value for each bam (on-target reads per Kb)/(off target reads per Kb)
            percontarget: multiprocessing.Array object to return percentaje of reads on target for each bam file
            tmpdir: string containing the path to a temporary directory where temporary files will be stored.
        Outputs: dictionary and json with information and bam results in the key 'results'
        ************************************************************************************************************************************************************"""
        read_on_results = {}
        global TMP


        # Check whether to use default tmpdir
        if (tmpdir is not None):
            TMP = tmpdir

        # Create a list of bam_file objects which includes self. Generates an array with the total number of reads in each bam.
        bamlist = [self] + bamlist
        tread = []


        # Calculate number of reads and duplicated reads on/off target per chromosome
        nread = []
        onperchr = []
        totalperchr = []
        onduplicates = []
        offduplicates = []
        percontarget = []
        perconperchr = []
        perconduplicates = []
        percoffduplicates = []
        enrichment = []

        #Adding data
        for bam in bamlist:
            nread_tmp, onperchr_tmp, totalperchr_tmp, onduplicates_tmp, offduplicates_tmp = bam.myReadsOnTarget(bed)
            tread.append(bam.nreads())
            nread.append(nread_tmp)
            onperchr.append(onperchr_tmp)
            totalperchr.append(totalperchr_tmp)
            onduplicates.append(onduplicates_tmp)
            offduplicates.append(offduplicates_tmp)

        bedobj = bed_file.bed_file(bed)
        targetsize = bedobj.size()

        retonduplicates= []
        retoffduplicates = []
        onduplicatesresult = []
        offduplicatesresult = []
        onoff_status = []
        duplicates_status = []

        # Generating list that will form the output dictionary and Json
        for i in range(len(bamlist)):

            # Calculate enrichment, same bed will be used for both bams
            if (tread[i] == nread[i]):
                enrichment.append(-1)
            else:
                enrichment.append((nread[i] * 1000.0 / targetsize) / ((tread[i] - nread[i]) * 1000.0 / (self.mappingsize() - targetsize)))

            percontarget.append(nread[i]*100.0/tread[i])
            retonduplicates.append(sum(onduplicates[i]) * 100.0 / tread[i])
            retoffduplicates.append(sum(offduplicates[i]) * 100.0 / tread[i])

            #Avoid 0 division
            perconperchr.append({key: (onperchr[i][key] * 100.0/totalperchr[i][key] if totalperchr[i][key] > 0 else 0) for key in onperchr[i]})


            # Select the largest one. Use the keys and fill with 0 until key lenght
            if len(onduplicates[i].tolist()) >= len(offduplicates[i].tolist()):
                onduplicatesresult.append({str(key + 1) + 'x': value for key, value in enumerate(onduplicates[i].tolist())})

                offduplicatesresult.append(dict(zip(list(onduplicatesresult[i].keys()), offduplicates[i].tolist() + [0] *
                                               (len(onduplicatesresult[i].keys()) - len(offduplicates[i].tolist())))))
            else:
                offduplicatesresult.append({str(key + 1) + 'x': value for key, value in enumerate(offduplicates[i].tolist())})
                onduplicatesresult.append(dict(zip(list(offduplicatesresult[i].keys()), onduplicates[i].tolist() + [0] *
                                               (len(offduplicatesresult[i].keys()) - len(onduplicates[i].tolist())))))

            # Compute percentage per number of replicates, each element, 1x is number of reads that appears 1 time /
            # total on reads.
            perconduplicates.append({key: (value / nread[i] * 100.0 if nread[i] > 0 else 0 )for key, value in onduplicatesresult[i].items()})
            percoffduplicates.append({key: (value / (tread[i]- nread[i]) * 100.0 if (tread[i] - nread[i]) > 0 else 0)
                                      for key, value in offduplicatesresult[i].items()})

            onoff_status.append(True if percontarget[i] >= warnthreshold else False)
            duplicates_status.append(True if retonduplicates[i] > retoffduplicates[i] else False)
        #Output generation, reads_on_target.json. Status variables doesn't go inside the json (only str is allowed not booleans)
        results = []
        for i in range(len(bamlist)):
            results.append(dict(
                [('bamfilename', bamlist[i].filename.decode('utf-8')),
                 ('enrichment', enrichment[i]),
                 ('onread', nread[i]),
                 ('percontotal', percontarget[i]),
                 ('totalread', tread[i]),
                 ('onperchr', onperchr[i]),
                 ('totalperchr', totalperchr[i]),
                 ('perconperchr', perconperchr[i]),
                 ('retonduplicates',retonduplicates),
                 ('retoffduplicates',retoffduplicates),
                 ('legend', (legend[i] if legend is not None else bamlist[i].filename.decode('utf-8').split('/')[-1])),
                 ('onduplicates', onduplicatesresult[i]),
                 ('offduplicates', offduplicatesresult[i]),
                 ('perconduplicates', perconduplicates[i]),
                 ('percoffduplicates', percoffduplicates[i]),
                 ('onoff_status', 'OK' if onoff_status[i] else 'Not OK'),
                 ('duplicates_status', 'OK' if duplicates_status[i] else 'Not OK')]))
        read_on_results['results'] = results
        read_on_results['bedfile'] = bed
        read_on_results['outdir'] = outdir
        read_on_results['maxduplicates'] = maxduplicates


        with open(outdir +'/read_on_results.json', 'w') as outfile:
            json.dump(read_on_results, outfile)
        return read_on_results,onoff_status,duplicates_status





    def reads_on_target(self, bed, outdir, bamlist=[], legend=None, maxduplicates = 5, executiongranted=None, onoff_status=None,
                        duplicates_status=None,
                        retonduplicates=None, retoffduplicates=None, enrichment=None, percontarget=None, tmpdir=None,
                        warnthreshold=80):
        """************************************************************************************************************************************************************
        Task: Print reads on traget and off target
        Inputs:
            bed: bed file with capture coordinates
            outdir: Output folder
            bamlist: list of bam_file objects representing bam files which are also wanted to be analyzed.
            legend: list of strings containing the label to tag each bam file in the xls and png files that will be created.

            onoff_status: multiprocessing.Value object to return whether the number of on-target reads is extremely low (False) or not (True)
            duplicates_status: multiprocessing.Value object to return whether the number of duplicated reads on-target is greater than the number of duplicated
                off-target (False) or not (True).
            enrichment: multiprocessing.Array objecto to return the enrichment value for each bam (on-target reads per Kb)/(off target reads per Kb)
            percontarget: multiprocessing.Array object to return percentaje of reads on target for each bam file
            tmpdir: string containing the path to a temporary directory where temporary files will be stored.
        Outputs: dictionary and json with information and bam results in the key 'results'
        ************************************************************************************************************************************************************"""
        read_on_results = {}
        global TMP


        # Check whether to use default tmpdir
        if (tmpdir is not None):
            TMP = tmpdir

        # Create a list of bam_file objects which includes self. Generates an array with the total number of reads in each bam.
        bamlist = [self] + bamlist
        tread = []


        # Calculate number of reads and duplicated reads on/off target per chromosome
        nread = []
        onperchr = []
        totalperchr = []
        onduplicates = []
        offduplicates = []
        percontarget = []
        perconperchr = []
        perconduplicates = []
        percoffduplicates = []
        enrichment = []

        #Adding data
        for bam in bamlist:
            nread_tmp, onperchr_tmp, totalperchr_tmp, onduplicates_tmp, offduplicates_tmp = bam.myReadsOnTarget(bed)
            tread.append(bam.nreads())
            nread.append(nread_tmp)
            onperchr.append(onperchr_tmp)
            totalperchr.append(totalperchr_tmp)
            onduplicates.append(onduplicates_tmp)
            offduplicates.append(offduplicates_tmp)

        bedobj = bed_file.bed_file(bed)
        targetsize = bedobj.size()

        retonduplicates= []
        retoffduplicates = []
        onduplicatesresult = []
        offduplicatesresult = []
        onoff_status = []
        duplicates_status = []

        # Generating list that will form the output dictionary and Json
        for i in range(len(bamlist)):

            # Calculate enrichment, same bed will be used for both bams
            if (tread[i] == nread[i]):
                enrichment.append(-1)
            else:
                enrichment.append((nread[i] * 1000.0 / targetsize) / ((tread[i] - nread[i]) * 1000.0 / (self.mappingsize() - targetsize)))

            percontarget.append(nread[i]*100.0/tread[i])
            retonduplicates.append(sum(onduplicates[i]) * 100.0 / tread[i])
            retoffduplicates.append(sum(offduplicates[i]) * 100.0 / tread[i])

            #Avoid 0 division
            perconperchr.append({key: (onperchr[i][key] * 100.0/totalperchr[i][key] if totalperchr[i][key] > 0 else 0) for key in onperchr[i]})


            # Select the largest one. Use the keys and fill with 0 until key lenght
            if len(onduplicates[i].tolist()) >= len(offduplicates[i].tolist()):
                onduplicatesresult.append({str(key + 1) + 'x': value for key, value in enumerate(onduplicates[i].tolist())})

                offduplicatesresult.append(dict(zip(list(onduplicatesresult[i].keys()), offduplicates[i].tolist() + [0] *
                                               (len(onduplicatesresult[i].keys()) - len(offduplicates[i].tolist())))))
            else:
                offduplicatesresult.append({str(key + 1) + 'x': value for key, value in enumerate(offduplicates[i].tolist())})
                onduplicatesresult.append(dict(zip(list(offduplicatesresult[i].keys()), onduplicates[i].tolist() + [0] *
                                               (len(offduplicatesresult[i].keys()) - len(onduplicates[i].tolist())))))

            # Compute percentage per number of replicates, each element, 1x is number of reads that appears 1 time /
            # total on reads.
            perconduplicates.append({key: (value / nread[i] * 100.0 if nread[i] > 0 else 0 )for key, value in onduplicatesresult[i].items()})
            percoffduplicates.append({key: (value / (tread[i]- nread[i]) * 100.0 if (tread[i] - nread[i]) > 0 else 0)
                                      for key, value in offduplicatesresult[i].items()})

            onoff_status.append(True if percontarget[i] >= warnthreshold else False)
            duplicates_status.append(True if retonduplicates[i] > retoffduplicates[i] else False)
        #Output generation, reads_on_target.json. Status variables doesn't go inside the json (only str is allowed not booleans)
        results = []
        for i in range(len(bamlist)):
            results.append(dict(
                [('bamfilename', bamlist[i].filename.decode('utf-8')),
                 ('enrichment', enrichment[i]),
                 ('onread', nread[i]),
                 ('percontotal', percontarget[i]),
                 ('totalread', tread[i]),
                 ('onperchr', onperchr[i]),
                 ('totalperchr', totalperchr[i]),
                 ('perconperchr', perconperchr[i]),
                 ('retonduplicates',retonduplicates),
                 ('retoffduplicates',retoffduplicates),
                 ('legend', (legend[i] if legend is not None else bamlist[i].filename.decode('utf-8').split('/')[-1])),
                 ('onduplicates', onduplicatesresult[i]),
                 ('offduplicates', offduplicatesresult[i]),
                 ('perconduplicates', perconduplicates[i]),
                 ('percoffduplicates', percoffduplicates[i]),
                 ('onoff_status', 'OK' if onoff_status[i] else 'Not OK'),
                 ('duplicates_status', 'OK' if duplicates_status[i] else 'Not OK')]))
        read_on_results['results'] = results
        read_on_results['bedfile'] = bed
        read_on_results['outdir'] = outdir
        read_on_results['maxduplicates'] = maxduplicates


        with open(outdir +'/read_on_results.json', 'w') as outfile:
            json.dump(read_on_results, outfile)
        return read_on_results,onoff_status,duplicates_status


# bam = bam_file('/home/agarcia/PycharmProjects/ngscat/example2.bam')
#
# json_file = bam.reads_on_target('/home/agarcia/PycharmProjects/ngscat/seqcap.example2.bed',
#                                       '/home/agarcia/PycharmProjects/ngscat/')
#
#
#
# onduplicates