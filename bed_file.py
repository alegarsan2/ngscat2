import region

# try:
#     from matplotlib import pyplot
# except ImportError:
#     print('WARNING: module pyplot was not loaded.')

try:
    import progressbar
except ImportError:
    print('WARNING: module progressbar was not loaded.')

import string
#import sets

import os

import numpy

import pysam

import bam_file

CHR_LENGTHS = '/usr/local/reference_genomes/human/human_g1k_v37.1-22XYM.genome'
TMP = '/tmp/'
BEDTOOLS = '/usr/local/bedtools/bin/'


class bed_file:

    def __init__(self, _filename):
        self.filename = _filename
        self.chrs = None
        self.nregions = None


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

    def count_lines(self, filename):
        return len (open(self.filename, errors= 'ignore').readlines())

    def load_custom(self, base):
        """************************************************************************************************************************************************************
        JPFLORIDO
        Task: loads data into self.chrs and self.nregions according to the base indicated as argument
        Inputs: BED file (self) and base (1 or 0). If base is not 1 or 0, bed file is loaded as it is....
        Output: self.chrs and self.nregions are modified.
            self.chrs: dictionary. Each key represents a chromosome. Values are lists of tuples (start,end) indicating
                       each of the regions in the chromosome.
            self.nregions = self.count_lines(self.filename)
            bed file loaded with coordinates in base "base"

        ************************************************************************************************************************************************************"""

        self.chrs = {}
        self.nregions = self.count_lines(self.filename)

        if (base == 0):  # Base 0 -> end position -1
            initOffset = 0
            endOffset = -1
        elif (base == 1):  # Base 1-> start position +1
            initOffset = 1
            endOffset = 0
        else:  # Any other value...
            initOffset = 0  # As it is...
            endOffset = 0

        # Each line is parsed to load the region it contains
        fd = open(self.filename)
        for i, line in enumerate(fd):
            parts = line.split('\t')
            # No coordinates for the chr in this line were already loaded
            if (parts[0] not in self.chrs):
                self.chrs[parts[0]] = []

            # Current coordinates are added to current chromosome. >>>> WARNING: check coordinates accoridng to "base" argument <<<<<<<<<<<<
            self.chrs[parts[0]].append((int(parts[1]) + initOffset, int(parts[2]) + endOffset))
    #			pbar.update(i+1)
    #		fd.close()
    #		pbar.finish()

    def checkformat(self):
        """************************************************************************************************************************************************************
        Task: checks the format of the bed file. The only requirements checked are that each line presents at least 3 tab separated columns, the
            two on the right must present integer values indicating the start/end position respectively. Right value must be greater than the
            left value.
        Outputs:
            err: string containing the detected error. Empty string in case of a correct format.
        ************************************************************************************************************************************************************"""

        fd = open(self.filename)

        line = fd.readline()
        fields = line.split('\t')
        lc = 1
        err = ''

        # Checks that the two columns on the right contain integer values
        try:
            # Parses each line and checks that there are at least 3 fields, the two on the right containing integer values and being the right one
            # greater than the left one
            while line != '' and len(fields) > 2 and int(fields[1]) < int(fields[2]):
                lc += 1
                line = fd.readline()
                fields = line.split('\t')
        except ValueError:
            err += 'Incorrect start/end values at line ' + str(lc) + '\n'
            err += 'Start/End coordinates must be indicated with integer values. The right value must be greater than the left value.\n'
            err += 'Line found: ' + line
            fd.close()

            return err

        # If it get to this point means that either the file ended or there is a line with less than 3 fields
        if line != '':
            err += 'Incorrect line format at line ' + str(lc) + '\n'
            err += 'At least three columns are expected in each line\n'
            err += 'The right value must be greater than the left value.\n'
            err += 'Line found: ' + line
            fd.close()

        return err

    def my_sort_bed(self, newbedfilename=None, tmpdir=None):
            """************************************************************************************************************************************************************
            JPFLORIDO
            Task: Sort a BED file by chromosome and start position. This function does not use bedtools
            Input: BED file (self)
            Output:
                An object with the sorted BED file
            ************************************************************************************************************************************************************"""

            global TMP

            if tmpdir is not None:
                TMP = tmpdir

            if newbedfilename is None:
                pid = str(os.getpid())
                newbedfilename = TMP + '/' + pid + os.path.basename(self.filename.replace('.bed', '_sorted.bed'))

            # Load bed file (chromosome, start and end information are stored -> standard bed format)
            self.load_custom(-1)  # as it is....
            # Sort it
            # Passes through the regions in each chromosome
            fdw = open(newbedfilename, 'w')
            for chr in self.chrs:
                # Sorts regions in current chromosome by first coordinate (and second coordinate in case of same first coordinate)
                #			print 'Chr '+chr
                #			print '	Sorting regions...'
                self.chrs[chr].sort()

                # Store it
                curr = 0
                while curr < len(self.chrs[chr]):
                    fdw.write(
                        str(chr) + '\t' + str(self.chrs[chr][curr][0]) + '\t' + str(self.chrs[chr][curr][1]) + '\n')
                    curr += 1

            fdw.close()

            # Return new bed file
            return bed_file(newbedfilename)

    def size(self):
        """*************************************************************************************************************
        Task: calculates the number of bases covered by regions in this bed appropriately merging overlaps.
        Output:
        nbases: integer containing the number of bases covered by regions in this bed.
         Overlapping regions are merged to count each "overlapped" base just once.
        *************************************************************************************************************"""

        # Check whether self.chrs is already loaded
        if (self.chrs == None):
            self.load()
        nbases = 0

        return self.covered_bases()

    def non_overlapping_exons(self, baseCodification, outputFile=None, tmpdir=None):
        """*************************************************************************************************************
        JPFLORIDO
        Task: Get exons of a given bed file removing overlapping areas
        Inputs:
            self: bed file
            baseCodification: whether exons are in "real" base 0, 1 or as it is....
        Output:
            A set of tuples for each exon: chromosome, exon begin position, exon end position
        Requirements: 	WARNING: BED FILE MUST BE SORTED BEFORE
        *************************************************************************************************************"""

        global TMP

        if (tmpdir != None):
            TMP = tmpdir

        chromosomes = []
        start_positions = []
        end_positions = []

        if (self.chrs == None):
            if (baseCodification == 1):
                self.load_custom(1)  # "Real" Base 1
            elif (baseCodification == 0):
                self.load_custom(0)  # "Real" Base 0
            else:
                self.load_custom(-1)  # Base 0 as it is described for standard BED format

        # Passes through the regions in each chromosome
        for chr in self.chrs:
            # Sorts regions in current chromosome by first coordinate (and second coordinate in case of same first coordinate)
            #			print 'Chr '+chr
            #			print '	Sorting regions...'
            self.chrs[chr].sort()

            # A progress bar is initialized
            #			print 'Ignoring overlaps...'

            #			widgets = ['Ignoring overlaps: ', progressbar.Percentage(), ' ',
            #						progressbar.Bar(marker=progressbar.RotatingMarker()), ' ', progressbar.ETA()]
            #			pbar = progressbar.ProgressBar(widgets=widgets, maxval=len(self.chrs[chr])).start()

            # Passes through each tuple (start,end) and processes its overlaps
            curr = 0
            while (curr < len(self.chrs[chr])):

                # "i" is the index of the next region to compare with curr
                # start: integer containing the starting position of the next set of overlapped regions
                i = curr + 1
                start = self.chrs[chr][curr][0]
                entered = 0

                if (chr == '20'):
                    a = 1
                # Check whether starting position of next region falls within current region
                while (i < len(self.chrs[chr]) and self.chrs[chr][curr][0] <= self.chrs[chr][i][0] and
                       self.chrs[chr][i][0] <= self.chrs[chr][curr][1]):
                    if (self.chrs[chr][curr][0] != self.chrs[chr][i][0] or self.chrs[chr][curr][1] != self.chrs[chr][i][
                        1]):
                        a = 1
                    # curr region is always the one with the higher ending position. If the ending point of "i" region is higher than the one of "curr" region, update
                    # curr and make "i" be the next region
                    if (self.chrs[chr][curr][1] < self.chrs[chr][i][1]):
                        curr = i
                        i = curr + 1
                    # otherwise  process the next overlapping region and compare it with curr
                    else:
                        i += 1

                    entered += 1

                if (entered > 1):
                    a = 1
                # The number of bases is calculated as the ending position of the region with the highest ending coordinate within this set of overlapping regions, minus
                # the starting coordinate of this whole set of overlapped regions.
                # curr is updated to "i", which is the first region that does not overlap with curr according to the "while" conditions
                chromosomes.append(chr)
                start_positions.append(start)
                end_positions.append(self.chrs[chr][curr][
                                         1])  # The ending coordinate depends on the base codification (baseCodification argument)
                curr = i
                #				pbar.update(curr)

        #		pbar.finish()

        # Write new BED file
        pid = str(os.getpid())

        if (outputFile == None):
            outputFile = TMP + '/' + pid + os.path.basename(self.filename.replace('.bed', '_noOverlapping.bed'))

        fdw = open(outputFile, 'w')
        for i, currentChromosome in enumerate(chromosomes):
            fdw.write(str(chromosomes[i]) + '\t' + str(start_positions[i]) + '\t' + str(end_positions[i]) + '\n')
        fdw.close()

        return bed_file(outputFile)

    def covered_bases(self):
        """************************************************************************************************************************************************************
        Task: calculates the number of bases covered by regions in this bed appropriately merging overlaps. It is required that self.chrs is already loaded by
            calling self.load()
        Output:
            nbases: integer containing the number of bases covered by regions in this bed. Overlapping regions are merged to count each "overlapped" base just once.
        ************************************************************************************************************************************************************"""

        nbases = 0

        # Passes through the regions in each chromosome
        for chr in self.chrs:
            # Sorts regions in current chromosome by first coordinate (and second coordinate in case of same first coordinate)
            #			print 'Chr '+chr
            #			print '	Sorting regions...'
            self.chrs[chr].sort()

            # Passes through each tuple (start,end) and processes its overlaps
            curr = 0
            while (curr < len(self.chrs[chr])):

                # "i" is the index of the next region to compare with curr
                # start: integer containing the starting position of the next set of overlapped regions
                i = curr + 1
                start = self.chrs[chr][curr][0]
                entered = 0

                if (chr == '20'):
                    a = 1
                # Check whether starting position of next region falls within current region
                while (i < len(self.chrs[chr]) and self.chrs[chr][curr][0] <= self.chrs[chr][i][0] and
                       self.chrs[chr][i][0] <= self.chrs[chr][curr][1]):
                    if (self.chrs[chr][curr][0] is not self.chrs[chr][i][0] or self.chrs[chr][curr][1] is not
                            self.chrs[chr][i][1]):
                        a = 1
                    # curr region is always the one with the higher ending position. If the ending point of "i" region is higher than the one of "curr" region, update
                    # curr and make "i" be the next region
                    if (self.chrs[chr][curr][1] < self.chrs[chr][i][1]):
                        curr = i
                        i = curr + 1
                    # otherwise  process the next overlapping region and compare it with curr
                    else:
                        i += 1

                    entered += 1

                if (entered > 1):
                    a = 1
                # The number of bases is calculated as the ending position of the region with the highest ending coordinate within this set of overlapping regions, minus
                # the starting coordinate of this whole set of overlapped regions.
                # curr is updated to "i", which is the first region that does not overlap with curr according to the "while" conditions
                nbases += (self.chrs[chr][curr][1] - start + 1)
                curr = i

        return nbases

    def load(self):
        """************************************************************************************************************************************************************
        Task: loads data into self.chrs and self.nregions.
        Output: self.chrs and self.nregions are modified.
            self.chrs: dictionary. Each key represents a chromosome. Values are lists of tuples (start,end) indicating each of the regions in the chromosome.
            self.nregions = self.count_lines(self.filename)
        ************************************************************************************************************************************************************"""

        # self.chrs: dictionary. Each key represents a chromosome. Values are lists of tuples (start,end) indicating each of the regions in the chromosome.
        #	 >>>>>>> WARNING: ending coordinate is also transformed to base 0!!!! <<<<<<<
        # self.nregions: total number of regions in the bed file.
        self.chrs = {}
        self.nregions = self.count_lines(self.filename)

        # Each line is parsed to load the region it contains
        fd = open(self.filename)
        for i, line in enumerate(fd):
            parts = line.split('\t')
            # No coordinates for the chr in this line were already loaded
            if (parts[0] not in self.chrs):
                self.chrs[parts[0]] = []

            # Current coordinates are added to current chromosome. >>>> WARNING: ending coordinate is also transformed to base 0!!! <<<<<
            self.chrs[parts[0]].append((int(parts[1]), int(parts[2]) - 1))

        fd.close()
        print('Done')

