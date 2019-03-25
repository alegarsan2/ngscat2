import numpy
import os
from ngscat2.misc import bed_file
from ngscat2.misc import region

class bedgraph_file:

    def __init__(self, _filename):
        self.filename = _filename

    def get_batch(self, fd, size):
        """*************************************************************************************************************
        Task: returns a list of n lines
        Input:
            fd:file handler
            size:~size of batch in bytes
        Output:
            batch: list of lines
            fd:file handler
        *************************************************************************************************************"""
        batch = []
        batch = fd.readlines(size)
        return batch, fd

    def getOffTarget(self, offset, coverageThreshold, target, outfile, tmpdir=None):
        """*************************************************************************************************************
        Task: selects off-target(+offset) regions with a coverage >  coverageThreshold
        Inputs:
            offset: integer indicating the number of bases to extend the target.
            coverageThreshold: integer indicating the coverage threshold to select the region
            target: ROIs bed file
        Ouputs: a new bedgraph file will be created containing selected regions.
        *************************************************************************************************************"""

        # pid = str(os.getpid())
        # tmpbed = tmpdir + '/' + pid + '.extended.bed'

        bed = bed_file.bed_file(target)
        # extendedBed = bed.extendnoref(offset, tmpbed)
        extendedBed = bed.extendnoref(offset)
        sortedBed = extendedBed.my_sort_bed()
        nonOverlappingBed = sortedBed.non_overlapping_exons(-1)  # Base 0, it is a standard BED
        finalBed = nonOverlappingBed.my_sort_bed()  # BED file in base 0
        finalBed.load_custom(-1)  # Load chromosome and positions in base 0
        bed_region = finalBed.get_region()
        bed_index = 0  # index to control bed_region position

        fd = open(outfile + "/" + self.filename)
        header = fd.readline()
        reading = True  # boolean to control while loop
        chr_found = False
        batch_n = 1
        fdw = open(outfile+ "/" + self.filename.replace('.bed', '.off.bed'), 'w')

        while reading:
            batch, fd = self.get_batch(fd, 10000000)
            #            print batch_n
            batch_n = batch_n + 1

            if batch == []:
                reading = False
            else:
                for line in batch:
                    aline = line.replace('\n', '').split(' ')
                    # new region
                    r = region.region(aline[0], aline[1], aline[2], aline[3])
                    search_open = True

                    while search_open:
                        type_overlap = r.overlap_type(bed_region[bed_index])

                        if type_overlap == 0:  # bed region comes before bedgraph region
                            search_open = True

                            if bed_index + 1 < len(bed_region) and (chr_found == False or (
                                    chr_found == True and r.chrom == bed_region[bed_index].chrom)):
                                bed_index = bed_index + 1
                            elif r.value >= coverageThreshold:
                                search_open = False
                                for region_selected in r - bed_region[bed_index]:
                                    fdw.write(str(region_selected))
                            else:
                                search_open = False


                        elif type_overlap == -1:  # bed region comes after bedgraph region
                            search_open = False
                            chr_found = True
                            if r.value >= coverageThreshold:
                                for region_selected in r - bed_region[bed_index]:
                                    fdw.write(str(region_selected))

                        else:
                            search_open = False
                            chr_found = True
                            if r.value >= coverageThreshold:
                                for region_selected in r - bed_region[bed_index]:
                                    fdw.write(str(region_selected))
        fd.close()
