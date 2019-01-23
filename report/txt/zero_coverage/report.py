
class Report():
    def __init__(self, outdir):
        self.outdir = outdir

    def report(self, coveragefiles):
        '''Generate Nocoverage.txt of regions within targets that have Zero coverage.
        Input: Coveragefile object, outdir '''
        for coveragefile in coveragefiles:
            zerosProcessor = ZeroesProcessor(coveragefile, self.outdir + "/data/NoCoverage.txt")
            coveragefile.iterateOverRegions(zerosProcessor.process)
            zerosProcessor.file.close()

class ZeroesProcessor():
    def __init__(self, coveragefile, outdir):
        self.file = open(outdir, "w")
        self.coverages = coveragefile.coverages

    def process(self, chromosome, region):
        #Write the zone of zeroes within the region.
        zeroregion = []
        initzero = None
        rinit = region.start
        regionlen = region.end - region.start
        regionidx = region.covStartIndex
        regioni = self.coverages[region.covStartIndex:region.covEndIndex]
        for idx, basecoverage in enumerate(regioni):
            if idx != regionlen:
                if basecoverage == 0:
                    if initzero is None:
                        initzero = idx
                else:
                    if initzero is not None:
                        self.file.write("\t".join(map(str,[chromosome.name, initzero + region.start, idx + region.start])) + "\n")
                        initzero = None
            else:
                if basecoverage == 0 and initzero is not None:
                    self.file.write("\t".join(map(str,[chromosome.name, initzero + region.start, idx + region.start]))+ "\n")
                else:
                    self.file.write("\t".join(map(str,[chromosome.name, idx + region.start, idx + region.start]))+ "\n")

# from 'report/html/depth_stdev/' import Reporter
# from depth_stdev_computer
#
# reporter = Report('/outdir')
# depth_stdev_computer(coveragefiles, reporter.report)