import os
import subprocess

class Report():
    def __init__(self, outdir):
        self.outdir = outdir + '/data'

    def report(self, results):
        '''Generate Nocoverage.txt of regions within targets that have Zero coverage.
        Input: Coveragefile object, outdir '''
        coveragefiles = results[0]
        annotation = results[1]
        for coveragefile in coveragefiles:
            zerosProcessor = ZeroesProcessor(coveragefile, self.outdir + "/NoCoverage.txt")
            coveragefile.iterateOverRegions(zerosProcessor.process)
            zerosProcessor.file.close()

        '''Annotation of NoCoverageBed'''
        if annotation is not None:
            bedname = "NoCoverage.txt"

            subprocess.run("bedmap --echo --delim \'\\t\' --echo-map-id-uniq " + os.path.join(self.outdir, bedname) \
                            + " " + annotation + " > " + os.path.join(self.outdir, bedname.replace(".txt", ".annotated.txt")), shell = True)
            os.remove(os.path.join(self.outdir, bedname))


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
                        self.file.write("\t".join(map(str,[chromosome.name, initzero + region.start, idx-1 + region.start])) + "\n")
                        initzero = None
            else:
                if basecoverage == 0:
                    if initzero is not None: #Current 0 region
                        self.file.write("\t".join(map(str,[chromosome.name, initzero + region.start, idx + region.start]))+ "\n")
                    else: #Last base with 0 coverage
                        self.file.write("\t".join(map(str,[chromosome.name, idx + region.start, idx + region.start]))+ "\n")

# from 'report/html/depth_stdev/' import Reporter
# from depth_stdev_computer
#
# reporter = Report('/outdir')
# depth_stdev_computer(coveragefiles, reporter.report)