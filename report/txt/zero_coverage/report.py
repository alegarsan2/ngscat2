
from metric.zero_coverage.metric import RegionsWithZeroesProcessor
class Report():
    def __init__(self, outdir):
        self.outdir = outdir

    def report(self, coveragefiles):
        '''Generate Nocoverage.txt of regions within targets that have Zero coverage.
        Input: Coveragefile object, outdir '''
        for coveragefile in coveragefiles:
            zerosProcessor = RegionsWithZeroesProcessor(coveragefile, self.outdir + "NoCoverage.txt")
            coveragefile.iterateOverRegions(zerosProcessor.process)
            zerosProcessor.close()

# from 'report/html/depth_stdev/' import Reporter
# from depth_stdev_computer
#
# reporter = Report('/outdir')
# depth_stdev_computer(coveragefiles, reporter.report)