from ngscat2.misc.bedgraph_file import bedgraph_file
import os

class Report():
    def __init__(self, outdir):
        self.outdir = outdir

    def report(self, offset, coverageThreshold, target):
        '''Generate Nocoverage.txt of regions within targets that have Zero coverage.
        Input: Coveragefile object, outdir'''
        bedgraphlist = [filename for filename in os.listdir(self.outdir) if filename[0]== "."]
        for bedgraph in bedgraphlist:
            bedgraph_object = bedgraph_file(bedgraph)
            bedgraph_object.getOffTarget(offset, coverageThreshold, target, self.outdir)

