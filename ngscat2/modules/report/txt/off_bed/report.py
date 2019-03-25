from ngscat2.misc.bedgraph_file import bedgraph_file
import os
import subprocess
class Report():
    def __init__(self, outdir):
        self.outdir = outdir + "/data"

    def report(self, offset, coverageThreshold, target, annotation = None):
        '''Generate Nocoverage.txt of regions within targets that have Zero coverage.
        Input: Coveragefile object, outdir '''
        bedgraphlist = [filename for filename in os.listdir(self.outdir) if ".bed" in filename]

        for bedgraph in bedgraphlist:
            bedgraph_object = bedgraph_file(bedgraph)
            bedgraph_object.getOffTarget(offset, coverageThreshold, target, self.outdir)

        if annotation is not None:
            bedlist = [filename for filename in os.listdir(self.outdir) if ".bed" in filename]
            for bedname in bedlist:
                directory = os.path.join(self.outdir, bedname)
                directory2 = os.path.join(self.outdir, bedname.replace(".bed", "annotated.bed"))
                print("bedmap --echo --delim \'\\t\' --echo-map-id-uniq " + os.path.join(self.outdir, bedname) \
                            + " " + annotation + " > " + os.path.join(self.outdir, bedname.replace(".bed", ".annotated.bed")))

                subprocess.run("bedmap --echo --delim \'\\t\' --echo-map-id-uniq " + os.path.join(self.outdir, bedname) \
                            + " " + annotation + " > " + os.path.join(self.outdir, bedname.replace(".bed", ".annotated.bed")), shell = True)
                os.remove(os.path.join(self.outdir, bedname))
                # subprocess.run("bedmap --echo --delim \'\\t\' --echo-map-id-uniq " + os.path.join(self.outdir, bedname) \
                #                             + " " + annotation, shell=True)
