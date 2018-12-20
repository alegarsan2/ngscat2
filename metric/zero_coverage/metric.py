


class RegionsWithZeroesProcessor():
    #Enter attribute coveragelist of coverage object

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
        self.file.close()
