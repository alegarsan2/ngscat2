


class Region():
    def __init__(self):
        self.start = None
        self.end = None
        self.covStartIndex = None
        self.covEndIndex = None
        self.mean = None
        self.std = None
        self.zeropos = None

    def iterateOverCoverages(self, processCoverages):
        processCoverages(self, coverage)


class Chromosome():
    def __init__(self, name):
        self.name = name
        self.regions = []

class Coveragefile():
    def __init__(self,name):
        self.name = name
        self.chromosomes = []
        self.coverages = None
        self.bedfilename = None
    def iterateOverRegions(self, processRegions):
        for chromosome in self.chromosomes:
            for region in chromosome.regions:
                processRegions(chromosome, region)

    def getChromosome(self, chromosomeName):
        for chrom in self.chromosomes:
            if chrom.name == chromosomeName:
                return chrom



    # def iterateOverCoverages(self,processCoverages):
    #     for chromosome in self.chromosomes:
    #         for region in chromosome.regions:
    #             for coverage in region.coverages:
    #                 processCoverages(chromosome,region,coverage)



## Esta parte la tendria que meter antes de los calculos de region.
class RegionFixer:
    def __init__(self, coverages):
        self.coverages = coverages
    def process(self, chromosome, region):
        region.coverages = self.coverages[region.covStartIndex:region.covEndIndex]
        del region.covStartIndex
        del region.covEndIndex
# processor = RegionFixer(coverages)
# coverages.iterateOverRegions(processor.process)