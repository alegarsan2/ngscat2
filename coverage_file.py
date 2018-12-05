


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

        # ... build the object
        fixer = RegionFixer(self.coverages)
        self.iterateOverRegions(fixer.process)

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
class RegionFixer():
    def __init__(self, coveragefile):
        self.coverages = coveragefile

    def process(self, chromosome, region):
        region.basecoverage = self.coverages[region.covStartIndex:region.covEndIndex]
        del region.covStartIndex
        del region.covEndIndex

    def iterateOverRegions(self, processRegions):
        for chromosome in self.chromosomes:
            for region in chromosome.regions:
                processRegions(chromosome, region)

# class ZeroProcessor:
#     def __init__(self, coveragefilecoverages):
#         self.coverages = coveragefilecoverages
#
#     def process(self, chromosome, region):
#         zeroregion = []
#         initzero = None
#         rinit = region.start
#         regionlen = region.end - region.start
#         regionidx = region.covStartIndex
#         regioni = self.coverages[region.covStartIndex:region.covEndIndex]
#         for idx, basecoverage in enumerate(regioni):
#             if idx != regionlen:
#                 if basecoverage == 0:
#                     if initzero is None:
#                         initzero = idx
#                 else:
#                     if initzero is not None:
#                         zeroregion.append([chromosome.name, initzero + region.start, idx + region.start])
#                         initzero = None
#             else:
#                 if basecoverage == 0 and initzero is not None:
#                     zeroregion.append([chromosome.name, initzero + region.start, idx +  region.start])
#                 else:
#                     zeroregion.append([chromosome.name, idx + region.start, idx +  region.start])
#
#         return zeroregion
#         # quedarme con las regiones
#
#
# # zerosProcessor = ZeroProcessor(coveragefile.coverages)
# # coveragefile.iterateOverRegions(zerosProcessor.process)
#
#
# # processor = RegionFixer(coverages)
# # coverages.iterateOverRegions(processor.process)