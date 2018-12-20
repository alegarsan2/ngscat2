import array
import re
import Bio

class GcBiasProcessor():
    def __init__(self):
        self.meanlist = []
        self.meanlists = []
        self.gcperclist = []
        #coveragefile.iterateOverChromosomes(self.process)
        #self.chromnamelist = coveragefile.getChromosomeNames()

    def process(self, coveragefiles, reference, callback):
        self.meanlists = [[]*len(coveragefiles)]
        chromnamelist = coveragefiles[0].getChromosomeNames()
        storeSequence = False
        refchrom = ''
        chromosomeSeq = ''
        with open(reference, "r") as fastafile:
            #Reference parsing looking for a chromosome, if it is in list, save it
            for line in fastafile:
                if line[0] == '>':

                    lastchrom = refchrom
                    refchrom = re.split(' +', line)[0][1:].rstrip()  # get first column divided by space and then delete the >
                    print(refchrom)

                    if storeSequence is True:
                        chromosomes =[]
                        #Chromosome objects list generation, in order to parse the reference just once
                        for coveragefile in coveragefiles:
                            chromosomes.append(coveragefile.getChromosome(lastchrom))
                        self.calculate_GC_chromosome(chromosomes, chromosomeSeq)
                        chromosomeSeq = ''
                        storeSequence = False

                    else:
                        if refchrom in chromnamelist:
                            storeSequence = True
                        else:
                            storeSequence = False
                else:
                    if storeSequence is True:
                        chromosomeSeq = chromosomeSeq + line.rstrip()

            #Last chromosome
            if storeSequence is True:
                chromosomes = []
                for coveragefile in coveragefiles:
                    chromosomes.append(coveragefile.getChromosome(lastchrom))
                self.calculate_GC_chromosome(coveragefiles[0].getChromosome(lastchrom), chromosomeSeq)
                chromosomeSeq = ''
                storeSequence = False

            if len(self.gcperclist) == 0:
                print('ERROR: G+C content values can not be calculated. Probably the provided reference file '
                      'does not match with the target file. That is, sequences of regions in the target file ar'
                      'e probably not included within the reference file. Check if chromosome names in the bed are'
                      'in the same format as the reference')

        callback(self.gcperclist, self.meanlists, coveragefiles)



    def calculate_GC_chromosome(self, chromosomes, chromosomeSeq):
        a = len(chromosomeSeq)
        for indx, region in enumerate(chromosomes[0].regions):
            self.gcperclist.append(100*((sum([1 for base in chromosomeSeq[region.start:region.end] if base == 'C' or base == 'G']))
                                        /(region.end - region.start)))

            for jdx, chromosome in enumerate(chromosomes):
                self.meanlists[jdx].append(chromosome.regions[indx].mean)
