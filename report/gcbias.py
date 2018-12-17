import array
import re
import Bio

class GcBiasProcessor():
    def __init__(self,coveragefile, reference):
        self.coveragefile = coveragefile
        self.reference = reference
        self.meanlist = array.array()
        self.gcratiolist = array.array()
        coveragefile.iterateOverChromosomes(process)

    def process(self,chromosome):
        with open(reference, "r") as fastafile:
            for line in fastafile:
                if line[0] == '>':
                    refchrom = re.split(' +', line)[0][1:] #get first column divided by space and then delete the >



    def calculate_results(self):


    def output(self,outdir):
        #Generation

