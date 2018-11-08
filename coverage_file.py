
import numpy as np

class coverage_file():

    def __init__(self, filename):

        lineCount = 0
        with open(filename, "r") as f:
            for line in f.readlines():
                lineCount = lineCount + 1

        chromindex = []
        data = np.empty(shape=(lineCount, 4),dtype= int)
        idx = 0
        with open(filename, "r") as f:
            for idx,line in enumerate(f.readlines()):
                part = line.rstrip().split()
                if part[0] not in chromindex:
                    chromindex.append(part[0])

                partindex = chromindex.index(part[0])
                row = [partindex, part[1],part[2],part[3]]

                data[idx,:] = row

        self.data = data
        self.chromindex = chromindex


    def getInit(self):
        return(self.data[:,1])

    def getEnd(self):
        return(self.data[:,2])

    def getCov(self):
        return(self.data[:,3])

    def getChrom(self):
        chromlist = [self.chromindex[i] for i in self.data[:,0]]

        return(chromlist)
