
import numpy as np

class coverage_file():
    def __init__(self, filename):
        # with open(filename,"r") as f:
        #     chromindex = []
        #     data = np.empty(shape=(0,4))
        #
        #     for line in f.readlines():
        #         line = line.rstrip()
        #         part = line.split()
        #
        #         if part[0] not in chromindex:
        #             chromindex.append(part[0])
        #
        #         partindex = chromindex.index(part[0])
        #         row = [partindex, part[1],part[2],part[3]]
        #         print(row)
        #         data.append([row],axis= 0)
        #
        # self.data = data
        # self.chromindex = chromindex

        lineCount = 0
        with open(filename, "r") as f:
            for line in f.readlines():
                lineCount = lineCount + 1

        # chromindex = []
        # data = np.empty(shape=(lineCount, 4))
        # idx = 0
        # with open(filename, "r") as f:
        #     for line in f.readlines():
        #         part = line.rstrip().split()
        #         if part[0] not in chromindex:
        #             chromindex.append(part[0])
        #         partindex = chromindex.index(part[0])
        #         row = [partindex, part[1],part[2],part[3]]
        #         # print(row)
        #         for jdx in range(0,4):
        #             data[idx, jdx] = row[jdx]
