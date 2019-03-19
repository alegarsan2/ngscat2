

class DepthPerPositionProcessor():
    def __init__(self, npoints = 50, warnregionsize=100, warnthreshold=6):
        self.coveragefiles = []
        self.npoints = npoints
        self.warnregionsize = warnregionsize
        self.warnthreshold = warnthreshold

    def process(self, coveragefiles):
        return coveragefiles, self.npoints, self.warnregionsize, self.warnthreshold
