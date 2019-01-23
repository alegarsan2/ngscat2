

class DepthPerPositionProcessor():
    def __init__(self, npoints = 50):
        self.coveragefiles = []
        self.npoints = npoints

    def process(self, coveragefiles):
        return coveragefiles, self.npoints