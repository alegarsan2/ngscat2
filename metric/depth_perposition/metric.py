

class DepthPerPositionProcessor():
    def __init__(self,npoints):
        self.coveragefiles = []
        self.npoints = npoints
    def process(self, coveragefiles, callback):
        callback(coveragefiles, self.npoints)