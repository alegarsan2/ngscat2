

class DepthPerPositionProcessor():
    def __init__(self, npoints = 50):
        self.coveragefiles = []
        self.npoints = npoints
    def process(self, coveragefiles, callback):
        callback(coveragefiles, self.npoints)