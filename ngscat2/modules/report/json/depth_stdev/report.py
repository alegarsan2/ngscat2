import json

class Report():
    def __init__(self, outdir):
        self.outdir = outdir

    def report(self, stdlist, region_stddistribution_result):
        '''Inputs: List of stdreport objects
             Output: std_wexons.json'''

        with open(self.outdir + '/data/std_wexons.json', 'w') as outfile:
            json.dump(region_stddistribution_result, outfile)
