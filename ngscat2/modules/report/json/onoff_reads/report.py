import json

class Report():
    def __init__(self, outdir):
        self.outdir = outdir

    def report(self, reads_on_results):

        with open(self.outdir + '/data/reads_on_results.json', 'w') as outfile:
            json.dump(reads_on_results, outfile)


