import numpy as np

class StdevIterProcessor():
    #Std exon report object, inputs: coveragefiles list of coverages objects, and warnthreshold(if mean of normalized
    #distribution is less than x)
    def __init__(self):
        self.stdlist = []

    def process(self, chromosome, region):
        #The normalization is std of a region divided with the mean of the region.
        if region.mean > 0:
            self.stdlist.append(region.std/region.mean)


class StdevProcessor:
    def __init__(self, warnthreshold= 0.3):
        self.warnthreshold = warnthreshold
        self.stdlists = []
    def process(self, coveragefiles, callback):
        region_stddistribution_result = {}
        region_stddistribution_result['warnthreshold'] = self.warnthreshold
        region_stddistribution_result['results'] = []
        for coverage in coveragefiles:
            iterator = StdevIterProcessor()
            coverage.iterateOverRegions(iterator.process)
            self.stdlists.append(iterator.stdlist)
            result = self.calculate_results(coverage, iterator.stdlist, self.warnthreshold)
            region_stddistribution_result['results'].append(result)

        callback(self.stdlists, region_stddistribution_result)

    def calculate_results(self, coverage, stdlist, warnthreshold):
        histlist = []
        percentile = []
        median = None
        maximum = None
        minimum = None
        mean = None

        number_std, bin_edges = np.histogram(stdlist, bins= np.arange(0, 1, 0.007))
        bin_edges = bin_edges.tolist()
        width = []
        xaxis = []
        for i in range(len(bin_edges)-1):
            xaxis.append((bin_edges[i]+ bin_edges[i +1])/2)
            width.append(bin_edges[1]- bin_edges[0])

        histlist = dict([('numberstd', number_std.tolist()),
                      ('binedges', bin_edges),
                      ('stdpos', [round(x,2) for x in xaxis]),
                      ('width', width)])

        percentile = (dict(
                [('Q1', np.percentile(stdlist, 25)),
                ('Q2', np.percentile(stdlist, 50)),
                ('Q3', np.percentile(stdlist, 75))
            ]))
        median = np.median(stdlist)
        maximum = np.max(stdlist)
        minimum = np.min(stdlist)
        mean = np.mean(stdlist)

        return dict(
                [('bamfilename', coverage.name.decode('utf-8')),
                 #('legend', (legend[i] if legend is not None else bamlist[i].filename.decode('utf-8').split('/')[-1])),
                 ('histdata', histlist),
                 ('percentile',percentile),
                 ('max', float(maximum)),
                 ('min', float(minimum)),
                 ('mean', mean),
                 ('median', median),
                 ('status', 'ok' if mean >= warnthreshold else 'warning')])
