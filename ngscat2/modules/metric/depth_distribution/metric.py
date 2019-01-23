import numpy as np

class DepthDistrProcessor():
    def __init__(self, bins = 40, warnthreshold=40):
        self.bins = bins
        self.warnthreshold = warnthreshold
        self.results = []

    def process(self, coveragefiles):
        results = []
        histlist = []
        percentile = []
        maximum = []
        minimum = []
        mean = []
        median= []
        zerocov = []

        for i, coveragefile in enumerate(coveragefiles):

            #Non zero indexs
            indnonzero = np.nonzero(coveragefile.coverages)
            #FIXME no existe normalizacion de tamaño de bins, solo toma segun el tamaño de la muestra(como es sin 0) QUitar 0??
            number_read, bin_edges = np.histogram(coveragefile.coverages[indnonzero], bins= self.bins)
            bin_edges = bin_edges.tolist()

            width = []
            xaxis = []
            for i in range(len(bin_edges)-1):
                xaxis.append((bin_edges[i]+ bin_edges[i +1])/2)
                width.append(bin_edges[1]- bin_edges[0])

            histlist.append(dict([('numberread', number_read.tolist()),
                                  ('binedges', bin_edges),
                                  ('coveragepos', [round(x,2) for x in xaxis]),
                                  ('width', width)]))

            percentile.append(dict(
                [('Q1', np.percentile(coveragefile.coverages, 25)),
                ('Q2', np.percentile(coveragefile.coverages, 50)),
                ('Q3', np.percentile(coveragefile.coverages, 75))
                ]))
            median.append(np.median(coveragefile.coverages))
            maximum.append(np.max(coveragefile.coverages))
            minimum.append(np.min(coveragefile.coverages))
            mean.append(coveragefile.coverages.mean())
            zerocov.append(len(coveragefile.coverages) - len(indnonzero))


        for i in range(len(coveragefiles)):
            results.append(dict(
                [('bamfilename', coveragefiles[i].name.decode('utf-8')),
                 ('legend', (coveragefiles[i].name.decode('utf-8').split('/')[-1])),
                 ('histdata', histlist[i]),
                 ('percentile', percentile[i]),
                 ('max', float(maximum[i])),
                 ('min', float(minimum[i])),
                 ('mean', mean[i]),
                 ('median', median[i]),
                 ('zerocov', float(zerocov[i])),
                 ('status', 'ok' if mean[i] >= self.warnthreshold else 'warning'),
                 ('warnthreshold', self.warnthreshold)]
            ))


        # callback(coveragefiles, results)
        return coveragefiles, results