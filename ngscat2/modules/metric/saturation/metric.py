import numpy as np

class DepthDistrProcessor():
    def __init__(self,bins, warnthreshold=40):
        self.bins = bins
        self.warnthreshold = warnthreshold
        self.results = []

    def process(self, coveragefiles, callback):

        # ntotal_positions = [0] * len(coveragefiles)
        # covered_positions_per_depth = [[0 for x in range(len(coveragethreshold))] for y in range(len(coveragefiles))]
        # covered_position = []  # list containing dictionary number of position with more depth than threshold
        # perc_covered_position = []  # percentage of covered position
        # perc_total_covered = []
        # results = []
        # target_coverage_status = []
        # # calculation of number of bed position and different coverages within thresholds
        # for i, coveragefile in enumerate(coveragefiles):
        #     ntotal = 0
        #
        #     for current_coverage in coveragefile.coverages:
        #         ntotal = ntotal + 1
        #         for j, cov in enumerate(coveragethreshold):
        #             if current_coverage >= cov:
        #                 covered_positions_per_depth[i][j] += 1
        #     ntotal_positions[i] = ntotal
        #
        #     covered_position.append(
        #         {'>=' + str(cov) + 'x': covered_positions_per_depth[i][indx] for indx, cov in
        #          enumerate(coveragethreshold)})
        #     perc_covered_position.append(
        #         {key: (value * 100 / ntotal_positions[i]) for key, value in covered_position[i].items()})
        #     perc_total_covered.append(perc_covered_position[i]['>=1x'])
        #     target_coverage_status.append(True if perc_total_covered[i] >= warnthreshold else False)
        #
        # for i in range(len(coveragefiles)):
        #     results.append(dict(
        #         [('bamfilename', coveragefiles[i].name.decode('utf-8').split("/")[-1]),
        #          ('ntotalposition', ntotal_positions[i]),
        #          ('perctotalcovered', perc_total_covered[i]),
        #          ('coveredposition', covered_position[i]),
        #          ('perccoveredposition', perc_covered_position[i]),
        #          ('targetstatus', 'OK' if target_coverage_status[i] else 'Not OK')])
        #     )
        # a = print()
        # callback(coveragefiles, results)

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

        for i in range(len(coveragefile)):
            results.append(dict(
                [('bamfilename', coveragefiles[i].filename.decode('utf-8')),
                 ('legend', (coveragefiles[i].decode('utf-8').split('/')[-1])),
                 ('histdata', histlist[i]),
                 ('percentile', percentile[i]),
                 ('max', float(maximum[i])),
                 ('min', float(minimum[i])),
                 ('mean', mean[i]),
                 ('median', median[i]),
                 ('zerocov', float(zerocov[i])),
                 ('status', 'OK' if mean[i] >= self.warnthreshold else 'Not OK')]
            ))


        callback(coveragefiles, results)