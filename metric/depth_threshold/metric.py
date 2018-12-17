

class DepthThresProcessor():
    def __init__(self, coveragethreshold, warnthreshold):
        self.results = []
        self.coveragethreshold = coveragethreshold
        self.warnthreshold = warnthreshold

    def process(self, coveragefiles, callback):
        ntotal_positions = [0] * len(coveragefiles)
        covered_positions_per_depth = [[0 for x in range(len(self.coveragethreshold))] for y in range(len(coveragefiles))]
        covered_position = []  # list containing dictionary number of position with more depth than threshold
        perc_covered_position = []  # percentage of covered position
        perc_total_covered = []
        results = []
        target_coverage_status = []

        # calculation of number of bed position and different coverages within thresholds
        for i, coveragefile in enumerate(coveragefiles):
            ntotal = 0

            for current_coverage in coveragefile.coverages:
                ntotal = ntotal + 1
                for j, cov in enumerate(self.coveragethreshold):
                    if current_coverage >= cov:
                        covered_positions_per_depth[i][j] += 1
            ntotal_positions[i] = ntotal

            covered_position.append(
                {'>=' + str(cov) + 'x': covered_positions_per_depth[i][indx] for indx, cov in
                 enumerate(self.coveragethreshold)})
            perc_covered_position.append(
                {key: (value * 100 / ntotal_positions[i]) for key, value in covered_position[i].items()})
            perc_total_covered.append(perc_covered_position[i]['>=1x'])
            target_coverage_status.append(True if perc_total_covered[i] >= self.warnthreshold else False)

        for i in range(len(coveragefiles)):
            results.append(dict(
                [('bamfilename', coveragefiles[i].name.decode('utf-8').split("/")[-1]),
                 ('ntotalposition', ntotal_positions[i]),
                 ('perctotalcovered', perc_total_covered[i]),
                 ('coveredposition', covered_position[i]),
                 ('perccoveredposition', perc_covered_position[i]),
                 ('targetstatus', 'OK' if target_coverage_status[i] else 'Not OK')])
            )

        callback(coveragefiles, results)

