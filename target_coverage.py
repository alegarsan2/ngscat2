def target_coverage(coveragefile, coveragelist):
    """************************************************************************************************************************************************************
    Task: draws statistics about the percentage of covered exons and transcripts at different coverage levels. A transcript is considered to be covered when
        at least the 90% of its positions present a coverage greater than the threshold.
    Inputs:
        filelist: list of strings indicating those files to be processed. For a file format example see
            /home/javi/MGP/capture_methods/data/coverage/GU_20120719_FC1_6L1S_AL_01_3376_BC1_AA_F3.filtered.singleHits.realigned.recalibrated.bam.coverage
        coveragelist: list of values with coverage thresholds to use.
        graph_legend: list of descriptions describing each of the files that will be processed. These descriptions will form the legend of the bar plot.
        dirout: string containing the full path to the directory where data will be saved.
    Output: a summary .xls file and two bar plots depicting coverage vs. %covered-positions and coverage vs. #covered transcripts. Figures will be saved as
        <dirout>/coverage_summary.xls, <dirout>/covered_positions.png and <dirout>/covered_transcripts.png
    ************************************************************************************************************************************************************"""

    ntotal_positions = count_lines(coveragefile)

    # covered_positions_per_depth: list of integers. There will be a position for each coverage threshold. Each value will be the count of positions
    #     covered for the corresponding threshold.
    # ccds_counts: dictionary. Keys are transcript ids. values are lists of two elements. The first element of this list will contain the length of the
    #     transcript. The second element will be a list of integers with as many positions as coverage thresholds, being each value the count of positions
    #     covered for the corresponding threshold.
    covered_positions_per_depth = [0 for i in range(len(coveragelist))]

    # A progress bar is initialized
    print
    'Counting covered positions...'
    #    widgets = ['Counting: ', progressbar.Percentage(), ' ',
    #                progressbar.Bar(marker=progressbar.RotatingMarker()), ' ', progressbar.ETA()]
    #    pbar = progressbar.ProgressBar(widgets=widgets, maxval=ntotal_positions).start()

    # Each line contains the coverage for a given position
    fd = file(coveragefile)
    for k, line in enumerate(fd):
        parts = line.split('\t')

        # Check whether the coverage is over each threshold
        for i, cov in enumerate(coveragelist):
            current_coverage = string.atof(parts[-1])
            # In case coverage is over the threshold, add 1 to the global number of covered positions and to the counts of the current transcript
            if (current_coverage >= cov):
                covered_positions_per_depth[i] += 1

            #        pbar.update(k+1)

    #    pbar.finish()
    fd.close()