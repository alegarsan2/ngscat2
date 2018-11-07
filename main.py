def main():
    # bam = bam_file('/home/agarcia/PycharmProjects/ngscat/example2.bam')
    #
    # [bamlist, nread, tread, onperchr, totalperchr, enrichment, percontarget, outdir, legend, onduplicates,
    # offduplicates] = bam.reads_on_target('/home/agarcia/PycharmProjects/ngscat/seqcap.example2.bed',
    #                                      '/home/agarcia/PycharmProjects/ngscat/')

    # npos, depth =target_coverage.precalculated_target_coverage("/home/agarcia/PycharmProjects/ngscat/coveragefile", [1,3,5,200])

    # x = "hola"
    # print("hola")

    import pysam

    samfile = pysam.AlignmentFile("/home/agarcia/PycharmProjects/ngscat/example2.bam")
    for pileupcolumn in samfile.pileup("chr2", 43333, 49000):
        print("\ncoverage at base %s = %s" %
              (pileupcolumn.pos, pileupcolumn.n))
        for pileupread in pileupcolumn.pileups:
            if not pileupread.is_del and not pileupread.is_refskip:
                # query position is None if is_del or is_refskip is set.
                print('\tbase in read %s = %s' %
                      (pileupread.alignment.query_name,
                       pileupread.alignment.query_sequence[pileupread.query_position]))


if __name__ == '__main__':
    main()