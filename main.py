import optparse
import sys
import os

from multiprocessing import Manager
from multiprocessing import Pool
from multiprocessing import cpu_count


def parse_arguments():
    usage = """	
       	************************************************************************************************************************************************************
       	Task: Assesses capture performance in terms of sensibility, specificity and uniformity of the coverage.
       	Output: An html report will be created at the path indicated with the --out option.
       	************************************************************************************************************************************************************
       	usage: %prog --bams <filename> --bed <filename> --out <path> --extendtarget <nbases> --reference <filename> --saturation <{y,n}> --depthlist <list> --tmp <path> --threads <integer>"""

    parser = optparse.OptionParser(usage)
    parser.add_option("--bams", dest="bams",
                      help="""Required. Comma separated list of bam files (2 maximum). E.g.: --bams /home/user/bam1.sorted.bam,/home/user/bam2.sorted.bam""")
    parser.add_option("--bed", dest="bed",
                      help="""Required. Full path to the bed file containing the target regions.""")
    parser.add_option("--out", dest="out", help="""Required. Full path to the directory where results will be saved.""")
    parser.add_option("--extendtarget", dest="extend",
                      help="""Optional. Integer indicating the number of bases to extend each target region up and down-stream. Default=None.""",
                      default=None)
    parser.add_option("--reference", dest="reference",
                      help="""Optional. String indicating the path to a .fasta file containing the reference chromosomes. Default=None.""",
                      default=None)
    parser.add_option("--saturation", dest="saturation",
                      help="""Optional. {y,n} to indicate whether saturation curve should be calculated. Default=n.""",
                      default='n')
    parser.add_option("--depthlist", dest="depthlist",
                      help="""Optional. Will only be used in case --saturation is "y". Comma separated list of real numbers (do not leave spaces between) indicating the number of millions of reads to simulate for the saturation curve. E.g.: 1,5,10 would indicate 1*10^6, 5*10^6 and 10*10^6. Default=auto.""",
                      default='auto')
    parser.add_option("--coveragethrs", dest="coveragethresholds",
                      help="""Optional. Comma separated list of real numbers (do not leave spaces between) indicating coverage thresholds to be used when calculating percentages of covered bases (first graph in the report). Default=1,5,10,20,30.""",
                      default='1,5,10,20,30')
    parser.add_option("--onefeature", dest="feature",
                      help="""Optional. Use this option if just one of the graphs/statistics should be calculated. String indicating one of the following features:  {'percbases','saturation','specificity','coveragefreq', 'coveragedistr', 'coveragestd', 'gcbias','coveragecorr'}.""",
                      default=None)
    parser.add_option("--tmp", dest="tmp",
                      help="""Optional. String indicating the full path to a temporary directory where temporary files will be created. Default=/tmp/.""",
                      default='/tmp/')
    parser.add_option("--threads", dest="nthreads",
                      help="""Optional. Integer indicating the number of concurrent threads to launch. Default=2.""",
                      default=2)
    return args, None


def check_parameters(options, parser):
    availablefeatures = ['percbases', 'saturation', 'specificity', 'coveragefreq', 'coveragedistr', 'coveragestd',
                         'gcbias', 'coveragecorr']
    textchars = ''.join(map(chr, [7, 8, 9, 10, 12, 13, 27] + range(0x20, 0x100)))
    is_binary_string = lambda bytes: bool(bytes.translate(None, textchars))

    # Check number of arguments
    if len(sys.argv) < 7:
        parser.print_help()
        print('ERROR: --bams, --bed and --out parameters are required.')
        sys.exit(1)

    # Check number of arguments
    if len(sys.argv) > 21:
        parser.print_help()
        print('ERROR: too many parameters. Please, check that there are no spaces between commas within the "depthlist" or "coveragethrs" arguments.')
        sys.exit(1)
    ## Bam
    try:
        bamlist = options.bams.split(',')
        if (len(bamlist) > 2):
            print('ERROR: please make sure that no more than two bam files are provided. Please, input a comma separated list. E.g.: --bams /home/user/bam1.sorted.bam,/home/user/bam2.sorted.bam')
            sys.exit(1)
    except AttributeError:
        print('ERROR: at least one bam file is required. Please, input a comma separated list. E.g.: --bams /home/user/bam1.sorted.bam,/home/user/bam2.sorted.bam')
        sys.exit(1)

    for bam in bamlist:
        if (not (os.path.isfile(bam) or os.path.islink(bam))):
            print('ERROR: ' + bam + ' does not exist.')
            sys.exit(1)

        if (not bam[-4:] == '.bam'):
            print('ERROR: ' + bam + ' must have .bam extension. Please, make sure that the bam file is appropriately formatted.')
            sys.exit(1)

        if (not is_binary_string(open(bam).read(3))):
            print('ERROR: ' + bam + ' must be a binary file. Please, make sure that the bam file is appropriately formatted.')
            sys.exit(1)
    ## Bed
    try:
        if (not (os.path.isfile(options.bed) or os.path.islink(options.bed))):
            print('ERROR: ' + options.bed + ' does not exist.')
            sys.exit(1)
    except AttributeError:
        print('ERROR: the --bed file is a required parameter. Please, provide one bed file indicating target regions to analyze.')
        sys.exit(1)

    err = bed_file.bed_file(options.bed).checkformat()
    if (err is not ''):
        print('ERROR: incorrect bed file format.')
        print('	' + err)
        sys.exit(1)

    ## out parameter
    try:
        if (not (os.path.isdir(os.path.dirname(options.out)) or os.path.islink(os.path.dirname(options.out)))):
            print('ERROR: ' + os.path.dirname(options.out) + ' does not exist.')
            sys.exit(1)
    except AttributeError:
        print('ERROR: the --out parameter is required. Please, provide full path to an existing directory where results can be saved.')
        sys.exit(1)

    if ((os.path.isdir(options.out) or os.path.islink(options.out)) and (
            os.path.isdir(options.out + '/data') or os.path.islink(options.out + '/data')) and len(
            glob.glob(options.out + '/data/*_Ontarget_Coverage.png')) > 0):
        print('WARNING: ' + options.out + ' directory seems to contain previous NGScat results. Saving results of current execution in this directory may cause incorrect report generation.')
        print('Continue with current setting? (y/n)')

        proceed = input().lower()
        while (proceed is not 'y' and proceed is not 'n'):
            proceed = input().lower()
        if (proceed is 'n'):
            sys.exit(1)

    ## reference
    if (options.reference is not None and (not (os.path.isfile(options.reference) or os.path.islink(options.reference)))):
        print('ERROR: ' + options.reference + ' does not exist.')
        sys.exit(1)

    ## saturation
    if (options.saturation is not 'y' and options.saturation is not 'n'):
        print('ERROR: incorrect value for --saturation parameters. Please indicate "y" or "n".')
        sys.exit(1)

    ## number of threads
    try:
        nthreads = int(options.nthreads)
    except ValueError:
        print('ERROR: invalid value for --nthreads option. Please, provide an integer value. Note that the application will launch as many processess as it needs between 1 and nthreads.')
        sys.exit(1)

    ## depth list
    if (options.depthlist is not 'auto'):
        try:
            depthlist = map(float, options.depthlist.split(','))
        except ValueError:
            print('ERROR: invalid values for --depthlist option. Please, provide a comma separated list of values without leaving spaces, e.g.: 1,2,10,20')
            sys.exit(1)

    ## coveragethresholds
    try:
        coveragetrhesholds = map(float, options.coveragethresholds.split(','))
    except ValueError:
        print('ERROR: invalid values for --coveragethrs option. Please, provide a comma separated list of values without leaving spaces, e.g.: 1,2,10,20')
        sys.exit(1)

    ## FIXME problably it is not going to be used

    if (options.feature is not None and options.feature.lower() not in availablefeatures):
        print('ERROR: ' + options.feature + " not available. Please, check that the selected feature is one of the following: 'percbases','saturation','specificity','coveragefreq', 'coveragedistr', 'coveragestd', 'gcbias'")
        sys.exit(1)

    if (not (os.path.isdir(options.tmp) or os.path.islink(options.tmp))):
        print('ERROR: ' + options.tmp + ' does not exist.')
        sys.exit(1)

    return True

def generate_report(args):

    # crear el pool, con nº hilos core-1
    # generar coverfile
    # crear métricas
    # crear main reports
    # crear subreports
    # lanzar todas las métricas
    # esperar a que acabe todo
    # generar informe final
#Creating Space of data to share between threads

mgr = Manager()
ns = mgr.Namespace()


#Bamfile object generation, if not sorted do it and
ns.covlist = coveragefiles

nthreads = if nthreads in ['nthreads'] is not None else cpu_count() - 1

#Maximum number of workers
mainpool = Pool(processes= cpu_count() - 1)


# reporter =
mainpool.apply_async(OnOffReadsProcessor().process, args=(ns.covlist,'/home/agarcia/PycharmProjects/ngscat/talidomida_v2_primary_targets.bed'), callback = reporter.report)






class CompoundReporter:
    def __init__(self, reporters):
        self.reporters = reporters

    def report(self, *args):
        for processor in self.reporters:
            processor(*args)

htmlReport = HtmlReport(...)
jsonReport = JsonReport(...)

# sensitivity

sensitivityHtmlReport = SensitivityHtmlReport(...)
sensitivityJsonReport = SensitivityJsonReport(...)
sensitivityReport = CompoundProcessor(sensitivityHtmlReport.process, sensitivityJsonReport.process)
mainpool.apply_async(SensitivityMetric().process, args=(ns.covlist,'/home/agarcia/PycharmProjects/ngscat/talidomida_v2_primary_targets.bed'), callback = sensitivityReport.process)

# specificity
sensitivityHtmlReport = SensitivityHtmlReport(...)
sensitivityJsonReport = SensitivityJsonReport(...)
sensitivityReport = CompoundProcessor(sensitivityHtmlReport.process, sensitivityJsonReport.process)
mainpool.apply_async(SpecificityMetric().process, args=(ns.covlist,'/home/agarcia/PycharmProjects/ngscat/talidomida_v2_primary_targets.bed'), callback = sensitivityReport.process)

if shouldExecuteMetricX(args):
    ...


...
pool.join()
...

htmlReport.output()
jsonReport.output()








def generate_report(args):

    # crear el pool, con nº hilos core-1
    # generar coverfile
    # crear métricas
    # crear main reports
    # crear subreports
    # lanzar todas las métricas
    # esperar a que acabe todo
    # generar informe final







def main():

    args, error = parse_arguments()
    if error:
        print(error)
        print_usage()
    else:
        generate_report(args)

    ################################################

    #### Options and arguments #####################

    ################################################

    usage = """	
    	************************************************************************************************************************************************************
    	Task: Assesses capture performance in terms of sensibility, specificity and uniformity of the coverage.
    	Output: An html report will be created at the path indicated with the --out option.
    	************************************************************************************************************************************************************
    	usage: %prog --bams <filename> --bed <filename> --out <path> --extendtarget <nbases> --reference <filename> --saturation <{y,n}> --depthlist <list> --tmp <path> --threads <integer>"""


    parser = optparse.OptionParser(usage)
    parser.add_option("--bams", dest="bams",
                      help="""Required. Comma separated list of bam files (2 maximum). E.g.: --bams /home/user/bam1.sorted.bam,/home/user/bam2.sorted.bam""")
    parser.add_option("--bed", dest="bed",
                      help="""Required. Full path to the bed file containing the target regions.""")
    parser.add_option("--out", dest="out", help="""Required. Full path to the directory where results will be saved.""")
    parser.add_option("--extendtarget", dest="extend",
                      help="""Optional. Integer indicating the number of bases to extend each target region up and down-stream. Default=None.""",
                      default=None)
    parser.add_option("--reference", dest="reference",
                      help="""Optional. String indicating the path to a .fasta file containing the reference chromosomes. Default=None.""",
                      default=None)
    parser.add_option("--saturation", dest="saturation",
                      help="""Optional. {y,n} to indicate whether saturation curve should be calculated. Default=n.""",
                      default='n')
    parser.add_option("--depthlist", dest="depthlist",
                      help="""Optional. Will only be used in case --saturation is "y". Comma separated list of real numbers (do not leave spaces between) indicating the number of millions of reads to simulate for the saturation curve. E.g.: 1,5,10 would indicate 1*10^6, 5*10^6 and 10*10^6. Default=auto.""",
                      default='auto')
    parser.add_option("--coveragethrs", dest="coveragethresholds",
                      help="""Optional. Comma separated list of real numbers (do not leave spaces between) indicating coverage thresholds to be used when calculating percentages of covered bases (first graph in the report). Default=1,5,10,20,30.""",
                      default='1,5,10,20,30')
    parser.add_option("--onefeature", dest="feature",
                      help="""Optional. Use this option if just one of the graphs/statistics should be calculated. String indicating one of the following features:  {'percbases','saturation','specificity','coveragefreq', 'coveragedistr', 'coveragestd', 'gcbias','coveragecorr'}.""",
                      default=None)
    parser.add_option("--tmp", dest="tmp",
                      help="""Optional. String indicating the full path to a temporary directory where temporary files will be created. Default=/tmp/.""",
                      default='/tmp/')
    parser.add_option("--threads", dest="nthreads",
                      help="""Optional. Integer indicating the number of concurrent threads to launch. Default=2.""",
                      default=2)

    (options, args) = parser.parse_args()


if __name__ == '__main__':
    main()