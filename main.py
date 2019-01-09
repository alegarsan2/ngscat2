import optparse
import sys
import os
import glob
import json
from multiprocessing import Manager
from multiprocessing import Pool
from multiprocessing import cpu_count
from multiprocessing import reduction

from bam_file import bam_file
import bed_file
from report.html.HtmlReport import HtmlReport

from metric.onoff_reads.metric import OnOffReadsProcessor
from report.html.onoff_reads.report import Report as OnOffHtml
from report.xls.onoff_reads.report import Report as OnOffXls
from report.json.onoff_reads.report import Report as OnOffJson


from metric.depth_threshold.metric import DepthThresProcessor
from report.html.depth_threshold.report import Report as ThresholdHtml
from report.json.depth_threshold.report import Report as ThresholdJson
from report.xls.depth_threshold.report import Report as ThresholdXls


from metric.depth_distribution.metric import DepthDistrProcessor
from report.html.depth_distribution.report import Report as DistributionHtml
from report.json.depth_distribution.report import Report as DistributionJson
from report.xls.depth_distribution.report import Report as DistributionXls

from metric.depth_perposition.metric import DepthPerPositionProcessor
from report.html.depth_perposition.report import Report as PerpositionHtml

from metric.depth_stdev.metric import StdevProcessor
from report.html.depth_stdev.report import Report as StdHtml
from report.xls.depth_stdev.report import Report as StdXls
from report.json.depth_stdev.report import Report as StdJson

from metric.zero_coverage.metric import RegionsWithZeroesProcessor
from report.txt.zero_coverage.report import Report as ZeroCoverageTxt

from metric.gc_bias.metric import GcBiasProcessor
from report.html.gc_bias.report import Report as GcBiasHtml


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
    parser.add_option("--out", dest="outdir", help="""Required. Full path to the directory where results will be saved.""")
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
                      help="""Optional. Integer indicating the number of concurrent threads to launch. Default=cpu_count() - 1.""",
                      default= cpu_count() - 1)
    #This is used in we want to return an arguments dictionary instead parser object
    #option_dict = {}
    #option_dict = vars(options)

    options, args = parser.parse_args()





    return options, args,parser


def check_parameters(options, parser):
    availablefeatures = ['percbases', 'saturation', 'specificity', 'coveragefreq', 'coveragedistr', 'coveragestd',
                         'gcbias', 'coveragecorr']
    #textchars = ''.join(map(chr, [7, 8, 9, 10, 12, 13, 27] + range(0x20, 0x100)))
    #is_binary_string = lambda bytes: bool(bytes.translate(None, textchars))

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

        # if (not is_binary_string(open(bam).read(3))):
        #     print('ERROR: ' + bam + ' must be a binary file. Please, make sure that the bam file is appropriately formatted.')
        #     sys.exit(1)
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
        if (not (os.path.isdir(os.path.dirname(options.outdir)) or os.path.islink(os.path.dirname(options.outdir)))):
            print('ERROR: ' + os.path.dirname(options.outdir) + ' does not exist.')
            sys.exit(1)
    except AttributeError:
        print('ERROR: the --out parameter is required. Please, provide full path to an existing directory where results can be saved.')
        sys.exit(1)

    if ((os.path.isdir(options.outdir) or os.path.islink(options.outdir)) and (
            os.path.isdir(options.outdir + '/data') or os.path.islink(options.outdir + '/data')) and len(
            glob.glob(options.outdir + '/data/*_Ontarget_Coverage.png')) > 0):
        print('WARNING: ' + options.outdir + ' directory seems to contain previous NGScat results. Saving results of current execution in this directory may cause incorrect report generation.')
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


    #Class wrapper of different kinds of results for instance (onoffplots, onoffjson, onoffxls)
class CompoundReporter:
    def __init__(self, reporters):
        self.reporters = reporters

    def report(self, *args):
        for report in self.reporters:
            report(*args)

def generate_report(options):

    # crear el pool, con nº hilos core-1
    # generar coverfile
    # crear métricas
    # crear main reports
    # crear subreports
    # lanzar todas las métricas
    # esperar a que acabe todo
    # generar informe final

    #Bamfile object generation, if not sorted do it and(sequentally made, maybe no improvement due I/O limitations)
    #Sorted and .bai Checking
    bamlist = []
    for bamdir in options.bams.split(','):
        bam = bam_file(bamdir)
        if not bam.issorted():
            bamsorted = bam.sort_bam()
        else:
            bamsorted = bam
        bamlist.append(bamsorted)

    #Generating coverage_file objects from bamlist
    #TODO arreglar bedgraphs
    coveragefiles = []
    for bam in bamlist:
        cover,_ = bam.myCoverageBed(options.bed)
        coveragefiles.append(cover)

    # Creating Space of data to share between threads
    mgr = Manager()
    ns = mgr.Namespace()
    ns.coveragefiles = coveragefiles
    #Creating the pool and designating number of workers
    mainpool = Pool(processes= options.nthreads)


    #Main reporter generator
    mainReporter = HtmlReport(options.outdir, options)


    #Metric reporter generator, main reporter will be passed as an argument

    #Sensitivity: Depth_Threshold OK

    thresholdhtml = ThresholdHtml(mainReporter).report
    thresholdjson = ThresholdJson(options.outdir).report
    thresholdxls = ThresholdXls(options.outdir).report

    thresholdreporter = CompoundReporter([thresholdhtml, thresholdjson, thresholdxls])

    mainpool.apply_async(DepthThresProcessor().process, args=(ns.coveragefiles, thresholdreporter.report))

    # Sensitivity:(Optional) Saturation


    #Specificity: OnOffReport
    onoffhtml = OnOffHtml(mainReporter).report
    onoffjson = OnOffJson(options.outdir).report
    onoffxls = OnOffXls(options.outdir).report

    # onoffreporter = CompoundReporter([onoffhtml,onoffjson,onoffxls])
    # mainpool.apply_async(OnOffReadsProcessor().process, args=(bamlist, options.bed, onoffreporter.report)).get()

    onoffreporter = CompoundReporter([onoffhtml, onoffjson, onoffxls])
    onoffprocessor = OnOffReadsProcessor(bamlist)
    # # mainpool.apply_async(onoffprocessor.process, args=(options.bed, onoffreporter.report)).get()
    # mainpool.apply(onoffprocessor.process, args=(bamlist,options.bed, onoffreporter.report)).get()
    # # onoffprocessor.process(bamlist, options.bed, onoffreporter.report)

    mainpool.apply(onoffprocessor.process, args=(options.bed, onoffreporter.report))



    #Uniformity: Depth_distribution Coverage Distribution
    distributionhtml = DistributionHtml(mainReporter).report
    distributionjson = DistributionJson(options.outdir).report
    distributionxls = DistributionXls(options.outdir).report

    distributionreporter = CompoundReporter([distributionhtml, distributionjson, distributionxls])

    mainpool.apply_async(DepthDistrProcessor().process, args=(ns.coveragefiles, distributionreporter.report))

    #Uniformity: depth_perposition Coverage_per_position
    perpositionreporter = PerpositionHtml(mainReporter)
    mainpool.apply_async(DepthPerPositionProcessor().process, args=(ns.coveragefiles, perpositionreporter.report))

    #Uniformity: Standart deviation within regions
    stdevhtml = StdHtml(mainReporter).report
    stdevjson = StdJson(options.outdir).report
    stdevxls = StdXls(options.outdir).report

    stdevreporter = CompoundReporter([stdevhtml,stdevjson,stdevxls])
    mainpool.apply_async(StdevProcessor().process, args=(ns.coveragefiles, stdevreporter.report)).get()

    #Uniformity: Nocoverage txt
    zerocoveragereporter = ZeroCoverageTxt(options.outdir).report
    mainpool.apply_async(RegionsWithZeroesProcessor, args= (ns.coveragefiles,zerocoveragereporter))


    #Uniformity(optional) GC Bias reference required
    if options.reference is not None:
        gcbiasreporter = GcBiasHtml(mainReporter)
        mainpool.apply_async(GcBiasProcessor().process, args=(ns.coveragefiles, options.reference, gcbiasreporter.report))


    #Waits until all threads are finished
    mainpool.close()
    # Collects all the data generated
    mainpool.join()

    #Html generation
    #mainReporter.report()


#TODO gnerate json with config arguments.
def generate_json(outdir):

    config = dict([('warnbasescovered', 90),
                   ('warnsaturation',1e-5 ),
                   ('warnontarget', 80),
                   ('warnmeancoverage',40),
                   ('warncoverageregion',100),
                   ('warncoveragethreshold',6),
                   ('warncoveragecorrelation',0.9),
                   ('warnstd',0.3),
                   ('bins', 40)
                   ])
    with open(outdir + '/config.json', 'w') as config_json:
        json.dump(config, config_json)


class ConfigArgs:
    def __init__(self, configdir):
        with open(configdir) as json_file:
            self.config = json.load(json_file)

    def getconfig(self):
        return self.config

def main():

    options, args, parser = parse_arguments()
    if check_parameters(options, parser):
        generate_report(options)
    else:
        print('Error')

    ################################################

    #### Options and arguments #####################

    ################################################




if __name__ == '__main__':
    main()