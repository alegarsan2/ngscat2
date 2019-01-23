 #!/usr/bin/env python3

import optparse
import sys
import os
import json
import shutil

from ngscat2.utils.getpaths import get_project_root
from multiprocessing import Manager
from multiprocessing import Pool
from multiprocessing import cpu_count

from ngscat2.utils.bam_file import bam_file
from ngscat2.utils import bed_file
from ngscat2.modules.report.html.HtmlReport import HtmlReport

from ngscat2.modules.metric.onoff_reads.metric import OnOffReadsProcessor
from ngscat2.modules.report.html.onoff_reads.report import Report as OnOffHtml
from ngscat2.modules.report.xls.onoff_reads.report import Report as OnOffXls
from ngscat2.modules.report.json.onoff_reads.report import Report as OnOffJson


from ngscat2.modules.metric.depth_threshold.metric import DepthThresProcessor
from ngscat2.modules.report.html.depth_threshold.report import Report as ThresholdHtml
from ngscat2.modules.report.json.depth_threshold.report import Report as ThresholdJson
from ngscat2.modules.report.xls.depth_threshold.report import Report as ThresholdXls


from ngscat2.modules.metric.depth_distribution.metric import DepthDistrProcessor
from ngscat2.modules.report.html.depth_distribution.report import Report as DistributionHtml
from ngscat2.modules.report.json.depth_distribution.report import Report as DistributionJson
from ngscat2.modules.report.xls.depth_distribution.report import Report as DistributionXls

from ngscat2.modules.metric.depth_perposition.metric import DepthPerPositionProcessor
from ngscat2.modules.report.html.depth_perposition.report import Report as PerpositionHtml

from ngscat2.modules.metric.depth_stdev.metric import StdevProcessor
from ngscat2.modules.report.html.depth_stdev.report import Report as StdHtml
from ngscat2.modules.report.xls.depth_stdev.report import Report as StdXls
from ngscat2.modules.report.json.depth_stdev.report import Report as StdJson

from ngscat2.modules.metric.zero_coverage.metric import RegionsWithZeroesProcessor
from ngscat2.modules.report.txt.zero_coverage.report import Report as ZeroCoverageTxt

from ngscat2.modules.metric.gc_bias.metric import GcBiasProcessor
from ngscat2.modules.report.html.gc_bias.report import Report as GcBiasHtml


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
    # parser.add_option("--extendtarget", dest="extend",
    #                   help="""Optional. Integer indicating the number of bases to extend each target region up and down-stream. Default=None.""",
    #                   default=None)
    parser.add_option("--reference", dest="reference",
                      help="""Optional. String indicating the path to a .fasta file containing the reference chromosomes. Default=None.""",
                      default=None)
    # parser.add_option("--saturation", dest="saturation",
    #                   help="""Optional. {y,n} to indicate whether saturation curve should be calculated. Default=n.""",
    #                   default='n')
    # parser.add_option("--depthlist", dest="depthlist",
    #                   help="""Optional. Will only be used in case --saturation is "y". Comma separated list of real numbers (do not leave spaces between) indicating the number of millions of reads to simulate for the saturation curve. E.g.: 1,5,10 would indicate 1*10^6, 5*10^6 and 10*10^6. Default=auto.""",
    #                   default='auto')
    parser.add_option("--coveragethrs", dest="coveragethresholds",
                      help="""Optional. Comma separated list of real numbers (do not leave spaces between) indicating coverage thresholds to be used when calculating percentages of covered bases (first graph in the report). Default=1,5,10,20,30.""",
                      default='1,5,10,20,30')

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
    if err is not '':
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


    #FIXME, lines down here commented for debugging
    # if ((os.path.isdir(options.outdir) or os.path.islink(options.outdir)) and (
    #         os.path.isdir(options.outdir + '/data') or os.path.islink(options.outdir + '/data')) and len(
    #         glob.glob(options.outdir + '/data/*_Ontarget_Coverage.html')) > 0):
    #     print('WARNING: ' + options.outdir + ' directory seems to contain previous NGScat results. Saving results of current execution in this directory may cause incorrect report generation.')
    #     print('Continue with current setting? (y/n)')
    #     proceed = input().lower()
    #     while proceed is not 'y' and proceed is not 'n':
    #         proceed = input().lower()
    #     if proceed is 'n':
    #         sys.exit(1)

    ## reference
    if (options.reference is not None and (not (os.path.isfile(options.reference) or os.path.islink(options.reference)))):
        print('ERROR: ' + options.reference + ' does not exist.')
        sys.exit(1)

    ## saturation
    # if (options.saturation is not 'y' and options.saturation is not 'n'):
    #     print('ERROR: incorrect value for --saturation parameters. Please indicate "y" or "n".')
    #     sys.exit(1)

    ## number of threads
    try:
        nthreads = int(options.nthreads)
    except ValueError:
        print('ERROR: invalid value for --nthreads option. Please, provide an integer value. Note that the application will launch as many processess as it needs between 1 and nthreads.')
        sys.exit(1)

    ## depth list
    # if (options.depthlist is not 'auto'):
    #     try:
    #         depthlist = map(float, options.depthlist.split(','))
    #     except ValueError:
    #         print('ERROR: invalid values for --depthlist option. Please, provide a comma separated list of values without leaving spaces, e.g.: 1,2,10,20')
    #         sys.exit(1)

    ## coveragethresholds
    try:
        coveragetrhesholds = map(float, options.coveragethresholds.split(','))
    except ValueError:
        print('ERROR: invalid values for --coveragethrs option. Please, provide a comma separated list of values without leaving spaces, e.g.: 1,2,10,20')
        sys.exit(1)

    ## FIXME problably it is not going to be used



    if (not (os.path.isdir(options.tmp) or os.path.islink(options.tmp))):
        print('ERROR: ' + options.tmp + ' does not exist.')
        sys.exit(1)

    options.coveragethresholds = tuple(int(x) for x in options.coveragethresholds.split(','))
    return True


#Class wrapper of different kinds of results for instance (onoffplots, onoffjson, onoffxls)
class CompoundReporter:
    def __init__(self, reporters):
        self.reporters = reporters

    def report(self, args):
        for report in self.reporters:
            report(*args)

def generate_report(options, config):
    #Directory generation, i
    if not os.path.exists(options.outdir):
        os.makedirs(options.outdir)
        os.mkdir(options.outdir + '/data')
        os.mkdir(options.outdir + '/img')
    else:
        shutil.rmtree(options.outdir)  # removes all the subdirectories!
        os.makedirs(options.outdir)
        os.mkdir(options.outdir + '/data')
        os.mkdir(options.outdir + '/img')

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
    bamlistdir = [bam.filename.decode('utf-8') for bam in bamlist]

    # Creating the pool and designating number of workers
    mainpool = Pool(processes=options.nthreads)
    # Main reporter generator
    mainReporter = HtmlReport(options.outdir, options)

    # #Specificity: OnOffReport

    # FIXME pasar las rutas de los bams ordenados a  esta funcion, e instanciar de nuevo clase bamfile
    onoffhtml = OnOffHtml(mainReporter).report
    onoffjson = OnOffJson(options.outdir).report
    onoffxls = OnOffXls(options.outdir).report

    # onoffreporter = CompoundReporter([onoffhtml,onoffjson,onoffxls])
    # mainpool.apply_async(OnOffReadsProcessor().process, args=(bamlist, options.bed, onoffreporter.report)).get()

    onoffreporter = CompoundReporter([onoffhtml, onoffjson, onoffxls])

    mainpool.apply_async(
        OnOffReadsProcessor(config.getconfig()['maxduplicates'], config.getconfig()['warnontarget']).process,
        args=(bamlistdir, options.bed), callback=onoffreporter.report)


    #Generating coverage_file objects from bamlist
    #TODO arreglar bedgraphs
    coveragefiles = []
    for bam in bamlist:
        cover,_ = bam.myCoverageBed(options.bed)
        coveragefiles.append(cover)

    del bamlist

    # Creating Space of data to share between threads
    mgr = Manager()
    ns = mgr.Namespace()
    ns.coveragefiles = coveragefiles

    # #Creating the pool and designating number of workers
    # mainpool = Pool(processes= options.nthreads)
    # #Main reporter generator
    # mainReporter = HtmlReport(options.outdir, options)

    #Metric reporter generator, main reporter will be passed as an argument

    #Sensitivity: Depth_Threshold OK

    thresholdhtml = ThresholdHtml(mainReporter).report
    thresholdjson = ThresholdJson(options.outdir).report
    thresholdxls = ThresholdXls(options.outdir).report

    thresholdreporter = CompoundReporter([thresholdhtml, thresholdjson, thresholdxls])

    mainpool.apply_async(DepthThresProcessor(options.coveragethresholds, config.getconfig()['warncoveragethreshold']).process,
                         args=(ns.coveragefiles,), callback=thresholdreporter.report).get()

    # Sensitivity:(Optional) Saturation


    # #Specificity: OnOffReport

    # #FIXME pasar las rutas de los bams ordenados a  esta funcion, e instanciar de nuevo clase bamfile
    # onoffhtml = OnOffHtml(mainReporter).report
    # onoffjson = OnOffJson(options.outdir).report
    # onoffxls = OnOffXls(options.outdir).report
    #
    # # onoffreporter = CompoundReporter([onoffhtml,onoffjson,onoffxls])
    # # mainpool.apply_async(OnOffReadsProcessor().process, args=(bamlist, options.bed, onoffreporter.report)).get()
    #
    # onoffreporter = CompoundReporter([onoffhtml, onoffjson, onoffxls])
    #
    # mainpool.apply_async(OnOffReadsProcessor(config.getconfig()['maxduplicates'], config.getconfig()['warnontarget']).process,
    #                      args=(bamlistdir, options.bed), callback=onoffreporter.report)



    #Uniformity: Depth_distribution Coverage Distribution
    distributionhtml = DistributionHtml(mainReporter).report
    distributionjson = DistributionJson(options.outdir).report
    distributionxls = DistributionXls(options.outdir).report

    distributionreporter = CompoundReporter([distributionhtml, distributionjson, distributionxls])

    mainpool.apply_async(DepthDistrProcessor(config.getconfig()['distributionbins'], config.getconfig()['warnmeancoverage']).process,
                         args=(ns.coveragefiles,), callback=distributionreporter.report)

    #Uniformity: depth_perposition Coverage_per_position
    perpositionHtml = PerpositionHtml(mainReporter).report
    perpositionreporter = CompoundReporter([perpositionHtml])

    mainpool.apply_async(DepthPerPositionProcessor(config.getconfig()['npointsperchrom']).process,
                         args=(ns.coveragefiles,), callback= perpositionreporter.report)

    # #Uniformity: Standart deviation within regions
    stdevhtml = StdHtml(mainReporter).report
    stdevjson = StdJson(options.outdir).report
    stdevxls = StdXls(options.outdir).report

    stdevreporter = CompoundReporter([stdevhtml,stdevjson,stdevxls])
    mainpool.apply_async(StdevProcessor(config.getconfig()['warnstd']).process,
                         args=(ns.coveragefiles,), callback= stdevreporter.report)

    #Uniformity: Nocoverage txt
    zerocoveragereporter = ZeroCoverageTxt(options.outdir).report
    mainpool.apply_async(RegionsWithZeroesProcessor, args= (ns.coveragefiles, zerocoveragereporter))



    # #Uniformity(optional) GC Bias reference required
    if options.reference is not None:
        gcbiashtml = GcBiasHtml(mainReporter).report
        gcbiasreporter = CompoundReporter([gcbiashtml])

        mainpool.apply_async(GcBiasProcessor().process, args=(ns.coveragefiles, options.reference), callback= gcbiasreporter.report)


    #Waits until all threads are finished
    mainpool.close()
    # Collects all the data generated
    mainpool.join()

    #Html generation
    mainReporter.report(options)


# #TODO gnerate json with config arguments.
# def generate_json(outdir):
#
#     config = dict([('warnbasescovered', 90),
#                    ('warnsaturation',1e-5 ),
#                    ('warnontarget', 80),
#                    ('maxduplicates', 5),
#                    ('warnmeancoverage',40),
#                    ('warncoveragethreshold',6),
#                    ('warncoveragecorrelation',0.9),
#                    ('warnstd',0.3),
#                    ('npointsperchrom', 50),
#                    ('distributionbins', 40)
#                    ])
#     with open(outdir + 'config.json', 'w') as config_json:
#         json.dump(config, config_json)

class ConfigArgs:
    def __init__(self, configdir):
        with open(configdir , 'r') as json_file:
            self.config = json.load(json_file)

    def getconfig(self):
        return self.config

def main():

    options, args, parser = parse_arguments()
    if check_parameters(options, parser):

        #TODO load the config json
        #generate_json(os.path.dirname(sys.argv[0]))
        config_folder = get_project_root()
        config_filepath = os.path.join(config_folder, "config.json")
        config =  ConfigArgs(config_filepath)
        #configargs = ConfigArgs(os.path.dirname(sys.argv[0]) + 'config.json')
        generate_report(options, config)

        print('a')
    else:
        print('Error')

    ################################################

    #### Options and arguments #####################

    ################################################




if __name__ == '__main__':
    main()
