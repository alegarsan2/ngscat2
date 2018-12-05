import numpy as np
from plotly import graph_objs as go

import math
def getChromosomeNames(coverage):
    chromosomeNames = []
    for chromosome in coverage.chromosomes:
        chromosomeNames.append(chromosome.name)
    return chromosomeNames


def computeWindowSize(coverage,chromosomeName, npoints):
    ''' Compute window size
    Inputs: chromosome name , npoints: number of npoints of the representation '''
    regionlens = []
    medianlen = 0
    totallen = 0
    npointsratio = 0

    chrom = coverage.getChromosome(chromosomeName)
    for region in chrom.regions:
        regionlens.append(region.covEndIndex - region.covStartIndex)

    medianlen = np.median(regionlens)
    totallen = sum(regionlens)

    # If median is used and the number of points is greater than (npoints) look for multiple of median until
    # npoints is reached.
    npointsratio = int(totallen / (medianlen if medianlen > 0 else 1))

    if npointsratio > npoints:
        windowsize = int(medianlen * (npointsratio // npoints))
    else:
        windowsize = medianlen

    return windowsize


# def renderCoveragePerChr(coverage, chromosomeName, windowsize):
#     '''Calculate traces and values per chromosome'''
#
# # FIXME Encontrar manera de sumar std.
#     y = []
#     error = []
#     text = []
#     regmean = []
#     regstd = []
#     regindx = []
#
#     chromosome = coverage.getChromosome(chromosomeName)
#     chromlen = len(chromosome.regions)
#
#     i = chromosome.regions[0].covStartIndex + windowsize
#
#     for idx, region in enumerate(chromosome.regions):
#         if region.covEndIndex < i and idx != chromlen:
#             # Check whether the region is
#             regmean.append(region.mean)
#             regstd.append(region.std)
#             regindx.append(idx)
#
#         else:  # save data points
#             regmean.append(region.mean)
#             regstd.append(region.std)
#             regindx.append(idx)
#
#             i = region.covEndIndex + windowsize
#
#             y.append(np.mean(regmean))
#             #Taking into account the independence of regions the global std is equal to sqrt of sum of variances.
#             error.append(math.sqrt(sum([std**2 for std in regstd])))
#
#             text.append(str(round(y[-1],2)) + '\t' + str(round(error[-1], 2)) +
#                         '<br>'+ 'Start\t\t\t End\t\t\t Mean\t\t\t\t  Std <br>' +
#                         "".join([str(chromosome.regions[x].start) + "\t " + str(chromosome.regions[x].end) + "\t " +
#                                  str(round(chromosome.regions[x].mean, 2)) + "\t " + str(round(chromosome.regions[x].std, 2)) +
#                                  "<br>" for x in regindx]))
#
#             # Current index will be the index of the end of last region.
#             # Reboot
#             regmean = []
#             regstd = []
#             regindx = []
#     colors = ['rgb(0,102,0)', 'rgb(255,178,178)', 'rgb(102,178,255)', 'rgb(178,102,255)']
#     trace = go.Scatter(
#         x=list(range(0, len(y))),
#         y=y,
#         error_y=dict(
#             type='data',
#             array=error,
#             visible=True,
#             thickness=0.5,
#             width=0.3,
#             color='#c2d6d6'),
#         hoverinfo='text',
#         text=text,
#         mode='lines+markers',
#         name=str(coverage.name),
#         #line=dict(color=colors[0]),
#     )
#
#     return trace


def renderCoveragePerChr(coverage, chromosomeName, windowsize):
    '''Calculate traces and values per chromosome'''

# FIXME Encontrar manera de sumar std.
    y = []
    error = []
    text = []
    regmean = []
    regstd = []
    regindx = []

    chromosome = coverage.getChromosome(chromosomeName)
    chromlen = len(chromosome.regions)

    i = chromosome.regions[0].covStartIndex + windowsize

    for idx, region in enumerate(chromosome.regions):
        if region.covEndIndex < i and idx != chromlen:
            # Check whether the region is
            regmean.append(region.mean)
            regstd.append(region.std)
            regindx.append(idx)

        else:  # save data points
            regmean.append(region.mean)
            regstd.append(region.std)
            regindx.append(idx)

            i = region.covEndIndex + windowsize

            #Calculating graph points data
            y.append(np.mean(coverage.coverages[chromosome.regions[regindx[0]].covStartIndex:chromosome.regions[regindx[-1]].covEndIndex]))
            error.append(np.std(coverage.coverages[chromosome.regions[regindx[0]].covStartIndex:chromosome.regions[regindx[-1]].covEndIndex]))

            text.append(str(round(y[-1],2)) + '\t' + str(round(error[-1], 2)) +
                        '<br>'+ 'Start\t\t\t End\t\t\t Mean\t\t\t\t  Std <br>' +
                        "".join([str(chromosome.regions[x].start) + "\t " + str(chromosome.regions[x].end) + "\t " +
                                 str(round(chromosome.regions[x].mean, 2)) + "\t " + str(round(chromosome.regions[x].std, 2)) +
                                 "<br>" for x in regindx]))

            # Current index will be the index of the end of last region.
            # Reboot
            regmean = []
            regstd = []
            regindx = []
    colors = ['rgb(0,102,0)', 'rgb(255,178,178)', 'rgb(102,178,255)', 'rgb(178,102,255)']
    trace = go.Scatter(
        x=list(range(0, len(y))),
        y=y,
        error_y=dict(
            type='data',
            array=error,
            visible=True,
            thickness=0.5,
            width=0.3,
            color='#c2d6d6'),
        hoverinfo='text',
        text=text,
        mode='lines+markers',
        name=str(coverage.name),
        #line=dict(color=colors[0]),
    )

    return trace



class RegionsWithZeroesProcessor():

    def __init__(self,coveragefilecoverages ,path):
        self.file = open(path, "w")
        self.coverages = coveragefilecoverages

    def process(self, chromosome, region):
        zeroregion = []
        initzero = None
        rinit = region.start
        regionlen = region.end - region.start
        regionidx = region.covStartIndex
        regioni = self.coverages[region.covStartIndex:region.covEndIndex]
        for idx, basecoverage in enumerate(regioni):
            if idx != regionlen:
                if basecoverage == 0:
                    if initzero is None:
                        initzero = idx
                else:
                    if initzero is not None:
                        self.file.write("\t".join(map(str,[chromosome.name, initzero + region.start, idx + region.start])) + "\n")
                        initzero = None
            else:
                if basecoverage == 0 and initzero is not None:
                    self.file.write("\t".join(map(str,[chromosome.name, initzero + region.start, idx + region.start]))+ "\n")
                else:
                    self.file.write("\t".join(map(str,[chromosome.name, idx + region.start, idx + region.start]))+ "\n")
    def close(self):
        self.file.close()


def zeroCoverageRegions(coveragefile, outdir):
    '''Generate Nocoverage.txt of regions within targets that have Zero coverage.
    Input: Coveragefile object, outdir '''


    zerosProcessor = RegionsWithZeroesProcessor(coveragefile.coverages, outdir + "NoCoverage.txt")
    coveragefile.iterateOverRegions(zerosProcessor.process)
    zerosProcessor.close()


def processRegionStd(chromosome,region):
    if region.mean > 0:
        return region.std

# def region_std_distribution(coveragefiles, outdir, legend=None, bins='auto', warnthreshold = 40):
#     # Histo_CV sustitute.
#     region_stddistribution_result = {}
#     results = []
#     histlist = []
#     percentile = []
#     maximum = []
#     minimum = []
#     mean = []
#     median= []
#     zerocov = []
#
#
#
#     for i, coveragefile in enumerate(coveragefiles):
#         stdlist = []
#         #Non zero indexs
#         stdlist = coveragefile.iterateOverRegions(processRegionStd)
#         #FIXME no existe normalizacion de tamaño de bins, solo toma segun el tamaño de la muestra(como es sin 0) QUitar 0??
#         number_read, bin_edges = np.histogram(coveragefile.coverages[indnonzero], bins= bins)
#         bin_edges = bin_edges.tolist()
#
#         width = []
#         xaxis = []
#         for i in range(len(bin_edges)-1):
#             xaxis.append((bin_edges[i]+ bin_edges[i +1])/2)
#             width.append(bin_edges[1]- bin_edges[0])
#
#         histlist.append(dict([('numberread', number_read.tolist()),
#                               ('binedges', bin_edges),
#                               ('coveragepos', [round(x,2) for x in xaxis]),
#                               ('width', width)]))
#
#         percentile.append(dict(
#             [('Q1', np.percentile(coveragefile.coverages, 25)),
#             ('Q2', np.percentile(coveragefile.coverages, 50)),
#             ('Q3', np.percentile(coveragefile.coverages, 75))
#             ]))
#         median.append(np.median(coveragefile.coverages))
#         maximum.append(np.max(coveragefile.coverages))
#         minimum.append(np.min(coveragefile.coverages))
#         mean.append(coveragefile.coverages.mean())
#         zerocov.append(len(coveragefile.coverages) - len(indnonzero))
#         #TODO ver que parametros metemos en el JSON
#     for i in range(len(bamlist)):
#         results.append(dict(
#             [('bamfilename', bamlist[i].filename.decode('utf-8')),
#              ('legend', (legend[i] if legend is not None else bamlist[i].filename.decode('utf-8').split('/')[-1])),
#              ('histdata', histlist[i]),
#              ('percentile',percentile[i]),
#              ('max', float(maximum[i])),
#              ('min', float(minimum[i])),
#              ('mean', mean[i]),
#              ('median', median[i]),
#              ('zerocov', float(zerocov[i])),
#              ('status', 'OK' if mean[i] >= warnthreshold else 'Not OK')]
#         ))
#
#     region_stddistribution_result['results'] = results
#     region_stddistribution_result['bins'] = bins
#     region_stddistribution_result['warnthreshold'] = warnthreshold
#     region_stddistribution_result['outdir'] = outdir
#
#     with open(outdir + '/region_stddistribution_result.json', 'w') as outfile:
#         json.dump(region_stddistribution_result, outfile)
#
#     return region_stddistribution_result, coveragefiles


class StdReport():
    def __init__(self, coveragefile,warnthreshold):
        self.stdlist = coveragefile.iterateOverRegions(self.process)
        self.results = self.calculate_results(self,coveragefile,warnthreshold)
    def process(self,chromosome,region):
        if region.mean > 0:
            self.stdlist.append(region.std)
    def calculate_results(self,coveragefile,warnthreshold):
        #Plot here
    def std_distr_hist(self,outdir):
        #Plot His
    def std_distr_box(self,outdir):
        #Plot
    def std_distr_json(self,outdir):
        #Dump self.results
    def write(self,outdir):
        self.stdDistrBox(outdir)
        self.stdDistrBox(outdir)
        self.stdDistrJson(outdir)

