import plotly.graph_objs as go
import plotly
import numpy as np

class Report():
    def __init__(self, mainreporter):
        #En este caso no hace falta nada, unicamente el outdir por lo que hay que
        #pasarle el mainreporter al report.

        self.mainreporter = mainreporter
        self.plot_dir = []


    def report(self, coveragefiles, npoints):

    # npoint,  warnregionsize = 100, warnthreshold=6):

        chromosomenames = coveragefiles[0].getChromosomeNames()
        for chromosomeName in chromosomenames:
            traces = []

            windowsize = self.computeWindowSize(coveragefiles[0], chromosomeName, npoints)
            traces = []
            for coverage in coveragefiles:
                trace = self.renderCoveragePerChr(coverage, chromosomeName, windowsize)
                traces.append(trace)
            maxrange = max([max(trace.y) for trace in traces])
            layout_comp = go.Layout(
                title=chromosomeName,
                hovermode='closest',
                xaxis=dict(showticklabels=False, showgrid=False),
                yaxis=dict(title="Coverage", range=[0, maxrange + maxrange / 30]),
                margin=go.layout.Margin(
                    l=50,
                    r=10,
                    b=10,
                    t=50,
                    pad=4
                ),
            )
            #TODO como obtengo el output dir desde el mainreporter
            fig = go.Figure(data=traces, layout=layout_comp)
            plotly.offline.plot(fig, filename= self.mainreporter.outdir + '/data/' + chromosomeName + '_Ontarget_Coverage.html',
                                auto_open=False,
                                config=dict(displaylogo=False, modeBarButtonsToRemove=['sendDataToCloud'],
                                            showLink=False))
            self.plot_dir.append(self.mainreporter.outdir + '/data/' + chromosomeName + '_Ontarget_Coverage.html')

            self.mainreporter.addsection('covperposition', self)

    def computeWindowSize(self, coverage, chromosomename, npoints):
        ''' Compute window size
        Inputs: chromosome name , npoints: number of npoints of the representation '''
        regionlens = []
        medianlen = 0
        totallen = 0
        npointsratio = 0

        chrom = coverage.getChromosome(chromosomename)
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

    def renderCoveragePerChr(self, coverage, chromosomename, windowsize):
        '''Calculate traces and values per chromosome'''


        y = []
        error = []
        text = []
        regmean = []
        regstd = []
        regindx = []

        chromosome = coverage.getChromosome(chromosomename)
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
            name=str(coverage.name.decode('utf-8').split('/')[-1]),
            #line=dict(color=colors[0]),
        )

        return trace
