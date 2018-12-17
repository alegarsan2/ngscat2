import numpy as np
from metric import region_coverage
from plotly import graph_objs as go
import plotly
import json
import xlwt


class StdReport():
    #Std exon report object, inputs: coveragefiles list of coverages objects, and warnthreshold(if mean of normalized
    #distribution is less than x)
    def __init__(self, coveragefiles, warnthreshold):
        self.coverages = coveragefiles
        self.stdlists = []
        for coverage in coveragefiles:
            self.stdlists.append([])
            coverage.iterateOverRegions(self.process)
        self.calculate_results(coveragefiles, warnthreshold)

    def process(self, chromosome, region):
        #The normalization is std of a region divided with the mean of the region.
        if region.mean > 0:
            self.stdlists[-1].append(region.std/region.mean)

    def calculate_results(self, coveragefiles, warnthreshold = 0.3):

        histlist = []
        percentile = []
        region_stddistribution_result = {}
        median = None
        maximum = None
        minimum = None
        mean = None
        results = []


        for indx, coveragefile in enumerate(self.coverages):
            number_std, bin_edges = np.histogram(self.stdlists[indx], bins= np.arange(0, 1, 0.007))
            bin_edges = bin_edges.tolist()
            width = []
            xaxis = []
            for i in range(len(bin_edges)-1):
                xaxis.append((bin_edges[i]+ bin_edges[i +1])/2)
                width.append(bin_edges[1]- bin_edges[0])

            histlist = dict([('numberstd', number_std.tolist()),
                                  ('binedges', bin_edges),
                                  ('stdpos', [round(x,2) for x in xaxis]),
                                  ('width', width)])

            percentile = (dict(
                [('Q1', np.percentile(self.stdlists[indx], 25)),
                ('Q2', np.percentile(self.stdlists[indx], 50)),
                ('Q3', np.percentile(self.stdlists[indx], 75))
                ]))
            median = np.median(self.stdlists[indx])
            maximum = np.max(self.stdlists[indx])
            minimum = np.min(self.stdlists[indx])
            mean = np.mean(self.stdlists[indx])

            results.append(dict(
                        [('bamfilename', coveragefile.name.decode('utf-8')),
                         #('legend', (legend[i] if legend is not None else bamlist[i].filename.decode('utf-8').split('/')[-1])),
                         ('histdata', histlist),
                         ('percentile',percentile),
                         ('max', float(maximum)),
                         ('min', float(minimum)),
                         ('mean', mean),
                         ('median', median),
                         ('status', 'OK' if mean >= warnthreshold else 'Not OK')]
                    ))


        region_stddistribution_result['results'] = results
        region_stddistribution_result['warnthreshold'] = warnthreshold

        self.results = region_stddistribution_result

    def std_distr_hist(self,outdir):
        # Plot His
        traces = []
        print(len(self.results['results']))
        for indx, stdlist in enumerate(self.stdlists):

            trace = go.Histogram(
                x= stdlist,
                xbins= dict(end=1, size= 0.007, start= 0),
                #hoverinfo='text',
                #hoverlabel=dict(font=dict(color=['black'])),
                #text=['Count: ' + str(x[i]) + '<br>' + 'Coverage: ' + str(results['histdata']['stdpos'][i])
                      #for i, value in enumerate(results['histdata']['numberstd'])],
                # mode='lines',
                opacity=0.7,
                name=self.results['results'][indx]['bamfilename'],
                marker=dict(
                        line=dict(
                        color='rgb(0,0,0)',
                        width=.6)))
            traces.append(trace)
        #Old takes the hist already calculated
        # def std_distr_hist(self,outdir):
        # # Plot His
        # traces = []
        # print(len(self.results['results']))
        # for results in self.results['results']:
        #     a = results
        #     b = results['histdata']['numberstd']
        #     c= results['histdata']['numberstd']
        #     trace = go.Bar(
        #         y=results['histdata']['numberstd'],
        #         x=results['histdata']['stdpos'],
        #         width=results['histdata']['width'],
        #         hoverinfo='text',
        #         hoverlabel=dict(font=dict(color=['black'])),
        #         text=['Count: ' + str(value) + '<br>' + 'Coverage: ' + str(results['histdata']['stdpos'][i])
        #               for i, value in enumerate(results['histdata']['numberstd'])],
        #         # mode='lines',
        #         opacity=0.7,
        #         name=results['bamfilename'],
        #         marker=dict(
        #                 line=dict(
        #                 color='rgb(0,0,0)',
        #                 width=.6)))
        #     traces.append(trace)

        layout_comp = go.Layout(
            title='Histogram',
            #hovermode='closest',
            barmode= 'overlay',
            xaxis=dict(showticklabels=True, showgrid=True, title='Normalized standart deviation',zeroline = True),
            yaxis=dict(showticklabels = True, title='Frequency', zeroline = True),
            showlegend=  True,
            legend=dict(x=0.75,
                        y=1),
            margin=go.layout.Margin(
                l=50,
                r=10,
                b=30,
                t=50,
                pad=4
            ),
        )

        fig = go.Figure(data=traces, layout=layout_comp)
        plotly.offline.plot(fig, filename=outdir + 'std_histogram.html',
                            auto_open=True, config=dict(displaylogo=False, modeBarButtonsToRemove=['sendDataToCloud']))


    def std_distr_box(self,outdir):
        # Plot
        traces =[]
        for indx, stdlist in enumerate(self.stdlists):
            trace = go.Box(
                # Random subsampling in order to represent fasther the data. Final size 100000 points
                y=np.random.choice(stdlist, size=int(len(stdlist) / ((len(stdlist) // 100000) if
                                                                              len(stdlist) > 100000 else 1))),
                name=self.results['results'][indx]['bamfilename'],
                boxpoints='suspectedoutliers',
                jitter=0.01)
            traces.append(trace)

        layout_comp = go.Layout(
            title='',
            hovermode='closest',
            # barmode='group',
            xaxis=dict(showticklabels=True, showgrid=True, title=''),
            yaxis=dict(title='Depth',
                       autorange=True),
            margin=go.layout.Margin(
                l=50,
                r=10,
                b=10,
                t=50,
                pad=4
            ),
        )

        fig = go.Figure(data=traces, layout=layout_comp)
        plotly.offline.plot(fig, filename=outdir + 'std_boxplot.html',
                            auto_open=True,
                            config=dict(displaylogo=False, modeBarButtonsToRemove=['sendDataToCloud']))

    def region_coverage_json(self,outdir):
        '''Inputs: List of stdreport objects
            Output: std_wexons.json'''

        with open(outdir + '/std_wexons.json', 'w') as outfile:
            json.dump(self.results, outfile)

    def region_coverage_xls(self, outdir):
        '''Inputs: List of stdreport objects
            Output: std_wexons.xls'''
        wb = xlwt.Workbook()
        ws = wb.add_sheet('std_exon_results')
        header_style = xlwt.easyxf('font: bold on')

        ws.write(0, 0, 'Sample', header_style)
        ws.write(0, 1, 'Q1', header_style)
        ws.write(0, 2, 'Q2', header_style)
        ws.write(0, 3, 'Q3', header_style)
        ws.write(0, 4, 'Maximum', header_style)
        ws.write(0, 5, 'Minumum', header_style)
        ws.write(0, 6, 'Mean', header_style)
        ws.write(0, 7, 'Median', header_style)


        # TODO Añadir a target_coverage argumento bed para obtener bedfilename
        # ws.write(1, 0, 'Input bedfile:', header_style)
        # ws.write(1, 1, read_on_results['bedfile'])
        for idx, results in enumerate(self.results['results']):
            ws.write(1 + idx, 0, results['bamfilename'], header_style)
            ws.write(1 + idx, 1, results['percentile']['Q1'])
            ws.write(1 + idx, 2, results['percentile']['Q2'])
            ws.write(1 + idx, 3, results['percentile']['Q3'])
            ws.write(1 + idx, 4, results['max'])
            ws.write(1 + idx, 5, results['min'])
            ws.write(1 + idx, 6, results['mean'])
            ws.write(1 + idx, 7, results['median'])

        wb.save(outdir + '/std_wexons.xls')

    def output(self, outdir):

        self.std_distr_hist(outdir)
        self.std_distr_box(outdir)
        self.region_coverage_xls(outdir)
        self.region_coverage_json(outdir)

