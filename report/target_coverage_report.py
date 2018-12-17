import numpy as np
import plotly.graph_objs as go
import plotly
import xlwt

from metric.region_coverage import getChromosomeNames, computeWindowSize, renderCoveragePerChr


def target_coverage_plot(target_coverage_result):
    colors = ['rgb(0,102,0)', 'rgb(255,0,0)', 'rgb(102,178,255)', 'rgb(178,102,255)']
    data = []
    for i, coveragefile in enumerate(target_coverage_result['results']):
        trace = go.Bar(
            x=list(coveragefile['perccoveredposition'].keys()),
            y=list(coveragefile['perccoveredposition'].values()),
            hoverinfo='text',
            hoverlabel=dict(font=dict(color=['black'])),
            text=['Percentage:' + str(value) + '%' for value in list(coveragefile['perccoveredposition'].values())],
            # mode='lines',
            name=coveragefile['legend'],
            marker=dict(color=colors[i],
                        line=dict(
                            color='rgb(0,0,0)',
                            width=.6)
                        ),)

        data.append(trace)

    layout_comp = go.Layout(
        title='% on positions covered',
        hovermode='closest',
        barmode='group',
        xaxis=dict(showticklabels=True, showgrid=True, title='Coverage threshold'),
        yaxis=dict(title='% covered positions',range=[0, 100]),
        margin=go.layout.Margin(
            l=50,
            r=10,
            b=10,
            t=50,
            pad=4
        ),
    )


    fig = go.Figure(data=data, layout=layout_comp)
    plotly.offline.plot(fig, filename=target_coverage_result['outdir'] + 'covered_positions.html',
                        auto_open=True, config=dict(displaylogo=False, modeBarButtonsToRemove=['sendDataToCloud'])
                        )


def target_coverage_xls(target_coverage_result):
    # Initialize the workbook and sheet
    wb = xlwt.Workbook()

    # A sheet is created in the xls for each coveragefile file
    for i, coveragefile in enumerate(target_coverage_result['results']):
        ws = wb.add_sheet(coveragefile['legend'])

        # Create header font
        header_style = xlwt.easyxf('font: bold on')

        ws.write(0, 0, 'Input bamfile: ', header_style)
        ws.write(0, 1,coveragefile['bamfilename'])


        # TODO Añadir a target_coverage argumento bed para obtener bedfilename
        # ws.write(1, 0, 'Input bedfile:', header_style)
        #ws.write(1, 1, read_on_results['bedfile'])

        ws.write(2, 0, 'Outdir:', header_style)
        ws.write(2, 1, target_coverage_result['outdir'])

        ws.write(3, 0, 'Table', header_style)

        ws.write(5, 0, 'Number of on positions covered', header_style)
        ws.write(6, 0, '% of on positions covered', header_style)
        ws.write(4, 1 + len(coveragefile['coveredposition']), 'Total',header_style)
        ws.write(5, 1 + len(coveragefile['coveredposition']), coveragefile['ntotalposition'])

        for j, depth in enumerate(coveragefile['coveredposition']):
            ws.write(4, 1+j, depth, header_style)
            ws.write(5, 1+j, coveragefile['coveredposition'][depth])
            ws.write(6, 1+j, coveragefile['perccoveredposition'][depth])

    wb.save(target_coverage_result['outdir'] + '/coverage_summary.xls')

def target_distribution_histplot(target_distribution_result):
    colors = ['rgb(0,102,0)', 'rgb(255,0,0)', 'rgb(102,178,255)', 'rgb(178,102,255)']
    data = []
    for i, coveragefile in enumerate(target_distribution_result['results']):
        trace = go.Bar(
            y=coveragefile['histdata']['numberread'],
            x= coveragefile['histdata']['coveragepos'],
            width = coveragefile['histdata']['width'],
            hoverinfo='text',
            hoverlabel=dict(font=dict(color=['black'])),
            text=['Count: ' + str(value) + ('<br>') + 'Coverage: ' + str(coveragefile['histdata']['coveragepos'][i])
                  for i, value in enumerate(coveragefile['histdata']['numberread'])],
            # mode='lines',
            opacity= 0.7,
            name=coveragefile['legend'],
            marker=dict(color=colors[i],
                        line= dict(
                            color = 'rgb(0,0,0)',
                            width = .6)
                        ),)

        data.append(trace)

    layout_comp = go.Layout(
        title='Histogram',
        hovermode='closest',
        barmode='group',
        xaxis=dict(showticklabels=True, showgrid=True, title='Coverage'),
        yaxis=dict(title='Count')
    )

    fig = go.Figure(data=data, layout=layout_comp)
    plotly.offline.plot(fig, filename=target_distribution_result['outdir'] + 'target_hist.html',
                        auto_open=True, config=dict(displaylogo=False, modeBarButtonsToRemove=['sendDataToCloud']))


def target_distribution_boxplot(target_distribution_result, coveragelist):
    colors = ['rgb(0,102,0)', 'rgb(255,0,0)', 'rgb(102,178,255)', 'rgb(178,102,255)']
    data = []
    for i, coveragefile in enumerate(coveragelist):
        trace = go.Box(
            #Random subsampling in order to represent fasther the data. Final size 100000 points
            y=np.random.choice(coveragefile.coverages,
                                 size=int(len(coveragefile.coverages)/((len(coveragefile.coverages) // 100000)
                                                                        if len(coveragefile.coverages) > 100000 else 1))),
            name=target_distribution_result['results'][i]['legend'],
            marker=dict(
                color=colors[i],
            ),
            boxpoints='suspectedoutliers',
            jitter=0.01)
        data.append(trace)

    layout_comp = go.Layout(
            title='',
            hovermode='closest',
            #barmode='group',
            xaxis=dict(showticklabels=True, showgrid=True, title=''),
            yaxis=dict(title='Depth',
                       autorange = True),
            margin=go.layout.Margin(
                l=50,
                r=10,
                b=10,
                t=50,
                pad=4
            ),
        )


    fig = go.Figure(data=data, layout=layout_comp)
    plotly.offline.plot(fig, filename=target_distribution_result['outdir'] + 'target_boxplot.html',
                        auto_open=True, config= dict(displaylogo=False, modeBarButtonsToRemove=['sendDataToCloud']))
def target_distribution_xls(target_distribution_result):
    wb = xlwt.Workbook()

    # A sheet is created in the xls for each coveragefile file
    for i, coveragefile in enumerate(target_distribution_result['results']):
        ws = wb.add_sheet(coveragefile['legend'])

        # Create header font
        header_style = xlwt.easyxf('font: bold on')

        ws.write(0, 0, 'Input bamfile: ', header_style)
        ws.write(0, 1, coveragefile['bamfilename'])

        # TODO Añadir a target_coverage argumento bed para obtener bedfilename
        # ws.write(1, 0, 'Input bedfile:', header_style)
        # ws.write(1, 1, read_on_results['bedfile'])

        ws.write(2, 0, 'Outdir:', header_style)
        ws.write(2, 1, target_distribution_result['outdir'])

        ws.write(4, 0, 'Table', header_style)

        ws.write(5, 0, 'Number bases coverage 0', header_style)
        ws.write(5, 1, 'Q1', header_style)
        ws.write(5, 2, 'Q2', header_style)
        ws.write(5, 3, 'Q3', header_style)
        ws.write(5, 4, 'Maximum', header_style)
        ws.write(5, 5, 'Minimum', header_style)
        ws.write(5, 6, 'Median', header_style)
        ws.write(5, 7, 'Mean', header_style)

        ws.write(6, 0, coveragefile['zerocov'])
        ws.write(6, 1, coveragefile['percentile']['Q1'])
        ws.write(6, 2, coveragefile['percentile']['Q2'])
        ws.write(6, 3, coveragefile['percentile']['Q3'])
        ws.write(6, 4, coveragefile['max'])
        ws.write(6, 5, coveragefile['min'])
        ws.write(6, 6, coveragefile['median'])
        ws.write(6, 7, coveragefile['mean'])

    wb.save(target_distribution_result['outdir'] + '/percentile.xls')


def coverage_per_chr(coveragelist, npoints, outdir):
                     # npoint,  warnregionsize = 100, warnthreshold=6):

        chromosomeNames = getChromosomeNames(coveragelist[0])
        for chromosomeName in chromosomeNames:
            traces = []

            windowsize = computeWindowSize(coveragelist[0], chromosomeName, npoints)
            traces = []
            for coverage in coveragelist:
                trace = renderCoveragePerChr(coverage, chromosomeName, windowsize)
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

            fig = go.Figure(data=traces, layout=layout_comp)
            plotly.offline.plot(fig, filename=outdir + chromosomeName + '_Ontarget_Coverage.html',
                                auto_open=True,
                                config=dict(displaylogo=False, modeBarButtonsToRemove=['sendDataToCloud'],
                                            showLink=False))



#def region_std_distribution_histplot(coveragelist, ):



# for coverage in coveragelist:
#     c = coverage
#     i = 0
#
#     for chr in coverage.chromosomes:
#         medianlen = []
#         regionlens = []
#         medianlen = 0
#         y = []
#         regmean = []
#         regstd = []
#         regindx = []
#         error = []
#         text = []
#         traces = []
#         numofregion = len(chr.regions)
#         totallen = 0
#         npointsratio = 0  # Number of points using median as window size
#
#         # Take median lenght in order to establish the window size.
#         if numofregion > 10:
#             print(len(chr.regions))
#             for region in chr.regions:
#                 regionlens.append(region.covEndIndex - region.covStartIndex)
#             medianlen = np.median(regionlens)
#             totallen = sum(regionlens)
#
#         # If median is used and the number of points is greater than (npoints) look for multiple of median until
#         # npoints is reached.
#         npointsratio = int(totallen / (medianlen if medianlen > 0 else 1))
#         if npointsratio >= 10:
#             windowsize = int(medianlen * (npointsratio // 10))
#             i = windowsize
#             for idx, region in enumerate(chr.regions):
#
#                 if region.covEndIndex < i and idx != numofregion:
#                     # Check wether the region is
#                     regmean.append(region.mean)
#                     regstd.append(region.std)
#                     regindx.append(idx)
#
#                 else:  # save data points
#                     regmean.append(region.mean)
#                     regstd.append(region.std)
#                     regindx.append(idx)
#
#                     i = region.covEndIndex + windowsize
#
#                     y.append(np.mean(regmean))
#                     error.append(np.std(regmean))
#                     text.append('Start\t\t\t End\t\t\t Mean\t\t\t\t   Std <br>' +
#                                 "".join([str(chr.regions[x].start) + "\t " + str(chr.regions[x].end) + "\t " +
#                                          str(round(chr.regions[x].mean, 2)) + "\t " + str(round(chr.regions[x].std, 2)) +
#                                          "<br>" for x in regindx]))
#                     regmean = []
#                     regstd = []
#                     regindx = []
#
#             colors = ['rgb(0,102,0)', 'rgb(255,178,178)', 'rgb(102,178,255)', 'rgb(178,102,255)']
#             trace = go.Scatter(
#                 x=list(range(0, len(y))),
#                 y=y,
#                 error_y=dict(
#                     type='data',
#                     array=error,
#                     visible=True,
#                     thickness=1.5,
#                     width=1,
#                     color='#c2d6d6'),
#                 hoverinfo='text',
#                 text=text,
#                 mode='lines+markers',
#                 name=str(coverage.name),
#                 line=dict(color=colors[0]),
#             )
#             traces.append(trace)
#
#             layout_comp = go.Layout(
#                 title=str(chr.name),
#                 hovermode='closest',
#                 xaxis=dict(showticklabels=False, showgrid=False),
#                 yaxis=dict(title="Coverage", range=[0, max(y) + max(y) / 30]),
#                 margin=go.layout.Margin(
#                     l=50,
#                     r=10,
#                     b=10,
#                     t=50,
#                     pad=4
#                 ),
#             )
#
#             fig = go.Figure(data=traces, layout=layout_comp)
#             plotly.offline.plot(fig, filename=outdir + chr.name + '_Ontarget_Coverage.html',
#                                 auto_open=True,
#                                 config=dict(displaylogo=False, modeBarButtonsToRemove=['sendDataToCloud'], showLink=False))
#
#             print("A")
#
#                 # Current index will be the index of the end of last region.
#                 # Reboot
