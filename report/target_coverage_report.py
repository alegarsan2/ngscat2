
import plotly.graph_objs as go
import plotly
import xlwt

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
                        ), )


        data.append(trace)

    layout_comp = go.Layout(
        title='% on positions covered',
        hovermode='closest',
        barmode='group',
        xaxis=dict(showticklabels=True, showgrid=True, title='Coverage threshold'),
        yaxis=dict(title='% covered positions',range=[0, 100])
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

def target_distribution_plot(target_distribution_result):
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



        # for j, depth in enumerate(coveragefile['coveredposition']):
        #     ws.write(4, 1 + j, depth, header_style)
        #     ws.write(5, 1 + j, coveragefile['coveredposition'][depth])
        #     ws.write(6, 1 + j, coveragefile['perccoveredposition'][depth])

    wb.save(target_distribution_result['outdir'] + '/percentile.xls')