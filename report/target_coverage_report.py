
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
            marker=dict(color=colors[i]))

        data.append(trace)

    layout_comp = go.Layout(
        title='Reads on target',
        hovermode='closest',
        barmode='group',
        xaxis=dict(showticklabels=True, showgrid=True, title='Coverage threshold'),
        yaxis=dict(title='% covered positions'))

    fig = go.Figure(data=data, layout=layout_comp)
    plotly.offline.plot(fig, filename=target_coverage_result['outdir'] + 'covered_positions.html',
                        auto_open=True, config=dict(displaylogo=False, modeBarButtonsToRemove=['sendDataToCloud']))


def target_coverage_xls(target_coverage_result):
    # Initialize the workbook and sheet
    wb = xlwt.Workbook()

    # A sheet is created in the xls for each coveragefile file
    for i, coveragefile in enumerate(target_coverage_result['results']):

        # Create header font
        header_style = xlwt.easyxf('font: bold on')

        ws.write(0, 0, 'Input bamfile: ', header_style)
        ws.write(0, 1,coveragefile['bamfilename'])

        ws.write(1, 0, 'Input bedfile:', header_style)
        ws.write(1, 1, read_on_results['bedfile'])

        ws.write(2, 0, 'Enrichment:', header_style)
        ws.write(2, 1,coveragefile['enrichment'])

        ws.write(4, 1, 'Reads on target', header_style)
        ws.write(4, 2, 'Reads off target', header_style)
        ws.write(4, 3, '% reads on target', header_style)
        ws.write(4, 4, '% reads off target', header_style)
        ws.write(5, 0, 'Total', header_style)
        ws.write(5, 1,coveragefile['onread'])
        ws.write(5, 2,coveragefile['totalread'])
        ws.write(5, 3,coveragefile['percontotal'])
        ws.write(5, 4, 100.0 -coveragefile['percontotal'])

        for j, chr in enumerate(listcoveragefile['totalperchr'].keys())):
            ws.write(j + 6, 0, chr, header_style)
            ws.write(j + 6, 1,coveragefile['onperchr'][chr])
            ws.write(j + 6, 2,coveragefile['totalperchr'][chr] -coveragefile['onperchr'][chr])
            ws.write(j + 6, 3,coveragefile['perconperchr'][chr])
            ws.write(j + 6, 4, 100.0 -coveragefile['perconperchr'][chr])

    wb.save(read_on_results['outdir'] + '/reads_on_target.xls')