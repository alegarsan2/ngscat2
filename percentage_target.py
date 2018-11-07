import bam_file

TMP = '/tmp/'
import xlwt

import plotly.graph_objs as go
import plotly


def on_target_plot(read_on_results):
    """*****************************************************************************************************************
        Task:ALEGARSAN this method is dependant on "reads_on_target".
            Generates a graph interactive plot
        Inputs:
            read_on_results: Dict, with same data as json generated.
            Data used results_on_target['results'] = list of dictionaries with reads data of each

            bamlist: list of bam_file objects.
            nread: list of integers, each containing the number of on-target reads in the corresponding bam file.
            tread: list of integers, each containing the total number of reads in the corresponding bam file.
            onperchr: list of dictionaries. There must be one dictionary for each bam file in bamlist.
                Each dictionary contains data about the number of on-target
                reads per chromosome. Keys of the dictionary are contig ids (e.g. 'chr1').
                Elements are the number of on-target reads mapped in that contig.
            totalperchr: list of dictionaries. There must be one dictionary for each bam file in bamlist.
                Each dictionary contains data about the total number of
                reads mapped in each chromosome. Keys of the dictionary are contig ids (e.g. 'chr1').
                Elements are the total number of reads mapped in that contig.
            enrichment: list or multiprocessing.
                Array object containing the enrichment value for each bam (on-target reads per Kb)/(off target reads per Kb)
            percontarget: list or multiprocessing.
                Array object containing the percentaje of reads on target for each bam file
            outdir: string containing the full path to the directory were result files will be saved.
            legend: list of strings containing the label to tag each bam file in the xls and png files that will be created.
        Outputs: a new bar plot will be created named outdir/reads_on_target.png.
                This is a bar plot which indicates for each bam, the percentaje of on-target
                reads in each chromosome. In addition, a new outdir/reads_on_target.html

        ************************************************************************************************************"""
    colors = ['rgb(0,102,0)', 'rgb(255,0,0)', 'rgb(102,178,255)', 'rgb(178,102,255)']
    data = []
    for i, bam in enumerate(read_on_results['results']):
        trace = go.Bar(
            x=list(bam['perconperchr'].keys()),
            y=list(bam['perconperchr'].values()),
            hoverinfo='text',
            hoverlabel=dict(font=dict(color=['black'])),
            text=['Percentage:' + str(value) + '%' for value in list(bam['perconperchr'].values())],
            # mode='lines',
            name=bam['legend'],
            marker=dict(color=colors[i]))

        data.append(trace)

    layout_comp = go.Layout(
        title='Reads on target',
        hovermode='closest',
        barmode='group',
        xaxis=dict(showticklabels=True, showgrid=True, title='Chromosome/Contig'),
        yaxis=dict(title='% on-target reads'),
        shapes =[{
            'type': 'line',
            'x0': -0.5,
            'y0': bam['percontotal'],
            'x1': len(bam['perconperchr'])-0.5,
            'y1': bam['percontotal'],
            'line': {
                'color': 'rgb(50, 171, 96)',
                'width': 4,
                'dash': 'dot'}}])

    fig = go.Figure(data=data, layout=layout_comp)
    plotly.offline.plot(fig, filename=read_on_results['outdir'] + 'reads_on_target.html',
                        auto_open=True, config=dict(displaylogo=False, modeBarButtonsToRemove=['sendDataToCloud']))


def on_target_xls(read_on_results):
    # Initialize the workbook and sheet
    wb = xlwt.Workbook()

    # A sheet is created in the xls for each bam file
    for i, bam in enumerate(read_on_results['results']):
        ws = wb.add_sheet(bam['legend'])

        # Create header font
        header_style = xlwt.easyxf('font: bold on')

        ws.write(0, 0, 'Input bamfile: ', header_style)
        ws.write(0, 1, bam['bamfilename'])

        ws.write(1, 0, 'Input bedfile:', header_style)
        ws.write(1, 1, read_on_results['bedfile'])

        ws.write(2, 0, 'Enrichment:', header_style)
        ws.write(2, 1, bam['enrichment'])

        ws.write(4, 1, 'Reads on target', header_style)
        ws.write(4, 2, 'Reads off target', header_style)
        ws.write(4, 3, '% reads on target', header_style)
        ws.write(4, 4, '% reads off target', header_style)
        ws.write(5, 0, 'Total', header_style)
        ws.write(5, 1, bam['onread'])
        ws.write(5, 2, bam['totalread'])
        ws.write(5, 3, bam['percontotal'])
        ws.write(5, 4, 100.0 - bam['percontotal'])

        for j, chr in enumerate(list(bam['totalperchr'].keys())):
            ws.write(j + 6, 0, chr, header_style)
            ws.write(j + 6, 1, bam['onperchr'][chr])
            ws.write(j + 6, 2, bam['totalperchr'][chr] - bam['onperchr'][chr])
            ws.write(j + 6, 3, bam['perconperchr'][chr])
            ws.write(j + 6, 4, 100.0 - bam['perconperchr'][chr])

    wb.save(read_on_results['outdir'] + '/reads_on_target.xls')


def duplicates_plot(read_on_results):

    for i, bam in enumerate(read_on_results['results']):

        # Extract dictionary values
        x = []
        yon= []
        yoff= []

        onkeys= list(bam['perconduplicates'].keys())

        offkeys= list(bam['percoffduplicates'].keys())
        onvalues= list(bam['perconduplicates'].values())
        offvalues= list(bam['percoffduplicates'].values())
        maxduplicates = read_on_results['maxduplicates']

        #Select the maximum number of times duplicates to show
        if maxduplicates > len(onkeys) and len(offkeys):
            maxduplicates = max(len(onkeys),len(offkeys))


        x = onkeys[0:maxduplicates-1]
        x += ['>=' + str(len(x)+1) + 'X']
        yon = onvalues[0:maxduplicates-1]
        yon.append(sum(onvalues[maxduplicates-1:]))
        yoff = offvalues[0:maxduplicates - 1]
        yoff.append(sum(offvalues[maxduplicates-1:]))

        data = []
        traceon = go.Bar(
            x=x,
            y=yon,
            hoverinfo='text',
            hoverlabel=dict(font=dict(color=['black']*len(x))),
            text=['Percentage on: ' + str(value) + '%' for value in yon],
            # mode='lines',
            name='on target duplicates')
            #marker=dict(color=colors[i]))
        data.append(traceon)

        traceoff = go.Bar(
            x=x,
            y=yoff,
            hoverinfo='text',
            hoverlabel=dict(font=dict(color=['black'])),
            text=['Percentage off: ' + str(value) + '%' for value in yoff],
            name='off target duplicates')
            #marker=dict(color=colors[i]))
        data.append(traceoff)



    layout_comp = go.Layout(
        title='Duplicates',
        hovermode='closest',
        barmode='group',
        xaxis=dict(showticklabels=True, showgrid=True, title='# of duplicates'),
        yaxis=dict(title='% of reads'))

    fig = go.Figure(data=data, layout=layout_comp)
    plotly.offline.plot(fig, filename=read_on_results['outdir'] + 'duplicates_' + bam['legend'] + '.html',
                        auto_open=True, config=dict(displaylogo=False, modeBarButtonsToRemove=['sendDataToCloud']))




def duplicates_xls(read_on_results):
    # Initialize the workbook and sheet

    maxduplicates = read_on_results['maxduplicates']

    wb = xlwt.Workbook()

    # A sheet is created in the xls for each bam file
    for i, bam in enumerate(read_on_results['results']):
        ws = wb.add_sheet(bam['legend'])
        # Extract dictionary values
        keys = []
        yon = []
        yoff = []
        ypercon = []
        ypercoff = []

        onkeys = list(bam['perconduplicates'].keys())

        onvalues = list(bam['onduplicates'].values())
        offvalues = list(bam['offduplicates'].values())

        perconvalues = list(bam['perconduplicates'].values())
        percoffvalues = list(bam['percoffduplicates'].values())

        # Select the maximum number of times duplicates to show
        if maxduplicates > len(onkeys):
            maxduplicates = len(onkeys)

        keys = onkeys[0:maxduplicates-1]
        keys += ['>=' + str(len(keys)+1) + 'X']

        # Generation of duplicate reads data
        yon = onvalues[0:maxduplicates - 1]
        yon.append(sum(onvalues[maxduplicates - 1:]))
        yoff = offvalues[0:maxduplicates - 1]
        yoff.append(sum(offvalues[maxduplicates - 1:]))

        # Generation of percentage data
        ypercon = perconvalues[0:maxduplicates-1]
        ypercon.append(sum(perconvalues[maxduplicates-1:]))
        ypercoff = percoffvalues[0:maxduplicates - 1]
        ypercoff.append(sum(percoffvalues[maxduplicates-1:]))



        # Create header font
        header_style = xlwt.easyxf('font: bold on')

        ws.write(0, 0, 'Input bamfile: ', header_style)
        ws.write(0, 1, bam['bamfilename'])

        ws.write(1, 0, 'Input bedfile:', header_style)
        ws.write(1, 1, read_on_results['bedfile'])

        ws.write(4, 0, '# on', header_style)
        ws.write(5, 0, '# off', header_style)
        ws.write(6, 0, '# total', header_style)
        ws.write(7, 0, '% on', header_style)
        ws.write(8, 0, '% off', header_style)
        j= []
        for j in range(len(keys)):
            ws.write(3, j + 1, keys[j], header_style)
            ws.write(4, j + 1, yon[j])
            ws.write(5, j + 1, yoff[j])
            ws.write(6, j + 1, yon[j] + yoff[j])
            ws.write(7, j + 1, ypercon[j])
            ws.write(8, j + 1, ypercoff[j])

        ws.write(3, maxduplicates + 1, 'Total',header_style)
        ws.write(4, maxduplicates + 1, bam['onread'])
        ws.write(5, maxduplicates + 1, bam['totalread'] - bam['onread'])


    wb.save(read_on_results['outdir'] + 'duplicates.xls')