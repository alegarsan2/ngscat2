import xlwt

class Report():
    def __init__(self, outdir):
        self.outdir = outdir

    def report(self, read_on_resutls):
        self.on_target_xls(read_on_resutls)
        self.duplicates_xls(read_on_resutls)


    def on_target_xls(self,read_on_results):
        # Initialize the workbook and sheet
        wb = xlwt.Workbook()

        # A sheet is created in the xls for each bam file
        for i, bam in enumerate(read_on_results['results']):
            ws = wb.add_sheet(bam['legend'] if len(bam['legend']) < 31 else bam['legend'][-31:])

            # Create header font
            header_style = xlwt.easyxf('font: bold on')

            ws.write(0, 0, 'Input bamfile: ', header_style)
            ws.write(0, 1, bam['bamfilename'])

            ws.write(1, 0, 'Input bedfile:', header_style)
            ws.write(1, 1, read_on_results['bedfile'])

            ws.write(2, 0, 'Enrichment:', header_style)
            ws.write(2, 1, bam['enrichment'])
            ws.write(4, 0, 'Table', header_style)

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

        wb.save(self.outdir + '/reads_on_target.xls')

    def duplicates_xls(self, read_on_results):
        # Initialize the workbook and sheet
        maxduplicates = read_on_results['maxduplicates']
        wb = xlwt.Workbook()

        # A sheet is created in the xls for each bam file
        for i, bam in enumerate(read_on_results['results']):
            ws = wb.add_sheet(bam['legend'] if len(bam['legend']) < 31 else bam['legend'][-31:])
            # Extract dictionary values
            # keys = []
            # yon = []
            # yoff = []
            # ypercon = []
            # ypercoff = []

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
            ws.write(3, 0, 'Table', header_style)

            ws.write(4, 0, '# on', header_style)
            ws.write(5, 0, '# off', header_style)
            ws.write(6, 0, '# total', header_style)
            ws.write(7, 0, '% on', header_style)
            ws.write(8, 0, '% off', header_style)

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

        wb.save(self.outdir + 'duplicates.xls')