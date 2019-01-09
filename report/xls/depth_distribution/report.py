import xlwt

class Report():
    def __init__(self, outdir):
        self.outdir = outdir

    def report(self, coveragefiles ,target_distribution_result):

        wb = xlwt.Workbook()

        # A sheet is created in the xls for each result file
        for i, result in enumerate(target_distribution_result):
            ws = wb.add_sheet(coveragefiles[i].name.decode('utf-8').split('/')[-1]
                              if len(coveragefiles[i].name.decode('utf-8').split('/')[-1]) < 31
                                else coveragefiles[i].name.decode('utf-8').split('/')[-1][-31:])

            # Create header font
            header_style = xlwt.easyxf('font: bold on')

            ws.write(0, 0, 'Input bamfile: ', header_style)
            ws.write(0, 1, result['bamfilename'])

            # TODO Añadir a target_coverage argumento bed para obtener bedfilename
            # ws.write(1, 0, 'Input bedfile:', header_style)
            # ws.write(1, 1, read_on_results['bedfile'])

            ws.write(2, 0, 'Outdir:', header_style)
            ws.write(2, 1, self.outdir)

            ws.write(4, 0, 'Table', header_style)

            ws.write(5, 0, 'Number bases coverage 0', header_style)
            ws.write(5, 1, 'Q1', header_style)
            ws.write(5, 2, 'Q2', header_style)
            ws.write(5, 3, 'Q3', header_style)
            ws.write(5, 4, 'Maximum', header_style)
            ws.write(5, 5, 'Minimum', header_style)
            ws.write(5, 6, 'Median', header_style)
            ws.write(5, 7, 'Mean', header_style)

            ws.write(6, 0, result['zerocov'])
            ws.write(6, 1, result['percentile']['Q1'])
            ws.write(6, 2, result['percentile']['Q2'])
            ws.write(6, 3, result['percentile']['Q3'])
            ws.write(6, 4, result['max'])
            ws.write(6, 5, result['min'])
            ws.write(6, 6, result['median'])
            ws.write(6, 7, result['mean'])

        wb.save(self.outdir + '/data/percentile.xls')

