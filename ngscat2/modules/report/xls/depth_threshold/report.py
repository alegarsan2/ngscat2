import xlwt

class Report():
    def __init__(self, outdir):
        self.outdir = outdir

    def report(self, coveragefiles ,target_coverage_result):

        wb = xlwt.Workbook()

        # A sheet is created in the xls for each result file
        for i, result in enumerate(target_coverage_result):
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

            ws.write(3, 0, 'Table', header_style)

            ws.write(5, 0, 'Number of on positions covered', header_style)
            ws.write(6, 0, '% of on positions covered', header_style)
            ws.write(4, 1 + len(result['coveredposition']), 'Total', header_style)
            ws.write(5, 1 + len(result['coveredposition']), result['ntotalposition'])

            for j, depth in enumerate(result['coveredposition']):
                ws.write(4, 1 + j, depth, header_style)
                ws.write(5, 1 + j, result['coveredposition'][depth])
                ws.write(6, 1 + j, result['perccoveredposition'][depth])

        wb.save(self.outdir + '/data/coverage_summary.xls')

