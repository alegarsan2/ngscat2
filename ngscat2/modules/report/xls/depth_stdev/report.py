import xlwt

class Report():
    def __init__(self, outdir):
        self.outdir = outdir

    def report(self, stdlists, region_stddistribution_result):
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
        for idx, results in enumerate(region_stddistribution_result['results']):
            ws.write(1 + idx, 0, results['bamfilename'], header_style)
            ws.write(1 + idx, 1, results['percentile']['Q1'])
            ws.write(1 + idx, 2, results['percentile']['Q2'])
            ws.write(1 + idx, 3, results['percentile']['Q3'])
            ws.write(1 + idx, 4, results['max'])
            ws.write(1 + idx, 5, results['min'])
            ws.write(1 + idx, 6, results['mean'])
            ws.write(1 + idx, 7, results['median'])

        wb.save(self.outdir + '/data/std_wexons.xls')


