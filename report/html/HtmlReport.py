import shutil
import string
import time
import sys
import os
import glob

DATASRC = os.path.dirname(sys.argv[0])+'/html/'
IMGSRC = os.path.dirname(sys.argv[0])+'/img/'
#TMP = /tmp/

class HtmlReport:
    def __init__(self, outdir, options):
        self.outdir = outdir
        self.sections = {}
        self.options = options
    def addsection(self, name, reporter):
        self.sections[name] = reporter

    def report(self):
        shutil.copy(IMGSRC + 'xls_icon.png', self.outdir + '/img')
        shutil.copy(IMGSRC + 'txt_icon.png', self.outdir + '/img')
        shutil.copy(IMGSRC + 'ok.jpg', self.outdir + '/img')
        shutil.copy(IMGSRC + 'warning.jpg', self.outdir + '/img')
        shutil.copy(IMGSRC + 'coverage_histogram_example.png', self.outdir + '/img')
        shutil.copy(DATASRC + 'styles.css', self.outdir)


        #Opening the data html template in order to write the data

        fd = open(DATASRC + 'captureQC.html')


        #esta parte iría en el main
        #Write the input files information

        reportcontent = string.join(fd.readlines(), sep='').\
            replace('bamfilename',string.join(bamfilenames, sep=', ')).replace('bedfilename', bedfilename).replace('reportdate', time.ctime()).replace('reference',
                    str(reference)).replace('saturationcurve', saturationcurve).replace('nthreads', str(nthreads)).replace('tmpdir', TMP)
        fd.close()


        # escribe la plantilla
        # Substitulle patrones clave que corresponda


        chromosomeimages = ''
        ontarget_coverage_files = glob.glob(self.outdir + '/data/*_Ontarget_Coverage.png')
        ontarget_coverage_files.sort()



        #TODO Coverage per chromosome, add 'covperposition'.plot_dir[]
        chromosomeimages = ''
        for afile in self.sections['covperposition'].plot_dir:
            #chromosomeimages += '<a href="data/' + os.path.basename(afile) + '"><img style="width: 33%; float: left;" src="data/' + os.path.basename(afile) + '" /></a>'
            chromosomeimages += '<a href="data/' + os.path.basename(afile) + '"><img style="width: 33%; float: left;" src="data/' + os.path.basename(afile) + '" /></a>'

        reportcontent = reportcontent.replace('<CHROMOSOMEIMAGES>', chromosomeimages)

        # Esta parte corresponde con la parte de

        if 'threshold' in self.sections:
            reportcontent = reportcontent.replace('<TARGETBASESSTATUS>', self.sections['threshold'].summary()['targetstatus'][0])
            reportcontent = reportcontent.replace('<WARNBASESCOVERED>', self.sections['threshold'].summary()['warnthreshold'])

        if 'onoff' in self.sections:
            percentagestr = '\n<ul>'
            enrichmentstr = '\n<ul>'
            for i, bamfilename in enumerate(bamfilenames):
                percentagestr += '<li>' + bamfilename + ': %.1f' % (self.sections['onoff'].summary()['percontotal'][i]) + '%</li>\n'
                enrichmentstr += '<li>' + bamfilename + ': %.1f' % (self.sections['onoff'].summary()['enrichment'][i]) + '</li>\n'
            percentagestr += '</ul>'
            enrichmentstr += '</ul>'
            reportcontent = reportcontent.replace('<PERCENTAGEONTARGET>', percentagestr)
            reportcontent = reportcontent.replace('<ENRICHMENT>', enrichmentstr)
            reportcontent = reportcontent.replace('<WARNONTARGET>', self.sections['onoff'].summary()['warnthreshold'])
            reportcontent = reportcontent.replace('<ONOFFSTATUS>', self.sections['onoff'].summary()['onoff_status'][0])

            dupimages = ''
            for plotdir in self.sections['onoff'].plot_dir_duplicates:
                dupimages += '<img style="width: 50%; float: left;"' + plotdir + '" /></a>'

            reportcontent = reportcontent.replace('<DUPIMAGES>', dupimages)
            reportcontent = reportcontent.replace('<DUPSTATUS>',  self.sections['onoff'].summary()['duplicates_status'][0])

        if 'distribution' in self.sections:
            reportcontent = reportcontent.replace('<WARNMEANCOVERAGE>',  self.sections['distribution'].summary()['warnthreshold'][0])
            reportcontent = reportcontent.replace('<DISTRIBUTIONSTATUS>', self.sections['distribution'].summary()['status'][0])



        if 'covcor' in self.sections:
            fd = open(DATASRC + '/coveragecorr_content.html')
            coveragecorr_content = fd.read()
            fd.close()
            reportcontent = reportcontent.replace('<COVERAGECORRCONTENT>', coveragecorr_content)
            reportcontent = reportcontent.replace('<WARNCOVERAGECORRELATION>', self.sections['covcor'].summary()['warnthreshold'])
            reportcontent = reportcontent.replace('<COVERAGECORRSTATUS>', self.sections['covcor'].summary()['status'])

        else:
            reportcontent = reportcontent.replace('<COVERAGECORRCONTENT>', '\n')


        #FIXME esta parte no vá en el report final, numero de veces que está 0 warning

        # reportcontent = reportcontent.replace('<WARNCOVERAGEREGION>', str(config.warncoverageregion))
        # reportcontent = reportcontent.replace('<WARNCOVERAGETHRESHOLD>', str(config.warncoveragethreshold))
        # if (throughtarget_status.value):
        #     reportcontent = reportcontent.replace('<COVERAGETHROUGHTARGETSTATUS>', 'ok')
        # else:
        #     reportcontent = reportcontent.replace('<COVERAGETHROUGHTARGETSTATUS>', 'warning')


        if 'stdev' in self.sections:
            reportcontent = reportcontent.replace('<WARNSTD>', self.sections['stdev'].summary()['warnthreshold'])
            reportcontent = reportcontent.replace('<COVERAGESTDSTATUS>', self.sections['stdev'].summary()['status'])


        if 'saturation' in self.sections:
            fd = open(DATASRC + '/saturation_content.html')
            saturation_content = fd.read()
            fd.close()
            reportcontent = reportcontent.replace('<SATURATIONCONTENT>',
                            saturation_content).replace('<DEPTHLIST>',string.join(map(str, depthlist[:-1]), sep='x10<sup>6</sup>, ')
                                                        + 'x10<sup>6</sup> and ' + str(depthlist[-1]) + 'x10<sup>6</sup>').replace('depthlist', str(depthlist)[1:-1])

            reportcontent = reportcontent.replace('<WARNSATURATION>', self.sections['saturation'].summary()['warnthreshold'])

            reportcontent = reportcontent.replace('<COVERAGESATURATIONSTATUS>',  self.sections['saturation'].summary()['status'])
            reportcontent = reportcontent.replace('coveragethrs', string.join(map(str, self.sections['saturation'].summary()['covthreshold']), sep=', '))

        else:
            reportcontent = reportcontent.replace('<SATURATIONCONTENT>', '\n').replace('depthlist', 'None')



        if 'gc_bias' in self.sections:
            fd = open(DATASRC + '/gcbias_content.html')
            gcbias_content = fd.read()
            fd.close()
            reportcontent = reportcontent.replace('<GCBIASCONTENT>', gcbias_content)

            gcbiasimages = ''
            for afile in self.sections['gc_bias'].plot_dir:
                gcbiasimages += '<img style="width:40%" src="data/' + os.path.basename(afile) + '" />'
            reportcontent = reportcontent.replace('<GCBIASIMAGES>', gcbiasimages)

        else:
            reportcontent = reportcontent.replace('<GCBIASCONTENT>', '\n')

        fd = open(outdir + '/captureQC.html', 'w')
        fd.write(reportcontent)
        fd.close()

        print('Results written at ' + self.outdir)