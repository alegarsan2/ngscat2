import shutil
import string
import time
import sys
import os
import glob

DATASRC = os.path.dirname(sys.argv[0])+'html/'
IMGSRC = os.path.dirname(sys.argv[0])+'img/'
#TMP = /tmp/

class HtmlReport:
    def __init__(self, outdir, options):
        self.outdir = outdir
        self.sections = {}
        self.options = options

    def addsection(self, name, reporter):
        self.sections[name] = reporter
    #FIXME don't use this argument, try to by pass throught a report class
    def report(self, options):

        #Copy all data required for html report

        shutil.copy(IMGSRC + 'xls_icon.png', self.outdir + '/img/xls_icon.png')
        shutil.copy(IMGSRC + 'txt_icon.png', self.outdir + '/img/txt_icon.png')
        shutil.copy(IMGSRC + 'close_icon.png', self.outdir + '/img/close_icon.png')
        shutil.copy(IMGSRC + 'ok.jpg', self.outdir + '/img/ok.jpg')
        shutil.copy(IMGSRC + 'warning.jpg', self.outdir + '/img/warning.jpg')
        shutil.copy(IMGSRC + 'coverage_histogram_example.png', self.outdir + '/img/coverage_histogram_example.png')
        shutil.copy(DATASRC + 'styles.css', self.outdir)
        shutil.copy(DATASRC + '/captureQC.html', self.outdir)


        #Opening the data html template in order to write the data

        fd = open(self.outdir + '/captureQC.html')
        reportcontent = ('').join(fd.readlines()) \
            # .\
        # replace('bamfilename', string.join(bamfilenames, sep=', ')).replace('bedfilename', bedfilename).replace('reportdate', time.ctime()).replace('reference',
        #         str(reference)).replace('saturationcurve', saturationcurve).replace('nthreads', str(nthreads)).replace('tmpdir', TMP)
        fd.close()



        # *************** Input parameters ***********************************
        #FIXME pass arguments throught another data type for instance Dictionary instead Options object
        reportcontent = reportcontent.replace('bamfilename',options.bams)
        reportcontent = reportcontent.replace('bedfilename',options.bed)
        reportcontent = reportcontent.replace('reportdate', time.ctime())
        reportcontent = reportcontent.replace('coveragethrs', str(options.coveragethresholds))
        reportcontent = reportcontent.replace('reference',options.reference if options.reference is not None else "No reference")
        reportcontent = reportcontent.replace('saturationcurve', options.saturation)
        reportcontent = reportcontent.replace('nthreads', str(options.nthreads))
        reportcontent = reportcontent.replace('tmpdir',options.tmp)


        # *************** Result summary ***********************************
        summaryrows = ''
        for i, bam in enumerate(self.sections['onoff'].getsummary()['bamfilename']):
            summaryrows += '<tr>\n'
            summaryrows += '<td class="table-cell"> ' + bam + '</td>'
            summaryrows += '<td class="table-cell">%.0f ' % self.sections['onoff'].getsummary()['totalread'][i] + ' </td>'
            summaryrows += '<td class="table-cell">%.1f' % self.sections['threshold'].getsummary()['perctotalcovered'][i] + '% </td>'

            if 'saturation' in self.sections:
                summaryrows += '<td class="table-cell">%.1e</td>\n' % saturationslopes[i]

            summaryrows += '<td class="table-cell">%.1f' % self.sections['onoff'].getsummary()['percontotal'][i] + '% </td>\n'

            summaryrows += ('<td class="table-cell">ON-%.1f%%' % self.sections['onoff'].getsummary()['perconduplicates'][i]['2x']) + '; OFF: %.1f' % (
                self.sections['onoff'].getsummary()['percoffduplicates'][i]['2x']) + '% </td>'

            summaryrows += '<td class="table-cell">%.1fx' % self.sections['distribution'].getsummary()['mean'][i] + '</td>\n'
            #Not added in this version
            summaryrows += '<td class="table-cell">%d consecutive bases<br>with coverage <= <WARNCOVERAGETHRESHOLD></td>\n' % (
            20)

            if 'coveragecorr' in self.sections:
                summaryrows += '<td class="table-cell">%.2f</td>\n' % corr.value

            summaryrows += '<td class="table-cell">%.2f</td>\n' % self.sections['stdev'].getsummary()['mean'][i]
            summaryrows += '</tr>\n'

        summarystatus = '<td class="table-header">Overall status</td>\n'
        summarystatus += '<td class="table-header"></td>\n'
        summarystatus += '<td class="table-header"><a href="#targetbases"><img src="img/<TARGETBASESSTATUS>.jpg" height=23px /></a></td>\n'

        if 'saturation' in self.sections:
            summarystatus += '<td class="table-header"><a href="#coveragesaturation"><img src="img/<COVERAGESATURATIONSTATUS>.jpg" height=23px /></a></td>\n'
        summarystatus += '<td class="table-header"><a href="#onoff"><img src="img/<ONOFFSTATUS>.jpg" height=23px /></a></td>\n'
        summarystatus += '<td class="table-header"><a href="#dup"><img src="img/<DUPSTATUS>.jpg" height=23px /></a></td>\n'
        summarystatus += '<td class="table-header"><a href="#distribution"><img src="img/<DISTRIBUTIONSTATUS>.jpg" height=23px /></a></td>\n'
        summarystatus += '<td class="table-header"><a href="#coveragethroughtarget"><img src="img/<COVERAGETHROUGHTARGETSTATUS>.jpg" height=23px /></a></td>\n'

        if 'coveragecorr' in self.sections:
            summarystatus += '<td class="table-header"><a href="#coveragecorr"><img src="img/<COVERAGECORRSTATUS>.jpg" height=23px /></a></td>\n'

        summarystatus += '<td class="table-header"><a href="#coveragestd"><img src="img/<COVERAGESTDSTATUS>.jpg" height=23px /></a></td>\n'

        reportcontent = reportcontent.replace('<SUMMARYROWS>', summaryrows)
        reportcontent = reportcontent.replace('<SUMMARYSTATUS>', summarystatus)

        if 'saturation' in self.sections:
            reportcontent = reportcontent.replace('<SUMMARYSATURATION>',
                                                  '<td class="table-header"><a href="#coveragesaturation">Coverage saturation<br>(slope at the end of the curve)</a></td>')
        else:
            reportcontent = reportcontent.replace('<SUMMARYSATURATION>', '')

        if 'coveragecorr' in self.sections:
            reportcontent = reportcontent.replace('<SUMMARYCOVCORRELATION>',
                                                  '<td class="table-header"><a href="#coveragecorr">Coverage correlation<br>per ROI</a></td>')
        else:
            reportcontent = reportcontent.replace('<SUMMARYCOVCORRELATION>', '')

        #reportcontent = reportcontent.replace('<SUMMARYCOVERAGETHRS>', str(coveragethresholds[0]))
        #reportcontent = reportcontent.replace('<SUMMARYTARGETSIZE>', str(bed_file.bed_file(bedfilename).size()))








        #TODO Coverage per chromosome, add 'covperposition'.plot_dir[]
        chromosomeimages = ''
        for afile in self.sections['covperposition'].plot_dir:
            #chromosomeimages += '<a href="data/' + os.path.basename(afile) + '"><img style="width: 33%; float: left;" src="data/' + os.path.basename(afile) + '" /></a>'
            chromosomeimages +=  "<iframe class =""coverageItem"" src=data/" + afile.split("/")[-1]+ " frameborder=""0""> </iframe>"

        reportcontent = reportcontent.replace('<CHROMOSOMEIMAGES>', chromosomeimages)

        # Esta parte corresponde con la parte de

        if 'threshold' in self.sections:
            reportcontent = reportcontent.replace('<TARGETBASESSTATUS>', self.sections['threshold'].getsummary()['targetstatus'][0])
            reportcontent = reportcontent.replace('<WARNBASESCOVERED>', str(self.sections['threshold'].getsummary()['warnthreshold']))

        if 'onoff' in self.sections:
            percentagestr = '\n<ul>'
            enrichmentstr = '\n<ul>'
            for i, bamfilename in enumerate(self.sections['onoff'].getsummary()['bamfilename']):
                percentagestr += '<li>' + bamfilename + ': %.1f' % self.sections['onoff'].getsummary()['percontotal'][i] + '%</li>\n'
                enrichmentstr += '<li>' + bamfilename + ': %.1f' % self.sections['onoff'].getsummary()['enrichment'][i] + '</li>\n'
            percentagestr += '</ul>'
            enrichmentstr += '</ul>'
            reportcontent = reportcontent.replace('<PERCENTAGEONTARGET>', percentagestr)
            reportcontent = reportcontent.replace('<ENRICHMENT>', enrichmentstr)
            reportcontent = reportcontent.replace('<WARNONTARGET>', str(self.sections['onoff'].getsummary()['warnthreshold']))
            reportcontent = reportcontent.replace('<ONOFFSTATUS>', self.sections['onoff'].getsummary()['onoffstatus'][0])

            dupimages = ''
            for plotdir in self.sections['onoff'].plot_dir_duplicates:
                dupimages += '<iframe src='+ "data/"+ plotdir.split("/")[-1] + ' height = "400" width = "60%" style = "float:left" frameborder="0"> </iframe>'

            reportcontent = reportcontent.replace('<DUPIMAGES>', dupimages)
            reportcontent = reportcontent.replace('<DUPSTATUS>',  self.sections['onoff'].getsummary()['duplicates_status'][0])

        if 'distribution' in self.sections:
            reportcontent = reportcontent.replace('<WARNMEANCOVERAGE>',  str(self.sections['distribution'].getsummary()['warnthreshold'])[0])
            reportcontent = reportcontent.replace('<DISTRIBUTIONSTATUS>', self.sections['distribution'].getsummary()['status'][0])



        if 'covcor' in self.sections:
            fd = open(DATASRC + '/coveragecorr_content.html')
            coveragecorr_content = fd.read()
            fd.close()
            reportcontent = reportcontent.replace('<COVERAGECORRCONTENT>', coveragecorr_content)
            reportcontent = reportcontent.replace('<WARNCOVERAGECORRELATION>', self.sections['covcor'].getsummary()['warnthreshold'])
            reportcontent = reportcontent.replace('<COVERAGECORRSTATUS>', self.sections['covcor'].getsummary()['status'])

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
            reportcontent = reportcontent.replace('<WARNSTD>', str(self.sections['stdev'].getsummary()['warnthreshold']))
            reportcontent = reportcontent.replace('<COVERAGESTDSTATUS>', self.sections['stdev'].getsummary()['status'][0])


        if 'saturation' in self.sections:
            fd = open(DATASRC + '/saturation_content.html')
            saturation_content = fd.read()
            fd.close()
            # reportcontent = reportcontent.replace('<SATURATIONCONTENT>',
            #                 saturation_content).replace('<DEPTHLIST>',string.join(map(str, depthlist[:-1]), sep='x10<sup>6</sup>, ')
            #                                             + 'x10<sup>6</sup> and ' + str(depthlist[-1]) + 'x10<sup>6</sup>').replace('depthlist', str(depthlist)[1:-1])

            reportcontent = reportcontent.replace('<WARNSATURATION>', self.sections['saturation'].getsummary()['warnthreshold'])

            reportcontent = reportcontent.replace('<COVERAGESATURATIONSTATUS>',  self.sections['saturation'].getsummary()['status'])
            reportcontent = reportcontent.replace('coveragethrs', (',').join(map(str, self.sections['saturation'].getsummary()['covthreshold'])))

        else:
            reportcontent = reportcontent.replace('<SATURATIONCONTENT>', '\n').replace('depthlist', 'None')



        if 'gc_bias' in self.sections:
            fd = open(DATASRC + '/gcbias_content.html')
            gcbias_content = fd.read()
            fd.close()
            reportcontent = reportcontent.replace('<GCBIASCONTENT>', gcbias_content)

            gcbiasimages = ''
            for afile in self.sections['gc_bias'].plot_dir:
                gcbiasimages += '<iframe src='+ "data/"+ afile.split("/")[-1] + ' height = "400" width="60%" style = "float:left"  frameborder=0> </iframe>'
            reportcontent = reportcontent.replace('<GCBIASIMAGES>', gcbiasimages)

        else:
            reportcontent = reportcontent.replace('<GCBIASCONTENT>', '\n')

        fd = open(self.outdir + '/captureQC.html', 'w')
        fd.write(reportcontent)
        fd.close()

        print('Results written at ' + self.outdir)