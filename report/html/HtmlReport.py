import shutil
import string
import time

# DATASRC = os.path.dirname(sys.argv[0])+'/html/'
# IMGSRC = os.path.dirname(sys.argv[0])+'/img/'
#TMP = /tmp/

class HtmlReport:
    def __init__(self, outdir):
        self.outdir = outdir
        self.sections = []

    def addsection(self, name, reporter):
        self.sections[name] = reporter

    def report(self):
        shutil.copy(IMGSRC + '/xls_icon.png', outdir + '/img')
        shutil.copy(IMGSRC + '/txt_icon.png', outdir + '/img')
        shutil.copy(IMGSRC + '/ok.jpg', outdir + '/img')
        shutil.copy(IMGSRC + '/warning.jpg', outdir + '/img')
        shutil.copy(IMGSRC + '/coverage_histogram_example.png', outdir + '/img')
        shutil.copy(DATASRC +'/styles.css', outdir)


        #Opening the data html template in order to write the data

        fd = open(DATASRC + '/captureQC.html')



        #esta parte iría en el main
        #Write the input files information

        reportcontent = string.join(fd.readlines(), sep='').\
            replace('bamfilename',string.join(bamfilenames, sep=', ')).replace('bedfilename', bedfilename).replace('reportdate', time.ctime()).replace('reference',
                    str(reference)).replace('saturationcurve', saturationcurve).replace('nthreads', str(nthreads)).replace('tmpdir', TMP)
        fd.close()




        # escribe la plantilla
        # Substitulle patrones clave que corresponda
        reportcontent.replace('%sensitivity'), self.sections['distribution'].summary()['lo que sea']
        # substituye las imagenes con ifram. Html
        reportcontent.replace('%sensitivity-loquesea%', self.sections['distribution'].getSensitivityPlotUri())
        #

