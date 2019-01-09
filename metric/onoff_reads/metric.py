import gc
import numpy as np
import bed_file
# due to pysam implementation and restrictions, we should make this
# global so we can share the memory by forking
gBamList = []
class OnOffReadsProcessor():
    def __init__(self, bamlist, maxduplicates = 5, warnthreshold = 80):
        global gBamlist
        self.results = []
        self.maxduplicates = maxduplicates
        self.warnthreshold = warnthreshold
        gBamlist = bamlist
        print('A')
    def process(self, beddir, callback):

        """************************************************************************************************************************************************************
        Task: Print reads on traget and off target
        Inputs:
            bed: bed file with capture coordinates
            outdir: Output folder
            bamlist: list of bam_file objects representing bam files which are also wanted to be analyzed.
            legend: list of strings containing the label to tag each bam file in the xls and png files that will be created.

            onoff_status: multiprocessing.Value object to return whether the number of on-target reads is extremely low (False) or not (True)
            duplicates_status: multiprocessing.Value object to return whether the number of duplicated reads on-target is greater than the number of duplicated
                off-target (False) or not (True).
            enrichment: multiprocessing.Array objecto to return the enrichment value for each bam (on-target reads per Kb)/(off target reads per Kb)
            percontarget: multiprocessing.Array object to return percentaje of reads on target for each bam file
            tmpdir: string containing the path to a temporary directory where temporary files will be stored.
        Outputs: dictionary and json with information and bam results in the key 'results'
        ************************************************************************************************************************************************************"""
        read_on_results = {}
        # global TMP

        # Calculate number of reads and duplicated reads on/off target per chromosome
        tread = []
        nread = []
        onperchr = []
        totalperchr = []
        onduplicates = []
        offduplicates = []
        percontarget = []
        perconperchr = []
        perconduplicates = []
        percoffduplicates = []
        enrichment = []

        bamlist = gBamlist
        print('A')
        #Adding data
        for bam in bamlist:
            nread_tmp, onperchr_tmp, totalperchr_tmp, onduplicates_tmp, offduplicates_tmp = bam.myReadsOnTarget(beddir)
            tread.append(bam.nreads())
            nread.append(nread_tmp)
            onperchr.append(onperchr_tmp)
            totalperchr.append(totalperchr_tmp)
            onduplicates.append(onduplicates_tmp)
            offduplicates.append(offduplicates_tmp)

        bedobj = bed_file.bed_file(beddir)
        targetsize = bedobj.size()

        retonduplicates= []
        retoffduplicates = []
        onduplicatesresult = []
        offduplicatesresult = []
        onoff_status = []
        duplicates_status = []
        legend = None
        # Generating list that will form the output dictionary and Json
        for i in range(len(bamlist)):

            # Calculate enrichment, same bed will be used for both bams
            if (tread[i] == nread[i]):
                enrichment.append(-1)
            else:
                enrichment.append((nread[i] * 1000.0 / targetsize) / ((tread[i] - nread[i]) * 1000.0 / (bamlist[0].mappingsize() - targetsize)))

            percontarget.append(nread[i]*100.0/tread[i])
            retonduplicates.append(sum(onduplicates[i]) * 100.0 / tread[i])
            retoffduplicates.append(sum(offduplicates[i]) * 100.0 / tread[i])

            #Avoid 0 division
            perconperchr.append({key: (onperchr[i][key] * 100.0/totalperchr[i][key] if totalperchr[i][key] > 0 else 0) for key in onperchr[i]})


            # Select the largest one. Use the keys and fill with 0 until key lenght
            if len(onduplicates[i]) >= len(offduplicates[i]):
                onduplicatesresult.append({str(key + 1) + 'x': value for key, value in enumerate(onduplicates[i].tolist())})

                offduplicatesresult.append(dict(zip(list(onduplicatesresult[i].keys()), offduplicates[i].tolist() + [0] *
                                               (len(onduplicatesresult[i].keys()) - len(offduplicates[i].tolist())))))
            else:
                offduplicatesresult.append({str(key + 1) + 'x': value for key, value in enumerate(offduplicates[i].tolist())})
                onduplicatesresult.append(dict(zip(list(offduplicatesresult[i].keys()), onduplicates[i].tolist() + [0] *
                                               (len(offduplicatesresult[i].keys()) - len(onduplicates[i].tolist())))))

            # Compute percentage per number of replicates, each element, 1x is number of reads that appears 1 time /
            # total on reads.
            perconduplicates.append({key: (value / nread[i] * 100.0 if nread[i] > 0 else 0 )for key, value in onduplicatesresult[i].items()})
            percoffduplicates.append({key: (value / (tread[i]- nread[i]) * 100.0 if (tread[i] - nread[i]) > 0 else 0)
                                      for key, value in offduplicatesresult[i].items()})

            onoff_status.append(True if percontarget[i] >= self.warnthreshold else False)
            duplicates_status.append(True if retonduplicates[i] > retoffduplicates[i] else False)

        #Output generation, reads_on_target.json. Status variables doesn't go inside the json (only str is allowed not booleans)
        results = []
        for i in range(len(bamlist)):
            results.append(dict(
                [('bamfilename', bamlist[i].filename.decode('utf-8')),
                 ('enrichment', enrichment[i]),
                 ('onread', nread[i]),
                 ('percontotal', percontarget[i]),
                 ('totalread', tread[i]),
                 ('onperchr', onperchr[i]),
                 ('totalperchr', totalperchr[i]),
                 ('perconperchr', perconperchr[i]),
                 ('retonduplicates',retonduplicates),
                 ('retoffduplicates',retoffduplicates),
                 ('legend', (legend[i] if legend is not None else bamlist[i].filename.decode('utf-8').split('/')[-1])),
                 ('onduplicates', onduplicatesresult[i]),
                 ('offduplicates', offduplicatesresult[i]),
                 ('perconduplicates', perconduplicates[i]),
                 ('percoffduplicates', percoffduplicates[i]),
                 ('onoff_status', 'ok' if onoff_status[i] else 'warning'),
                 ('duplicates_status', 'ok' if duplicates_status[i] else 'warning'),
                 ('warnthreshold', self.warnthreshold),
                 ('maxduplicates', self.maxduplicates)]))
        read_on_results['results'] = results
        read_on_results['bedfile'] = beddir
        read_on_results['maxduplicates'] = self.maxduplicates


        callback(read_on_results)

        # with open(self.outdir +'/read_on_results.json', 'w') as outfile:
        #     json.dump(read_on_results, outfile)
        # return read_on_results, onoff_status, duplicates_status

