#!/usr/bin/python

# import glob
# import time
import pysam
# import sets
import sys
import os
import subprocess
import optparse
import string
import numpy
# import progressbar
import xlwt
import multiprocessing
import shutil
import math

TMP = '/tmp/'

# DATASRC = '/home/javi/MGP/ngscat/data/'
# IMGSRC = '/home/javi/MGP/ngscat/data/'

DATASRC = os.path.dirname(sys.argv[0]) + '/html/'  # Ver ruta de esto, no lo tengo muy claro
IMGSRC = os.path.dirname(sys.argv[0]) + '/img/'  # esto tambien


def run(command):
    """*****************************************************************************************************************
    Task: launches a system call
    Inputs:
        command: string containing the system call.
    *****************************************************************************************************************"""

    print('CMD: ' + command)
    # Checks whether an error occurred during the execution of the system call
    try:
        subprocess.check_call(command)
    except:
        print('Some error ocurred while executing: ')
        print('', command)


def launch_coveragebed(bamfilenames, bedfilename, legend, outdir, executiongranted):
    global TMP

    coveragefiles = []
    Pcoveragebeds = []
    pid = str(os.getpid())

    for i, bamfilename in enumerate(bamfilenames):
        coveragefile = TMP + '/' + os.path.basename(bamfilename).replace('.bam', '.' + pid + '.coverage')
        coveragebedgraph = outdir + '/data/' + legend[i].replace('.bam', '.bed')

        print('Coveragefile = ' + coveragefile)
        bam = bam_file.bam_file(bamfilename, 'rb')

        print('Launching coverageBed...')
        Pcoveragebed = multiprocessing.Process(target=bam.myCoverageBed,
                                               args=(bedfilename, None, coveragefile, executiongranted,
                                                     TMP, coveragebedgraph,))
        Pcoveragebed.start()

        print('Done.')

        coveragefiles.append(coveragefile)
        Pcoveragebeds.append(Pcoveragebed)

    return [Pcoveragebeds, coveragefiles]


def launch_onoff_reads(bamfilenames, bedfilename, legend, outdir, executiongranted):
    global TMP

    #Mirar aqui para sacar las variables globales, o dejar el diccionario en memoria. A priori el diccionario
    #que pasa no se pisa por lo tanto podria meterlo en el diccionario para no tener on off, duplicates...
    # onoff_status = multiprocessing.Value('b', False)
    # duplicates_status = multiprocessing.Value('b', False)
    # enrichment = multiprocessing.Array('f', len(bamfilenames))
    # percontarget = multiprocessing.Array('f', len(bamfilenames))
    # onduplicates = multiprocessing.Array('f', len(bamfilenames))
    # offduplicates = multiprocessing.Array('f', len(bamfilenames))

    bam = bam_file.bam_file(bamfilenames[0], 'rb')
    print('Launching on/off target enrichment calculation...')
    Ponoff_reads = multiprocessing.Process(target=bam.reads_on_target,
                                           args=(bedfilename, outdir,
                                            [bam_file.bam_file(bamfilenames[i]) for i in range(1, len(bamfilenames))],
                                            legend, executiongranted, onoff_status, duplicates_status, onduplicates,
                                            offduplicates, enrichment, percontarget, TMP, config.warnontarget,))
    Ponoff_reads.start()
    bam.close()

    return Ponoff_reads, onoff_status, onduplicates, offduplicates, duplicates_status, enrichment, percontarget



def launch_covered_positions(coveragefiles, coveragethresholds, outdir, legend, executiongranted):
    status = multiprocessing.Value('b', False)
    coveredbases = multiprocessing.Array('f', len(coveragefiles))


    Pcoveredpositions = multiprocessing.Process(target=target_coverage.target_coverage_lite,
                                                args=(coveragefiles, coveragethresholds, outdir, legend, None,
                                                      executiongranted, status, coveredbases, config.warnbasescovered,))

    print('Launching covered positions calculation...')
    Pcoveredpositions.start()

    return Pcoveredpositions, status, coveredbases

