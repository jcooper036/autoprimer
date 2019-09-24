#! /usr/bin/env python3
"""
Loops over the input location looking for .fasta files, submits those to be run as genes. This is usually called by the autobot.
"""
import AutoPrimer as ntp
import glob



def find_files(input_loc=""):
    """
    Searches the input location for fasta files. Yields file names so they can be submitted. 
    """

    if input_loc[-1] != '/':
        input_loc += '/'

    for f in glob.glob(input_loc + '*/*.fasta'):
        yield f