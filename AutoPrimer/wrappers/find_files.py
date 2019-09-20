#! /usr/bin/env python3
"""
Loops over the input location looking for .fasta files, submits those to be run as genes. This is usually called by the autobot.
"""
import AutoPrimer as ntp



def find_files(input_loc=""):
    """
    Searches the input location for fasta files. Yields file names so they can be submitted. 
    """


    pass
    # find all .fasta files in the input location (and sub-directories)
    # loop over those files, submit them
    
    # for file in fasta_files:
        # ntp.submission(file)

    # compile into master list

    # consider yielding these so that the bot can call submit