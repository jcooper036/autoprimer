#! /usr/bin/env python3

import AutoPrimer as ntp
import sys

def submission(file):
    """
    Input: Fasta file
    Main wrapper that coordinates finding primers for a gene
    """

    fasta = npt.read_fasta(file)
    print(fasta)
    # find primers with primer3
    

    # filter those primers on certain criteria

    # organize and return information

if __name__ == "__main__":
    submission(sys.argv[1])