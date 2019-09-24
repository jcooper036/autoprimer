#! /usr/bin/env python3

"""
For testing AutoPrimer functions and workflow.
"""

import AutoPrimer as ntp
import sys

settings = {
    
}

folder = sys.argv[1]

# iterator that returns file paths
for file in ntp.find_files(input_loc=folder):
    ntp.Submission(file)

# # read in info

# # gather the sample and crispr targets
# infile = sys.argv[1]
# ntp.submission(sys.argv[1])

# gene = ntp.Gene(infile.split('/')[-1].split('.fasta')[0])
# fasta = ntp.read_fasta(infile)
# gene = ntp.parse_fasta(fasta, gene)
# # terminal output for clarity
# print('\n\n######################\nFinding primers for the CRISPR sites in {}'.format(gene))


# find primers with primer3

# filter those primers on certain criteria

# organize and return information