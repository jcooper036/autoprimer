#! /usr/bin/env python3

import AutoPrimer as ntp
import sys

class Submission(object):
    """
    Main wrapper that coordinates finding primers for a gene
    """

    def __init__(self, file):
        self.file = file
        self.gene = ntp.Gene(self.file.split('/')[-1].split('.fasta')[0])
        self.step = "Init"
        # call the main set of functions
        self.main()
        
    def main(self):
        self.step = 'Read Fasta'
        self.handle_fasta()
        self.step = 'Find Primers'
        self.find_primers()
        self.step = 'Sorting Primers'
        self.sort_primers()

    def handle_fasta(self):
        assert self.file
        self.fasta = ntp.read_fasta(self.file)
        self.gene = ntp.parse_fasta(self.fasta, self.gene)
        # terminal output for clarity
        print('\n\n######################\nFinding primers for the CRISPR sites in {}'.format(self.gene))

    def find_primers(self):
        """
        Call primer3 to find primers
        """
        assert self.gene.cds
        assert self.gene.crisprs
        self.gene.find_primers()
    
    def sort_primers(self):
        """
        Sorts through potential primers to find good ones
        """
        pass

    def output(self):
        """
        Writes a csv and returns information to autoprimer.
        """
        pass

if __name__ == "__main__":
    Submission(sys.argv[1])