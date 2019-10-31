#! /usr/bin/env python3

import AutoPrimer as ntp
import sys

class Submission(object):
    """
    Submissions are done on a gene level. They take a folder that contains a fasta file for that gene which contains the cds and the crispr sequences. See the examples for details.

    They are meant so that parrallelization can be achieved on the gene level.
    """

    def __init__(self, file, run=True):
        self.file = file
        self.gene = ntp.Gene(self.file.split('/')[-1].split('.fasta')[0])
        self.out = None
        self.step = "Init"
        # call the main set of functions
        if run:
            # allows for manual control when debugging
            self.main()
        
    def main(self):
        self.step = 'Read Fasta'
        self.handle_fasta()
        self.step = 'Find Primers'
        self.find_primers()
        self.step = 'Sorting Primers'
        self.sort_primers()
        self.step = 'Expanding Serach'
        self.expand_search()
        self.step = 'Writing Output'
        self.output()

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
        self.gene.sort_primers()

    def expand_search(self):
        """
        Asks the gene which CRISPRs need their search expanded.
        Then kicks off iterative rounds of the expanded search.
        """
        expansion = self.gene.completeness_check()
        expanded_search_parameters = [(600,500), (700,600), (800,700), (900,800), (1000,900)]
        while expansion and expanded_search_parameters:
            
            # terminal print for process clarity
            print(f'Expanding search for {self.gene} to {expanded_search_parameters[0][1]}-{expanded_search_parameters[0][0]} from target site')
            
            # change the search parameters
            self.gene.end_buffer, self.gene.inside_buffer = expanded_search_parameters.pop(0)
            
            # search the primers
            self.find_primers()
            
            # sort the primers
            self.sort_primers()

            # check for expansion again
            expansion = self.gene.completeness_check()

    def output(self):
        """
        Writes a csv and returns information to autoprimer.
        """
        self.out = self.gene.sort_output()
        header = 'GENE,CRISPR,PRIMER,SEQUENCE,SIDE,START,TM,GC%'
        fname = self.file.split('.fasta')[0] + '_autoprimer_beta_out.csv'
        ntp.write_csv(fname, self.out, header=header)

if __name__ == "__main__":
    Submission(sys.argv[1])