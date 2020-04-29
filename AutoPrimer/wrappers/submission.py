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
        self.step = 'Writing Output'
        self.output()
        self.step = 'Evaluating Output'
        if self.evaluate_new_fasta():
            self.step = 'Rerunning against homopolymers'
            self.find_secondary_primers()
            self.sort_secondary_primers()
            self.secondary_output()

    def secondary_output(self):
        """
        Writes a csv and returns information to autoprimer.
        """
        self.out2, self.all_out2 = self.gene.sort_output_2()
        header = 'GENE,CRISPR,PRIMER,SEQUENCE,SIDE,START,TM,GC%,DISTANCE_TO_CRISPR,MAX_AT_HOMOPOLYMER,MAX_GC_HOMOPOLYMER,FLAG'
        fname = self.file.split('.fasta')[0] + '_autoprimer_out_2.csv'
        ntp.write_csv(fname, self.out2, header=header)
        fname = self.file.split('.fasta')[0] + '_autoprimer_out_raw_2.csv'
        ntp.write_csv(fname, self.all_out2, header=header)

    def sort_secondary_primers(self):
        """
        Sorts through potential primers to find good ones
        """
        self.gene.sort_primers_2()

    def find_secondary_primers(self):
        """
        Call primer3 to find primers with restricted sequence
        """
        assert self.gene.cds
        assert self.gene.crisprs
        self.gene.find_primers_2()

    def evaluate_new_fasta(self):
        """
        evaluates need design additional primers
        """
        new_fasta = False
        for cr in self.gene.crisprs:
            target = cr.start
            fnum = rnum = F_seq = R_seq = F_PCR = R_PCR = 0

            for pr in cr.best_fprimers:
                if not cr.best_fprimers[pr]['pr'].fail_case:
                    AT_runs = cr.best_fprimers[pr]['pr'].maxAT
                    GC_runs = cr.best_fprimers[pr]['pr'].maxGC
                    all_runs = AT_runs+';'+GC_runs
                    F_seq = ntp.eval_homopolymers(all_runs, target, 7)
                    F_PCR = ntp.eval_homopolymers(all_runs, target, 11)
                    fnum += 1
            for pr in cr.best_rprimers:
                if not cr.best_rprimers[pr]['pr'].fail_case:
                    AT_runs = cr.best_rprimers[pr]['pr'].maxAT
                    GC_runs = cr.best_rprimers[pr]['pr'].maxGC
                    all_runs = AT_runs+';'+GC_runs
                    R_seq = ntp.eval_homopolymers(all_runs, target, 7)
                    R_PCR = ntp.eval_homopolymers(all_runs, target, 11)
                    rnum += 1

            # Make a decision on what to do based on homopolymers:
            cr.new_fasta_L = False
            cr.new_fasta_R = False
            if F_PCR != 10000:
                #print('F_PCR',F_PCR)
                cr.new_fasta_L = True
                if (R_PCR != 10000) or (R_seq!= 10000):
                    #print('R_PCR or R_seq',R_PCR,R_seq)
                    cr.new_fasta_R = True
            elif F_seq != 10000:
                #print('F_seq', F_seq)
                if (R_PCR != 10000) or (R_seq!= 10000):
                    #print('R_PCR or R_seq',R_PCR,R_seq)
                    cr.new_fasta_L = True
                    cr.new_fasta_R = True
            elif not fnum and (R_seq != 10000):
                #print('not fnum and R_seq')
                cr.new_fasta_R = True
            elif R_PCR != 10000:
                #print('R_PCR',R_PCR)
                cr.new_fasta_R = True

            # define new boundaries to be searched
            if cr.new_fasta_L or cr.new_fasta_R:
                new_fasta = True
                cr.new_R = min([R_seq, R_PCR, int(self.gene.end_buffer)])
                cr.new_L = min([F_seq, F_PCR, int(self.gene.end_buffer)])
                print('new boundaries for ',cr.name,cr.new_L,cr.new_R)
        return new_fasta


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

    def write_errors(self):
        """
        Writes any errors to the gene folder in a new file
        """
        self.gene.write_errors(self.file)

    def output(self):
        """
        Writes a csv and returns information to autoprimer.
        """
        self.out, self.all_out = self.gene.sort_output()
        header = 'GENE,CRISPR,PRIMER,SEQUENCE,SIDE,START,TM,GC%,DISTANCE_TO_CRISPR,MAX_AT_HOMOPOLYMER,MAX_GC_HOMOPOLYMER,FLAG'
        fname = self.file.split('.fasta')[0] + '_autoprimer_out.csv'
        ntp.write_csv(fname, self.out, header=header)
        fname = self.file.split('.fasta')[0] + '_autoprimer_out_raw.csv'
        ntp.write_csv(fname, self.all_out, header=header)

if __name__ == "__main__":
    Submission(sys.argv[1])

    # potential expansion of search, not used right now
    # def expand_search(self):
    #     """
    #     Asks the gene which CRISPRs need their search expanded.
    #     Then kicks off iterative rounds of the expanded search.
    #     """
    #     expansion = self.gene.completeness_check()
    #     expanded_search_parameters = [(600,500), (700,600), (800,700), (900,800), (1000,900)]
    #     while expansion and expanded_search_parameters:

    #         # terminal print for process clarity
    #         print(f'Expanding search for {self.gene} to {expanded_search_parameters[0][1]}-{expanded_search_parameters[0][0]} from target site')

    #         # change the search parameters
    #         self.gene.end_buffer, self.gene.inside_buffer = expanded_search_parameters.pop(0)

    #         # search the primers
    #         self.find_primers()

    #         # sort the primers
    #         self.sort_primers()

    #         # check for expansion again
    #         expansion = self.gene.completeness_check()
