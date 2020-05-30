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
            self.step = 'Finding primers excluding homopolymeric regions'
            self.find_primers(method = 'second')
            self.sort_primers(method = 'second')
            self.output(method = 'second')
        self.step = 'Evaluating Output'
        if self.evaluate_new_criteria():
            print('Finding primers using relaxed criteria')
            self.step = 'Finding primers with relaxed criteria'
            self.find_primers(method = 'third')
            self.sort_primers(method = 'third')
            self.output(method = 'third')

    def evaluate_new_criteria(self):
        """
        evaluate need to find new primers using relaxed thermodynamic criteria
        """
        relaxed_criteria = False
        for cr in self.gene.crisprs:
            cr.loosened_f = False
            cr.loosened_r = False
            if (len(cr.best_fprimers) < 1) or (cr.new_fasta_L and len(cr.best_fprimers2) < 1):
                cr.loosened_f = True
                relaxed_criteria = True
            if (len(cr.best_rprimers) < 1) or (cr.new_fasta_R and len(cr.best_rprimers2) < 1):
                cr.loosened_r = True
                relaxed_criteria = True

        return relaxed_criteria

    def evaluate_new_fasta(self):
        """
        evaluates need design additional primers
        """
        min_internal_buffer = 45
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
            # If sequencing is blocked, move the outside buffer
            # If PCR is blocked, move inside and outside buffers
            # move inside buffer right away, outside with separate logic later
            cr.new_fasta_L = False
            cr.new_fasta_R = False
            if F_PCR != 10000:
                cr.new_fasta_L = True
                cr.in_buffer_L = min_internal_buffer
                if R_PCR != 10000:
                    cr.new_fasta_R = True
                    cr.in_buffer_R = min_internal_buffer
                elif R_seq!= 10000:
                    cr.new_fasta_R = True
            elif F_seq != 10000:
                #print('F_seq', F_seq)
                if R_PCR != 10000:
                    cr.new_fasta_R = True
                    cr.new_fasta_L = True
                    cr.in_buffer_R = min_internal_buffer
                elif R_seq!= 10000:
                    cr.new_fasta_R = True
                    cr.new_fasta_L = True
            elif not fnum and (R_seq != 10000):
                cr.new_fasta_R = True
            elif R_PCR != 10000:
                cr.new_fasta_R = True
                cr.in_buffer_R = min_internal_buffer

            # define new outside buffers here
            if cr.new_fasta_L or cr.new_fasta_R:
                new_fasta = True
                cr.out_buffer_R = min([R_seq, R_PCR, int(cr.out_buffer_R)])
                cr.out_buffer_L = min([F_seq, F_PCR, int(cr.out_buffer_L)])
                #print('new boundaries for ',cr.name,cr.out_buffer_L,cr.out_buffer_R)
        return new_fasta

    def handle_fasta(self):
        assert self.file
        self.fasta = ntp.read_fasta(self.file)
        self.gene = ntp.parse_fasta(self.fasta, self.gene)
        # terminal output for clarity
        print('\n\n######################\nFinding primers for the CRISPR sites in {}'.format(self.gene))

    def find_primers(self, method = 'first'):
        """
        Call primer3 to find primers
        """
        assert self.gene.cds
        assert self.gene.crisprs
        if method == 'first': self.gene.find_primers()
        elif (method == 'second') or (method == 'third'):
            self.gene.find_primers(method = method)
        else: print ('method not recognized')

    def sort_primers(self, method = 'first'):
        """
        Sorts through potential primers to find good ones
        """
        if method == 'first': self.gene.sort_primers()
        elif (method == 'second') or (method == 'third'):
            self.gene.sort_primers(method = method)
        else: print ('method not recognized')

    def write_errors(self):
        """
        Writes any errors to the gene folder in a new file
        """
        self.gene.write_errors(self.file)

    def output(self, method = 'first'):
        """
        Writes a csv and returns information to autoprimer.
        """
        header = 'GENE,CRISPR,PRIMER,SEQUENCE,SIDE,START,TM,GC%,DISTANCE_TO_CRISPR,MAX_AT_HOMOPOLYMER,MAX_GC_HOMOPOLYMER,FLAG'
        if method == 'first':
            self.out, self.all_out = self.gene.sort_output()
            out1 = self.out
            out2 = self.all_out
            fname1 = self.file.split('.fasta')[0] + '_autoprimer_out.csv'
            fname2 = self.file.split('.fasta')[0] + '_autoprimer_out_raw.csv'
        elif method == 'second':
            self.out2, self.all_out2 = self.gene.sort_output(method = method)
            out1 = self.out2
            out2 = self.all_out2
            fname1 = self.file.split('.fasta')[0] + '_autoprimer_out_2.csv'
            fname2 = self.file.split('.fasta')[0] + '_autoprimer_out_raw_2.csv'
        elif method == 'third':
            self.out3, self.all_out3 = self.gene.sort_output(method = method)
            out1 = self.out3
            out2 = self.all_out3
            fname1 = self.file.split('.fasta')[0] + '_autoprimer_out_3.csv'
            fname2 = self.file.split('.fasta')[0] + '_autoprimer_out_raw_3.csv'
        else: print('output method not recognized')
        ntp.write_csv(fname1, out1, header=header)
        ntp.write_csv(fname2, out2, header=header)

if __name__ == "__main__":
    Submission(sys.argv[1])
