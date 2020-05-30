#! /usr/bin/env python3

import AutoPrimer as ntp

class Crispr(object):
    """
    Holds information about a CRISPR site
    - sequence
    - naming information
    - primers that surround that site
    - gene it goes with
    """
    def __init__(self, name):
        self.name = name
        self.fprimers = {}
        self.rprimers = {}
        self.fprimers2 = {}
        self.rprimers2 = {}
        self.fprimers3 = {}
        self.rprimers3 = {}
        self.best_fprimers = {}
        self.best_rprimers = {}
        self.best_fprimers2 = {}
        self.best_rprimers2 = {}
        self.best_fprimers3 = {}
        self.best_rprimers3 = {}
        self.seq = None
        self.start = None
        self.stop = None
        self.primerTag = None
        self.gene = None
        self.fprimercount = 0
        self.rprimercount = 0
        self.fprimercount2 = 0
        self.rprimercount2 = 0
        self.fprimercount3 = 0
        self.rprimercount3 = 0
        self.complete = False
        self.max_primer_number = 1
        self.out_buffer_L = 500
        self.out_buffer_R = 500
        self.in_buffer_L = 150
        self.in_buffer_R = 150


    def sort_primers(self, method = 'first'):
        """
        Iterates over the left and right primers to evaluate quality.
        """
        assert self.gene

        print(f'\nAnalyzing primers for {self.name}')
        
        if method == 'first':
            fprimers = self.fprimers
            rprimers = self.rprimers
            best_fprimers = self.best_fprimers
            best_rprimers = self.best_rprimers
            fprimercount = self.fprimercount
            rprimercount = self.rprimercount

        elif method == 'second':
            fprimers = self.fprimers2
            rprimers = self.rprimers2
            best_fprimers = self.best_fprimers2
            best_rprimers = self.best_rprimers2
            fprimercount = self.fprimercount2
            rprimercount = self.rprimercount2

        elif method == 'third':
            fprimers = self.fprimers3
            rprimers = self.rprimers3
            best_fprimers = self.best_fprimers3
            best_rprimers = self.best_rprimers3
            fprimercount = self.fprimercount3
            rprimercount = self.rprimercount3
        else: print ('sort_primers method unknown')

        for pr in fprimers:
            start, stop = ntp.find_match(self.gene.cds, fprimers[pr]['pr'].seq)
            fprimers[pr]['pr'].start = start
            fprimers[pr]['pr'].stop = stop

        for pr in rprimers:
            start, stop = ntp.find_match(self.gene.cds, rprimers[pr]['pr'].seq)
            rprimers[pr]['pr'].start = start
            rprimers[pr]['pr'].stop = stop

        # For some reason, the call to evaluate_primers must be isolated to work
        if method == 'first':
            self.best_fprimers = ntp.evaluate_primers(self.fprimers)
            self.best_rprimers = ntp.evaluate_primers(self.rprimers)

        if method == 'second':
            self.best_fprimers2 = ntp.evaluate_primers(self.fprimers2)
            self.best_rprimers2 = ntp.evaluate_primers(self.rprimers2)

        if method == 'third':
            self.best_fprimers3 = ntp.evaluate_primers(self.fprimers3)
            self.best_rprimers3 = ntp.evaluate_primers(self.rprimers3)

        fprimercount = len(best_fprimers)
        rprimercount = len(best_rprimers)
