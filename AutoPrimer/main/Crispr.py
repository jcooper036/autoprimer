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
        self.best_fprimers = {}
        self.best_rprimers = {}
        self.seq = None
        self.start = None
        self.stop = None
        self.primerTag = None
        self.gene = None
        self.fprimercount = 0
        self.rprimercount = 0
        self.complete = False
        self.max_primer_number = 1
    
    def sort_primers(self):
        """
        Iterates over the left and right primers to evaluate quality.
        """

        assert self.gene

        print(f'\nAnalyzing primers for {self.name}')

        # get the real starting postitions (primer3 is bad at this)
        for pr in self.fprimers:
            start, stop = ntp.find_match(self.gene.cds, self.fprimers[pr]['pr'].seq)
            self.fprimers[pr]['pr'].start = start
            self.fprimers[pr]['pr'].stop = stop
        
        for pr in self.rprimers:
            start, stop = ntp.find_match(self.gene.cds, self.rprimers[pr]['pr'].seq)
            self.rprimers[pr]['pr'].start = start
            self.rprimers[pr]['pr'].stop = stop
        
        # for the forward
        # print('Searching for FORWARD primers:')
        self.best_fprimers = ntp.evaluate_primers(self.fprimers)
        self.fprimercount = len(self.best_fprimers)

        # print('Searching for REVERSE primers:')
        self.best_rprimers = ntp.evaluate_primers(self.rprimers)
        self.rprimercount = len(self.best_rprimers)