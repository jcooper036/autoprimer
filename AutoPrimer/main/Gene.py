#! /usr/bin/env python3

import AutoPrimer as ntp

class Gene(object):
    """
    Holds information about a gene
    - what CRISPRS a gene has
    - the cds of that gene
    """
    def __init__(self, name):
        self.name = name
        self.crisprs = []
        self.cds = None
    
    def __repr__(self):
        return self.name
    
    def find_primers(self):
        """
        Iterate over CRISPRS to find possible primers
        """
        for cr in self.crisprs:
            cr.start, cr.stop = ntp.find_match(self.cds, cr.seq)
            leftseg, rightseg = ntp.find_cds_segs(self.cds, cr)
            left_raw = ntp.pick_primers(leftseg, 'left')
            right_raw = ntp.pick_primers(leftseg, 'right')
            cr.fprimers = ntp.raw_to_primers(left_raw, 'left')
            cr.rprimers = ntp.raw_to_primers(right_raw, 'right')
    
    def sort_primers(self):
        """
        Iterate over CRISPRS to sort the primers
        """
        for cr in self.crisprs:
            cr.sort_primers()
    
    def sort_output(self):
        """
        Writes a csv of the output to the gene folder.
        Also returns the output so the submission can use it
        """
        self.out = []
        for cr in self.crisprs:
            for pr in cr.best_fprimers:
                primer_name = 'F' + str(cr.best_fprimers[pr]['num'])
                primer_info = cr.best_fprimers[pr]['pr'].output()
                outstring = f'{self.name},{cr.name},{primer_name},{primer_info}'
                self.out.append(outstring)
            for pr in cr.best_rprimers:
                primer_name = 'F' + str(cr.best_rprimers[pr]['num'])
                primer_info = cr.best_rprimers[pr]['pr'].output()
                outstring = f'{self.name},{cr.name},{primer_name},{primer_info}'
                self.out.append(outstring)                
        return self.out
