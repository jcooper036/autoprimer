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
            primer_raw = ntp.pick_primers(leftseg, 'left')
            print(cr.name, primer_raw)
