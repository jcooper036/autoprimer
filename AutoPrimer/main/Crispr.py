#! /usr/bin/env python3

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
        self.fprimers = []
        self.rprimers = []
        self.seq = None
        self.start = None
        self.stop = None
        self.primerTag = None
        self.gene = None
