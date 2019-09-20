#! /usr/bin/env python3

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