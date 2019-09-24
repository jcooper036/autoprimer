#!/usr/bin/env python3

class Primer(object):

    def __init__(self, side, pen, seq, start, length, tm, gc, anyTH, endTH, hairpin, stab):
        self.side = side
        self.pen = pen
        self.seq = seq
        self.start = start
        self.stop = None
        self.length = length
        self.tm = tm
        self.gc = gc
        self.anyTH = anyTH
        self.endTH = endTH
        self.hairpin = hairpin
        self.stab = stab
        self.crispr = None
        self.fail_case = None
    
    def output(self):
        out = "{},{},{},{},{}".format(self.seq, self.side, self.start, self.tm, self.gc)
        return out