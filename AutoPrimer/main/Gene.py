#! /usr/bin/env python3

import AutoPrimer as ntp
#import find_max_polymer as ex
import copy

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
        self.end_buffer = 500
        self.inside_buffer = 150
        self.error_log = []

    def __repr__(self):
        return self.name

    def find_primers(self):
        """
        Iterate over CRISPRS to find possible primers
        """
        for cr in self.crisprs:
            if not cr.complete:
                cr.start, cr.stop = ntp.find_match(self.cds, cr.seq)
                if cr.start == 'ERROR':
                    continue
                leftseg, rightseg = ntp.find_cds_segs(self.cds, cr, end_buffer=self.end_buffer, inside_buffer=self.inside_buffer)

                # need to make sure that both values are true (will fail if too close to the front or the back)
                if leftseg and rightseg:
                    left_raw = ntp.pick_primers(leftseg, 'left')
                    right_raw = ntp.pick_primers(rightseg, 'right')

                    temp = ntp.raw_to_primers(left_raw, 'left')
                    for pr in temp:
                        cr.fprimers[pr] = copy.deepcopy(temp[pr])

                    temp = ntp.raw_to_primers(right_raw, 'right')
                    for pr in temp:
                        cr.rprimers[pr] = copy.deepcopy(temp[pr])

                else:
                    self.error_log.append(f'{leftseg} - {cr.name}')
                    print(f'{leftseg} - {cr.name}')

    def sort_primers(self):
        """
        Iterate over CRISPRS to sort the primers
        """
        for cr in self.crisprs:
            cr.sort_primers()

    def completeness_check(self):
        for cr in self.crisprs:
            if (cr.fprimercount == cr.max_primer_number) and (cr.rprimercount == cr.max_primer_number):
                cr.complete = True
        return all(cr.complete for cr in self.crisprs)

    def sort_output(self):
        """
        Writes a csv of the output to the gene folder.
        Also returns the output so the submission can use it
        """
        self.out = []
        for cr in self.crisprs:
            cr_mid = abs(int((cr.start+cr.stop)/2))
            for pr in cr.best_fprimers:
                primer_name = 'F' + str(cr.best_fprimers[pr]['num'])
                primer_info = cr.best_fprimers[pr]['pr'].output()
                pr_start = int(cr.best_fprimers[pr]['pr'].start) # new here
                distance, max_AT, max_GC = ntp.find_max_polymer(self.cds, cr_mid, pr_start) # new here
                outstring = f'{self.name},{cr.name},{primer_name},{primer_info},{distance},{max_AT},{max_GC}' # new here
                self.out.append(outstring)
            for pr in cr.best_rprimers:
                primer_name = 'R' + str(cr.best_rprimers[pr]['num'])
                primer_info = cr.best_rprimers[pr]['pr'].output()
                pr_start = int(cr.best_rprimers[pr]['pr'].start) # new here
                distance, max_AT, max_GC = ntp.find_max_polymer(self.cds, cr_mid, pr_start) # new here
                outstring = f'{self.name},{cr.name},{primer_name},{primer_info},{distance},{max_AT},{max_GC}' # new here
                self.out.append(outstring)
        return self.out

    def write_errors(self, file):
        """
        Input:
            self
            file name, str that ends in a suffix (should be .fasta)
        """
        if self.error_log:
            file = file.split('.')[0] + '_errors.txt'
            with open(file, 'w') as f:
                for error in self.error_log:
                    f.write(f'{error}+\n')
