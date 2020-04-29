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

    def sort_output_2(self):
        """
        Writes a csv of the output to the gene folder.
        Also returns the output so the submission can use it
        """

        self.out2 = []
        self.all_out2 = []
        for cr in self.crisprs:
            if not cr.new_fasta_L and not cr.new_fasta_L: continue
            cr_mid = abs(int((cr.start+cr.stop)/2))
            cr_name = cr.name

            for pr in cr.best_fprimers2:
                primer = cr.best_fprimers2[pr]
                primer_name = 'F' + str(primer['num'])
                good, outstring = self.eval_primer(primer, primer_name, cr_name, cr.start, cr.stop)
                self.all_out2.append(outstring)
                if good:
                    self.out2.append(outstring)

            for pr in cr.best_rprimers2:
                primer = cr.best_rprimers2[pr]
                primer_name = 'R' + str(primer['num'])
                good, outstring = self.eval_primer(primer, primer_name, cr_name, cr.start, cr.stop)
                self.all_out2.append(outstring)
                if good:
                    self.out2.append(outstring)

        return self.out2, self.all_out2

    def sort_primers_2(self):
        """
        Iterate over CRISPRS to sort the primers
        """
        for cr in self.crisprs:
            if cr.new_fasta_L or cr.new_fasta_L:
                cr.sort_primers_2()

    def find_primers_2(self):
        """
        Iterate over CRISPRS to find possible primers with restricted sequence
        """
        for cr in self.crisprs:
            if not cr.new_fasta_L and not cr.new_fasta_L: continue
            if cr.start == 'ERROR':
                continue
            leftseg, rightseg = ntp.find_cds_segs(self.cds, cr, left_buffer=cr.new_L, right_buffer=cr.new_R, inside_buffer=self.inside_buffer)#cr.new_L cr.new_R,
            print (cr.name,'left buffer',cr.new_L,'right buffer',cr.new_R)
            print('len leftseg',len(leftseg),'len rightseg',len(rightseg))

            # qualify segments individually if they have new endpoints
            if leftseg and (cr.new_L != self.end_buffer):
                print('finding leftseg')
                left_raw = ntp.pick_primers(leftseg, 'left')
                temp = ntp.raw_to_primers(left_raw, 'left')
                for pr in temp:
                    cr.fprimers2[pr] = copy.deepcopy(temp[pr])

            if rightseg and (cr.new_R != self.end_buffer):
                print('finding rightseg')
                right_raw = ntp.pick_primers(rightseg, 'right')
                temp = ntp.raw_to_primers(right_raw, 'right')
                for pr in temp:
                    cr.rprimers2[pr] = copy.deepcopy(temp[pr])

            else:
                self.error_log.append(f'{leftseg} - {cr.name}')
                print(f'{leftseg} - {cr.name}')

    def find_primers(self):
        """
        Iterate over CRISPRS to find possible primers
        """
        for cr in self.crisprs:
            if not cr.complete:
                cr.start, cr.stop = ntp.find_match(self.cds, cr.seq)
                if cr.start == 'ERROR':
                    continue
                leftseg, rightseg = ntp.find_cds_segs(self.cds, cr, left_buffer=self.end_buffer, right_buffer=self.end_buffer, inside_buffer=self.inside_buffer)

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

    def eval_primer(self, primer, primer_name, cr_name, cr_start, cr_stop):
        """
        evaluates a given primer for problems
        """
        good = False
        primer_info = primer['pr'].output()
        pr_start = int(primer['pr'].start) # new here
        distance, max_AT, max_GC = ntp.find_max_polymer(self.cds, cr_start, cr_stop, pr_start) # new here
        primer['pr'].distance = distance # adds attribute -how far between primer and crispr site
        primer['pr'].maxAT = max_AT # adds information on homopolymeric tract location and length
        primer['pr'].maxGC = max_GC # adds information on homopolymeric tract location and length
        flag = primer['pr'].fail_case
        outstring = f'{self.name},{cr_name},{primer_name},{primer_info},{distance},{max_AT},{max_GC},{flag}' # new here
        if not primer['pr'].fail_case: good = True
        return good, outstring
        #self.all_out.append(outstring)
        #if not primer['pr'].fail_case:
        #    self.out.append(outstring)

    def sort_output(self):
        """
        Writes a csv of the output to the gene folder.
        Also returns the output so the submission can use it
        """

        self.out = []
        self.all_out = []
        for cr in self.crisprs:
            cr_mid = abs(int((cr.start+cr.stop)/2))
            cr_name = cr.name

            for pr in cr.best_fprimers:
                primer = cr.best_fprimers[pr]
                primer_name = 'F' + str(primer['num'])
                good, outstring = self.eval_primer(primer, primer_name, cr_name, cr.start, cr.stop)
                self.all_out.append(outstring)
                if good:
                    self.out.append(outstring)

            for pr in cr.best_rprimers:
                primer = cr.best_rprimers[pr]
                primer_name = 'R' + str(primer['num'])
                good, outstring = self.eval_primer(primer, primer_name, cr_name, cr.start, cr.stop)
                self.all_out.append(outstring)
                if good:
                    self.out.append(outstring)

        return self.out, self.all_out

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
