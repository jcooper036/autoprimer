#! /usr/bin/env python3

import AutoPrimer as ntp
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

    def find_primers(self, method = 'first'):
        """
        Iterate over CRISPRS to find possible primers
        """
        settings = {
            'sequence_id' : 'blank',
            'task' : 'generic',
            'pick_internal_oligo' : '0',
            'explain_flag' : '1',
            'opt_size' : '20',
            'min_size' : '19',
            'max_size' : '25',
            'min_gc' : '45',
            'max_gc' : '65',
            'min_tm' : '60',
            'opt_tm' : '61',
            'max_tm' : '62',
            'gc_clamp' : '1',
            'max_poly_x' : '4',
            'primer_return_number' : '200'
        }

        for cr in self.crisprs:
            # conditions to skip iteration
            if (method == 'first') and cr.complete: continue
            if (method == 'second') and not cr.new_fasta_L and not cr.new_fasta_L: continue

            # method-specific variables and tasks
            pick_left = True
            pick_right = True

            if method == 'first':
                cr.start, cr.stop = ntp.find_match(self.cds, cr.seq)
                fprimers = cr.fprimers
                rprimers = cr.rprimers

            elif method == 'second':
                if cr.out_buffer_L == self.end_buffer: pick_left = False
                if cr.out_buffer_R == self.end_buffer: pick_right = False
                fprimers = cr.fprimers2
                rprimers = cr.rprimers2

            elif method == 'third':
                settings = {
                    'sequence_id' : 'blank',
                    'task' : 'generic',
                    'pick_internal_oligo' : '0',
                    'explain_flag' : '1',
                    'opt_size' : '20',
                    'min_size' : '18',
                    'max_size' : '30',
                    'min_gc' : '40',
                    'max_gc' : '70',
                    'min_tm' : '55',
                    'opt_tm' : '61',
                    'max_tm' : '62',
                    'gc_clamp' : '0',
                    'max_poly_x' : '5',
                    'primer_return_number' : '200'
                }
                if not cr.loosened_f: pick_left = False
                if not cr.loosened_r: pick_right = False
                fprimers = cr.fprimers3
                rprimers = cr.rprimers3
                #print ('find primers third')
            else: print ('find method not recognized')

            if cr.start == 'ERROR': continue
            leftseg, rightseg = ntp.find_cds_segs(self.cds, cr, left_buffer=cr.out_buffer_L, right_buffer=cr.out_buffer_R, l_in_buffer=cr.in_buffer_L, r_in_buffer=cr.in_buffer_R)

            # qualify segments individually
            if leftseg and pick_left:
                left_raw = ntp.pick_primers(leftseg, 'left', sett = settings)
                temp = ntp.raw_to_primers(left_raw, 'left')
                for pr in temp:
                    fprimers[pr] = copy.deepcopy(temp[pr])

            if rightseg and pick_right:
                right_raw = ntp.pick_primers(rightseg, 'right', sett = settings)
                temp = ntp.raw_to_primers(right_raw, 'right')
                for pr in temp:
                    rprimers[pr] = copy.deepcopy(temp[pr])
            #else:
            #    self.error_log.append(f'{leftseg} - {cr.name}')
            #    print(f'{leftseg} - {cr.name}')

    def sort_primers(self, method = 'first'):
        """
        Iterate over CRISPRS to sort the primers
        """
        if method == 'first':
            for cr in self.crisprs:
                cr.sort_primers(method = method)
        if method == 'second':
            for cr in self.crisprs:
                if cr.new_fasta_L or cr.new_fasta_R:
                    cr.sort_primers(method = method)
        if method == 'third':
            for cr in self.crisprs:
                if cr.loosened_f or cr.loosened_r:
                    cr.sort_primers(method = method)

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
        # Awkward place to add this as a fail case, neaten up in later versions
        all_runs = max_AT+';'+max_GC #new
        if ntp.eval_homopolymers(all_runs, cr_start, 11) < 10000: #new
            primer['pr'].fail_case = 'Poly-N > 11' #new
        flag = primer['pr'].fail_case
        outstring = f'{self.name},{cr_name},{primer_name},{primer_info},{distance},{max_AT},{max_GC},{flag}' # new here
        if not primer['pr'].fail_case: good = True
        return good, outstring

    def sort_output(self, method = 'first'):
        """
        Writes a csv of the output to the gene folder.
        Also returns the output so the submission can use it
        """
        keep_primers = 1
        if method == 'first':
            self.out = []
            self.all_out = []
        elif method == 'second':
            self.out2 = []
            self.all_out2 = []
        elif method == 'third':
            self.out3 = []
            self.all_out3 = []
        else: print ('sort method not recognized')

        for cr in self.crisprs:
            f_found = 0
            r_found = 0
            if method == 'first':
                out = self.out  #
                all_out = self.all_out #
                fprimers = cr.best_fprimers #
                rprimers = cr.best_rprimers #

            elif method == 'second':
                out = self.out2  #
                all_out = self.all_out2 #
                fprimers = cr.best_fprimers2 #
                rprimers = cr.best_rprimers2 #

            elif method == 'third':
                out = self.out3  #
                all_out = self.all_out3 #
                fprimers = cr.best_fprimers3 #
                rprimers = cr.best_rprimers3 #
            else: print ('sort method not recognized')

            if (method == 'second') and not cr.new_fasta_L and not cr.new_fasta_L: continue
            if (method == 'third') and not cr.loosened_f and not cr.loosened_r: continue

            cr_mid = abs(int((cr.start+cr.stop)/2))
            cr_name = cr.name

            for pr in fprimers:
                primer = fprimers[pr]
                primer_name = 'F' + str(primer['num'])
                good, outstring = self.eval_primer(primer, primer_name, cr_name, cr.start, cr.stop)
                all_out.append(outstring)
                if good and (f_found < keep_primers):
                    out.append(outstring)
                    f_found += 1

            for pr in rprimers:
                primer = rprimers[pr]
                primer_name = 'R' + str(primer['num'])
                good, outstring = self.eval_primer(primer, primer_name, cr_name, cr.start, cr.stop)
                all_out.append(outstring)
                if good and (r_found < keep_primers):
                    out.append(outstring)
                    r_found += 1

        return out, all_out

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
