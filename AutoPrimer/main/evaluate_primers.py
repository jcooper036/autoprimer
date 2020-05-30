#! /usr/bin/env python3

import AutoPrimer as ntp
import copy

def evaluate_primers(primers):
    """
    Input: Dicitonary primers = {'name' : {'pr' : PrimerObject}}
    Returns: Dictionary of all primers until it has found N primers that meet criteria
    """
    # settings
    th_max = 8
    start_buffer = 25
    max_blast = 8
    keep_primers = 1

    best_primers = {}
    tried_sequences = []
    start_postitions = []
    blastCount = 0
    keepCount = 0
    primerCount = 1

    for pr in primers:
        prm = primers[pr]['pr']
        add = True

        primers[pr]['num'] = primerCount
        primerCount += 1

        # make sure the primer hasn't already failed
        if prm.fail_case:
            add = False

        # make sure the primers hasn't been tried yet
        if prm.seq in tried_sequences:
            add = False
            prm.fail_case = 'Duplicate sequence'

        # check the anyTH and endTH
        if (float(prm.anyTH) + float(prm.endTH)) > th_max and add:
            add = False
            prm.fail_case = f'anyTH + endTH > {th_max}'

        # might remove otherwise good candidates before the BLAST thing
        # start positions need to be appart from each other
        #if start_postitions and add:
        #    if any(abs(prm.start - x) <= start_buffer for x in start_postitions):
        #        add = False
        #        prm.fail_case = f'Too close to another primer'

        # check blast
        primers[pr]['blast'] = 'Not run'
        if add and (blastCount < max_blast):
            # print('Blasting {}, try {} / {}'.format(prm.seq, blastCount+1, max_blast)) #@
            primers[pr]['blast'] = ntp.blast_primer(prm.seq)
            hits = ntp.check_blast_result(primers[pr]['blast'])
            if hits > 1:
                add = False
                prm.fail_case = f'Aligns to {hits} locations'
            blastCount += 1

        # if all is good, add the primer
        if add: keepCount += 1

        best_primers[pr] = copy.deepcopy(primers[pr])
        start_postitions.append(prm.start)

        # always add the sequnece to the tried sequences
        tried_sequences.append(prm.seq)

        # break conditions, which will cause best primers to be returned
        # if we found enough or tried blast too many times
        #if (keepCount == keep_primers) or (blastCount == max_blast):
        #    break

    return best_primers
