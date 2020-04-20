#! /usr/bin/env python3

import AutoPrimer as ntp
import copy

def evaluate_primers(primers):
    """
    Input: Dicitonary primers = {'name' : {'pr' : PrimerObject}}
    Returns: Dictionary of the best primers
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

    for pr in primers:
        prm = primers[pr]['pr']
        add = True

        # make sure the primer hasn't already failed
        if prm.fail_case:
            add = False

        # make sure the primers hasn't been tried yet
        if prm.seq in tried_sequences:
            add = False
            prm.fail_case = 'Duplicate sequence'
            # print(prm.fail_case) #@

        # check the anyTH and endTH
        if (float(prm.anyTH) + float(prm.endTH)) > th_max and add:
            add = False
            prm.fail_case = f'anyTH + endTH > {th_max}'

        # start positions need to be appart from each other
        if start_postitions and add:
            if any(abs(prm.start - x) <= start_buffer for x in start_postitions):
                add = False
                prm.fail_case = f'Too close to another primer'

        # check blast
        primers[pr]['blast'] = 'Not run'
        if add and (blastCount < max_blast):
            # print('Blasting {}, try {} / {}'.format(prm.seq, blastCount+1, max_blast)) #@
            primers[pr]['blast'] = ntp.blast_primer(prm.seq)
            if not ntp.check_blast_result(primers[pr]['blast']):
                add = False
                prm.fail_case = f'Too close to another primer'
            blastCount += 1

        # if all is good, add the primer
        primers[pr]['num'] = 'NA'
        if add:
            primers[pr]['num'] = keepCount + 1
            # print(f'Found primer {keepCount + 1} : {prm.seq}\n') #@
            best_primers[pr] = copy.deepcopy(primers[pr])
            start_postitions.append(prm.start)
            keepCount += 1

        # always add the sequnece to the tried sequences
        tried_sequences.append(prm.seq)

        # break conditions, which will cause best primers to be returned
        # if we found enough or tried blast too many times
        if (keepCount == keep_primers) or (blastCount == max_blast):
            # if blastCount == max_blast:
            #     print(f'Too many BLAST attempts, moving on.')
            break

    return best_primers
