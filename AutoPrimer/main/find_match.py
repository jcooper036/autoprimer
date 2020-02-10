#! /usr/bin/env python3

import AutoPrimer as ntp

def find_match(seq1, seq2):
    """
    Input: sequence 1 (primary sequence) and sequence 2 (the search)
    Returns: The start and end position of seq2 in seq1
    """
    seq1 = seq1.lower()
    seq2 = seq2.lower()
    found = False
    len2 = len(seq2)

    # try the sequence forward
    for idx, site in enumerate(seq1):
        if len(seq1) >= idx+len2:
            if seq1[idx:idx+len2] == seq2:
                found = True
                break
    if found:
        return (idx, idx+len2)

    # try the reverse complement
    seq2 = ntp.reverse_comp(seq2).lower()
    for idx, site in enumerate(seq1):
        if len(seq1) >= idx+len2:
            if seq1[idx:idx+len2] == seq2:
                found = True
                break
    if found:
        return (idx, idx+len2)
    
    # return('ERROR', 'ERROR')