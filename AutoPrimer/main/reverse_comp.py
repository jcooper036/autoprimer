#! /usr/bin/env python3

def reverse_comp(sequence):
    """Reverse complement a DNA string"""
    rc_seq = ''
    rcDict = {
        'A' : 'T',
        'T' : 'A',
        'C' : 'G',
        'G' : 'C'
    }
    for let in sequence:
        rc_seq += rcDict[let.upper()]
    rc_seq = rc_seq[::-1]
    return rc_seq