#! /usr/bin/env python3

import AutoPrimer as ntp

def raw_to_primers(raw, side):
    """
    Input : Raw output from primer3
    Returns : Dictionary of primer objects
        primers = {0 : {'pr' : PrimerObject}}
        Done this way so other metadata can be stored later about primer '0'
    """
    primers = {}
    raw = raw.split('\n')
    count = 0
    reset = True
    
    for line in raw:
        if reset:
            pen = None
            seq = None
            start = None
            length = None
            gc = None
            tm = None
            anyTH = None
            endTH = None
            hairpin = None
            stab = None
            reset = False

        if "PENALTY=" in line:
            pen = line.split('=')[1]
        if "SEQUENCE=" in line:
            seq = line.split('=')[1]
        if "," in line and 'EXPLAIN' not in line:
            # print(line)
            start = int(line.split('=')[1].split(',')[0])
            length = line.split('=')[1].split(',')[1]
        if "GC_PERCENT" in line:
            gc = line.split('=')[1]
        if "TM" in line:
            tm = line.split('=')[1]
        if "SELF_ANY" in line:
            anyTH = line.split('=')[1]
        if "SELF_END" in line:
            endTH = line.split('=')[1]
        if "HAIRPIN" in line:
            hairpin = line.split('=')[1]
        if "END_STABILITY" in line:
            stab = line.split('=')[1]

        # only save a primer if it has all the right elements
        if all([pen, seq, start, length, gc, tm, anyTH, endTH, hairpin, stab]):
            reset = True
            primers[count] = {
                'pr' : ntp.Primer(side, pen, seq, start, length, tm, gc, anyTH, endTH, hairpin, stab),
                'side' : side
            }
            count += 1

    return primers