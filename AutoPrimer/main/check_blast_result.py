#! /usr/bin/env python3

def check_blast_result(out):
    """
    Input : -6 formated blast hit, split on new lines
    Returns : Bool True if only 1 hit above the cutoff
    cutoff of 36.2 will include things that are 100% matching over 90% of sequence
    """
    cutoff = 36.2

    examine = []
    for output in out:
        output = output.split('\t')
        if output[-1]:
            if float(output[-1]) > cutoff:
                examine.append(output)
                
    return len(examine)
    #if len(examine) == 1:
    #    return True
    #else:
    #    return False
