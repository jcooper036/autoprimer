#! /usr/bin/env python3

def read_fasta(file):
    """
    Input: fasta formatted file
    Returns: dictionary of fasta file
    """
    
    fasta = {}
    with open(file, 'r') as f:
        # set some flags and variables
        header = False
        entry = ''
        head = ''
        
        # loop over the lines
        for line in f:
            line = line.rstrip()

            ## check if a new header
            if '>' in line:

                # add the recorded entry
                if head and entry:
                    fasta[head] = entry
                # set the header to true
                header = True

            # reset the header
            if header:
                entry = ''
                head = line.split('>')[1]
                header = False
            # otherwize keep adding to the entry
            else:
                entry += line
    
        # add in the final entry
        if head and entry:
            fasta[head] = entry

    return fasta
