#! /usr/bin/env python3

def write_csv(file, data, header=False):
    """
    Inputs : File name, data = ['c1a,c2a,c3a', 'c1b,c2b,c3b'], header='c1,c2,c3'
    Writes to a file
    """
    with open(file, 'w') as f:
        if header:
            f.write(f'{header}\n')
        for line in data:
            f.write(f'{line}\n')