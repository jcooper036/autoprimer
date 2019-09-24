#! /usr/bin/env python3

import AutoPrimer as ntp

def parse_fasta(fasta, gene):
    """
    Input: fasta dictionary
    Returns: CDS sequence, dictionary of the CRISPR targets with info about them.
    """
    for entry in fasta:
        if 'crispr' in entry:
            entrysp = entry.split(' ')
            for info in entrysp:
                if 'crispr' in info:
                    cr = ntp.Crispr(info.split(':')[1])
                if 'segment' in info:
                    segment = info.split(':')[1]
                    pass
            cr.seq = fasta[entry]
            cr.primerTag = f'{gene.name}_{segment}_{cr.name}'
            gene.crisprs.append(cr)
        else:
            gene.cds = fasta[entry]
    return gene