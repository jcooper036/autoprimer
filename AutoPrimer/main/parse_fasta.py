#! /usr/bin/env python3

import AutoPrimer as ntp

def parse_fasta(fasta, gene):
    """
    Input: fasta dictionary
    Returns: CDS sequence, dictionary of the CRISPR targets with info about them.
    """
    for entry in fasta:
        if 'crispr' in entry:
            entry = entry.split(' ')
            for info in entry:
                if 'crispr' in info:
                    cr = ntp.Crispr(info.split(':')[1])
                if 'segment' in info:
                    # exon = info.split(':')[1]
                    pass
            print(fasta)
            cr.seq = fasta[entry]
            crisprs[name] = {
                'seq' : crseq,
                'name' : name,
                # 'segment' : exon,
                # 'gene' : gene,
                'primers' : {},
            }
            cr.primerTag = '{}'.format(crisprs[name]['name']) #@
            # cr.primerTag = '{}_{}_{}'.format(crisprs[name]['gene'], crisprs[name]['segment'], crisprs[name]['name'])
            gene.crisprs.append(cr)
        else:
            gene.cds = fasta[entry]
    return gene