

def find_cds_segs(cds, cr, left_buffer=510, right_buffer = 510, l_in_buffer=150, r_in_buffer =150):
    """
    Input: CDS seqeunce, Crispr object
    Returns: Left and right CDS segment for primer searching
    """

    if (cr.start - l_in_buffer) < 0:
        return 'ERROR: CRISPR too close to CDS start', None

    leftseg = cds[cr.start-left_buffer : cr.start-l_in_buffer]

    if (cr.stop + r_in_buffer) > len(cds):
        return 'ERROR: CRISPR too close to CDS end', None

    rightseg = cds[cr.stop+r_in_buffer : cr.stop+right_buffer]

    return leftseg, rightseg
