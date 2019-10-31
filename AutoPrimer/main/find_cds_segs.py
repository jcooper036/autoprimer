

def find_cds_segs(cds, cr, end_buffer=510, inside_buffer=150):
    """
    Input: CDS seqeunce, Crispr object
    Returns: Left and right CDS segment for primer searching
    """
    
    if (cr.start - end_buffer) < 0:
        return 'ERROR: CRISPR too close to CDS start', None
        
    leftseg = cds[cr.start-end_buffer : cr.start-inside_buffer]
     
    if (cr.stop + end_buffer) > len(cds):
        return 'ERROR: CRISPR too close to CDS end', None
        
    rightseg = cds[cr.stop+inside_buffer : cr.stop+end_buffer]
    
    return leftseg, rightseg