
def find_max_polymer(cds, start, stop):
    """
    Input: CDS seqeunce, Crispr object, primer object
    Returns: length of longest polyhomomeric tract between the primer and crispr
    """

    distance = abs(start-stop)
    if stop > start: seq_in = cds[start:stop]
    elif start > stop: seq_in = cds[stop:start]
    elif start == stop: return 0,0


    sequence = seq_in.upper()
    AT_list = []
    GC_list = []
    run = 1
    for i in range(len(sequence)):
        letter = sequence[i]
        if sequence[i-1] == letter:
            run += 1
        else:
            run = 1
        if letter == 'A' or letter == 'T':AT_list.append(run)
        if letter == 'G' or letter == 'C':GC_list.append(run)

    return distance, max(AT_list), max(GC_list)
