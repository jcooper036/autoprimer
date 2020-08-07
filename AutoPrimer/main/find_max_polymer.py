
def find_max_polymer(cds, cr_start, cr_stop, stop):
    """
    Input: CDS seqeunce, Crispr object, primer object
    Returns: length of homopolymeric tracts between the primer and crispr
            and distance between tract and crispr
    """
    if stop > cr_stop:
        base = cr_stop
        seq_in = cds[cr_stop:stop]
        orientation = 'R'
    elif cr_start > stop:
        base = stop
        seq_in = cds[stop:cr_start]
        orientation = 'L'
    else:
        return 0,0
    distance = len(seq_in)
    sequence = seq_in.upper()

    AT_list = ''
    GC_list = ''
    run = 1
    for i in range(1,len(sequence)):
        letter = sequence[i]
        retro = sequence[i-1] # retro is the letter preceding the current letter
        if retro == letter:
            run += 1
        else: # if we're at the end of a run, decide to record or ignore
            if orientation == 'R': dist = i-run
            if orientation == 'L': dist = len(sequence) - i
            if run > 7 and retro in 'AT':
                AT_list += str(dist)+' '+str(run)+';'
            if run > 7 and retro in 'GC':
                GC_list += str(dist)+' '+str(run)+';'
            run = 1

    return distance, AT_list.strip(';'), GC_list.strip(';')


def eval_homopolymers(runs, target, limit):
    """
    evaluates lists of homopolymers and returns locations closest to target location
    """
    seq_start = 0
    seq_dist = 10000
    run_list = runs.split(';')
    for run in run_list:
        if run == '': continue
        dist, length = run.split()
        if (int(length) > limit) and (int(dist) < seq_dist):
            seq_dist = int(dist)

    return int(seq_dist)
