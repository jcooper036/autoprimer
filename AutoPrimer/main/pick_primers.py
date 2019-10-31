#! /usr/bin/env python3

import AutoPrimer as ntp

settings = {
    'sequence_id' : 'blank',
    'task' : 'generic',
    'pick_internal_oligo' : '0',
    'explain_flag' : '1',
    'opt_size' : '20',
    'min_size' : '19',
    'max_size' : '25',
    'min_gc' : '45',
    'max_gc' : '65',
    'min_tm' : '60',
    'opt_tm' : '61',
    'max_tm' : '62',
    'gc_clamp' : '1',
    'max_poly_x' : '3',
    'primer_return_number' : '200'
}

def pick_primers(template, side, sett=settings):
    """
    Input: Template sequence, side is left or right
    Returns: dictionary of primer objects, somewhere between 0 and 4
    """

    # settings variables
    sequence_template = template
    if side == 'left':
        pick_left_primer = '1'
        pick_right_primer = '0'
    elif side == 'right':
        pick_left_primer = '0'
        pick_right_primer = '1'        
    
    sequence_id = sett['sequence_id']
    task = sett['task']
    pick_internal_oligo = sett['pick_internal_oligo']
    opt_size = sett['opt_size']
    min_size = sett['min_size']
    max_size = sett['max_size']
    explain_flag = sett['explain_flag']
    min_gc = sett['min_gc']
    max_gc = sett['max_gc']
    min_tm = sett['min_tm']
    opt_tm = sett['opt_tm']
    max_tm = sett['max_tm']
    gc_clamp = sett['gc_clamp']
    poly_x = sett['max_poly_x']
    primer_return_num = sett['primer_return_number']

    # string concat for settings
    primer3_settings = 'SEQUENCE_ID=' + sequence_id + '\n' + 'SEQUENCE_TEMPLATE=' + sequence_template + '\n' + 'PRIMER_TASK=' + task + '\n' + 'PRIMER_PICK_LEFT_PRIMER=' + pick_left_primer + '\n' + 'PRIMER_PICK_INTERNAL_OLIGO=' + pick_internal_oligo + '\n' + 'PRIMER_PICK_RIGHT_PRIMER=' + pick_right_primer + '\n' + 'PRIMER_OPT_SIZE=' + opt_size + '\n' + 'PRIMER_MIN_SIZE=' + min_size + '\n' + 'PRIMER_MAX_SIZE=' + max_size + '\n'  + 'PRIMER_EXPLAIN_FLAG=' + explain_flag + '\n' + 'PRIMER_INTERNAL_MIN_TM=' + min_tm + '\n' + 'PRIMER_INTERNAL_OPT_TM=' + opt_tm + '\n' + 'PRIMER_INTERNAL_MAX_TM=' + max_tm + '\n' + 'PRIMER_INTER_MIN_GC=' + min_gc + '\n' + 'PRIMER_INTERNAL_MAX_GC=' + max_gc + '\n' + 'PRIMER_NUM_RETURN=' + primer_return_num + '\n' + 'PRIMER_GC_CLAMP=' + gc_clamp + '\n' + 'PRIMER_MAX_POLY_X=' + poly_x + '\n' '='

    # test for primer3
    return ntp.run_primer3(primer3_settings)