#! /usr/bin/env python3

import AutoPrimer as ntp

def pick_primers(template, side, lowTm=60, highTm=65):
    """
    Input: Template sequence, side is left or right
    Returns: dictionary of primer objects, somewhere between 0 and 4
    """

    # settings variables
    sequence_id = 'blank'
    sequence_template = template
    # print(template) #@
    task = 'generic'
    if side == 'left':
        pick_left_primer = '1'
        pick_right_primer = '0'
    elif side == 'right':
        pick_left_primer = '0'
        pick_right_primer = '1'        
    
    pick_internal_oligo = '0'
    opt_size = '20'
    min_size = '19'
    max_size = '25'
    product_size_range = '75-150'
    explain_flag = '1'
    min_gc = '35'
    max_gc = '75'
    min_tm = '59'
    opt_tm = '63'
    max_tm = '67'
    primer_return_num = '200'

    # + 'PRIMER_PRODUCT_SIZE_RANGE=' + product_size_range + '\n'

    # string concat for settings
    primer3_settings = 'SEQUENCE_ID=' + sequence_id + '\n' + 'SEQUENCE_TEMPLATE=' + sequence_template + '\n' + 'PRIMER_TASK=' + task + '\n' + 'PRIMER_PICK_LEFT_PRIMER=' + pick_left_primer + '\n' + 'PRIMER_PICK_INTERNAL_OLIGO=' + pick_internal_oligo + '\n' + 'PRIMER_PICK_RIGHT_PRIMER=' + pick_right_primer + '\n' + 'PRIMER_OPT_SIZE=' + opt_size + '\n' + 'PRIMER_MIN_SIZE=' + min_size + '\n' + 'PRIMER_MAX_SIZE=' + max_size + '\n'  + 'PRIMER_EXPLAIN_FLAG=' + explain_flag + '\n' + 'PRIMER_INTERNAL_MIN_TM=' + min_tm + '\n' + 'PRIMER_INTERNAL_OPT_TM=' + opt_tm + '\n' + 'PRIMER_INTERNAL_MAX_TM=' + max_tm + '\n' + 'PRIMER_INTER_MIN_GC=' + min_gc + '\n' + 'PRIMER_INTERNAL_MAX_GC=' + max_gc + '\n' + 'PRIMER_NUM_RETURN=' + primer_return_num + '\n='

    # test for primer3
    return ntp.run_primer3(primer3_settings)