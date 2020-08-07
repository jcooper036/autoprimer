#! /usr/bin/env python3

import AutoPrimer as ntp
import os

def run_primer3(primer3_settings):
    """Appends the necessary quotes to the primer3 settings, runs primer3, returns the shell output"""
    primer3_paths = ['~/primer3/src/primer3_core',
        '~/miniconda3/bin/primer3_core',
        '/Users/matlab/miniconda3/pkgs/primer3-2.5.0-pl526h6de7cb9_0/bin/primer3_core',
        '/Users/chris.johnson/miniconda3/envs/autoprimer/bin/primer3_core']
    primer3_path = ''
    for path in primer3_paths:
        if os.path.isfile(path): primer3_path = path
    if primer3_path == '': print ('primer3_core not correctly specified in run_primer3.py')
    primer3_settings = '"' + primer3_settings + '"'
    out = ntp.bash('printf ' + primer3_settings + ' | ' + primer3_path) # Chris J
    return out
