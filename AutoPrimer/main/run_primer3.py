#! /usr/bin/env python3

import AutoPrimer as ntp

def run_primer3(primer3_settings):
    """Appends the necessary quotes to the primer3 settings, runs primer3, returns the shell output"""
    primer3_settings = '"' + primer3_settings + '"'
    out = ntp.bash('printf ' + primer3_settings + ' | ~/primer3/src/primer3_core')
    return out