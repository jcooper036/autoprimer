#!/usr/bin/env python3
"""init for the AutoPrimer module"""

# autobot
from AutoPrimer.autobot.autobot import *

# wrappers
from AutoPrimer.wrappers.find_files import find_files
from AutoPrimer.wrappers.Submission import Submission

# main
from AutoPrimer.main.Primer import Primer
from AutoPrimer.main.Crispr import Crispr
from AutoPrimer.main.Gene import Gene
from AutoPrimer.main.parse_fasta import parse_fasta
from AutoPrimer.main.read_fasta import read_fasta
from AutoPrimer.main.read_settings import read_settings
from AutoPrimer.main.find_match import find_match
from AutoPrimer.main.find_cds_segs import find_cds_segs
from AutoPrimer.main.reverse_comp import reverse_comp
from AutoPrimer.main.bash import bash
from AutoPrimer.main.run_primer3 import run_primer3
from AutoPrimer.main.pick_primers import pick_primers