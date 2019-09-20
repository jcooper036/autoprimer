#!/usr/bin/env python3
"""init for the AutoPrimer module"""

# autobot
from AutoPrimer.autobot.autobot import *

# wrappers
from AutoPrimer.wrappers.find_files import find_files
from AutoPrimer.wrappers.submission import submission

# main
from AutoPrimer.main.Primer import Primer
from AutoPrimer.main.Crispr import Crispr
from AutoPrimer.main.Gene import Gene
from AutoPrimer.main.parse_fasta import parse_fasta
from AutoPrimer.main.read_fasta import read_fasta
from AutoPrimer.main.read_settings import read_settings