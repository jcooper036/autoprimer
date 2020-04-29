#! /usr/bin/env python3

import AutoPrimer as ntp
import os

def blast_primer(seq):
    """
    Input: primer sequence
    Returns True if primer is unique, false if the primer is not.
    """

    #@ this is a bit janky - ther are other possible options here
    GENOME = '/Users/chris.johnson/Documents/resources/genomes/GCF_000001405.39_GRCh38.p13_genomic.fna' #Chris J
    #GENOME = '/Users/$USER/resources/genomes/GRCh38_latest_genomic.fasta'
    # GENOME = '/Volumes/i_bio/Crispr_F0_Screens/checkprimer/genome/GRCh38_latest_genomic.fasta'
    # GENOME = '/mnt/i_bio/Crispr_F0_Screens/checkprimer/genome/GRCh38_latest_genomic.fasta'
    # GENOME = '/Users/jacob.cooper/resources/genomes/GRCh38_latest_genomic.fasta'

    # if not os.path.exists(GENOME + '.nhr'):
    #     print("UPDATE: making BLAST database for " + GENOME)
    #     command = 'makeblastdb -in ' + GENOME + ' -parse_seqids -dbtype nucl'
    #     ntp.bash(command)
    command = 'blastn -task blastn-short -outfmt "6" -query <(echo ' + seq + ') -db ' + GENOME + ' | sort -k12 -n -r | head -5'

    out = str(ntp.bash(command)).split('\n')

    return out
