#! /usr/bin/env python3
"""
gff output line
26973  NC_000001.11  BestRefSeq%2CGnomon  gene  11980181  12013515     .      +     .  ID=gene-MFN2;Dbxref=GeneID:9927,HGNC:HGNC:1687...  MFN2   9927  protein_coding
"""

from pyfaidx import Faidx


file = '/Users/jacob.cooper/resources/genomes/GRCh38_latest_genomic.fasta'
chromosome = 'NC_000001.11'
start = 11980181
end = 12013515

fa = Faidx(file)
print(fa.fetch(chromosome, start, end))

