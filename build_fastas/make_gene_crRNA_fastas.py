#! /usr/bin/env python3

import sys
import pandas as pd
import os

# ./make_gene_crRNA_fastas.py genome_mRNA.fasta all_genes_and_crisprs.csv

def read_fasta(file):
    """
    Input: fasta formatted file
    Returns: dictionary of fasta file
    """

    fasta = {}
    with open(file, 'r') as f:
        # set some flags and variables
        header = False
        entry = ''
        head = ''

        # loop over the lines
        for line in f:
            line = line.rstrip()

            ## check if a new header
            if '>' in line:

                # add the recorded entry
                if head and entry:
                    fasta[head] = entry
                # set the header to true
                header = True

            # reset the header
            if header:
                entry = ''
                head = line.split('>')[1]
                header = False
            # otherwize keep adding to the entry
            else:
                entry += line

        # add in the final entry
        if head and entry:
            fasta[head] = entry

    return fasta

def breakformat(s, space=60):
    '''
    Input : str
    Reutrn : str with a new line character added every space characters
    '''
    return '\n'.join(s[i:i+space] for i in range(0, len(s), space))

def write_fasta(d, file):
    with open(file, 'w') as f:
        for key, value in d.items():
            f.write(f'>{key}\n')
            f.write(breakformat(value))
            f.write('\n')

def main():
    # load in the mRNA fasta
    mRNA_FASTA = read_fasta(sys.argv[1])

    # load in the crRNA table
    crRNA_TABLE = pd.read_csv(sys.argv[2])

    working_path = sys.argv[3]

    if not os.path.isdir(f'{working_path}/genes'):
        os.mkdir(f'{working_path}/genes')

    # for each gene
    for gene in mRNA_FASTA:

        if len(mRNA_FASTA[gene]) > 10:

            # look for the crRNAs
            assert gene in crRNA_TABLE['gene'].values, f'{gene} is not in the crRNA table'
            df = crRNA_TABLE[crRNA_TABLE['gene'] == gene]
            temp = df.set_index('rec_id')['sequence'].to_dict()

            # update the names
            towrite = {}
            for seq in temp:
                newname = f'crispr:{seq} gene:{gene} segment:E900'
                towrite[newname] = temp[seq]

            towrite[gene] = mRNA_FASTA[gene]

            # write a new folder and new fasta file with that gene
            if not os.path.isdir(f'{working_path}/genes/{gene}'):
                os.mkdir(f'{working_path}/genes/{gene}')

            write_fasta(towrite, f'{working_path}/genes/{gene}/{gene}.fasta')
        else:
            print(gene)




if __name__ == "__main__":
    main()
