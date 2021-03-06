#! /usr/bin/env python3

import sys
import gff

human_gff = '/Users/jacob.cooper/resources/genomes/GRCh38_latest_genomic.gff'
human_genome = '/Users/jacob.cooper/resources/genomes/GRCh38_latest_genomic.fasta'

def read_list(file, lower=False):
    l = []
    with open(file, 'r') as f:
        for line in f:
            line = line.rstrip()
            if lower:
                line = line
            l.append(line)
    return l

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

def find_gene_sequences(gene_list, parser):
    """
    Inputs:
        List
        GFF object
    Returns
        Dictionary of [gene] = sequence
    """
    fasta = {}
    found_genes = set(parser.genes['gene'].values)
    for gene in gene_list:
        if gene in found_genes:
            seq = parser.retrieve_gene_sequence(gene.upper())
            seq = str(seq)
            fasta[gene.upper()] = seq
        else:
            print(f'WARNING: Did not find {gene} in data table')
            fasta[gene] = ''
    return fasta


###############################
######## main

def main():
    gene_list = read_list(sys.argv[1], lower=True)
    parser = gff.GFF(human_gff, human_genome)
    parser.build_genes()
    print('Loaded GFF')

    # make a fasta dictionary that has all of the genes
    fasta = find_gene_sequences(gene_list, parser)

    # write a new fasta file
    write_fasta(fasta, 'genome_mRNA.fasta')


if __name__ == "__main__":
    main()