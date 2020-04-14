#! /usr/bin/env python3

"""
Class for parsing gff objects, returning gene info, etc
Paired with the file that it annotates
GFF documentation :
    https://en.wikipedia.org/wiki/General_feature_format
    https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md

Column 1: "seqid"
    [a-zA-Z0-9.:^*$@!+_?-|]. IDs may not contain unescaped whitespace and must not begin with an unescaped ">".
Column 2: "source"
Column 3: "type"
Columns 4 & 5: "start" and "end"
Column 6: "score"
Column 7: "strand"
Column 8: "phase"
For features of type "CDS", the phase indicates where the next codon begins relative to the 5' end (where the 5' end of the CDS is relative to the strand of the CDS feature) of the current CDS feature. For clarification the 5' end for CDS features on the plus strand is the feature's start and and the 5' end for CDS features on the minus strand is the feature's end. The phase is one of the integers 0, 1, or 2, indicating the number of bases forward from the start of the current CDS feature the next codon begins. A phase of "0" indicates that a codon begins on the first nucleotide of the CDS feature (i.e. 0 bases forward), a phase of "1" indicates that the codon begins at the second nucleotide of this CDS feature and a phase of "2" indicates that the codon begins at the third nucleotide of this region. Note that ‘Phase’ in the context of a GFF3 CDS feature should not be confused with the similar concept of frame that is also a common concept in bioinformatics. Frame is generally calculated as a value for a given base relative to the start of the complete open reading frame (ORF) or the codon (e.g. modulo 3) while CDS phase describes the start of the next codon relative to a given CDS feature.

The phase is REQUIRED for all CDS features.

Column 9: "attributes"
A list of feature attributes in the format tag=value. Multiple tag=value pairs are separated by semicolons. URL escaping rules are used for tags or values containing the following characters: ",=;". Spaces are allowed in this field, but tabs must be replaced with the %09 URL escape. Attribute values do not need to be and should not be quoted. The quotes should be included as part of the value by parsers and not stripped.

"""

import pandas as pd
from pyfaidx import Faidx


# d = gff.GFF('/Users/jacob.cooper/resources/genomes/GRCh38_latest_genomic.gff', '/Users/jacob.cooper/resources/genomes/GRCh38_latest_genomic.fasta')

def gene_name_parse(s):
    """
    Input:
        s - str, format ID=gene-ND6;.... (; sep)
        ID=gene-TRNT;Dbxref=GeneID:4576,HGNC:HGNC:7499,MIM:590090;Name=TRNT;gbkey=Gene;gene=TRNT;gene_biotype=tRNA;gene_synonym=MTTT'
    Returns:
        str, gene name
    """

    g = s.split(';')
    for x in g:
        if 'gene=' in x:
            return x.split('gene=')[1]
    return None

def gene_id_parse(s):
    '''
    Input :
        s - str, any of gene, exon, intron feature
    '''
    gid = s.split(';')
    for x in gid:
        if 'Dbxref=GeneID:' in x:
            return x.split('Dbxref=GeneID:')[1].split(',')[0]
    return None

def biotype_parse(s):
    '''
    Input :
        s - str, any of gene, exon, intron feature
    '''
    gid = s.split(';')
    for x in gid:
        if 'gene_biotype=' in x:
            return x.split('gene_biotype=')[1]
    return None

class GFF(object):

    def __init__(self, file, target):
        self.file = file # this is the .gff file
        self.target = target # this is the file that it is an annotation of
        self.read_gff_file()
        self.faidx = Faidx(target)
    
    def read_gff_file(self):
        assert self.file
        l = []
        with open(self.file, 'r') as f:
            for line in f:
                if line[0] != '#':
                    line = line.rstrip()
                    l.append(line.split('\t'))
        self.df = pd.DataFrame(l, columns = 'seqid source type start end score strand phase attributes'.split())
    
    def build_genes(self):
        """
        For the entries that have type == gene

        """
        self.genes = self.df[self.df['type'].isin(['gene', 'exon', 'mRNA'])]
        self.genes['gene'] = self.genes['attributes'].apply(gene_name_parse)
        self.genes['geneID'] = self.genes['attributes'].apply(gene_id_parse)
        self.genes['biotype'] = self.genes['attributes'].apply(biotype_parse)
    
    def find_genes(self, genelist):
        """
        Input: list of genes to query
        Returns : self.genes info for those genes
        """
        
        return self.genes[self.genes['gene'].isin(genelist)]
    
    def retrieve_gene_sequence(self, gene):
        start = int(self.genes[(self.genes['gene'] == gene)&(self.genes.type=='gene')]['start'].values[0])
        end = int(self.genes[(self.genes['gene'] == gene)&(self.genes.type=='gene')]['end'].values[0])
        source = str(self.genes[(self.genes['gene'] == gene)&(self.genes.type=='gene')]['seqid'].values[0])
        strand = str(self.genes[(self.genes['gene'] == gene)&(self.genes.type=='gene')]['strand'].values[0])
        
        # add buffers for sequencing primers
        start -= 1000
        end += 1000

        seq = str(self.faidx.fetch(source, start, end))
        if strand == '-':
            seq = self.reverse_comp(seq)
        return seq
    
    def retrieve_gene_structures(self, gene):
        """
        Input
            - gene : str
        Returns
            - fasta_dict : 
                keys are gene;transcript;segment_type;start;stop;strand
                values are sequences
        """
        
    
    def reverse_comp(self, seq):
        """
        Returns the reverse complmiment of a DNA sequence
        """
        d = {
            'A' : 'T',
            'G' : 'C',
            'C' : 'G',
            'T' : 'A',
            'N' : 'N'
        }
        new_seq = ''
        for base in seq:
            new_seq += d[base.upper()]
        return new_seq[::-1]
    

        
        
    
