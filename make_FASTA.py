# make_FASTA.py
# takes a list of CRISPR guide RNA rec IDs and runs the sequence of things
# needed to produce the FASTA files for all that
# The list should be a guides.csv file in a folder that the user points to
# can use tkinter to make that super easy

# python make_FASTA.py ../path/to/folder

import cauldron_sdk as api
import pandas as pd
import sys
import csv
import subprocess
import os
from tkinter import *
from tkinter.filedialog import askopenfilename

class GetInfile:
    "GUI that asks the user to identify list of CRISPR guide REC IDs"
    def __init__(self, window):
        row_num = 0
        self.window = window
        self.label_text = 'none'
        self.window.title('make FASTAs')
        window_label = 'Find list of guide RNA REC IDs'
        self.label_Picks = Label(window,text=window_label).grid(row = row_num, columnspan = 2, sticky = W)
        row_num += 1
        self.sourceplate_entry = self.create_browse_widget(row_num)
        row_num += 1
        self.OK_button = Button(window, text = 'OK', command = self.submitFiles)
        self.OK_button.grid(row = row_num, column=0, sticky = E)

    def create_browse_widget(self,r):
        label_text = StringVar()
        label_text.set('none')
        new_button = Button(self.window, text = 'browse for file', command = lambda: label_text.set(askopenfilename()))
        new_label = Label(self.window, textvariable = label_text)
        new_button.grid(row = r, column = 0, sticky = E)
        new_label.grid(row = r, column = 1, sticky = W)
        return new_label

    def submitFiles(self):
        global lead_file
        flags = 0
        x = self.sourceplate_entry.cget('text')
        if x == 'none':
            flags = 1
        lead_file = x
        if flags == 1: print ('Please identify file')
        else:
            self.window.destroy()

def GetGuides(guide_list):
    """
    Query Cauldron with a list of guide REC IDs and return a dataframe with
    gene, rec_id, target_sequence
    """
    GnC = [['gene','rec_id','sequence']]
    for guide in guide_list:
        X = [i for i in api.guides.find(rec_ids = guide)][0]
        GnC.append([X['gene'], X['rec_id'], X['target_sequence']])
    return pd.DataFrame(GnC[1:], columns = GnC[0])

lead_file = ''

get_genes = 'build_fastas/get_gene_mrna.py'
make_fasta = 'build_fastas/make_gene_crRNA_fastas.py'
genes_out = 'build_fastas/temp/gene_list.csv'
df_out = 'build_fastas/temp/temp.csv'
collected_genes = 'build_fastas/temp/genome_mRNA.fasta'

# Have the user browse for a list of guide RNA rec ids
root = Tk()
my_gui = GetInfile(root)
root.mainloop()

working_path = os.path.dirname(lead_file)

# read that list in
guides = []
with open(lead_file,'r') as infile:
    readfile = csv.reader(infile, delimiter = ',')
    for row in readfile: guides.append(row)

# Get a dataframe with guide rna info
print('Fetching guide information from Cauldron')
GnC_df = GetGuides(guides)

# Run the subprocesses
print ('collecting genes from genome')
subprocess.call(['python',get_genes,df_out])

print ('building individual FASTA files')
subprocess.call(['python',make_fasta,collected_genes,df_out, working_path])
