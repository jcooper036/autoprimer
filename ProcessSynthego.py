# ProcessSynthego.py
# 6/3/2020
# Chris Johnson
#
# Script to execute some processing of Synthego output files
# User browses for a file, which is then processed and autosaved
# to use, type:
# python ProcessSynthego.py

import pandas as pd
import os
import cauldron_sdk as api
from tkinter import *
from tkinter.filedialog import askopenfilename

class GetInfile:
    "GUI that asks the user to browse for a file"
    def __init__(self, window):
        row_num = 0
        self.window = window
        self.label_text = 'none'
        self.window.title('Process Synthego File')
        window_label = 'Identify file to be processed'
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

def SortSynthego(df_in):
    df = df_in.copy()
    out_columns = ['cell_type','gene_name','rec_id','exp_name','sample_name','knockout_rate','mutation_rate','r_squared','approx_seed_date','Guide Sequences','Synthego Notes','Notes']
    col_rename = {'Label':'sample_name','KO-Score':'knockout_rate','ICE':'mutation_rate','R Squared':'r_squared','Notes':'Synthego Notes'}
    df.rename(columns = col_rename, inplace = True)
    df['Notes'] = ''
    df['cell_type'] = ''
    df = AddGeneRec(df)
    df['exp_name'] = ''
    df['approx_seed_date'] = ''
    df = df[out_columns]
    r_idx = df.columns.get_loc('r_squared')
    n_idx = df.columns.get_loc('Notes')
    for i in range(len(df)):
        if df.iloc[i,7] == 'None': continue
        if float(df.iloc[i,r_idx]) < 0.8: df.iloc[i,n_idx] = 'R squared < 0.8'
        if float(df.iloc[i,r_idx]) < 0.5: df.iloc[i,n_idx] = 'R squared < 0.5'
    df1 = df.query("Notes != 'R squared < 0.5' and mutation_rate != 'None'").sort_values('sample_name')
    df2 = df.dropna(subset = ['Synthego Notes'])
    df3 = df[df['Notes'] == 'R squared < 0.5']
    df4 = pd.concat([df2, df3]).sort_values('sample_name')
    df_blank = pd.DataFrame([['','','','','','','','','','','',''],
                             ['Failed','','','','','','','','','','','']], columns = out_columns)
    return pd.concat([df1,df_blank,df4])

def AddGeneRec(df_in):
    # Adds gene and rec_id for each element
    df = df_in.copy()
    df['gene_name'] = ''
    df['rec_id'] = ''
    # build dicts of sequences in df
    sequences = df['Guide Sequences'].unique().tolist()
    REC_IDS = {}
    GENES = {}
    for sequence in sequences:
        try:
            x = [i for i in api.guides.find(target_sequences = sequence)][0]
            REC_IDS[sequence] = x['rec_id']
            GENES[sequence] = x['gene']
        except:
            REC_IDS[sequence] = 'unknown'
            GENES[sequence] = 'unknown'
    # assign values to df
    g_idx = df.columns.get_loc('Guide Sequences')
    for i in range(len(df)):
        sequence = df.iloc[i, g_idx]
        gene = GENES[sequence]
        rec_id = REC_IDS[sequence]
        df.iloc[i,-1] = rec_id
        df.iloc[i,-2] = gene
    return df

lead_file = ''

root = Tk()
my_gui = GetInfile(root)
root.mainloop()

df = pd.read_csv(lead_file)
out = SortSynthego(df)
out.to_csv(lead_file[:-4]+'_processed.csv', index = False)
