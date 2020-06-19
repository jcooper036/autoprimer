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
    df_in.rename(columns = {'Notes':'Synthego Notes'}, inplace = True)
    df_in['Notes'] = ''
    df = df_in[['Label','ICE','R Squared','KO-Score','Guide Sequences','Synthego Notes','Notes']].copy()
    for i in range(len(df)):
        if df.iloc[i,2] == 'None': continue
        if float(df.iloc[i,2]) < 0.8: df.iloc[i,6] = 'R squared < 0.8'
        if float(df.iloc[i,2]) < 0.5: df.iloc[i,6] = 'R squared < 0.5'
    df1 = df.query("Notes != 'R squared < 0.5' and ICE != 'None'").sort_values('Label')
    df2 = df.dropna(subset = ['Synthego Notes'])
    df3 = df[df['Notes'] == 'R squared < 0.5']
    df4 = pd.concat([df2, df3]).sort_values('Label')
    df_blank = pd.DataFrame([['','','','','','',''],['Failed','','','','','','']], columns = ['Label','ICE','R Squared','KO-Score','Guide Sequences','Synthego Notes','Notes'])
    return pd.concat([df1,df_blank,df4])

lead_file = ''

root = Tk()
my_gui = GetInfile(root)
root.mainloop()

df = pd.read_csv(lead_file)
out = SortSynthego(df)
out.to_csv(lead_file[:-4]+'_processed.csv', index = False)
