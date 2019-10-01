#! /usr/bin/env python3

import sys
import AutoPrimer as ntp
import datetime

def submit_folder(folder):
    # data here ends up being a list of the Submission objects
    data = ntp.spawn_processes(folder)

    # gather data from the submissions
    out = []
    for sub in data:
        for ent in sub.out:
            out.append(ent)
    
    # write to a file
    header = 'GENE,CRISPR,PRIMER,SEQUENCE,SIDE,START,TM,GC%'
    now = datetime.datetime.now()
    fname = folder + 'autoprimer_combine_beta_' + now.strftime("%Y-%m-%d.%H-%M") + '.csv'
    ntp.write_csv(fname, out, header=header)
    return data

if __name__ == "__main__":
    submit_folder('/Volumes/i_bio/Crispr_F0_Screens/0-Genes_for_design/Genes_for_autoprimer/')
