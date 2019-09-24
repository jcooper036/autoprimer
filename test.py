#! /usr/bin/env python3

"""
For testing AutoPrimer functions and workflow.
"""

import AutoPrimer as ntp
import multiprocessing
import sys
import time

settings = {
    
}

folder = sys.argv[1]

# start the multiprocessing
jobs = []

files = ntp.find_files(input_loc=folder)
proces = int(len(files))
p = multiprocessing.Pool(processes=proces)
data = p.map(ntp.Submission, [file for file in files])
p.close()
print(data)

# for file in ntp.find_files(input_loc=folder):
#     d = multiprocessing.Process(target=ntp.Submission, args=(file,))
#     jobs.append(d)
#     d.start()