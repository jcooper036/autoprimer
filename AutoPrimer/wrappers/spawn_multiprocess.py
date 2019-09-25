#! /usr/bin/env python3

import AutoPrimer as ntp
import multiprocessing
import sys

def spawn_processes(folder):
    """
    Uses multiprocess to spawn jobs
    """
    # start the multiprocessing
    jobs = []
    files = ntp.find_files(input_loc=folder)
    print(files)
    proces = int(len(files))
    p = multiprocessing.Pool(processes=proces)
    data = p.map(ntp.Submission, [file for file in files])
    p.close()
    return data