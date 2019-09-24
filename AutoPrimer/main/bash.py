#! /usr/bin/env python3

import subprocess

def bash(command):
    """Makes running shell commands less verbose"""
    return subprocess.check_output(command, shell=True, executable='/bin/bash', universal_newlines=True)
    # return subprocess.check_output(command, shell=True, executable='/bin/bash', universal_newlines=True).decode('utf-8').strip()