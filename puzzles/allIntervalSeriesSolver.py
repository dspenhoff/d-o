#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
from subprocess import Popen, PIPE


def solveIt(n_string):



    # Runs the command: java QueensSolver n
    process = Popen(['java', 'AllIntervalSeriesSolver', str(n_string)], stdout=PIPE)
    # print process.stdout.read()
    (stdout, stderr) = process.communicate()


    return stdout.strip()


import sys

if __name__ == "__main__":
    if len(sys.argv) > 1:
        try:
            n = int(sys.argv[1].strip())
        except:
            print sys.argv[1].strip(), 'is not an integer'
        print 'Solving Size:', n
        print(solveIt(sys.argv[1].strip()))

    else:
        print('This test requires an instance size.  Please select the size of problem to solve. (i.e. python queensSolver.py 8)')
