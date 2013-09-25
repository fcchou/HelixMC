#!/usr/bin/env python

from sys import argv
import os
import subprocess
import numpy as np

test_path = os.path.join( os.path.abspath(os.path.dirname(__file__)),
    'cmdline_tests')

result_folder = 'test_results'
if not os.path.isdir(result_folder):
    os.mkdir(result_folder)
os.chdir(result_folder)

if len(argv) > 1:
    tests = argv[1:]
else:
    tests = []
    for test in os.listdir(test_path):
        abspath = os.path.join(test_path, test)
        if test[-3:] != '.py' and os.path.isfile(abspath):
            tests.append(test)

for test in tests:
    cmdline = open(os.path.join(test_path, test)).read()
    cmdline = cmdline.replace('\\', '')
    if not os.path.isdir(test):
        os.mkdir(test)
    os.chdir(test)
    out_str = subprocess.check_output(cmdline.split())
    with open('log', 'w') as out:
        for line in out_str.splitlines():
            if 'time' not in line:
                out.write(line)

    data = np.load('MC_data.npz')
    for key in data:
        if key != 'frame_terminal':
            np.savetxt(key, data[key], fmt='%.6f')
    os.chdir('../')

