#!/usr/bin/env python

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('input_file', help='binary file output by tsunami simulator')
args = parser.parse_args()

input_file = args.input_file

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

dims = 100, 100
field = np.reshape(np.fromfile(input_file, dtype='float32'), dims)

fig = plt.figure(figsize=(8,7))
ax = fig.add_subplot(111, aspect='equal')
plt.contourf(field[1:-2, 1:-2], np.arange(-0.2, 0.22, 0.02))
plt.title(input_file)
plt.savefig(input_file[:-2] + '.png')
plt.close(fig)
