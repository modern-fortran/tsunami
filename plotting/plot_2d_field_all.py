#!/usr/bin/env python

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np

dims = 101, 101

for i in range(1, 1001, 10):

    input_file = 'tsunami_h_' + '%4.4i' % i + '.dat'
    print('Plotting ' + input_file)

    field = np.reshape(np.fromfile(input_file, dtype='float32'), dims)

    fig = plt.figure(figsize=(8, 7))
    ax = fig.add_subplot(111, aspect='equal')
    plt.contourf(field[1:-2, 1:-2], np.arange(-0.2, 0.22, 0.02), cmap=cm.bwr)
    plt.title(input_file)
    plt.savefig(input_file[:-2] + '.png')
    plt.close(fig)
