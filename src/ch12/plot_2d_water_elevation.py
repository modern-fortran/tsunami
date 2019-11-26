#!/usr/bin/env python

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np

plt.rcParams.update({'font.size': 16})

dt = 0.02
dims = 201, 201
x = range(-100, 101)

for i in range(0,1100, 100):

    input_file = 'tsunami_h_' + '%4.4i' % i + '.dat'
    print('Plotting ' + input_file)

    field = np.reshape(np.fromfile(input_file, dtype='float32'), dims)
    ticks = np.arange(-0.1, 0.11, 0.01)

    field[field > 0.0999] = 0.0999
    field[field <-0.0999] =-0.0999

    fig = plt.figure(figsize=(8, 7))
    ax = fig.add_subplot(111, aspect='equal')
    cnt = plt.contourf(x, x, field, ticks, cmap=cm.Spectral)
    for c in cnt.collections:
        c.set_edgecolor('face')
    plt.colorbar(shrink=0.8)
    plt.xlabel('Distance [m]')
    plt.ylabel('Distance [m]')
    plt.title('Water height @ time = ' + '%3.1f' % (i * dt) + ' s')
    plt.savefig(input_file[:-2] + '.png')
    #plt.savefig(input_file[:-2] + '.svg')
    plt.close(fig)
