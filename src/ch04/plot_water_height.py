"""
plot_water_height.py

Plots the water elevation profile at one time step.
"""
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('input_file', help='text file output by tsunami simulator')
parser.add_argument('time_step', help='time step to plot')
args = parser.parse_args()

input_file = args.input_file
time_step = int(args.time_step)

import numpy as np
import matplotlib.pyplot as plt
import matplotlib

matplotlib.use('Agg')
matplotlib.rcParams.update({'font.size': 16})

# read data into a list
data = [line.rstrip().split() for line in open(input_file).readlines()]

time = [float(line[0]) for line in data]
h = np.array([[float(x) for x in line[1:]] for line in data])
x = np.arange(1, h.shape[1]+1)

fig = plt.figure(figsize=(8, 3))
ax = fig.add_axes((0.12, 0.2, 0.8, 0.7))
plt.ylim(-0.2, 1.2)
plt.xlim(1, 100)
plt.xticks(range(25, 125, 25))
plt.yticks(np.arange(-0.2, 1.4, 0.2))
x = np.arange(1, 101)
plt.plot(x, h[time_step], 'b-')
plt.fill_between(x, -0.5, h[time_step], color='b', alpha=0.4)
plt.grid()
plt.xlabel('Distance [m]')
plt.ylabel('Water elevation [m]')
plt.title(r'Water elevation [m], time step ' + str(time_step))
plt.savefig('water_height_' + '%4.4i' % time_step + '.svg')
plt.close(fig)
