"""
plot_fig_shallow1d.py

Reads output of tsunami and plots the results in an image file.
"""
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('input_file', help='text file output by tsunami simulator')
args = parser.parse_args()

input_file = args.input_file

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

dt = 0.02

# read data into a list
data = [line.rstrip().split() for line in open(input_file).readlines()]

time = [float(line[0]) for line in data]
u = np.array([[float(x) for x in line[1:]] for line in data])
x = np.arange(1, u.shape[1]+1)
#time_steps = [0, 160, 7420, 9990]
time_steps = [0, 80, 3152, 4862]

fig = plt.figure(figsize=(8,10))
axes = [plt.subplot2grid((4, 1), (row, 0), colspan=1, rowspan=1)
    for row in range(4)]

for ax in axes:
    n = axes.index(ax)
    ax.plot(x, u[time_steps[n], :], 'b-')
    ax.fill_between(x, -0.2, u[time_steps[n], :], color='b', alpha=0.4)
    #ax.grid(True)
    ax.set_xlim(1, 100)
    ax.set_ylim(-0.2, 1.2)
    ax.set_ylabel('Height [m]', fontsize=16)
    ax.set_xticks([25, 50, 75, 100])
    ax.set_yticks(np.arange(-0.2, 1.4, 0.2))
    ax.tick_params(axis='both', which='major', labelsize=16)

for ax in axes:
    n = axes.index(ax)
    ax.set_title(r'Water elevation [m], time = '+'%5.1f' % (time_steps[n]*dt)+' s', fontsize=16)

for ax in axes[:-1]:
    ax.set_xticklabels([])

axes[3].set_xlabel('', fontsize=16)
axes[-1].set_xlabel('Distance [m]', fontsize=16)

axes[1].arrow(15, 0.6, -10, 0, color='k', transform=axes[1].transData, width=0.04, head_length=4)
axes[1].arrow(35, 0.6, 10, 0, color='k', transform=axes[1].transData, width=0.04, head_length=4)

plt.savefig('fig_shallow1d.png', dpi=100)
plt.savefig('fig_shallow1d.svg')
plt.close(fig)
