import numpy as np
import matplotlib.pyplot as plt
import matplotlib
"""
plot_adv1d.py

Reads output of tsunami and plots the results in an image file.
"""
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('input_file', help='text file output by tsunami simulator')
args = parser.parse_args()

input_file = args.input_file

matplotlib.use('Agg')

# read data into a list
data = [line.rstrip().split() for line in open(input_file).readlines()]

time = [float(line[0]) for line in data]
u = np.array([[float(x) for x in line[1:]] for line in data])
x = np.arange(1, u.shape[1]+1)
time_steps = [0, 25, 50, 75]

fig = plt.figure(figsize=(8, 10))
axes = [plt.subplot2grid((4, 1), (row, 0), colspan=1, rowspan=1)
        for row in range(4)]

for ax in axes:
    n = axes.index(ax)
    ax.plot(x, u[time_steps[n], :], 'b-')
    ax.fill_between(x, 0, u[time_steps[n], :], color='b', alpha=0.4)
    ax.grid(True)
    ax.set_xlim(1, 100)
    ax.set_ylim(0, 1)
    ax.set_ylabel('Amplitude', fontsize=16)
    ax.set_xticks([25, 50, 75, 100])
    ax.set_yticks([0, 0.25, 0.5, 0.75, 1])
    ax.tick_params(axis='both', which='major', labelsize=16)

axes[0].set_title(r'Initial state, $u(x, t=0)$', fontsize=16)
for ax in axes[1:]:
    n = axes.index(ax)
    ax.set_title(r'$u(x, t=$'+str(time_steps[n])+r'$)$', fontsize=16)

for ax in axes[:-1]:
    ax.set_xticklabels([])

axes[3].set_xlabel('', fontsize=16)
axes[-1].set_xlabel('Spatial grid index', fontsize=16)

plt.savefig('fig_adv1d.png', dpi=100)
plt.close(fig)
