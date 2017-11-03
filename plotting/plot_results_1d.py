"""
plot_results_1d.py

Reads output of tsunami and plots 
the results in an image file.
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

# read data into a list
data = [line.rstrip().split() for line in open('adv_1d.txt').readlines()]

time = [float(line[0]) for line in data]
u = np.array([[float(x) for x in line[1:]] for line in data])
x = np.arange(1, u.shape[1]+1)
time_steps = [0, 25, 50, 75]

fig = plt.figure(figsize=(8,10))
axes = [plt.subplot2grid((4, 1), (row, 0), colspan=1, rowspan=1)
    for row in range(4)]

for ax in axes:
    n = axes.index(ax)
    ax.plot(x, u[time_steps[n], :])
    ax.grid(True)
    ax.set_xlim(1, 100)
    ax.set_ylim(0, 1)
    ax.set_ylabel('Amplitude', fontsize=16)
    ax.set_title('Solution at time step '+str(time_steps[n]), fontsize=16)
    ax.tick_params(axis='both', which='major', labelsize=16)

for ax in axes[:-1]:
    ax.set_xticklabels([])

axes[3].set_xlabel('', fontsize=16)
axes[-1].set_xlabel('Spatial grid index', fontsize=16)

plt.savefig('adv_1d.png')
plt.close(fig)
