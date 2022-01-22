import matplotlib.pyplot as plt
from matplotlib import rcParams
import numpy as np
## Delay Plotting
K_axis = range(1,10)
J10S = [320 for _ in K_axis]
J10P = [313,220,189,173,173,158,142,142,142]
J10G = [313,237,206,206,206,206,206,206,206]
J6S = [192 for _ in K_axis]
J6P = [189,158,142,142,142,142,142,142,142]
J6G = [189,159,159,159,159,159,159,159,159]

# rcParams['font.family'] = 'Times New Roman'
rcParams['font.family'] = 'cmr10'
rcParams['mathtext.default'] = 'regular'



fig, ax = plt.subplots(1,1,figsize=(8,8),dpi=200)
ax.set_xlabel("$K$")
ax.set_ylabel("$t_J$ [ms]")
ax.set_ylim(125,325)
# ax.set_xlim(1,9)
ax.plot(K_axis,J10S,'r--',marker='d',markersize=5,label='Serial DAG with ${J=10}$')
ax.plot(K_axis,J10P,'y--',marker='d',markersize=5,label='Parallel DAG with $J=10$')
ax.plot(K_axis,J10G,'c--',marker='d',markersize=5,label='General DAG with $J=10$')
ax.plot(K_axis,J6S,'r',marker='o',markersize=5,label='Serial DAG with $J=6$')
ax.plot(K_axis,J6P,'y',marker='o',markersize=5,label='Parallel DAG with $J=6$')
ax.plot(K_axis,J6G,'c',marker='o',markersize=5,label='General DAG with $J=6$')
ax.legend(loc=(0.6,0.5))
ax.grid(ls='-.')
# fig.savefig('./DelayPlot.svg')

## Energy Plotting

Serl = [102.4 for _ in range(9)]
Gnrl = [123.5,123.5,123.5,123.5,123.5,123.5,123.5, 136,123.5]
Prl = [123.5,130.5,137.5,137.5,137.5,144.5,151.6,144.5,151.6]

K_axis = np.arange(1,10)
fig1, ax1 = plt.subplots(1,1,figsize=(8,8),dpi=200)
ax1.set_xlabel("$K$")
ax1.set_ylabel("$E^{max}$ [mJ]")
ax1.set_ylim(0, 250)
ax1.set_xticks(range(1,10))
ax1.bar(K_axis-0.25, Serl, color = 'xkcd:gold', width = 0.25, label='Serial DAG')
ax1.bar(K_axis, Prl, color = 'xkcd:dull blue', width = 0.25, label ='Parrel DAG')   # https://xkcd.com/color/rgb/
ax1.bar(K_axis+0.25, Gnrl, color = 'r', width = 0.25, label = 'General DAG')
ax1.legend()
plt.show()
# fig.savefig('./EnergyPlot.svg')