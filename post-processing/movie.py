import numpy as np
import sys

import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter

data = np.load('S_checkpoint_pv_all.npz')
pos=data['pos']
print(pos.shape)

x = pos[:, 0, 0]
y = pos[:, 0, 1]
z = pos[:, 0, 2]

fig = plt.figure()
ax = fig.gca(projection='3d')

# Plot the surface.
surf = ax.plot_surface(X, Y, Z, linewidth=0, antialiased=False)

# Labels.
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')

if azim is not None:
    ax.azim = azim
if dist is not None:
    ax.dist = dist
if elev is not None:
    ax.elev = elev

print('ax.azim = {}'.format(ax.azim))
print('ax.dist = {}'.format(ax.dist))
print('ax.elev = {}'.format(ax.elev))

plt.savefig(
    'main_{}_{}_{}.png'.format(ax.azim, ax.dist, ax.elev),
    format='png',
    bbox_inches='tight'
)


