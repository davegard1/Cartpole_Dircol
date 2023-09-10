import matplotlib.pyplot as plt
import numpy as np
from load_data import load_data




data = load_data("cartpole_dircol.dat", delim=",")
data_u = load_data("cartpole_dircol_u.dat", delim=",")



    




fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2,2)

ax1.plot(data[:,0],data[:,1])
ax2.plot(data[:,0],data[:,2])
ax3.plot(data[:,0],data[:,3])
ax4.plot(data[:,0],data[:,4])

ax1.set_ylabel("x (m)")
ax2.set_ylabel("xdot (m/s)")
ax3.set_ylabel("th (rad)")
ax4.set_ylabel("thdot (rad/s)")

fig, ax5 = plt.subplots()

ax5.plot(data_u[:,0],data_u[:,1])
ax5.set_xlabel("t (s)")
ax5.set_ylabel("u (N)")


# extents = np.array([getattr(ax, 'get_{}lim'.format(dim))() for dim in 'xyz'])
# sz = extents[:,1] - extents[:,0]
# centers = np.mean(extents, axis=1)
# maxsize = max(abs(sz))
# r = maxsize/2
# for ctr, dim in zip(centers, 'xyz'):
    # getattr(ax, 'set_{}lim'.format(dim))(ctr - r, ctr + r)


plt.show()