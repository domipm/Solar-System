import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

#   Define constantes, arrays, and re-scaling

#   Sun and Earth Mass
mS = 1.99*10**30
mT = 5.97*10**24

c = 1.496*10**11
#   Gravitational constant
G = 6.67*10**(-11)
#   Constant Sun position (heliocentric approx.)
xS = 0
yS = 0

#   Initial Earth's position
xT = 149.6*10**6
yT = 0
#   Rewrite position to meter, rescale
xT = (xT * 1000) / c
yT = (yT * 1000) / c

#   Initial orbital velocity
vxT = 0
vyT = 29.8
#   Rewrite velocities to meter/second and rescale
vxT = (vxT * 1000 / c) * ( (G*mS/(c**3))**(-1/2) )
vyT = (vyT * 1000 / c) * ( (G*mS/(c**3))**(-1/2) )

#   Rescale mass
mS = mS / mS
mT = mT / mS

#   x,y-axis accelerations
axT = - mS * (xT - xS) / (abs(xT - xS)**3)
ayT = 0

#   Initial h,t parameters
h = 0.1
t = h

#   W function
wX = vxT + h/2 * axT
wY = vyT + h/2 * ayT

#   Open output file
open("earth_pos.txt", "w").close()

#   Append data to file
with open("earth_pos.txt", "w") as f:

    np.savetxt(f, np.c_[xT, yT])

    for k in range(0, 1000):
        #   New position
        xT = xT + h * vxT + h**2/2 * axT
        yT = yT + h * vyT + h**2/2 * ayT
        #   Update W function
        wX = vxT + h/2 * axT
        wY = vyT + h/2 * ayT
        #   New accelerations
        axT =  - mS * (xT - xS) / ((np.sqrt( (xT - xS)**2 + (yT - yS)**2 ) )**3)
        ayT =  - mS * (yT - yS) / ((np.sqrt( (xT - xS)**2 + (yT - yS)**2 ) )**3)
        #   New velocity
        vxT = wX + h/2 * axT
        vyT = wY + h/2 * ayT

        np.savetxt(f, np.c_[xT, yT])

#   Load file with positions
data = np.loadtxt("earth_pos.txt")

#   Generate animated plot
fig, ax = plt.subplots()

def animate(i):

    ax.clear()
    ax.set_title("Earth Orbit")
    ax.set_xlabel("x [AU]")
    ax.set_ylabel("y [AU]")
    ax.set_xlim([-1.5,1.5])
    ax.set_ylim([-1.5,1.5])
    ax.plot(0,0,'.',markersize=25,color="Yellow")
    ax.plot(data[i,0],data[i,1],".",color="Green",markersize=5)
    ax.plot(data[:,0],data[:,1],color="Blue",alpha=0.15)

anim = FuncAnimation(fig, animate, frames=60, interval=5, repeat=True)
anim.save("earth_orbit.gif", dpi=300, writer='imagemagick', fps=10)

plt.show()