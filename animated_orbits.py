import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib import rcParams

data = np.loadtxt("out.txt")    # Heliocentric coordinates

# Planets, display size and color
colors = ["Gold", "Red", "Black", "Blue", "Green", "Orange", "Olive", "darkslategrey", "teal", "slategray"]
planets = ["Sun", "Mercury", "Venuns", "Earth", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune", "Pluto"]
size = [15, 1.3, 4.5, 5, 2.5, 10, 9, 7.5, 7, 1]

PL = int(len(data[0,:])/2)  # Planet array length
IT = int(len(data[:,0]))    # Iteration length

# Heliocentric coordinates arrays
x = np.zeros((IT, PL))
y = np.zeros((IT, PL))

for i in range(0,IT):
    for j in range(0,PL):
        x[i][j] = data[i,(2*j)]   # (x) vector for iteration (i) planet (j)
        y[i][j] = data[i,(2*j+1)] # (y) vector for iteration (i) planet (j)

# Initialize subplots
fig, ax = plt.subplots()

# Plotting parameters
zlim = 6        # Limit of coordinate axes
trail = 950     # Iterations after which delete trail
interval = 10   # Interval (?)
skip = 10       # Number of iterations to skip in each frame

# Empty arrays of points and lines of each planet
xdata = np.empty(PL, dtype=object)
ydata = np.empty(PL, dtype=object)

xldata = np.empty(PL, dtype=object)
yldata = np.empty(PL, dtype=object)

# Point and lines objects for each planet
points = [plt.plot([],[],linestyle='',color=colors[i],marker='.',markersize=size[i],label=planets[i],animated=True)[0] for i in range(PL)]
lines = [plt.plot([],[],linestyle='-',color=colors[i],alpha=0.5,marker='',animated=True)[0] for i in range(PL)]

patches = points + lines

# Initialization function
def init():

    ax.set_xlim(-zlim,zlim)
    ax.set_ylim(-zlim,zlim)
    ax.set_title("Solar System")
    ax.set_xlabel("x [AU]")
    ax.set_ylabel("y [AU]")
    #ax.legend(loc='center left', bbox_to_anchor=(1,0.5))
    ax.legend()

    return patches

# Animation function
def update(n):

    # Restart plot at beginning of each loop
    if (n == 0):
        for i in range(PL):
            xldata[i] = np.empty(PL, dtype=object)
            yldata[i] = np.empty(PL, dtype=object)

    for i in range(PL):

        # Set data of point and line objects
        xdata[i] = x[n][i]
        ydata[i] = y[n][i]
        points[i].set_data(xdata[i], ydata[i])

        xldata[i] = np.append(xldata[i],x[n][i])
        yldata[i] = np.append(yldata[i],y[n][i])
        lines[i].set_data(xldata[i], yldata[i])

        # Set trail to a finite length
        if (n > trail):
            xldata[i] = np.delete(xldata[i], 0)
            yldata[i] = np.delete(yldata[i], 0)

    return patches

# Generate and save animation
anim = animation.FuncAnimation(fig, update, frames=np.arange(0,IT,skip),interval=interval,init_func=init,blit=True)
anim.save("animated_orbits.gif", writer='imagemagick', fps=30, dpi=300)

#plt.show()