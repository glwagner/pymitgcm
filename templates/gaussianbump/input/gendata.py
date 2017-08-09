import numpy as np

# Parameters
Lx = 4.0e5
Lz = 4500.0

nx = 80
ny = 42
nz = 8

dx = Lx / nx
dy = dx
dz = Lz / nz

# Physical parameters
g     = 9.81
alpha = 2.0e-4
f     = 1.0e-4
N     = 8.33e-4 #20.0*f

# Topography parameters
Lbump = 2.5e4
dbump = 0.5*Lz

Tz = N**2.0 / (g*alpha)

# Grid
x = np.arange(dx, Lx+dx, dx)
y = np.arange(dy, nx*dy+dy, dy)
z = np.arange(-0.5*dz, -Lz, -dz)

Y, X = np.meshgrid(y, x)

# Reference temperature profile
Tref = Tz*z - (Tz*z).mean()
print(Tref)

# Topography
bump = -Lz + dbump*np.exp( -(X**2.0 + Y**2.0) / (2*Lbump**2.0) )

print("dz(T) = {:.2e}".format(Tz))
