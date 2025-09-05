import numpy as np

# 1D potentials
def harmonic_oscillator_1d(x, omega=1.0):
    return 0.5 * (omega * x)**2

def particle_in_box_1d(x, L=10.0):
    V = np.zeros_like(x)
    V[x < -L/2] = 1e6
    V[x > L/2] = 1e6
    return V

def double_well_1d(x, a=1.0, b=1.0):
    return a * x**4 - b * x**2

def finite_square_well_1d(x, V0=5.0, L=5.0):
    V = np.zeros_like(x)
    V[x < -L/2] = V0
    V[x > L/2] = V0
    return V

def potential_barrier_1d(x, V0=2.0, width=1.0, center=0.0):
    V = np.zeros_like(x)
    V[(x > center - width/2) & (x < center + width/2)] = V0
    return V

# 2D potentials
def harmonic_oscillator_2d(X, Y, omega_x=1.0, omega_y=1.0):
    return 0.5 * (omega_x**2 * X**2 + omega_y**2 * Y**2)

def particle_in_box_2d(X, Y, Lx=10.0, Ly=10.0):
    V = np.zeros_like(X)
    V[(X < -Lx/2) | (X > Lx/2) | (Y < -Ly/2) | (Y > Ly/2)] = 1e6
    return V

def double_well_2d(X, Y, a=1.0, b=1.0):
    return a * (X**4 + Y**4) - b * (X**2 + Y**2)