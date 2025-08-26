import numpy as np

def particle_in_box(x, L=10):
    """Infinite potential well from -L/2 to L/2"""
    V = np.zeros_like(x)
    V[x < -L/2] = 1e10
    V[x > L/2] = 1e10
    return V

def harmonic_oscillator(x, k=1.0):
    """Harmonic oscillator potential: V = 0.5 * k * x^2"""
    return 0.5 * k * x**2

def double_well(x, a=5.0, b=0.5):
    """Double well potential: V = b*(x^2 - a^2)^2"""
    return b * (x**2 - a**2)**2
