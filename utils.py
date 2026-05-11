import numpy as np
import cmath

# Our infinity element in C_inf
inf = complex(float('inf'), float('inf'))

def mobius_transform(z, a, b, c, d):
    """Basic Mobius Transform on z"""
    with np.errstate(divide='ignore', invalid='ignore'): # Suppresses division by zero warnings
        return (a * z + b) / (c * z + d)

def is_identity(a,b,c,d):
    return b==0 and c==0 and a==d

def stereographic(x1, x2, x3):
    """Projects R^3 -> C_inf"""
    with np.errstate(divide='ignore', invalid='ignore'): # Suppresses division by zero warnings
        z = (x1 + 1j * x2) / (1.0 - x3)
    return z

def inverse_stereographic(z):
    """Projects C_inf -> R^3"""
    z = np.asarray(z)
    is_inf = np.isinf(z)
    
    mag_sq = np.abs(z)**2
    denom = mag_sq + 1

    with np.errstate(invalid='ignore'):
        x1 = np.where(is_inf, 0, 2 * np.real(z) / denom)
        x2 = np.where(is_inf, 0, 2 * np.imag(z) / denom)
        x3 = np.where(is_inf, 1, (mag_sq - 1) / denom)

    return np.column_stack((x1, x2, x3))

def find_fixed_points(a, b, c, d):
    """Finds the fixed points of a given Mobius transform"""
    
    # Equation: cz^2 + (d-a)z - b = 0
    A, B, C = c, (d - a), -b
    
    if A == 0:
        p1 = complex('inf')
        p2 = b / (d - a) if (d - a) != 0 else complex('inf')
        return [p1, p2]

    # Quadratic formula
    discriminant = B**2 - 4*A*C
    sqrt_disc = cmath.sqrt(discriminant)
    
    z1 = (-B + sqrt_disc) / (2 * A)
    z2 = (-B - sqrt_disc) / (2 * A)
    return [z1, z2]
