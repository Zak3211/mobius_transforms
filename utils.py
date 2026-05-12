import numpy as np
from numpy import linalg as la

def σ_C(z,w):
    """chordal metric on C"""
    Z = inverse_stereographic(z)
    W = inverse_stereographic(w)
    return la.norm(Z-W)

def stereographic(x1, x2, x3):
    """Projects R^3 -> C_inf"""
    with np.errstate(divide='ignore', invalid='ignore'): # Suppresses division by zero warnings
        z = (x1 + 1j * x2) / (1.0 - x3)
    return z

def inverse_stereographic(z):
    """Projects C_inf -> R^3"""
    is_inf = np.isinf(z)
    
    mag_sq = np.abs(z)**2
    denom = mag_sq + 1

    with np.errstate(invalid='ignore'):
        x1 = np.where(is_inf, 0, 2 * np.real(z) / denom)
        x2 = np.where(is_inf, 0, 2 * np.imag(z) / denom)
        x3 = np.where(is_inf, 1, (mag_sq - 1) / denom)

    return np.column_stack((x1, x2, x3))
