import numpy as np
from numpy import linalg as la
import cmath
from scipy.optimize import differential_evolution
import warnings

# Our infinity element in C_inf
inf = complex(float('inf'), float('inf'))

class Mobius:
    """Mobius Transform"""
    def __init__(self, a:complex, b:complex, c:complex, d:complex) -> None:
        sqrt_det = cmath.sqrt(a*d-b*c)

        if sqrt_det == 0:
            raise ValueError("This Mobius transformations is not invertable, ad-bc=0")
        
        self.a = a/sqrt_det
        self.b = b/sqrt_det
        self.c = c/sqrt_det
        self.d = d/sqrt_det

    def get_coefficients(self):
        return np.array([self.a,self.b,self.c,self.d])

    def is_identity(self,similarity=0):
        """Returns true if the matrix is the identity, similaity allows for a tolarence on how similar the two mobius transforms are."""
        if similarity == 0:
            return self.b == 0 and self.c == 0 and self.a == self.d
        else:
            d = self.uniform_metric(Mobius(1,0,0,1))
            if d <= similarity:
                return True
            else:
                return False

    def get_SL_2_C(self):
        """Returns the SL_2_C_homomorphism as [[a,b],[c,d]]"""
        return np.array([[self.a, self.b],[self.c, self.d]])
    
    def __matmul__(self, other: 'Mobius'):
        return GL_2_C_to_Mobius(self.get_SL_2_C()@other.get_SL_2_C())

    def apply_mobius(self,z:complex):
        """Basic Mobius Transform on z"""
        with np.errstate(divide='ignore', invalid='ignore'): # Suppresses division by zero warnings
            return (self.a * z + self.b) / (self.c * z + self.d)

    def get_conjugacy_class(self):
        """Returns conjugacy class of Mobius Transform as (k,is_identity(a,b,c,d)), note:
            k ≠ 1 iff two fixed points, with subclassificiations of:
                k ∈ (0,2)U(2,∞) iff hyperbolic transform,
                |k| = 1 iff eliptic transform, and
                otherwise iff loxodromic transform;
            k = 1 and non-identity iff parabolic and 1 fixed point; and 
            k = 1 and identity"""
        Tr=np.trace(self.get_SL_2_C()) #trace
        k=(Tr-cmath.sqrt(Tr**2-4))/(Tr+cmath.sqrt(Tr**2-4))
        return np.array([k,self.is_identity()])

    def uniform_metric(self,other:'Mobius'):
        """"Returns the uniform metric between two Mobius Transforms.
            Defined as d(f,g)=sup_{x∈ℂ}(σ(f(z),g(z))),
            where σ(z,w) is the length of the cord between the projection of z and w on the Riemann-sphere"""

        if (self@other).is_identity() == True:
            return 0
        else:
            def neg_dist(params):
                θ, φ = params
                if θ == 0:
                    z = inf
                else:
                    z = 1/np.tan(θ/2)*np.exp(1j*φ)

                return -σ_C(self.apply_mobius(z),other.apply_mobius(z))

            bounds = [(0, np.pi), (0, 2*np.pi)]
            result = differential_evolution(neg_dist, bounds, seed=42, polish=True)
            return -result.fun
    
    def get_fixed_points(self):
        eigenvalues, eigenvectors=la.eig(self.get_SL_2_C())
        with warnings.catch_warnings():
            warnings.filterwarnings('error', category=RuntimeWarning)
            try:
                z_1 = eigenvectors[0,0]/eigenvectors[1,0]
            except RuntimeWarning:
                z_1 = inf
        with warnings.catch_warnings():
            warnings.filterwarnings('error', category=RuntimeWarning)
            try:
                z_2 = eigenvectors[0,1]/eigenvectors[1,1]
            except RuntimeWarning:
                z_2 = inf
        return np.array([z_1,z_2])


def GL_2_C_to_Mobius(M):
    """GL(2,C) to Mobius transformation (obvious homomorphism)"""
    return Mobius(M[0,0],M[0,1],M[1,0],M[0,0])

def σ_C(z,w):
    #chordal metric on C
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
