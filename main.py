from mobius_plotter import mobius_plot

a = complex(1, 0)
b = complex(1, 0)
c = complex(0, 0)
d = complex(1, 0)

mobius_plot(a, b, c, d)

"""
Example parameters.

Elliptic - f(z)=iz
a = complex(0, 1)
b = complex(0, 0)
c = complex(0, 0)
d = complex(1, 0)

Parabolic - f(z)=z+1
a = complex(1, 0)
b = complex(1, 0)
c = complex(0, 0)
d = complex(1, 0)

Hyperbolic - f(z)=2z
a = complex(2, 0)
b = complex(0, 0)
c = complex(0, 0)
d = complex(1, 0)

Loxodromic - f(z)=(1+1.5i)z
a = complex(1, 1.5)
b = complex(0, 0)
c = complex(0, 0)
d = complex(1, 0)
"""