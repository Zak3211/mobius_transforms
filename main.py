from mobius import Mobius

a = complex(1, 0)
b = complex(1, 0)
c = complex(0, 0)
d = complex(1, 0)

f = Mobius(a, b, c, d)
f.plot(vectors_scaled=False)