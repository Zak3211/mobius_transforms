# Möbius Transformation Visualizer

A Python-based tool for visualizing **Möbius Transformations** on the Riemann Sphere ($\mathbb{S}^2$). This project uses stereographic projection to map complex functions onto a 3D spherical surface, representing the transformation as a conformal vector field and a warped grid.

<p align="center">
  <img width="668" height="463" src="https://github.com/user-attachments/assets/222bbf90-a51a-48d7-9503-4f8719d3c107" />
</p>

## The Mathematics
A Möbius transformation is a bijective conformal map of the Riemann sphere onto itself, defined by the following rational function:

$$f(z) = \frac{az + b}{cz + d}$$

where $z, a, b, c, d \in \mathbb{C}$ and $ad - bc \neq 0$.

## Installation
This project requires **NumPy** for numerical computations and **PyVista** for 3D rendering. Install the dependencies using your terminal:

```powershell
py -m pip install numpy pyvista
```

## Usage
1. Open `main.py`.
2. Modify the values for the complex parameters `a`, `b`, `c`, and `d`.
3. Run the script to generate the interactive 3D visualization.

## Example Parameters
Möbius transformations are classified based on their trace and geometric behavior:

| Type | Function | Parameter Example |
| :--- | :--- | :--- |
| **Elliptic** (Rotation) | $f(z) = iz$ | `a=complex(0,1), b=0, c=0, d=1` |
| **Parabolic** (Translation) | $f(z) = z + 1$ | `a=1, b=1, c=0, d=1` |
| **Hyperbolic** (Dilation) | $f(z) = 2z$ | `a=2, b=0, c=0, d=1` |
| **Loxodromic** (Spiral) | $f(z) = (1+1.5i)z$ | `a=complex(1,1.5), b=0, c=0, d=1` |

## Features
*   **Vector Field Flow:** Arrows indicate the direction and magnitude of the transformation at every point.
*   **Grid Warping:** A transformed wireframe shows how the sphere's geometry is "stretched" or "compressed."
*   **Fixed Point Detection:** Automatic calculation and plotting of points where $f(z) = z$ using the quadratic formula.
*   **Riemann Sphere Mapping:** Full support for the "point at infinity" via robust stereographic projection handling.
