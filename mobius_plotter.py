from utils import stereographic, inverse_stereographic, Mobius
import pyvista as pv
import numpy as np

def mobius_plot(a,b,c,d, r=0.01):

    f = Mobius(a,b,c,d)
    
    fixed_z = f.get_fixed_points()
    fixed_coords = inverse_stereographic(fixed_z)

    sphere = pv.Icosphere(radius=1.0,nsub=3)
    z_points = stereographic(*sphere.points.T)
    w_points = f.apply_mobius(z_points)

    transformed_coords = inverse_stereographic(w_points)
    vectors = transformed_coords - sphere.points

    if not f.is_identity() and len(fixed_coords) > 0:
        distances = np.linalg.norm(sphere.points[:, None, :] - fixed_coords[None, :, :], axis=2)
        min_distance = np.min(distances, axis=1)
        near_fixed = min_distance < r
        vectors[near_fixed] = 0.0

    sphere["vectors"] = vectors

    plotter = pv.Plotter()
    arrows = sphere.glyph(orient="vectors", scale=True, factor=0.2)
    plotter.add_mesh(arrows, color="grey")
    plotter.add_title("Möbius Transformation on S2")

    # Adding text
    def fmt(val):
        return f"{val.real:g} + {val.imag:g}j" if val.imag != 0 else f"{val.real:g}"
    info_text = (
        f"Möbius Parameters:\n"
        f"a = {fmt(a)}\n"
        f"b = {fmt(b)}\n"
        f"c = {fmt(c)}\n"
        f"d = {fmt(d)}"
    )
    plotter.add_text(
        info_text, 
        position='lower_left', 
        font_size=12, 
        color='black', 
        shadow=True, 
        name='params_label' # naming it allows you to update it later
    )

    # Adding fixed points
    if f.is_identity() == False:
        for i, coord in enumerate(fixed_coords):
            point_mesh = pv.Sphere(radius=0.01, center=coord)
            plotter.add_mesh(point_mesh, color="red", label=f"Fixed Point {i+1}" if i==0 else "")

    
    plotter.show()
