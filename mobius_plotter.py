from utils import mobius_transform, stereographic, inverse_stereographic, find_fixed_points, is_identity
import pyvista as pv
import numpy as np

def mobius_plot(a,b,c,d):
    
    
    fixed_z = find_fixed_points(a, b, c, d)
    fixed_coords = inverse_stereographic(np.array(fixed_z))

    sphere = pv.Icosphere(radius=1.0,nsub=1)

    # Logic to remove vectors on the fixed points
    temp_sphere=np.True_

    for i in range(np.shape(sphere.points)[0]):
        if (~np.all((np.round(sphere.points,6)[i]==np.round(fixed_coords[0,:],6))) and ~np.all((np.round(sphere.points,6)[i]==np.round(fixed_coords[1,:],6)))) and temp_sphere.dtype!='bool_':
            temp_sphere=np.vstack((temp_sphere,sphere.points[i]))
        elif (~np.all((np.round(sphere.points,6)[i]==np.round(fixed_coords[0,:],6))) and ~np.all((np.round(sphere.points,6)[i]==np.round(fixed_coords[1,:],6)))) and temp_sphere.dtype=='bool_':
            temp_sphere=sphere.points[i]

    sphere.points=temp_sphere

    z_points = stereographic(*sphere.points.T)
    w_points = mobius_transform(z_points, a, b, c, d)

    transformed_coords = inverse_stereographic(w_points)
    sphere["vectors"] = transformed_coords - sphere.points



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
    if is_identity(a,b,c,d) == False:
        for i, coord in enumerate(fixed_coords):
            print('hi')
            point_mesh = pv.Sphere(radius=0.01, center=coord)
            plotter.add_mesh(point_mesh, color="red", label=f"Fixed Point {i+1}" if i==0 else "")

    
    plotter.show()
