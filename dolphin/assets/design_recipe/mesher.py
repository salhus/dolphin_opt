import pygmsh
import numpy as np


def auv_mesh(Length, radius, resolution):
    ""


    with pygmsh.occ.Geometry() as geom:
    
        body    = geom.add_cylinder([0,0,0], [0,0,-Length], 
                                    radius, 
                                    angle=2 * np.pi,
                                    mesh_size=resolution)
        body.dim_tag = [3,1]

        head    = geom.add_ball([0, 0, -Length], radius, 0,mesh_size = resolution)
        geom.rotate(head, [0,0,-Length,],1*np.pi,[0,1,0]) 
        head.dim_tag = [3,2]
    
        geom.boolean_union([head,body])
        auv = geom.generate_mesh()
        auv.write('auv_mesh.stl')

def flap(corner, thickness, width, depth, name, resolution):
    ""
    with pygmsh.occ.Geometry() as geom:
    
        extends = [thickness, width, depth]  # Size in x, y, z directions
        panel = geom.add_box(corner, extends, mesh_size=resolution)
        flap  = geom.generate_mesh()
        flap.write(f"{name}.stl")


if __name__ == "__main__":
    
    auv_length  = 10    # AUV Length
    auv_radius  = 1     # AUV Radius
    depth       = 4     # Flap Depth (- z axis)
    width       = 4     # Flap Width (+ y-axis)
    thickness   = 0.125 # Flap Thickness (+- x-axis)
    resolution  = 0.1
    auv_mesh(auv_length, auv_radius, resolution)

    # flap(corner,thickness, width, depth, name)
    flap([0, auv_radius,0], 
         thickness, width, -depth,
         'flap_left', 5*resolution) 
    flap([0,-(auv_radius+width),0], 
         thickness, width, -depth,
         'flap_right', 5*resolution)