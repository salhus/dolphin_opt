import pygmsh
import numpy as np


def auv_mesh(Length, radius):
    ""


    with pygmsh.occ.Geometry() as geom:
    
        body    = geom.add_cylinder([0,0,0], [0,0,-Length], radius, angle=2 * np.pi)
        body.dim_tag = [3,1]

        head    = geom.add_ball([0, 0, -Length], radius, 0,mesh_size = 0.1)
        geom.rotate(head, [0,0,-Length,],1*np.pi,[0,1,0]) 
        head.dim_tag = [3,2]
    
        geom.boolean_union([head,body])
        auv = geom.generate_mesh()
        auv.write('auv_mesh.stl')

def flap(corner, width, depth, thickness, name):
    ""
    with pygmsh.occ.Geometry() as geom:
    
        extends = [width, depth, thickness]  # Size in x, y, z directions
        panel = geom.add_box(corner, extends, mesh_size=0.1)
        flap  = geom.generate_mesh()
        flap.write(f"{name}.stl")


if __name__ == "__main__":
    auv_mesh(10, 1)
    flap([0, 0.5,0], .125,1,-2,'flap_left')
    flap([0,-1.5,0], .125,1,-2,'flap_right')