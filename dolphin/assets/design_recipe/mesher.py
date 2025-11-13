import pygmsh
import numpy as np
import os

def auv_mesh(Length, radius, resolution, name , file_location):
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
        auv_mesh = geom.generate_mesh()
        # Combine file_location and name for the output file
        output_path = os.path.join(file_location, f"{name}.stl")
        auv_mesh.write(output_path)

def flap_mesh(corner, thickness, width, depth, resolution, name , file_location):
    ""
    with pygmsh.occ.Geometry() as geom:
    
        extends = [thickness, width, depth]  # Size in x, y, z directions
        panel = geom.add_box(corner, extends, mesh_size=resolution)
        flap_mesh = geom.generate_mesh()
        # Combine file_location and name for the output file
        output_path = os.path.join(file_location, f"{name}.stl")
        flap_mesh.write(output_path)


