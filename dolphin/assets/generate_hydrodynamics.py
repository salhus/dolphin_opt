from design_recipe.mesher import auv_mesh, flap_mesh
import os

# Get the absolute path to the existing design_recipe folder
design_recipe_folder = os.path.abspath(os.path.join(os.path.dirname(__file__), 'design_recipe'))

auv_length  = 10    # AUV Length
auv_radius  = 1     # AUV Radius
depth       = 0.40*auv_length     # Flap Depth (- z axis)
width       = 0.40*auv_length      # Flap Width (+ y-axis)
thickness   = 0.125*auv_radius # Flap Thickness (+- x-axis)
resolution  = 0.1

# auv_mesh(Length, radius, resolution, name , file_location)

auv_mesh(auv_length, auv_radius, resolution, 'auv_mesh' , design_recipe_folder)

# flap_mesh(corner, thickness, width, depth, resolution, name , file_location)

flap_mesh([-0.5*thickness, auv_radius, 0], 
          thickness, width, -depth, 
          5.0*resolution, 'flap_left' , 
          design_recipe_folder)

flap_mesh([-0.5*thickness, -(auv_radius+width), 0], 
          thickness, width, -depth, 
          5.0*resolution, 'flap_right' , 
          design_recipe_folder)