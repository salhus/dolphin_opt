import numpy as np
from helpers import compute_volume
import pygmsh
import time
import os
import capytaine as cpt
from capytaine import (FloatingBody,
                       BEMSolver,
                       RadiationProblem,
                       DiffractionProblem,
                       assemble_dataset)
from capytaine.io.xarray import separate_complex_values
from capytaine.post_pro import impedance
from capytaine.post_pro import rao
import sys
import vtk
import matplotlib.pyplot
import shutil
from design_recipe.mesher import auv_mesh, flap_mesh

# Setup directories

print('Current working directory:', os.getcwd())

# Get the absolute path to the existing design_recipe folder
design_recipe_folder = os.path.abspath(
    os.path.join(os.path.dirname(__file__), 'design_recipe'))


def generate_hydrodynamics_poseidon(
    auv_length  = 10,    # AUV Length
    auv_radius  = 1,     # AUV Radius
    depth       = 4,     # Flap Depth (- z axis)
    width       = 4,      # Flap Width (+ y-axis)
    thickness   = 0.1, # Flap Thickness (+- x-axis)
    resolution  = 0.1):

    h5_filename = f"poseidon_L{auv_length}_R{auv_radius}_flap{width}x{depth}"
    # auv_mesh(Length, radius, resolution, name , file_location)

    auv_mesh(auv_length, 
             auv_radius, 
             resolution, 'auv_mesh' , 
             design_recipe_folder)

    # flap_mesh(corner, thickness, width, depth, resolution, name , file_location)

    flap_mesh([-0.5*thickness, auv_radius, 0], 
            thickness, width, -depth, 
            5.0*resolution, 'flap_left' , 
            design_recipe_folder)

    flap_mesh([-0.5*thickness, -(auv_radius+width), 0], 
            thickness, width, -depth, 
            5.0*resolution, 'flap_right' , 
            design_recipe_folder)
    

    resDir = os.path.abspath(
             os.path.join(os.path.dirname(__file__), 'hydroData'))
    
    print('Results directory:', resDir)

    auv_mesh_file  = os.path.abspath(
             os.path.join(os.path.dirname(__file__), 
                          'design_recipe',
                          'auv_mesh.stl'))
    
    flap_left_file = os.path.abspath(
             os.path.join(os.path.dirname(__file__), 
                          'design_recipe',
                          'flap_left.stl'))

    flap_right_file = os.path.abspath(
             os.path.join(os.path.dirname(__file__), 
                          'design_recipe',
                          'flap_right.stl'))
    

    #######################
    # Initialize bodies and wave frequency range

    bodies               = []
    rhoWater             = 1025.0
    freq_range           = np.linspace(0.1, 10.1, 500)
    # bem_headings         = np.linspace(0,np.pi*2,36)             
    #  # wave headings 0,np.pi*2,36   

    #######################
    ## Import AUV Body
    auv_mesh_body = cpt.load_mesh(auv_mesh_file, file_format='stl')
    auv_mesh_body.keep_immersed_part()
    lid_auv_mesh_body = auv_mesh_body.generate_lid()
    auv_obj = FloatingBody(auv_mesh_body, lid_mesh = lid_auv_mesh_body)
    auv_obj.center_of_mass = auv_obj.center_of_buoyancy - np.array([0, 0, -0.25 * auv_length])
    auv_obj.rotation_center= auv_obj.center_of_mass
    auv_obj.add_all_rigid_body_dofs()
    auv_obj.compute_rigid_body_inertia()
    auv_obj.compute_hydrostatics()

    # auv_obj.show()

    #######################
    ## Import Flap Bodies

    flap_left_mesh_body = cpt.load_mesh(flap_left_file, file_format='stl')
    flap_left_mesh_body.keep_immersed_part()
    lid_flap_left_mesh_body = flap_left_mesh_body.generate_lid()
    flap_left_obj = FloatingBody(flap_left_mesh_body, 
                                 lid_mesh = lid_flap_left_mesh_body)
    flap_left_obj.center_of_mass = flap_left_obj.center_of_buoyancy - np.array([0, 0, -0.25 * depth])
    flap_left_obj.rotation_center= flap_left_obj.center_of_mass
    flap_left_obj.add_all_rigid_body_dofs()
    flap_left_obj.compute_rigid_body_inertia()
    flap_left_obj.compute_hydrostatics()

    flap_right_mesh_body = cpt.load_mesh(flap_right_file, file_format='stl')
    flap_right_mesh_body.keep_immersed_part()
    lid_flap_right_mesh_body = flap_right_mesh_body.generate_lid()
    flap_right_obj = FloatingBody(flap_right_mesh_body, 
                                 lid_mesh = lid_flap_right_mesh_body)
    flap_right_obj.center_of_mass = flap_right_obj.center_of_buoyancy - np.array([0, 0, -0.25 * depth])
    flap_right_obj.rotation_center= flap_right_obj.center_of_mass
    flap_right_obj.add_all_rigid_body_dofs()
    flap_right_obj.compute_rigid_body_inertia()
    flap_right_obj.compute_hydrostatics()
    
    #######################

    bodies = cpt.FloatingBody.join_bodies(auv_obj, 
                                          flap_left_obj, 
                                          flap_right_obj)

    bodies.show()


    from helpers import export_hydrostatics
    export_hydrostatics(resDir, [auv_obj, 
                                 flap_left_obj, 
                                 flap_right_obj]) 


    problems = [RadiationProblem(body=bodies, rho=rhoWater, radiating_dof=dof,
                                omega=omega,  free_surface=0, water_depth=250) 
                                for dof in bodies.dofs for omega in freq_range]
    problems += [DiffractionProblem(body=bodies, rho=rhoWater, omega=omega,
                                    free_surface=0, water_depth=250.0)
                                    for omega in freq_range]

    # solve all radiation and diffraction problems
    timeStart = time.time()
    solver = BEMSolver()

    results = [solver.solve(pb) for pb in sorted(problems)]

    data = assemble_dataset(results)

    from helpers import convert_categoricals_to_strings

    data_cleaned = convert_categoricals_to_strings(separate_complex_values(data))

    data_cleaned.to_netcdf(
        os.path.join(resDir, f"{h5_filename}.nc"),
        encoding={
            'radiating_dof': {'dtype': 'U'},
            'influenced_dof': {'dtype': 'U'}
        }
    )


    import xarray as xr
    Zi = impedance(data)
    Zi.name = 'impedance'
    RAO      =  rao(data)
    RAO.name =  'RAO'
    data = xr.merge([data, Zi])
    data = xr.merge([data, RAO])

    print(data.keys)

    print('Capytaine run(s) complete')

    import scipy.io as sio

    output_file = os.path.join(resDir, h5_filename)
    sio.savemat(output_file + '.mat',
    mdict={'added_mass':data['added_mass'].values,
        'radiation_damping':data['radiation_damping'].values,
        'Hexc':(data['Froude_Krylov_force'] + data['diffraction_force']).squeeze().values,
        'omega':data['omega'].values,
        'mass':data['inertia_matrix'].values,
        'hydrostatic_stiffness':data['hydrostatic_stiffness'].values,
        'impedance':data['impedance'],
        'rao':data['RAO'],
        },
                do_compression=True,
    )

    print('Capytaine runs were saved as a .mat file')

    from helpers import run_bemio

    run_bemio(h5_filename, resDir)

if __name__ == "__main__":
    generate_hydrodynamics_poseidon(
    auv_length  = 10,    # AUV Length
    auv_radius  = 1,     # AUV Radius
    depth       = 3,     # Flap Depth (- z axis)
    width       = 2,      # Flap Width (+ y-axis)
    thickness   = 0.1, # Flap Thickness (+- x-axis)
    resolution  = 0.1)

    
    



