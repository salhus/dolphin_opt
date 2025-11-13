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


# Setup directories

print('Current working directory:', os.getcwd())

resDir                        = './hydroData'
wec_mesh_file                 = './geometry/flap_N1p24.stl'
base_mesh_file                = './geometry/dolphin_N5p99.stl'

#######################
# Initialize bodies and wave frequency range

bodies               = []
rhoWater             = 1025.0
freq_range           = np.linspace(0.1, 10.1, 500)
# bem_headings         = np.linspace(0,np.pi*2,36)              # wave headings 0,np.pi*2,36


#######################
# Body 1
wec_mesh = cpt.load_mesh(wec_mesh_file, file_format='stl')
wec_mesh.translate([0, 0, -1.24])
wec_mesh.keep_immersed_part()
lid_mesh_wec = wec_mesh.generate_lid()
wec_obj  = FloatingBody(wec_mesh, lid_mesh = lid_mesh_wec)
wec_obj.center_of_mass = ([0,0,-1.24])#wec_mesh.center_of_buoyancy
wec_obj.rotation_center = wec_obj.center_of_mass
wec_obj.add_all_rigid_body_dofs()
wec_obj.compute_rigid_body_inertia()
wec_obj.compute_hydrostatics()

from helpers import export_mesh_to_stl
export_mesh_to_stl(wec_obj, './geometry', 'flap')

#######################
# Body 2
Length            = 8
pitch             = 0*0.25*np.pi
axle_depth        = 0.1875*Length
from capytaine import Axis
dolphin_mesh = cpt.load_mesh(base_mesh_file, file_format='stl')
my_axis = Axis(vector=[0, 1, 0], point=[0, 0, 0])
dolphin_mesh.rotate(axis=my_axis, angle=pitch)
dolphin_mesh.translate([0, 0, -5.99])
dolphin_mesh.keep_immersed_part()
lid_mesh_dolphin = dolphin_mesh.generate_lid()
dolphin_obj  = FloatingBody(dolphin_mesh, lid_mesh = lid_mesh_dolphin)
dolphin_obj.center_of_mass = ([0,0,-6])#dolphin_mesh.center_of_buoyancy
dolphin_obj.rotation_center = dolphin_obj.center_of_mass
dolphin_obj.add_all_rigid_body_dofs()
dolphin_obj.compute_rigid_body_inertia()
dolphin_obj.compute_hydrostatics()

from helpers import export_mesh_to_stl
export_mesh_to_stl(dolphin_obj, './geometry', 'dolphin')

bodies = cpt.FloatingBody.join_bodies(wec_obj, dolphin_obj)
# bodies.show()

# bodies.show_matplotlib()

from helpers import export_hydrostatics
export_hydrostatics(resDir, [wec_obj, dolphin_obj]) 


problems = [RadiationProblem(body=bodies, rho=rhoWater, radiating_dof=dof,
                             omega=omega,  free_surface=0, water_depth=250) 
                             for dof in bodies.dofs for omega in freq_range]
problems += [DiffractionProblem(body=bodies, rho=rhoWater, omega=omega,
                                free_surface=0, water_depth=250.0)
                                for omega in freq_range]

# solve all radiation and diffraction problems
timeStart = time.time()
solver = BEMSolver()

# save capytaine data




results = [solver.solve(pb) for pb in sorted(problems)]

data = assemble_dataset(results)

from helpers import convert_categoricals_to_strings

data_cleaned = convert_categoricals_to_strings(separate_complex_values(data))

data_cleaned.to_netcdf(
    os.path.join(resDir, 'poseidon.nc'),
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

output_file = os.path.join(resDir, 'poseidon')
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

run_bemio()