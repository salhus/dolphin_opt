import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'assets'))

from assets.generate_hydrodynamics import generate_hydrodynamics_poseidon


generate_hydrodynamics_poseidon(
    auv_length  = 10,    # AUV Length
    auv_radius  = 1,     # AUV Radius
    depth       = 3,     # Flap Depth (- z axis)
    width       = 2,      # Flap Width (+ y-axis)
    thickness   = 0.1, # Flap Thickness (+- x-axis)
    resolution  = 0.1)
