import time
import os
import capytaine as cpt
from capytaine import (FloatingBody,
                       BEMSolver,
                       RadiationProblem,
                       DiffractionProblem,
                       assemble_dataset)
import math
import numpy as np
from capytaine import FloatingBody
import meshio
from stl import mesh  # Make sure numpy-stl is installed

def export_hydrostatics(hydrostatics_directory, bodies):
    """Determine filenames (following Nemoh convention) and call the .dat file writer"""

    # if os.path.isdir(hydrostatics_directory):
    #     LOG.warning(f"""Exporting problem in already existing directory: {hydrostatics_directory}
    #          You might be overwriting existing files!""")
    # else:
    #     os.makedirs(hydrostatics_directory)
        
    def hydrostatics_writer(hydrostatics_file_path, kh_file_path, body):
        """Write the Hydrostatics.dat and KH.dat files"""
        with open(hydrostatics_file_path, 'w') as hf:
            for j in range(3):
                line =  f'XF = {body.center_of_buoyancy[j]:7.3f} - XG = {body.center_of_mass[j]:7.3f} \n'
                hf.write(line)
            line = f'Displacement = {body.volume:E}'
            hf.write(line)
            hf.close()
        np.savetxt(kh_file_path, body.hydrostatic_stiffness.values)

    if isinstance(bodies, FloatingBody):
        bodies = [bodies]
    
    hydrostatics_file_name = "Hydrostatics.dat"
    kh_file_name = "KH.dat"
    
    body_count = len(bodies)
    if body_count == 1:
        body = bodies[0]
        hydrostatics_file_path = os.path.join(hydrostatics_directory, hydrostatics_file_name)
        kh_file_path = os.path.join(hydrostatics_directory, kh_file_name)
        hydrostatics_writer(hydrostatics_file_path, kh_file_path, body)
    else:
        for (i, body) in enumerate(bodies):
            hydrostatics_file_path = os.path.join(hydrostatics_directory, f"Hydrostatics_{i}.dat")
            kh_file_path = os.path.join(hydrostatics_directory, f"KH_{i}.dat")
            hydrostatics_writer(hydrostatics_file_path, kh_file_path, body)

def prune_nodes(points, cells):
    # Only points/cells that actually used
    uvertices, uidx = np.unique(cells, return_inverse=True)
    cells = uidx.reshape(cells.shape)
    points = points[uvertices]
    return points, cells


def get_triangle_volumes(pts, cells):
    # Works in any dimension; taken from voropy
    local_idx = np.array([[1, 2], [2, 0], [0, 1]]).T
    idx_hierarchy = cells.T[local_idx]

    half_edge_coords = pts[idx_hierarchy[1]] - pts[idx_hierarchy[0]]
    ei_dot_ej = np.einsum(
        "ijk, ijk->ij", half_edge_coords[[1, 2, 0]], half_edge_coords[[2, 0, 1]]
    )

    vols = 0.5 * np.sqrt(
        +ei_dot_ej[2] * ei_dot_ej[0]
        + ei_dot_ej[0] * ei_dot_ej[1]
        + ei_dot_ej[1] * ei_dot_ej[2]
    )
    return vols


def get_simplex_volumes(pts, cells):
    """Signed volume of a simplex in nD. Note that signing only makes sense for
    n-simplices in R^n.
    """
    n = pts.shape[1]
    assert cells.shape[1] == n + 1

    p = pts[cells]
    p = np.concatenate([p, np.ones(list(p.shape[:2]) + [1])], axis=-1)
    return np.abs(np.linalg.det(p) / math.factorial(n))


def compute_volume(mesh):
    if "tetra" in mesh.cells_dict:
        vol = math.fsum(
            get_simplex_volumes(*prune_nodes(mesh.points, mesh.cells_dict["tetra"]))
        )
    elif "triangle" in mesh.cells_dict or "quad" in mesh.cells_dict:
        vol = 0.0
        if "triangle" in mesh.cells_dict:
            # triangles
            vol += math.fsum(
                get_triangle_volumes(
                    *prune_nodes(mesh.points, mesh.cells_dict["triangle"])
                )
            )
        if "quad" in mesh.cells_dict:
            # quad: treat as two triangles
            quads = mesh.cells_dict["quad"].T
            split_cells = np.column_stack(
                [[quads[0], quads[1], quads[2]], [quads[0], quads[2], quads[3]]]
            ).T
            vol += math.fsum(
                get_triangle_volumes(*prune_nodes(mesh.points, split_cells))
            )
    else:
        assert "line" in mesh.cells_dict
        segs = np.diff(mesh.points[mesh.cells_dict["line"]], axis=1).squeeze()
        vol = np.sum(np.sqrt(np.einsum("...j, ...j", segs, segs)))

    return vol


def plot(filename, points, triangles):
    from matplotlib import pyplot as plt

    pts = points[:, :2]
    for e in triangles:
        for idx in [[0, 1], [1, 2], [2, 0]]:
            X = pts[e[idx]]
            plt.plot(X[:, 0], X[:, 1], "-k")
    plt.gca().set_aspect("equal", "datalim")
    plt.axis("off")

    # plt.show()
    plt.savefig(filename, transparent=True)


def export_mesh_to_stl(wec_obj, destination_folder, filename):

    """
    Export the mesh from a WEC object to an STL file.

    Parameters:
    - wec_obj: object containing the mesh attribute
    - destination_folder: folder path where STL will be saved
    - filename: filename for STL (without extension)
    """

    mesh_obj = wec_obj.mesh

    # Check for vertex attribute (try _vertices, then vertices)
    if hasattr(mesh_obj, '_vertices'):
        vertices = mesh_obj._vertices
    elif hasattr(mesh_obj, 'vertices'):
        vertices = mesh_obj.vertices
    else:
        raise ValueError("mesh_obj must have '_vertices' or 'vertices' attribute")

    # Check for faces attribute (triangles)
    if hasattr(mesh_obj, '_faces'):
        faces = mesh_obj._faces
    elif hasattr(mesh_obj, 'faces'):
        faces = mesh_obj.faces
    else:
        raise ValueError("mesh_obj must have '_faces' or 'faces' attribute")

    # Convert vertices and faces to numpy arrays (if not already)
    vertices = np.array(vertices)
    faces = np.array(faces)

    # Create an empty mesh for numpy-stl
    # Each face is a triangle defined by 3 vertex indices
    stl_mesh = mesh.Mesh(np.zeros(faces.shape[0], dtype=mesh.Mesh.dtype))

    for i, face in enumerate(faces):
        for j in range(3):
            stl_mesh.vectors[i][j] = vertices[face[j]]

    # Create destination folder if it does not exist
    if not os.path.exists(destination_folder):
        os.makedirs(destination_folder)

    # Build full path with .stl extension
    filepath = os.path.join(destination_folder, f"{filename}.stl")

    # Save STL file
    stl_mesh.save(filepath)
    print(f"STL file saved to {filepath}")


def export_mesh_to_obj(wec_obj, destination_folder, filename):


    """
    Export the mesh from a WEC object to an obj file.

    Parameters:
    - wec_obj: object containing the mesh attribute
    - destination_folder: folder path where obj will be saved
    - filename: filename for obj (without extension)
    """

    mesh_obj = wec_obj.mesh

    # Check for vertex attribute (try _vertices, then vertices)
    if hasattr(mesh_obj, '_vertices'):
        vertices = mesh_obj._vertices
    elif hasattr(mesh_obj, 'vertices'):
        vertices = mesh_obj.vertices
    else:
        raise ValueError("mesh_obj must have '_vertices' or 'vertices' attribute")

    # Check for faces attribute (triangles)
    if hasattr(mesh_obj, '_faces'):
        faces = mesh_obj._faces
    elif hasattr(mesh_obj, 'faces'):
        faces = mesh_obj.faces
    else:
        raise ValueError("mesh_obj must have '_faces' or 'faces' attribute")

    # Convert vertices and faces to numpy arrays (if not already)
    vertices = np.array(vertices)
    faces = np.array(faces)

    # Create an empty mesh for numpy-obj
    # Each face is a triangle defined by 3 vertex indices
    obj_mesh = mesh.Mesh(np.zeros(faces.shape[0], dtype=mesh.Mesh.dtype))

    for i, face in enumerate(faces):
        for j in range(3):
            obj_mesh.vectors[i][j] = vertices[face[j]]

    # Create destination folder if it does not exist
    if not os.path.exists(destination_folder):
        os.makedirs(destination_folder)

    # Build full path with .obj extension
    filepath = os.path.join(destination_folder, f"{filename}.obj")

    # Save obj file
    obj_mesh.save(filepath)
    print(f"obj file saved to {filepath}")


def run_bemio():
    """
    Run the BEMIO MATLAB script to compute hydrodynamic coefficients.
    This function assumes that the MATLAB script is set up to generate the necessary plots.
    """
    # Ensure the MATLAB engine is started and the script is run
    try:
        import matlab.engine
        eng = matlab.engine.start_matlab()
        eng.run('./hydroData/bemio.m', nargout=0)
        eng.quit()
    except Exception as e:
        print(f"An error occurred while running BEMIO: {e}")
    print("BEMIO calculations completed successfully.")
# This function runs the BEMIO MATLAB script to compute hydrodynamic coefficients.


