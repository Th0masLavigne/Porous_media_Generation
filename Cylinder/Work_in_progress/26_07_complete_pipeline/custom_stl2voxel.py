# Based on the toolbox developed for stltofoxel package
# See licence file https://github.com/cpederkoff/stl-to-voxel/blob/master/LICENCE.md

def calculate_mesh_limits(meshes):
    import numpy as np
    mesh_min = meshes[0].min(axis=(0, 1))
    mesh_max = meshes[0].max(axis=(0, 1))
    for mesh in meshes[1:]:
        mesh_min = np.minimum(mesh_min, mesh.min(axis=(0, 1)))
        mesh_max = np.maximum(mesh_max, mesh.max(axis=(0, 1)))
    return mesh_min, mesh_max

def calculate_resolution(meshes, i_res=5):
    import numpy as np
    mesh_min, mesh_max = calculate_mesh_limits(meshes)
    ROI = mesh_max-mesh_min
    idx_L = np.argmax(ROI)
    list_potential_resolution = []
    for ii in range(len(ROI)):
        if ii != idx_L:
            list_potential_resolution.append(int(ROI[idx_L]*2**i_res/ROI[ii]))
        else:
            list_potential_resolution.append(2**(i_res+1))
    resolution = np.min(list_potential_resolution)
    return int(resolution), idx_L


def UL_convert_files(input_file_paths, output_file_path, output_file_path_voxel,i_res=5, parallel=False):
    import stltovoxel
    from stl import mesh
    import numpy as np
    meshes = []
    for input_file_path in input_file_paths:
        mesh_obj = mesh.Mesh.from_file(input_file_path)
        org_mesh = np.hstack((mesh_obj.v0[:, np.newaxis], mesh_obj.v1[:, np.newaxis], mesh_obj.v2[:, np.newaxis]))
        meshes.append(org_mesh)
    output = open(output_file_path, 'w')
    output_voxel = open(output_file_path_voxel, 'w')
    resolution, idx_L = calculate_resolution(meshes, i_res)
    voxels, scale, shift = stltovoxel.convert.convert_meshes(meshes, resolution-1, parallel)
    # padding
    pad = np.zeros((int(2**(i_res+1)),int(2**i_res),int(2**i_res)))
    pad[:voxels.shape[0],:voxels.shape[1],:voxels.shape[2]]=voxels
    for z in range(pad.shape[0]):
        for y in range(pad.shape[1]):
            for x in range(pad.shape[2]):
                    point = (np.array([x, y, z]) / scale) # + shift
                    output.write('%s\t%s\t%s\t' % tuple(point))
                    output_voxel.write('%s\t' % pad[z][y][x])
    output.close()
    output_voxel.close()