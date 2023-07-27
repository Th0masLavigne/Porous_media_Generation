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
    vx_size  = ROI[idx_L]/resolution
    res_vec = [2**(i_res),2**(i_res), 2**(i_res)]
    res_vec[idx_L]=2**(i_res+1)
    return int(resolution), idx_L, vx_size, np.array(res_vec)


def UL_convert_files(input_file_paths, output_file_path, output_file_path_voxel, output_file_path_meta, seedn, i_res=5, parallel=False):
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
    output_meta = open(output_file_path_meta, 'w')
    resolution, idx_L,vx_size, res_vec = calculate_resolution(meshes, i_res)
    voxels, scale, shift = stltovoxel.convert.convert_meshes(meshes, resolution-1, parallel)
    # padding
    pad = np.zeros((int(res_vec[2]),int(res_vec[1]),int(res_vec[0])),dtype=np.float32)
    offset = (np.asarray(pad.shape)-np.asarray(voxels.shape))//2
    print(offset)
    pad[offset[0]:offset[0]+voxels.shape[0],offset[1]:offset[1]+voxels.shape[1],offset[2]:offset[2]+voxels.shape[2]]=voxels
    for z in range(pad.shape[0]):
        for y in range(pad.shape[1]):
            for x in range(pad.shape[2]):
                    point = (np.array([x, y, z]) / scale) # + shift
                    output.write('%s\t%s\t%s\t' % tuple(point))
                    output_voxel.write('%s\t' % pad[z][y][x])
    output_meta.write('seeding number: ')
    output_meta.write('%s\n' % seedn)
    output_meta.write('Voxel size = L_max/resolution : ')
    output_meta.write('%s\n' % vx_size)
    output_meta.write('Resolution = ')
    output_meta.write('%s x %s x %s \n' % tuple(res_vec))
    output_meta.write('domain dims = ')
    output_meta.write('%s x %s x %s' % tuple(res_vec*vx_size))
    output.close()
    output_voxel.close()
    output_meta.close()
    from tifffile import imsave
    imsave('voxelized.tif', pad)
    

def UL_convert_files_cyl(input_file_paths, output_file_path, output_file_path_voxel, output_file_path_meta, Diameter, internal_Diameter, seedn, i_res=5, parallel=False):
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
    output_meta = open(output_file_path_meta, 'w')
    resolution, idx_L,vx_size, res_vec = calculate_resolution(meshes, i_res)
    voxels, scale, shift = stltovoxel.convert.convert_meshes(meshes, resolution-1, parallel)
    # padding
    pad = np.zeros((int(2**(i_res+1)),int(2**i_res),int(2**i_res)),dtype=np.float32)
    offset = (np.asarray(pad.shape)-np.asarray(voxels.shape))//2
    rest_offset = (np.asarray(pad.shape)-np.asarray(voxels.shape))%2
    pad[offset[0]:offset[0]+voxels.shape[0],offset[1]:offset[1]+voxels.shape[1],offset[2]:offset[2]+voxels.shape[2]]=voxels
    # stl to voxel adds an empty slice at the beginning and two at the end of the stack 
    pad[0,:,:]=pad[1,:,:]
    pad[pad.shape[0]-1,:,:]=pad[pad.shape[0]-3,:,:]
    pad[pad.shape[0]-2,:,:]=pad[pad.shape[0]-3,:,:]
    for z in range(pad.shape[0]):
        for y in range(pad.shape[1]):
            for x in range(pad.shape[2]):
                    point = (np.array([x, y, z]) / scale) # + shift
                    # tag external of cylinder as ssolid to avoid flow
                    center_cyl = np.array([Diameter/2, Diameter/2])
                    xp_m_offset = np.array([(x-offset[2])/scale[2], (y-offset[1])/scale[1]]) # + shift
                    if np.linalg.norm(xp_m_offset-center_cyl)-internal_Diameter/2 >= 0:
                        pad[z][y][x]=1
                    output.write('%s\t%s\t%s\t' % tuple(point))
                    output_voxel.write('%s\t' % pad[z][y][x])
    output_meta.write('seeding number: ')
    output_meta.write('%s\n' % seedn)
    output_meta.write('Voxel size = L_max/resolution : ')
    output_meta.write('%s\n' % vx_size)
    output_meta.write('Resolution = ')
    output_meta.write('%s x %s x %s \n' % tuple(res_vec))
    output_meta.write('domain dims = ')
    output_meta.write('%s x %s x %s' % tuple(res_vec*vx_size))
    output.close()
    output_voxel.close()
    output_meta.close()
    from tifffile import imsave
    imsave('voxelized_cyl.tif', pad)