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
    voxels, scale, shift = stltovoxel.convert.convert_meshes(meshes, resolution=resolution-1, parallel=parallel)
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


def calculate_resolution_2(meshes, resolution_expected):
    # Expectation, the idx_L is already known when asking the resolution
    import numpy as np
    mesh_min, mesh_max = calculate_mesh_limits(meshes)
    ROI = mesh_max-mesh_min
    idx_L = np.argmax(resolution_expected)
    list_potential_resolution = []
    for ii in range(len(ROI)):
        if ii != idx_L:
            list_potential_resolution.append(int(ROI[idx_L]*resolution_expected[ii]/ROI[ii]))
        else:
            list_potential_resolution.append(resolution_expected[idx_L])
    resolution = np.min(list_potential_resolution)
    vx_size  = ROI[idx_L]/resolution
    res_vec = resolution_expected
    return int(resolution), idx_L, vx_size, np.array(res_vec)


def UL_convert_files_cyl_2(input_file_paths, output_file_path, output_file_path_voxel, output_file_path_meta, Diameter, internal_Diameter, seedn, resolution_expected=[64,64,128], parallel=False):
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
    resolution, idx_L,vx_size, res_vec = calculate_resolution_2(meshes, resolution_expected)
    voxels, scale, shift = stltovoxel.convert.convert_meshes(meshes, resolution=resolution-1, parallel=parallel)
    # padding
    pad = np.zeros((int(res_vec[2]),int(res_vec[1]),int(res_vec[0])),dtype=np.float32)
    offset = (np.asarray(pad.shape)-np.asarray(voxels.shape))//2
    rest_offset = (np.asarray(pad.shape)-np.asarray(voxels.shape))%2
    pad[offset[0]:offset[0]+voxels.shape[0],offset[1]:offset[1]+voxels.shape[1],offset[2]:offset[2]+voxels.shape[2]]=voxels
    # stl to voxel adds an empty slice at the beginning and two at the end of the stack 
    if offset[0]==0:
        pad[0,:,:]=pad[1,:,:]
        pad[pad.shape[0]-1,:,:]=pad[pad.shape[0]-3,:,:]
        pad[pad.shape[0]-2,:,:]=pad[pad.shape[0]-3,:,:]
    elif rest_offset[0]==0:
        for jj in range(offset[0]+1):
            pad[jj,:,:]=pad[offset[0]+1,:,:]
        for jj in range(offset[0]+voxels.shape[0]-2,pad.shape[0]):
            pad[jj,:,:]=pad[offset[0]+voxels.shape[0]-3,:,:]
    elif rest_offset[0]!=0:
        for jj in range(offset[0]):
            pad[jj,:,:]=pad[offset[0]+1,:,:]
        for jj in range(offset[0]+voxels.shape[0]-2,pad.shape[0]):
            pad[jj,:,:]=pad[offset[0]+voxels.shape[0]-3,:,:]
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
    imsave('voxelized_cyl_2.tif', pad)



def cut_before_inlet(voxels, voxel_size, chamber_height, expected_inlet_vx):
    # Expectation, the idx_L is already known when asking the resolution
    import numpy as np
    for ii in range(2,voxels.shape[0]-1):
        if np.sum(voxels[ii,:,:])<np.sum(voxels[ii+1,:,:]):
            begin_porous = ii
            break
    voxels_cut = voxels[begin_porous-(expected_inlet_vx-1):,:,:]
    res_vec = [voxels_cut.shape[2],voxels_cut.shape[1],voxels_cut.shape[0]]
    return voxels_cut, res_vec


def UL_convert_files_cyl_3(input_file_paths, output_file_path, output_file_path_voxel, output_file_path_meta, Diameter, internal_Diameter, seedn, chamber_height, expected_inlet_vx=6, voxel_size__=1e-3, parallel=False):
    import stltovoxel
    from stl import mesh
    import numpy as np
    meshes = []
    for input_file_path in input_file_paths:
        mesh_obj = mesh.Mesh.from_file(input_file_path)
        org_mesh = np.hstack((mesh_obj.v0[:, np.newaxis], mesh_obj.v1[:, np.newaxis], mesh_obj.v2[:, np.newaxis]))
        meshes.append(org_mesh)
    # 
    output = open(output_file_path, 'w')
    output_voxel = open(output_file_path_voxel, 'w')
    output_meta = open(output_file_path_meta, 'w')
    # 
    mesh_min, mesh_max = calculate_mesh_limits(meshes)
    ROI = mesh_max-mesh_min
    idx_L = np.argmax(ROI)
    # 
    resolution = int(round(ROI[idx_L]/voxel_size__))
    voxels, scale, shift = stltovoxel.convert.convert_meshes(meshes, resolution=resolution-1, parallel=parallel)
    # padding
    voxels_cut[0,:,:]=voxels_cut[1,:,:]
    voxels_cut[voxels_cut.shape[0]-2:,:,:]=voxels_cut[voxels_cut.shape[0]-3,:,:]
    for z in range(voxels_cut.shape[0]):
        for y in range(voxels_cut.shape[1]):
            for x in range(voxels_cut.shape[2]):
                    point = (np.array([x, y, z]) / scale) # + shift
                    # tag external of cylinder as ssolid to avoid flow
                    center_cyl = np.array([Diameter/2, Diameter/2])
                    if np.linalg.norm([x/scale[2],y/scale[1]]-center_cyl)-internal_Diameter/2 >= 0:
                        voxels_cut[z][y][x]=1
                    output.write('%s\t%s\t%s\t' % tuple(point))
                    output_voxel.write('%s\t' % voxels_cut[z][y][x])
    # 
    output_meta.write('seeding number: ')
    output_meta.write('%s\n' % seedn)
    output_meta.write('Voxel size = vx/m : ')
    output_meta.write('%s\n' % np.mean(scale))
    output_meta.write('Resolution = ')
    output_meta.write('%s x %s x %s \n' % tuple(res_vec))
    output_meta.write('domain dims = ')
    output_meta.write('%s x %s x %s' % tuple(res_vec/scale))
    output.close()
    output_voxel.close()
    output_meta.close()
    from tifffile import imsave
    tif_matrix = np.asarray(voxels_cut,dtype=np.float32)
    imsave('voxelized_cyl_3.tif', tif_matrix)

def enforce_size_x_y(voxels):
    # Expectation, the idx_L is already known when asking the resolution
    voxels_cut = voxels[:,:min(voxels.shape[0],voxels.shape[1]),:min(voxels.shape[0],voxels.shape[1])]
    res_vec = [voxels_cut.shape[2],voxels_cut.shape[1],voxels_cut.shape[0]]
    return voxels_cut, res_vec

def UL_convert_files_cyl_4(input_file_paths, output_file_path, output_file_path_voxel, output_file_path_meta, Diameter, external_Diameter, seedn, chamber_height, expected_inlet_vx=6, voxel_size__=1e-3, parallel=False):
    import stltovoxel
    from stl import mesh
    import numpy as np
    meshes = []
    for input_file_path in input_file_paths:
        mesh_obj = mesh.Mesh.from_file(input_file_path)
        org_mesh = np.hstack((mesh_obj.v0[:, np.newaxis], mesh_obj.v1[:, np.newaxis], mesh_obj.v2[:, np.newaxis]))
        meshes.append(org_mesh)
    # 
    output = open(output_file_path, 'w')
    output_voxel = open(output_file_path_voxel, 'w')
    output_meta = open(output_file_path_meta, 'w')
    # 
    mesh_min, mesh_max = calculate_mesh_limits(meshes)
    ROI = mesh_max-mesh_min
    idx_L = np.argmax(ROI)
    # 
    resolution = int(round(ROI[idx_L]/voxel_size__))
    voxels, scale, shift = stltovoxel.convert.convert_meshes(meshes, resolution=resolution-1, parallel=parallel)
    voxels[0,:,:]=voxels[1,:,:]
    voxels[voxels.shape[0]-2:,:,:]=voxels[voxels.shape[0]-3,:,:]
    # padding
    voxels_cut, res_vec = cut_before_inlet(voxels, scale, chamber_height, expected_inlet_vx)
    in_x, in_y, in_z, out_x, out_y, out_z = [], [], [], [], [], []
    for z in range(voxels_cut.shape[0]):
        for y in range(voxels_cut.shape[1]):
            for x in range(voxels_cut.shape[2]):
                    point = (np.array([x, y, z]) / scale) # + shift
                    # tag external of cylinder as ssolid to avoid flow
                    center_cyl = np.array([Diameter/2, Diameter/2])
                    if not (np.linalg.norm([x/scale[0]+ shift[0],y/scale[1]+ shift[1]]-center_cyl)-(external_Diameter)/2 >= voxel_size__/2 and voxels_cut[z][y][x]==0):
                        output.write('%s\t%s\t%s\t' % tuple(point))
                        output_voxel.write('%s\t' % voxels_cut[z][y][x])
                    else:
                        voxels_cut[z][y][x]=3
    # 
    voxels_cut, res_vec = enforce_size_x_y(voxels_cut)
    output_meta.write('seeding number: ')
    output_meta.write('%s\n' % seedn)
    output_meta.write('Voxel size = vx/m : ')
    output_meta.write('%s\n' % np.mean(scale))
    output_meta.write('Resolution = ')
    output_meta.write('%s x %s x %s \n' % tuple(res_vec))
    output_meta.write('domain dims = ')
    output_meta.write('%s x %s x %s' % tuple(res_vec/scale))
    output.close()
    output_voxel.close()
    output_meta.close()
    from tifffile import imsave
    tif_matrix = np.asarray(voxels_cut,dtype=np.float32)
    imsave('voxelized_cyl_4.tif', tif_matrix)



def UL_convert_files_cyl_4_debug(input_file_paths, output_file_path, output_file_path_voxel, output_file_path_meta, Diameter, external_Diameter, seedn, chamber_height, expected_inlet_vx=6, voxel_size__=1e-3, parallel=False):
    import stltovoxel
    from stl import mesh
    import numpy as np
    meshes = []
    for input_file_path in input_file_paths:
        mesh_obj = mesh.Mesh.from_file(input_file_path)
        org_mesh = np.hstack((mesh_obj.v0[:, np.newaxis], mesh_obj.v1[:, np.newaxis], mesh_obj.v2[:, np.newaxis]))
        meshes.append(org_mesh)
    # 
    output = open(output_file_path, 'w')
    output_voxel = open(output_file_path_voxel, 'w')
    output_meta = open(output_file_path_meta, 'w')
    # 
    mesh_min, mesh_max = calculate_mesh_limits(meshes)
    ROI = mesh_max-mesh_min
    idx_L = np.argmax(ROI)
    # 
    resolution = int(round(ROI[idx_L]/voxel_size__))
    voxels, scale, shift = stltovoxel.convert.convert_meshes(meshes, resolution=resolution-1, parallel=parallel)
    voxels[0,:,:]=voxels[1,:,:]
    voxels[voxels.shape[0]-2:,:,:]=voxels[voxels.shape[0]-3,:,:]
    # padding
    voxels_cut, res_vec = cut_before_inlet(voxels, scale, chamber_height, expected_inlet_vx)
    in_x, in_y, in_z, out_x, out_y, out_z = [], [], [], [], [], []
    for z in range(voxels_cut.shape[0]):
        for y in range(voxels_cut.shape[1]):
            for x in range(voxels_cut.shape[2]):
                    point = (np.array([x, y, z]) / scale) # + shift
                    # tag external of cylinder as ssolid to avoid flow
                    center_cyl = np.array([Diameter/2, Diameter/2])
                    if not (np.linalg.norm([x/scale[0]+ shift[0],y/scale[1]+ shift[1]]-center_cyl)-(external_Diameter)/2 >= 0 and voxels_cut[z][y][x]==0):
                        if voxels_cut[z][y][x]!=0:
                            in_x.append(x / scale[0])
                            in_y.append(y/ scale[1])
                            in_z.append(z / scale[2])
                        else:
                            out_x.append(x / scale[0])
                            out_y.append(y/ scale[1])
                            out_z.append(z / scale[2])
                        output.write('%s\t%s\t%s\t' % tuple(point))
                        output_voxel.write('%s\t' % voxels_cut[z][y][x])
    # 
    voxels_cut, res_vec = enforce_size_x_y(voxels_cut)
    import matplotlib.pyplot as plt
    import numpy as np
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.set_aspect('equal', 'box')
    ax.scatter(in_x, in_y, in_z, color='red')
    ax.scatter(out_x, out_y, out_z,color='black')
    plt.show()
    exit()

def UL_convert_files_cyl_5(input_file_paths, output_file_path, output_file_path_voxel, output_file_path_meta, Diameter, external_Diameter, seedn, chamber_height, expected_inlet_vx=6, voxel_size__=1e-3, parallel=False):
    import stltovoxel
    from stl import mesh
    import numpy as np
    meshes = []
    for input_file_path in input_file_paths:
        mesh_obj = mesh.Mesh.from_file(input_file_path)
        org_mesh = np.hstack((mesh_obj.v0[:, np.newaxis], mesh_obj.v1[:, np.newaxis], mesh_obj.v2[:, np.newaxis]))
        meshes.append(org_mesh)
    # 
    output = open(output_file_path, 'w')
    output_voxel = open(output_file_path_voxel, 'w')
    output_meta = open(output_file_path_meta, 'w')
    # 
    mesh_min, mesh_max = calculate_mesh_limits(meshes)
    ROI = mesh_max-mesh_min
    idx_L = np.argmax(ROI)
    # 
    resolution = int(round(ROI[idx_L]/voxel_size__))
    voxels, scale, shift = stltovoxel.convert.convert_meshes(meshes, resolution=resolution-1, parallel=parallel)
    voxels[0,:,:]=voxels[1,:,:]
    voxels[voxels.shape[0]-2:,:,:]=voxels[voxels.shape[0]-3,:,:]
    # padding
    voxels_cut, res_vec = cut_before_inlet(voxels, scale, chamber_height, expected_inlet_vx)
    in_x, in_y, in_z, out_x, out_y, out_z = [], [], [], [], [], []
    for z in range(voxels_cut.shape[0]):
        for y in range(voxels_cut.shape[1]):
            for x in range(voxels_cut.shape[2]):
                    point = (np.array([x, y, z]) / scale) # + shift
                    # tag external of cylinder as ssolid to avoid flow
                    center_cyl = np.array([Diameter/2, Diameter/2])
                    if (np.linalg.norm([x/scale[0]+ shift[0],y/scale[1]+ shift[1]]-center_cyl)-(external_Diameter)/2 >= 0 and voxels_cut[z][y][x]==0):
                        voxels_cut[z][y][x]=3
                    output.write('%s\t%s\t%s\t' % tuple(point))
                    output_voxel.write('%s\t' % voxels_cut[z][y][x])
    # 
    voxels_cut, res_vec = enforce_size_x_y(voxels_cut)
    output_meta.write('seeding number: ')
    output_meta.write('%s\n' % seedn)
    output_meta.write('Voxel size = vx/m : ')
    output_meta.write('%s\n' % np.mean(scale))
    output_meta.write('Resolution = ')
    output_meta.write('%s x %s x %s \n' % tuple(res_vec))
    output_meta.write('domain dims = ')
    output_meta.write('%s x %s x %s' % tuple(res_vec/scale))
    output.close()
    output_voxel.close()
    output_meta.close()
    from tifffile import imsave
    tif_matrix = np.asarray(voxels_cut,dtype=np.float32)
    imsave('voxelized_cyl_5.tif', tif_matrix)