# import open3d as o3d
# import numpy as np

# points = np.vstack((point_cloud.x, point_cloud.y, point_cloud.z)).transpose()

# pcd = o3d.geometry.PointCloud()
# pcd.points = o3d.utility.Vector3dVector(points)
# pcd.colors = o3d.utility.Vector3dVector(colors/65535)
# pcd.normals = o3d.utility.Vector3dVector(normals)o3d.visualization.draw_geometries([pcd])

import bonsai # my custom wrapper library
import tipsy # my custom .tipsy file utility 

step = 0.0625/1
data_prefix='data/plummer_snap_mpi'
figure_prefix='fig/plummer'
nStars = 65000 # when we run MPI this is number of stars per process
T_ = 0.0625/1

# import os
# path = os.getcwd()
# print(path)
# res = os.listdir('./data/')
# print(res)

star_obj = tipsy.read_tipsy(data_prefix,None, pointsize = 2,lim=5,nRed=int(nStars/2)) #no data is returned when we supply 'figure_prefix'
# print(star_obj)
a =12312