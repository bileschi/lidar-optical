#!/usr/bin/python
#stdlib
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from random import random
import math, sys, getopt
from scipy import spatial
import copy
import pdb
#local
import associations
from camera import Camera
from optimization import estimate_projection_params
from projection_3d_pt import project_3d, proj_3d_jacobian
from rotations_3d import rot_mat
#local illustrations
from simple_illustrations import illustrate_toy_problem

def rand_sample_cube(n_pts=1, x_min=0, x_max=1, y_min=0, y_max=1, z_min=0, z_max=1):
	" picks n_pts randomly in 3d. unit box.  pts returned as list-of-xyz-tuples"
	space_pts = []
	for i in range(0, n_pts):
		x = random() * (x_max - x_min) + x_min
		y = random() * (y_max - y_min) + y_min
		z = random() * (z_max - z_min) + z_min
		space_pts.append((x, y, z))
	return space_pts

def select_initial_camera_pose():
	cam = Camera(t_z=-10.0, k=350, c_x=320, c_y=240)
	return cam

def gen_3d_toy_problem():
	space_pts = rand_sample_cube(n_pts=10,
		x_min=-10, x_max=10,
		y_min=-10, y_max=10,
		z_min=10, z_max=50)
	true_params = select_initial_camera_pose()
	guess_params = true_params.mutate(r_x=0.1)
	img_pts = project_3d(space_pts, true_params)
	return [space_pts, img_pts, true_params, guess_params]

if __name__ == "__main__":
	[space_pts, img_pts, true_cam, guess_cam] = gen_3d_toy_problem()
	true_P_as_dict = true_cam.as_dict() 
	guess_P_as_dict = guess_cam.as_dict()
	print "true cam"
	print true_cam
	print "guess cam"
	print guess_cam
	illustrate_toy_problem(space_pts, true_cam, guess_cam)
	# Build KDtree (if using approximate knn as association)
	kdtree = spatial.KDTree([pt.tolist()[0] for pt in img_pts])

	# perform optimization procedure to estimate true projection params
	if False:
		est_params = estimate_projection_params(
				img_pts = img_pts, 
				space_pts = space_pts,
				projection_fcn = project_3d,
				jacobian_fcn = proj_3d_jacobian,
				# associate_fcn = associations.all_to_all,
				# associate_fcn = associations.all_to_nearest,
				# associate_fcn = associations.cheating,
				associate_fcn = lambda img_pts, proj_pts: associations.knn(
					proj_pts=proj_pts, kdtree=kdtree, k=2, eps=np.inf),
				guess_params = guess_P_as_dict,
				iterations = 25,
				# valid illustrate includes 'projection', 'association', 'jacobian'
				illustrate = set(),
				# illustrate = set(['projection', 'association']),
				# illustrate = set(['projection', 'association', 'jacobian']),
				verbose_on = True)

	# illustrate scaling
	# (2) estimate true camera parameters 
	#     from space points and true image points 
	# (3) report results









