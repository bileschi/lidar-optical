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
from camera_illustration import move_camera, get_default_camera, draw_camera, draw_3d_axis
from optimization import estimate_projection_params
from projection_3d_pt import project_3d, proj_3d_jacobian, cam_P_to_dict
from rotations_3d import rot_mat
from simple_illustrations import illustrate_points, illustrate_forces, illustrate_assoc, illustrate_jacobian, illustrate_toy_problem

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
	true_cam_params = {}
	# external params
	true_cam_params['t'] = np.mat([0.0, 0.0, -10.0])
	true_cam_params['R'] = rot_mat(0.0, 0.0, 0.0)
	# internal params
	true_cam_params['k'] = 350
	true_cam_params['cx'] = 320
	true_cam_params['cy'] = 240
	return true_cam_params

def mutate_camera_params(in_params, 
	x_rot=0.0, y_rot=0.0, z_rot=0.0,
	t_x=0.0, t_y=0.0, t_z=0.0,
	k_mult=1.0):
	out_params = {}
	out_params['t'] = in_params['t'].copy() + np.mat([t_x, t_y, t_z])
	R = in_params['R'].copy()
	out_params['R'] = rot_mat(x_rot, y_rot, z_rot).dot(R)
	out_params['k'] = in_params['k'] * k_mult
	out_params['cx'] = in_params['cx']
	out_params['cy'] = in_params['cy']
	return out_params

def gen_3d_toy_problem():
	space_pts = rand_sample_cube(n_pts=10,
		x_min=-10, x_max=10,
		y_min=-10, y_max=10,
		z_min=10, z_max=50)
	true_params = select_initial_camera_pose()
	img_pts = project_3d(space_pts, true_params)
	guess_params = mutate_camera_params(true_params, x_rot=0.1)
	return [space_pts, img_pts, true_params, guess_params]

if __name__ == "__main__":
	[space_pts, img_pts, true_params, guess_params] = gen_3d_toy_problem()
	true_P_as_dict = cam_P_to_dict(true_params)
	guess_P_as_dict = cam_P_to_dict(guess_params)
	illustrate_toy_problem(space_pts, true_params, guess_params)
	# Build KDtree (if using approximate knn as association)
	kdtree = spatial.KDTree([pt.tolist()[0] for pt in img_pts])

	# perform optimization procedure to estimate true projection params
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









