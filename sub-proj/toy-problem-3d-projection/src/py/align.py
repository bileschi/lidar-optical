#!/usr/bin/python
import numpy as np
from random import random
from rotations_3d import rot_mat

import math, sys, getopt
from scipy import spatial
import associations
from simple_illustrations import illustrate_points, illustrate_forces, illustrate_assoc, illustrate_jacobian
import copy
import pdb

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
	true_cam_params['k'] = 35
	true_cam_params['cx'] = 320
	true_cam_params['cy'] = 120
	return true_cam_params

def mutate_camera_params(in_params):
	out_params = {}
	out_params['t'] = in_params['t'].copy()
	R = in_params['R'].copy()
	out_params['R'] = rot_mat(0, 0, 0.05).dot(R)
	out_params['k'] = in_params['k']
	out_params['cx'] = in_params['cx']
	out_params['cy'] = in_params['cy']
	return out_params
	
def image_points(space_points, cam_params):
	imgs = []
	RI = cam_params['R'].I
	k = cam_params['k']
	cx = cam_params['cx']
	cy = cam_params['cy']
	for pt in space_points:
		tpt = RI.dot((pt - cam_params['t']).T)
		denom = tpt.item(2)
		if(denom == 0):
			imgs.append(np.mat([np.nan, np.nan]))
		elif denom < 0:
			continue # point is behind camera
		else:
			px = k * tpt.item(0) / denom + cx
			py = k * tpt.item(1) / denom + cy
			if (px < 0) or (py < 0) or (px > (cx * 2)) or (py > (cy * 2)):
				continue # point outside camera frame
			imgs.append(np.mat([px, py]))
	return imgs


def gen_3d_toy_problem():
	space_pts = rand_sample_cube(n_pts=10,
		x_min=-10, x_max=10,
		y_min=-10, y_max=10,
		z_min=10, z_max=20)
	true_params = select_initial_camera_pose()
	img_pts = image_points(space_pts, true_params)
	guess_params = mutate_camera_params(true_params)
	return [space_pts, img_pts, true_params, guess_params]

if __name__ == "__main__":
	# (1) problem generation
	[space_pts, img_pts, true_params, guess_params] = gen_3d_toy_problem()
	# (2) estimate true camera parameters 
	#     from space points and true image points 
	# (3) report results
	pdb.set_trace()









