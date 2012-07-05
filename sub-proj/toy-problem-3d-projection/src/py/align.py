#!/usr/bin/python
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from rotations_3d import rot_mat
from camera_illustration import move_camera, get_default_camera, draw_camera, draw_3d_axis
from random import random
import math, sys, getopt
from scipy import spatial
import associations
from simple_illustrations import illustrate_points, illustrate_forces, illustrate_assoc, illustrate_jacobian, illustrate_toy_problem
import copy
import pdb
from projection_3d_pt import image_points

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
	img_pts = image_points(space_pts, true_params)
	guess_params = mutate_camera_params(true_params)
	return [space_pts, img_pts, true_params, guess_params]

if __name__ == "__main__":
	[space_pts, img_pts, true_params, guess_params] = gen_3d_toy_problem()
	# illustrate rotations
	if False:
		for i in range(0,10):
			guess_params = mutate_camera_params(guess_params, z_rot=0.1)
			illustrate_toy_problem(space_pts, true_params, guess_params)
			plt.pause(.2)
		guess_params = true_params
		for i in range(0,10):
			guess_params = mutate_camera_params(guess_params, y_rot=0.1)
			illustrate_toy_problem(space_pts, true_params, guess_params)
			plt.pause(.2)
		guess_params = true_params
		for i in range(0,10):
			guess_params = mutate_camera_params(guess_params, x_rot=0.1)
			illustrate_toy_problem(space_pts, true_params, guess_params)
			plt.pause(.2)
		guess_params = true_params
	# illustrate translations
	for i in range(0,10):
		guess_params = mutate_camera_params(guess_params, t_z=2)
		illustrate_toy_problem(space_pts, true_params, guess_params)
		plt.pause(.2)
	guess_params = true_params
	for i in range(0,10):
		guess_params = mutate_camera_params(guess_params, t_y=2)
		illustrate_toy_problem(space_pts, true_params, guess_params)
		plt.pause(.2)
	guess_params = true_params
	for i in range(0,10):
		guess_params = mutate_camera_params(guess_params, t_x=2)
		illustrate_toy_problem(space_pts, true_params, guess_params)
		plt.pause(.2)
	guess_params = true_params
	guess_params = mutate_camera_params(guess_params, k_mult=2)
	for i in range(0,10):
		guess_params = mutate_camera_params(guess_params, k_mult=0.8)
		illustrate_toy_problem(space_pts, true_params, guess_params)
		plt.pause(.2)
	# illustrate scaling
	# (2) estimate true camera parameters 
	#     from space points and true image points 
	# (3) report results









