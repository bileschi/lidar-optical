#!/usr/bin/python
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from rotations_3d import rot_mat
from color_lists import rainbow_colors
from camera_illustration import move_camera, get_default_camera, draw_camera, draw_3d_axis


import numpy as np
from random import random
from rotations_3d import rot_mat
from color_lists import rainbow_colors
import matplotlib.pyplot as plt
from camera_illustration import move_camera, get_default_camera, draw_camera, draw_3d_axis
from mpl_toolkits.mplot3d import Axes3D

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
	space_pts.append((0, 0, 100))
	#todo remove center test point

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

def mutate_camera_params(in_params):
	out_params = {}
	out_params['t'] = in_params['t'].copy()
	R = in_params['R'].copy()
	out_params['R'] = rot_mat(0.10, 0.10, 0.25).dot(R)
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
		z_min=10, z_max=30)
	true_params = select_initial_camera_pose()
	img_pts = image_points(space_pts, true_params)
	guess_params = mutate_camera_params(true_params)
	return [space_pts, img_pts, true_params, guess_params]

def draw_line(img_pt, proj_pt, color='g'):
	plt.plot([img_pt[0], proj_pt[0]], [img_pt[1], proj_pt[1]], color = color)

def illustrate_toy_problem(space_pts, true_params, guess_params):
	true_color = [0, 0, 0]
	guess_color = [1, 0, 0]
	frame_color = [0, 0, 0]
	# draw 3d frame.
	fig = plt.figure(2)
	plt.clf()
	ax1 = fig.add_subplot(211, projection='3d')
	pts_colors = rainbow_colors(len(space_pts))
	for p, c in zip(space_pts, pts_colors):
		ax1.scatter(p[0], p[1], p[2], c=c, marker='^')
	draw_3d_axis(ax1)
	true_cam = move_camera(get_default_camera(), true_params['t'], true_params['R'])
	draw_camera(ax1, true_cam, color=true_color)
	guess_cam = move_camera(get_default_camera(), guess_params['t'], guess_params['R'])
	draw_camera(ax1, guess_cam, color=guess_color)
	# draw image as seen by cameras
	true_img = image_points(space_pts, true_params)
	guess_img = image_points(space_pts, guess_params)
	ax2 = fig.add_subplot(212)
	for p in true_img:
		ax2.scatter(p.item(0), p.item(1), c=true_color, marker='^')
	for p in guess_img:
		ax2.scatter(p.item(0), p.item(1), c=guess_color, marker='^')
	for p1, p2 in zip(true_img, guess_img):
		plt.plot([p1.item(0), p2.item(0)], [p1.item(1), p2.item(1)], color='g')
	ax2.set_xlabel('U')
	ax2.set_ylabel('V')
	# draw box around image area
	plt.plot([0, 640], [0, 0], color=frame_color)
	plt.plot([0, 640], [480, 480], color=frame_color)
	plt.plot([0, 0], [0, 480], color=frame_color)
	plt.plot([640, 640], [0, 480], color=frame_color)
	plt.plot([320, 320], [0, 480], color=[.2, .2, .2], ls=':')
	plt.plot([0, 640], [240, 240], color=[.2, .2, .2], ls=':')
	ax2.grid(True)
	ax2.set_xbound(-20, 660)
	ax2.set_ybound(-20, 500)
	plt.show()
	plt.pause(0.1)	

if __name__ == "__main__":
	# (1) problem generation
	[space_pts, img_pts, true_params, guess_params] = gen_3d_toy_problem()
	illustrate_toy_problem(space_pts, true_params, guess_params)
	# (2) estimate true camera parameters 
	#     from space points and true image points 
	# (3) report results









