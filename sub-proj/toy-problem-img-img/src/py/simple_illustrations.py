#!/usr/bin/python
#stdlib
import matplotlib.pyplot as plt
import numpy as np
from pylab import quiver, quiverkey
#local
from color_lists import rainbow_colors
from camera_illustration import draw_3d_axis, move_camera, draw_camera, get_default_camera
from projection_3d_pt import project_3d

def illustrate_points(img_pts, proj_pts, draw_associations = True):
	"""scatter plot of two point sets"""
	plt.figure(0)
	plt.clf()
	plt.plot([x for [x,y] in img_pts], [y for [x,y] in img_pts], 'ro')
	plt.plot([x for [x,y] in proj_pts], [y for [x,y] in proj_pts], 'bx')
	plt.axis([-2, 2, -2, 2])
	plt.grid(True)
	plt.title('image points in red.  projection points in blue.')
	plt.xlabel('x')
	plt.ylabel('y')
	if draw_associations:
		for (img_idx, img_pt) in enumerate(img_pts):
			draw_line(img_pts[img_idx], proj_pts[img_idx], color='#dddddd')
	plt.show()
	plt.pause(0.01)


def illustrate_forces(force_list):
	"""draws bar plot of force in x and y direction for a list of forces
	force_list is like [(f_x1, f_y1), (f_x2, f_y2), ... ]"""
	x_force = [x for (x,y) in force_list]
	y_force = [y for (x,y) in force_list]
	plt.figure(1)
	plt.clf()
	ax1 = plt.subplot(2, 1, 1)
	plt.bar(range(len(x_force)), x_force)
	plt.ylim([-20, 20])
	ax1.grid(True)
	plt.title('x direction')
	ax2 = plt.subplot(2, 1, 2)
	plt.bar(range(len(y_force)), y_force)
	plt.ylim([-20, 20])
	ax2.grid(True)
	plt.title('y direction')
	plt.show()
	plt.pause(0.01)

def illustrate_assoc(img_pts, proj_pts, assocs):
	"""draws lines connecting associated img_pt / proj_pt pairs"""
	for (img_idx, proj_idx) in assocs['pair_indicies']:
		draw_line(img_pts[img_idx], proj_pts[proj_idx], color='g')
	plt.show()
	plt.pause(0.01)


def illustrate_jacobian(
	space_pts = [(x/4.0,y/4.0) for x in range(-12, 12, 4) for y in range(-12, 12, 4)],
	proj_fcn = None,
	jacobian_fcn = None, 
	current_params = None,
	figure_idx = 1,
	clear_figure_first = True,
	param_idx = 0
	):
	""" for a set of space_pts, calculate their current image: img_pts.  
	Also calculate the derivative of the projection function 
	vs the parameter current_params[param_idx].
	This results in a vector representing how the img_pt will move if the 
	parameter is increased.   

	At each img_pt, draw an arrow indicating the direction and magnitude of the
	movement resulting from a change in current_params[param_idx]
	"""
	plt.figure(figure_idx)
	if clear_figure_first:
		plt.clf()
	img_pts = [proj_fcn(s_pt, current_params) for s_pt in space_pts]	
	X = np.array([img_pt[0] for img_pt in img_pts])
	Y = np.array([img_pt[1] for img_pt in img_pts])
	dX = np.zeros(X.shape)
	dY = np.zeros(Y.shape)
	for i, space_pt in enumerate(space_pts):
		J = jacobian_fcn(space_pt, current_params = current_params)
		dxdp = J[0,param_idx]
		dydp = J[1,param_idx]
		dX[i] = dxdp
		dY[i] = dydp
	Q = quiver( X, Y, dX, dY, units='width')
	plt.axis([-4, 4, -4, 4])
	plt.show()
	plt.pause(0.01)

def draw_line(img_pt, proj_pt, color='g'):
	plt.plot([img_pt[0], proj_pt[0]], [img_pt[1], proj_pt[1]], color = color)

def illustrate_toy_problem(space_pts, true_params, guess_params, 
	title1 = '3d space', title2 = 'image space'):
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
	true_img = project_3d(space_pts, true_params)
	guess_img = project_3d(space_pts, guess_params)
	ax1.set_title(title1)
	ax2 = fig.add_subplot(212)
	for p in true_img:
		if p is None:
			continue
		ax2.scatter(p.item(0), p.item(1), c=true_color, marker='^')
	for p in guess_img:
		if p is None:
			continue
		ax2.scatter(p.item(0), p.item(1), c=guess_color, marker='^')
	for p1, p2 in zip(true_img, guess_img):
		if (p1 is None) or (p2 is None):
			continue
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
	ax2.set_title(title2)
	plt.show()
	plt.pause(0.1)  


