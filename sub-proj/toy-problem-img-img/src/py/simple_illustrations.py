#!/usr/bin/python
import matplotlib.pyplot as plt
import numpy as np
from pylab import quiver, quiverkey

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
