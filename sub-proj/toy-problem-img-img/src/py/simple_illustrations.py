#!/usr/bin/python
import matplotlib.pyplot as plt

def illustrate_points(img_pts, proj_pts):
	"""scatter plot of two point sets"""
	plt.figure(0)
	plt.clf()
	plt.plot([x for [x,y] in img_pts], [y for [x,y] in img_pts], 'ro')
	plt.plot([x for [x,y] in proj_pts], [y for [x,y] in proj_pts], 'bx')
	plt.axis([-7, 7, -5, 5])
	plt.grid(True)
	plt.title('image points in red.  projection points in blue.')
	plt.xlabel('x')
	plt.ylabel('y')
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


def draw_line(img_pt, proj_pt, color='g'):
	plt.plot([img_pt[0], proj_pt[0]], [img_pt[1], proj_pt[1]], color = color)
