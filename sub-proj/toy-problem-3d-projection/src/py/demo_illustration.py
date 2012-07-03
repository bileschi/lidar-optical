#!/usr/bin/python
import matplotlib.pyplot as plt
import numpy as np
from pylab import quiver, quiverkey
from mpl_toolkits.mplot3d import Axes3D
import colorsys
import pdb

""" illustrates the projection of 3d points into an image.  8 3d points in cube
formation are illustrated in 3d, along with a movable camera.  The image of these
points in the camera is shown in a second figure"""

def get_box():
	""" a 2x2x2 box at (x,y,z) = (0,0,0) """
	box = []
	box.append((-1,+1,+1.0))
	box.append((-1,-1,+1.0))
	box.append((+1,+1,+1.0))
	box.append((+1,-1,+1.0))
	box.append((-1,+1,-1.0))
	box.append((-1,-1,-1.0))
	box.append((+1,+1,-1.0))
	box.append((+1,-1,-1.0))
	return box

def get_default_camera(size=2, k=1):
	""" lines describing a camera oriented down the z axis"""
	cam = []
	k = k * size/2.0
	cam.append([(0,0,0), (+k,+k,size)]) # edges
	cam.append([(0,0,0), (+k,-k,size)])
	cam.append([(0,0,0), (-k,-k,size)])
	cam.append([(0,0,0), (-k,+k,size)])
	cam.append([(+k,+k,size), (+k,-k,size)])
	cam.append([(+k,-k,size), (-k,-k,size)])
	cam.append([(-k,-k,size), (-k,+k,size)])
	cam.append([(-k,+k,size), (+k,+k,size)])
	cam.append([(0,0,0), (0,+k,size)]) # top marker
	return cam

def move_camera(in_cam, t=np.mat([0,0,0]), R=np.mat(np.diag([1,1,1]))):
	out_cam = []
	for line in in_cam:
		transformed_line = map(lambda pt: R.dot(pt) + t, line) # rotate then transalte
		out_cam.append(transformed_line)
	return out_cam

def rainbow_colors(n = 8):
	return [colorsys.hsv_to_rgb(float(h)/n, 1.0, 1.0) for h in range(n)]

def draw_all(
	cam_t = np.mat([0,0,-10.0]),
	cam_R = np.mat(np.diag([1,1,1]))):
	# draw 3d frame.
	fig = plt.figure(2)
	plt.clf()
	ax1 = fig.add_subplot(211, projection='3d')
	box = get_box()
	box_colors = rainbow_colors(len(box))
	for p, c in zip(box, box_colors):
		ax1.scatter(p[0], p[1], p[2], c=c, marker='^')
	draw_axis(ax1)
	ax1.set_xlabel('X')
	ax1.set_ylabel('Y')
	ax1.set_zlabel('Z')
	cam = move_camera(get_default_camera(), cam_t, cam_R)
	draw_camera(ax1, cam)
	# draw image as seen by camera
	box_img = projection(box, cam_t, cam_R)
	ax2 = fig.add_subplot(212)
	for p, c in zip(box_img, box_colors):
		ax2.scatter(p.item(0), p.item(1), c=c, marker='^')
	ax2.set_xlabel('U')
	ax2.set_ylabel('V')
	ax2.set_xbound(-2, 2)
	ax2.set_ybound(-2, 2)
	plt.show()
	plt.pause(0.1)

def projection(pts, t, R):
	imgs = []
	RI = R.I
	for pt in pts:
		tpt = RI.dot((pt - t).T)
		denom = tpt.item(2)
		if(denom == 0):
			imgs.append(np.mat([np.nan, np.nan]))
		else:
			imgs.append(np.mat([tpt.item(0) / denom, tpt.item(1) / denom]))
	return imgs

def draw_camera(ax, cam, color=[0,0,0]):
	for line in cam:
		ax.plot([line[0].item(0), line[1].item(0)], 
			[line[0].item(1), line[1].item(1)],
			[line[0].item(2), line[1].item(2)],
			c=color)

def draw_axis(ax, x_max = 12, y_max = 12, z_max = 12):
	for direction, color in zip((1, -1), [[0.3, 0.3, 0.3], [0.7, 0.3, 0.3]]):
		for point in np.diag(direction * np.array([x_max, y_max, z_max])):
			ax.plot([point[0]], [point[1]], [point[2]], 'w')
			ax.plot([0, point[0]], [0, point[1]], [0, point[2]], c=color)

def randrange(n, vmin, vmax):
    return (vmax-vmin)*np.random.rand(n) + vmin



def x_rot(rads):
  return np.matrix([
  [1,0,0], 
  [0, np.cos(rads), -np.sin(rads)],
  [0, np.sin(rads), np.cos(rads)]])

def y_rot(rads):
  return np.matrix([
  [np.cos(rads), 0, np.sin(rads)],
  [0, 1, 0], 
  [-np.sin(rads), 0, np.cos(rads)]])

def z_rot(rads):
  return np.matrix([
  [np.cos(rads), -np.sin(rads), 0],
  [np.sin(rads), np.cos(rads), 0],
  [0,0,1]]) 

def rot_mat(x_rads = 0 , y_rads = 0, z_rads = 0):
	x_mat = x_rot(x_rads)
	y_mat = y_rot(y_rads)
	z_mat = z_rot(z_rads)
	return x_mat.dot(y_mat.dot(z_mat))

def param_sweep():
	# cam_R = np.mat(np.diag([1,1,1]))
	# for z in range(-10,-2):
	# 	cam_t = np.mat([0.0, 0.0, z])
	# 	draw_all(cam_t = cam_t, cam_R = cam_R)
	cam_t = np.mat([0.0, 0.0, -5])
	for rot_z_deg in range(0,180, 15):
		cam_R = rot_mat(z_rads = np.deg2rad(rot_z_deg))
		cam_t = np.mat([0.0, 0.0, -5.0])
		draw_all(cam_t = cam_t, cam_R = cam_R)
	for rot_x_deg in range(0, 30, 5):
		cam_R = rot_mat(x_rads = np.deg2rad(rot_x_deg), 
			z_rads = np.deg2rad(rot_z_deg))
		cam_t = np.mat([0.0, 0.0, -5.0])
		draw_all(cam_t = cam_t, cam_R = cam_R)

##########
# tests #
##########

def test_raibow_colors():
	r8 = rainbow_colors(8)
	assert(len(r8) == 8)
	assert(all([len(c) == 3 for c in r8]))


