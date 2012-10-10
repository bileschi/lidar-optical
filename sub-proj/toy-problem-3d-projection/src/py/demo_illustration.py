#!/usr/bin/python
#stdlib
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
#local
from rotations_3d import rot_mat
from color_lists import rainbow_colors
from camera_illustration import draw_camera, draw_3d_axis
from projection_3d_pt import project_3d
from camera import Camera

""" illustrates the projection of 3d points into an image.  8 3d points in cube
formation are illustrated in 3d, along with a movable camera.  The image of these
points in the camera is shown in a second figure"""

def get_sample_points():
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

def draw_all(
	pts3d,
	cam_t = (0,0,-10.0),
	cam_R = np.mat(np.diag([1,1,1]))):
	# draw 3d frame.
	fig = plt.figure(2)
	plt.clf()
	ax1 = fig.add_subplot(211, projection='3d')
	pts_colors = rainbow_colors(len(pts3d))
	for p, c in zip(pts3d, pts_colors):
		ax1.scatter(p[0], p[1], p[2], c=c, marker='^')
	draw_3d_axis(ax1)
	cam = Camera(t_x=cam_t[0], t_y=cam_t[1], t_z=cam_t[2], R=cam_R, k=500)
	draw_camera(ax1, cam)
	# draw image as seen by camera
	pts3d_img = project_3d(pts3d, cam)
	print pts3d_img
	ax2 = fig.add_subplot(212)
	for p, c in zip(pts3d_img, pts_colors):
		if p is not None:
			ax2.scatter(p.item(0), p.item(1), c=c, marker='^')
	ax2.set_xlabel('U')
	ax2.set_ylabel('V')
	ax2.set_xbound(0, 640)
	ax2.set_ybound(0, 480)
	plt.show()
	plt.pause(0.1)

def param_sweep():
	# cam_R = np.mat(np.diag([1,1,1]))
	# for z in range(-10,-2):
	# 	cam_t = np.mat([0.0, 0.0, z])
	# 	draw_all(cam_t = cam_t, cam_R = cam_R)
	pts3d = get_sample_points()
	cam_t = np.mat([0.0, 0.0, -5])
	for rot_z_deg in range(0,180, 15):
		cam_R = rot_mat(z_rads = np.deg2rad(rot_z_deg))
		cam_t = (0.0, 0.0, -5.0)
		draw_all(pts3d, cam_t = cam_t, cam_R = cam_R)
	for rot_x_deg in range(0, 30, 5):
		cam_R = rot_mat(x_rads = np.deg2rad(rot_x_deg), 
			z_rads = np.deg2rad(rot_z_deg))
		cam_t = (0.0, 0.0, -5.0)
		draw_all(pts3d, cam_t = cam_t, cam_R = cam_R)

if __name__ == "__main__":
	param_sweep()

