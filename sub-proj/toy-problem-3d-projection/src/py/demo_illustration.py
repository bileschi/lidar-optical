#!/usr/bin/python
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from rotations_3d import rot_mat
from color_lists import rainbow_colors
from camera_illustration import move_camera, get_default_camera

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

def draw_camera(ax, cam, color=[0,0,0]):
	for line in cam:
		ax.plot([line[0].item(0), line[1].item(0)], 
			[line[0].item(1), line[1].item(1)],
			[line[0].item(2), line[1].item(2)],
			c=color)

def draw_3d_axis(ax, x_max = 12, y_max = 12, z_max = 12):
	for direction, color in zip((1, -1), [[0.3, 0.3, 0.3], [0.7, 0.3, 0.3]]):
		for point in np.diag(direction * np.array([x_max, y_max, z_max])):
			ax.plot([point[0]], [point[1]], [point[2]], 'w')
			ax.plot([0, point[0]], [0, point[1]], [0, point[2]], c=color)
	ax.set_xlabel('X')
	ax.set_ylabel('Y')
	ax.set_zlabel('Z')

def draw_all(
	pts3d,
	cam_t = np.mat([0,0,-10.0]),
	cam_R = np.mat(np.diag([1,1,1]))):
	# draw 3d frame.
	fig = plt.figure(2)
	plt.clf()
	ax1 = fig.add_subplot(211, projection='3d')
	pts_colors = rainbow_colors(len(pts3d))
	for p, c in zip(pts3d, pts_colors):
		ax1.scatter(p[0], p[1], p[2], c=c, marker='^')
	draw_3d_axis(ax1)
	cam = move_camera(get_default_camera(), cam_t, cam_R)
	draw_camera(ax1, cam)
	# draw image as seen by camera
	pts3d_img = projection(pts3d, cam_t, cam_R)
	ax2 = fig.add_subplot(212)
	for p, c in zip(pts3d_img, pts_colors):
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

def param_sweep():
	# cam_R = np.mat(np.diag([1,1,1]))
	# for z in range(-10,-2):
	# 	cam_t = np.mat([0.0, 0.0, z])
	# 	draw_all(cam_t = cam_t, cam_R = cam_R)
	pts3d = get_sample_points()
	cam_t = np.mat([0.0, 0.0, -5])
	for rot_z_deg in range(0,180, 15):
		cam_R = rot_mat(z_rads = np.deg2rad(rot_z_deg))
		cam_t = np.mat([0.0, 0.0, -5.0])
		draw_all(pts3d, cam_t = cam_t, cam_R = cam_R)
	for rot_x_deg in range(0, 30, 5):
		cam_R = rot_mat(x_rads = np.deg2rad(rot_x_deg), 
			z_rads = np.deg2rad(rot_z_deg))
		cam_t = np.mat([0.0, 0.0, -5.0])
		draw_all(pts3d, cam_t = cam_t, cam_R = cam_R)

if __name__ == "__main__":
	param_sweep()

