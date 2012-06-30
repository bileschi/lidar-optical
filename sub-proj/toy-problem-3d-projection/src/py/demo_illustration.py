#!/usr/bin/python
import matplotlib.pyplot as plt
import numpy as np
from pylab import quiver, quiverkey
from mpl_toolkits.mplot3d import Axes3D
import colorsys

""" illustrates the projection of 3d points into an image.  8 3d points in cube
formation are illustrated in 3d, along with a movable camera.  The image of these
points in the camera is shown in a second figure"""

def get_box():
	""" a 2x2x2 box at z = 10 """
	box = []
	box.append((-1,1,10))
	box.append((-1,-1,10))
	box.append((1,1,10))
	box.append((1,-1,10))
	box.append((-1,1,12))
	box.append((-1,-1,12))
	box.append((1,1,12))
	box.append((1,-1,12))
	return box

def get_default_camera(size=1, k=1):
	""" lines describing a camera oriented down the z axis"""
	cam = []
	k = k * size/2.0
	cam.append([(0,0,0), (+k,+k,size)]) # edges
	cam.append([(0,0,0), (+k,-k,size)])
	cam.append([(0,0,0), (-k,-k,size)])
	cam.append([(0,0,0), (-k,+k,size)])
	cam.append([(+k,+k,size), (+k,-k,size)])
	cam.append([(-k,+k,size), (-k,-k,size)])
	cam.append([(-k,-k,size), (-k,+k,size)])
	cam.append([(-k,+k,size), (+k,+k,1)])
	cam.append([(0,0,0), (0,+k,1)]) # top marker
	return cam

def move_camera(in_cam, t=np.mat([0,0,0]), R=np.mat(np.diag([1,1,1]))):
	out_cam = []
	for line in in_cam:
		transformed_line = map(lambda pt: R.dot(pt) + t, line) # rotate then transalte
		out_cam.append(transformed_line)
	return out_cam

def rainbow_colors(n = 8):
	return [colorsys.hsv_to_rgb(float(h)/n, 1.0, 1.0) for h in range(n)]

def draw_all():
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
	cam_t = np.mat([0,0,0])
	cam_R = np.mat(np.diag([1,1,1]))
	cam = move_camera(get_default_camera(), cam_t, cam_R)
	draw_camera(ax1, cam)
	for offset in [1,3,5,9]:
		cam_t2 = np.mat([offset, 0, 0])
		cam2 = move_camera(get_default_camera(), cam_t2, cam_R)
		draw_camera(ax1, cam2)

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

##########
# tests #
##########

def test_raibow_colors():
	r8 = rainbow_colors(8)
	assert(len(r8) == 8)
	assert(all([len(c) == 3 for c in r8]))


