#stdlib
import numpy as np
from numpy.testing import assert_allclose
from random import random
from math import sqrt, pi
#local
from copy import copy
from projection_3d_pt import *
from rotations_3d import rot_mat

def default_space_points():
	space_points = [(0, 0, 10)]
	space_points.append((1, 0, 10))
	space_points.append((0, 1, 10))
	return space_points

def rand_sample_cube(n_pts=1, x_min=0, x_max=1, y_min=0, y_max=1, z_min=0, z_max=1):
	" picks n_pts randomly in 3d. unit box.  pts returned as list-of-xyz-tuples"
	space_pts = []
	for i in range(0, n_pts):
		x = random() * float(x_max - x_min) + x_min
		y = random() * float(y_max - y_min) + y_min
		z = random() * float(z_max - z_min) + z_min
		space_pts.append((x, y, z))
	return space_pts

def rand_point_sphere(n_pts=50):
	space_pts = []
	for i in range(0, n_pts):
		(x, y, z) = (random()-.5, random()-.5, random()-.5)
		r = sqrt(x*x + y*y + z*z)
		space_pts.append((x/r, y/r, z/r))
	return space_pts

def assert_allclose_or_bothnone(a, b, rtol=1e-8):
	if a is not None or b is not None:
		assert_allclose(a, b, rtol)

###########
## Tests ##
###########

def test_rand_point_sphere():
	for (x,y,z) in rand_point_sphere(10):
		assert_allclose([1], (x*x + y*y + z*z))
	[(x1,y1,z1), (x2,y2,z2)] = rand_point_sphere(2)
	assert(x1 != x2)
	assert(y1 != y2)
	assert(z1 != z2)

def test_translation_polarity():
	"""In twin scenarios, randomly translate the camera in one, and the points
	in the other.  The images from both cameras should be the same.
	"""
	# pick a random offset from the unit cube
	rand_offset = rand_sample_cube()
	[(t_x, t_y, t_z)] = rand_offset
	(r_x, r_y, r_z) = (.05, .05, .05) # small rot to check commutativity
	# scenaro 1, translate the camera.	
	space_pts_1 = [(0, 0, 1), (0, 0, 10), (0, 0, 100)]
	cam_1 = Camera() # creates a default camera.
	cam_1 = cam_1.mutate(t_x=t_x, t_y=t_y, t_z=t_z)
	cam_1 = cam_1.mutate(r_x=r_x, r_y=r_y, r_z=r_z)
	# scenareo 2, translate the points in the opposite dir.
	space_pts_2 = [(x-t_x, y-t_y, z-t_z) for (x,y,z) in space_pts_1]
	cam_2 = Camera() # creates a default camera.
	cam_2 = cam_2.mutate(r_x=r_x, r_y=r_y, r_z=r_z)
	# the images should be the same.
	img_pts_1 = project_3d(space_pts_1, cam_1)
	img_pts_2 = project_3d(space_pts_2, cam_2)
	for i1, i2 in zip(img_pts_1, img_pts_2):
		assert_allclose_or_bothnone(i1, i2, rtol=1e-8)


def test_rotation_polarity():
	"""In twin scenarios, randomly rotate the camera in one, and the points
	in the other.  The images from both cameras should be the same.
	"""
	# pick a random rotation
	rand_rotation = rand_sample_cube(n_pts=1, 
		x_min = -pi, x_max=pi,
		y_min = -pi, y_max=pi,
		z_min = -pi, z_max=pi)
	[(r_x, r_y, r_z)] = rand_rotation
	R = rot_mat(r_x, r_y, r_z)
	cam_1 = Camera().mutate(k_mult = .1)
	cam_2 = Camera().mutate(k_mult = .1)
	space_pts_1 = rand_point_sphere(10)
	# scenario 1, rotate the camera.	
	cam_1 = cam_1.mutate(r_x=r_x, r_y=r_y, r_z=r_z)
	# scenario 2, rotate the points in the opposite dir.
	space_pts_2 = [R.I.dot(sp_1) for sp_1 in space_pts_1]
	# the images should be the same.
	img_pts_1 = project_3d(space_pts_1, cam_1)
	img_pts_2 = project_3d(space_pts_2, cam_2)
	for i1, i2 in zip(img_pts_1, img_pts_2):
		assert_allclose_or_bothnone(i1, i2, rtol=1e-8)


def test_project():
	"test projection vs. pen and paper"
	space_pts = default_space_points()
	cam = Camera()
	img_pts = project_3d(space_pts, cam)
	expected_results = []
	expected_results.append(np.mat([320.0, 240.0]))
	expected_results.append(np.mat([420.0, 240.0]))
	expected_results.append(np.mat([320.0, 340.0]))
	for (i, e) in zip(img_pts, expected_results):
		assert_allclose(i, e, rtol=1e-2)

def test_homg():
	"tests for geometric homoginization of sequence and mat types"
	x1 = [1.0, 2, 3]
	x2 = np.mat([1.0, 2, 3])
	for x in [x1, x2]:
		x_h = to_homg(x)
		assert(x_h.shape == (1, 4))
		assert(np.all(x_h == np.mat([1.0, 2, 3, 1])))

def test_unhomg():
	x_h1 = np.mat([1.0, 2, 3, 1])
	x_h2 = np.mat([2.0, 4, 6, 2])
	x_h3 = np.mat([-1.0, -2, -3, -1])
	for x_h in [x_h1, x_h2, x_h3]:
		x = to_unhomg(x_h)
		assert(x.shape == (1, 3))
		assert(np.all(x == np.mat([1.0, 2, 3])))

def test_cam_as_list():
	l = Camera().as_list()
	assert_allclose(l, [1000.0, 0.0, 320.0, 0.0, 0.0, 1000.0, 240.0, 0.0, 0.0, 0.0, 1.0, 0.0])

def test_cam_to_dict():
	l = Camera().as_list()
	d = Camera().as_dict()
	for (li, (k, di)) in zip(l, sorted(d.items())):
		assert_allclose(li, di)

def test_jacobian():
	"""each element of jacobian J should represent how a small change in params 
	moves the image.  To test, we modify the projection matrix element (12 of them)
	to see if projecting by the mutation moves the image in the direction / amount
	predicted by the Jacobian. """
	# return None
	delta = 0.001
	cam = Camera()
	P_as_list = cam.as_list()
	space_pt = (0, 0, 10)
	J = proj_3d_jacobian(space_pt, P_as_list)
	img_center = project_3d([space_pt], P_as_list)[0]
	print img_center
	for i_param in range(0,12):
		print "Param: %d" % i_param
		P_mod = copy(P_as_list)
		print "  %f -> %f" % (P_mod[i_param], P_mod[i_param] + delta)
		P_mod[i_param] = P_mod[i_param] + delta
		img_mod = project_3d([space_pt], cam=P_mod)[0]
		print "img_mod: (%.3f,%.3f)" % (img_mod[0], img_mod[1])
		print "img_center: (%.3f,%.3f)" % (img_center[0], img_center[1])
		dx = img_mod[0] - img_center[0]
		dy = img_mod[1] - img_center[1]
		dx_predict = delta * J.item((0, i_param))
		dy_predict = delta * J.item((1, i_param))
		print "prediction: dx=%f, dy=%f" % (dx_predict, dy_predict)
		print "    actual: dx=%f, dy=%f\n" % (dx, dy)
		assert_allclose(dx, dx_predict, rtol=1e-2)
		assert_allclose(dy, dy_predict, rtol=1e-2)











