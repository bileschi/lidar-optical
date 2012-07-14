#stdlib
import numpy as np
from numpy.testing import assert_allclose
from random import random
#local
from projection_3d_pt import project_3d, project_3d_legacy, to_homg, to_unhomg, cam_as_list, proj_mat_from_cam_params, proj_3d_jacobian
from rotations_3d import rot_mat
from copy import copy

def default_camera():
	true_cam_params = {}
	# external params
	true_cam_params['t'] = np.mat([0.0, 0.0, 0.0])
	true_cam_params['R'] = rot_mat(0.0, 0.0, 0.0)
	# internal params
	true_cam_params['k'] = 1000
	true_cam_params['cx'] = 320
	true_cam_params['cy'] = 240
	return true_cam_params

def default_space_points():
	space_points = [(0, 0, 10)]
	space_points.append((1, 0, 10))
	space_points.append((0, 1, 10))
	return space_points

def rand_sample_cube(n_pts=1, x_min=0, x_max=1, y_min=0, y_max=1, z_min=0, z_max=1):
	" picks n_pts randomly in 3d. unit box.  pts returned as list-of-xyz-tuples"
	space_pts = []
	for i in range(0, n_pts):
		x = ((0.0 + i) / n_pts) * (x_max - x_min) + x_min
		y = ((0.0 + i) / n_pts) * (y_max - y_min) + y_min
		z = ((0.0 + i) / n_pts) * (z_max - z_min) + z_min
		space_pts.append((x, y, z))
	return space_pts

###########
## Tests ##
###########

def test_translation_polarity_consistency():
	space_pts = [(0, 0, 10)]
	cam = default_camera()
	cam['t'] = np.mat([0.0, 1, 0])
	img_pts_legacy = project_3d_legacy(space_pts, cam)
	img_pts = project_3d(space_pts, cam)
	for i1, i2 in zip(img_pts, img_pts_legacy):
		assert_allclose(i1, i2, rtol=1e-8)


def test_projection_match_legacy():
	space_pts = rand_sample_cube(n_pts=15, 
		x_min=-10, x_max=10,
		y_min=-10, y_max=10,
		z_min=-10, z_max=10)
	cam = default_camera()
	cam['R'] = rot_mat(random(), random(), random())
	cam['t'] = np.mat([random(), random(), random()])
	cam['k'] = random()
	img_pts_1 = project_3d_legacy(space_pts, cam)
	img_pts_2 = project_3d(space_pts, cam)
	for (i1, i2) in zip(img_pts_1, img_pts_2):
		if (i1 is None) and (i2 is None):
			continue
		if (i1 is None) != (i2 is None):
			assert(False)
		assert_allclose(i1, i2, rtol=1e-8)

def test_project():
	space_pts = default_space_points()
	cam = default_camera()
	img_pts = project_3d(space_pts, cam)
	expected_results = []
	expected_results.append(np.mat([320.0, 240.0]))
	expected_results.append(np.mat([420.0, 240.0]))
	expected_results.append(np.mat([320.0, 340.0]))
	for (i, e) in zip(img_pts, expected_results):

		assert_allclose(i, e, rtol=1e-2)

def test_homg():
	"tests for homoginization of sequence and mat types"
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
	l = cam_as_list(default_camera())
	assert_allclose(l, [1000.0, 0.0, 320.0, 0.0, 0.0, 1000.0, 240.0, 0.0, 0.0, 0.0, 1.0, 0.0])

def test_jacobian():
	"""each element of jacobian J should represent how a small change in params 
	moves the image """
	delta = 0.001
	cam = default_camera()
	P = proj_mat_from_cam_params(cam_params = cam)
	P_as_list = P.reshape((1,12)).tolist()[0]
	space_pt = (0, 0, 10)
	J = proj_3d_jacobian(space_pt, P_as_list)
	img_center = project_3d([space_pt], cam_params=P_as_list)[0]
	print img_center
	for i_param in range(0,12):
		P_mod = copy(P_as_list)
		P_mod[i_param] = P_mod[i_param] + delta
		img_mod = project_3d([space_pt], cam_params=P_mod)[0]
		dx = to_unhomg(img_mod).item(0) - to_unhomg(img_center).item(0)
		dy = to_unhomg(img_mod).item(1) - to_unhomg(img_center).item(1)
		dx_predict = delta* J.item((0, i_param))
		dy_predict = delta* J.item((1, i_param))
		print i_param
		print to_unhomg(img_center)
		print to_unhomg(img_mod)
		print "prediction: dx=%f, dy=%f" % (dx_predict, dy_predict)
		print "    actual: dx=%f, dy=%f\n" % (dx, dy)
		assert_allclose(dx, dx_predict, rtol=1e-2)
		assert_allclose(dy, dy_predict, rtol=1e-2)











