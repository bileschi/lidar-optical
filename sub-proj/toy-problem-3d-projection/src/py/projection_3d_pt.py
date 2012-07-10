#!/usr/bin/python
import numpy as np

"""
Tools to project 3d points into a 2d image.
Data is expected to be in np.mat format.  Some functions allow for array-like input.
"""

def project_3d(space_points, cam_params):
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

def norm(mat):
	v = mat.dot(mat.transpose()).item(0) ** .5

def to_unhomg(p3):
	return np.take(p3, range(0, p3.size-1), mode='wrap') / np.take(p3, [-1], mode='wrap')
	 
def to_homg(pt):
	return np.append(np.mat(pt), np.mat(1), 1)

def project_homg_pt(pt_3d = np.mat([[1], [1], [1], [1]]), 
	projection_mat = np.mat([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0]])):
	""" pt_3d must be a 1d matrix of length 4 representing the homogenous coordinates
	of a point in 3d space."""
	pt_2d = projection_mat * pt_3d
	return pt_2d


###########
## Tests ##
###########

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



def whatever():
	y = to_unhomg(x)
	print y
	assert(y.len() == 3)
	assert(all([abs(x[i] - y[i]) < .01 for i in range(0,3)]))