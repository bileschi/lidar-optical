#!/usr/bin/python
import math, sys, getopt
from random import random
from associations import *
from simple_illustrations import illustrate_points, illustrate_forces, illustrate_assoc, illustrate_jacobian
import numpy as np
import copy
import pdb

"""
align_rot_4_param:

Four free parameters, one for each element in 2x2 projection matrix.  

True parameter settings are a rotation matrix.
"""


def estimate_projection_params(
	img_pts,
	space_pts,
	projection_fcn = None,
	jacobian_fcn = None,
	associate_fcn = None,
	guess_params = None,
	iterations = 20,
	illustrate = set(),
	verbose_on = True):
	""" Determines a set of params to the projection function
	which align the projection of the space_pts onto the img_points, 
	according to some metric.  Caller must provide a guess as to the projection_params.
	Method returns optimized params.
	"""
	proj_params = [v for (k, v) in sorted(guess_params.items())]
	param_names = [k for (k, v) in sorted(guess_params.items())]
	n_params = len(param_names)
	n_pts2 = len(space_pts)
	for i in range(iterations):
		# project the space points into the image domain
		proj_pts = [projection_fcn(s_pt, proj_params) for s_pt in space_pts]
		# Draw current state of alignment
		if ('projection' in illustrate):
			illustrate_points(img_pts = img_pts, proj_pts = proj_pts)
		if (verbose_on):
			print "iteration = " + `i+1`
		# Associate projected space-points to nearby image-points
		# creating associated pairs.  Each pair infers an offset in image space.
		# This offset infers a confidence in the association and a parameter
		# manipulation to align the points.
		assocs = associate_fcn(img_pts = img_pts, proj_pts = proj_pts)
		if ('association' in illustrate):
			illustrate_assoc(img_pts = img_pts, proj_pts = proj_pts, assocs = assocs)
		if ('jacobian' in illustrate):
			for i_param in range(0, n_params):
				illustrate_jacobian(
					space_pts = [(x/4.0,y/4.0) for x in range(-12, 12, 3) for y in range(-12,12,3)],
					proj_fcn = projection_fcn,
					jacobian_fcn = jacobian_fcn, 
					current_params = proj_params,
					figure_idx = i_param + 10,
					clear_figure_first = True,
					param_idx = i_param)
		# From each association, compute a suggested parameter update.
		param_update_list = assocs_to_updates(assocs = assocs, 
			img_pts = img_pts,
			proj_pts = proj_pts,
			space_pts = space_pts,
			move_portion = 0.5,
			n_params = n_params,
			jacobian_fcn = jacobian_fcn,
			current_params = proj_params)
		new_proj_params = apply_param_update(
			proj_params, param_update_list, assocs['confidences'])
		proj_params = copy.deepcopy(new_proj_params)
	est_params = {}
	for i, k in enumerate(param_names):
		est_params[k] = proj_params[i]
	return est_params

def apply_param_update(proj_params, param_updates, confidences):
	"""returns a new estimate of proj_params by incorporating a concensus update
	from the param_update_list.  Uses (confidence) variable in assocs as input
	in aggregating the concensus"""
	n_params = np.size(param_updates, 1)
	weighted_updates = param_updates * (np.tile(confidences, (n_params, 1)).T)
	mean_update = param_updates.mean(axis=0)
	return proj_params + mean_update

def assocs_to_updates(assocs, img_pts, proj_pts, space_pts, 
	move_portion = 0.1, 
	grace_dist = .010,
	n_params = None, 
	jacobian_fcn = None,
	current_params = None):
	""" Input associations between img_pts and projected space points. Return
	for each asociation an additive parameter update which will map the projection
	of the space_point a fraction (move_portion) closer to the img_pt.

	returns updates as a np.array(n_assocs, n_params)
	"""
	n_assocs = len(assocs['pair_indicies'])
	param_updates = np.zeros((n_assocs, n_params))
	for [assoc_idx, (img_idx, space_idx)] in enumerate(assocs['pair_indicies']):
		proj_pt = proj_pts[space_idx]
		space_pt = space_pts[space_idx]
		img_pt = img_pts[img_idx]
		offset = assocs['offsets'][assoc_idx]
		dist = assocs['distances'][assoc_idx]
		conf = assocs['confidences'][assoc_idx]
		if dist > grace_dist: # no update if projection is close enough.
			# derivative of projection WRT params
			J = jacobian_fcn(space_pt, current_params = current_params)
			Jinv = J.I # deriv of params WRT projection
			param_updates[assoc_idx, :] = -move_portion * Jinv.dot(offset)
	return param_updates

def gen_rand_space_pts(n_pts = 5):
	" picks n_pts randomly in 2d. unit box.  pts returned as list-of-tuples"
	space_pts = []
	for i in range(0, n_pts):
		x = 2 * random() - 1
		y = 2 * random() - 1
		space_pts.append((x, y))

	return space_pts

def img_pt_subtract(pt1, pt2):
	return (pt1[0] - pt2[0], pt1[1] - pt2[1] )

def proj_rotate_jacobian(space_pt = None, current_params = None):
	"""
	Let i = f(s,p) where f is a function from space S to space I.
	f is the function calculating the projection of s, parameterized by p. 

	proj_rotate_jacobian calculates the partial derivatives of i with respect to 
	p.  In general, the derivative may depend on the space_point s.

	
	Input:
	space_pt: the tuple representing the point in space to project
	current_params: p.  a dict or sequence indicating the current param settings p

	Output:
	Output array[m,n] is partial derivative of i[m] W.R.T. p[n]
	"""
	(sx, sy) = space_pt
	m00 = current_params[0]
	m01 = current_params[1]
	m10 = current_params[2]
	m11 = current_params[3]
	return -np.mat([[ sx, sy, 0, 0],
				   [  0, 0, sx, sy]], np.double);

def project_rotate(space_pt, proj_params = {}):
	"""applies proj_params as a matrix multiply """
	(space_x, space_y) = space_pt
	M = np.zeros((2,2))
	try:
		M[0,0] = proj_params['mat_00']
		M[0,1] = proj_params['mat_01']
		M[1,0] = proj_params['mat_10']
		M[1,1] = proj_params['mat_11']
	except TypeError:
		M[0,0] = proj_params[0]
		M[0,1] = proj_params[1]
		M[1,0] = proj_params[2]
		M[1,1] = proj_params[3]
	except ValueError:
		M[0,0] = proj_params[0]
		M[0,1] = proj_params[1]
		M[1,0] = proj_params[2]
		M[1,1] = proj_params[3]
	i = M.dot(np.array(space_pt))
	img_pt = (i[0], i[1])
	return img_pt

##########
## MAIN ##
##########

if __name__ == "__main__":
	from time import time

	# Simulation parameters
	time_start = time()
	noise_std = .0025;
	n_pts = 25;
	real_rot = random()
	tru_proj_params = {
	  'mat_00': np.cos(real_rot),
	  'mat_01': np.sin(real_rot),
	  'mat_10': -np.sin(real_rot),
	  'mat_11': np.cos(real_rot)
	}
	guess_rot = real_rot + .5 
	guess_params = {}
	guess_params = {
	  'mat_00': np.cos(guess_rot),
	  'mat_01': np.sin(guess_rot),
	  'mat_10': -np.sin(guess_rot),
	  'mat_11': np.cos(guess_rot)
	}

	# print to console the problem to solve.
	print "guess params = [%.2f, %.2f, %.2f, %.2f]" % \
	   (guess_params['mat_00'], guess_params['mat_01'],
	    guess_params['mat_10'], guess_params['mat_11'])
	print "true  params = [%.2f, %.2f, %.2f, %.2f]" % \
	   (tru_proj_params['mat_00'], tru_proj_params['mat_01'],
	    tru_proj_params['mat_10'], tru_proj_params['mat_11'])

	# generate random space points
	space_pts = gen_rand_space_pts(n_pts)

	# img_pts = true projection of space_pts + some noise
	img_pts = []
	for space_pt in space_pts:
		(img_x, img_y) = project_rotate(space_pt, tru_proj_params)
		noise_x = random() * noise_std
		noise_y = random() * noise_std
		img_pts.append((img_x + noise_x, img_y + noise_y))

	# Build KDtree (if using approximate knn as association)
	kdtree = spatial.KDTree(img_pts)

	# perform optimization procedure to estimate true projection params
	est_params = estimate_projection_params(
			img_pts = img_pts, 
			space_pts = space_pts,
			projection_fcn = project_rotate,
			jacobian_fcn = proj_rotate_jacobian,
			# associate_fcn = associate_points_all_to_all,
			# associate_fcn = associate_points_all_to_nearest,
			# associate_fcn = associate_points_cheating,
			associate_fcn = lambda img_pts, proj_pts: associate_knn(proj_pts=proj_pts, kdtree=kdtree, k=2, eps=np.inf),
			guess_params = guess_params,
			iterations = 25,
			# valid illustrate includes 'projection', 'association', 'jacobian'
			illustrate = set(),
			# illustrate = set(['projection', 'association']),
			# illustrate = set(['projection', 'association', 'jacobian']),
			verbose_on = True)

	# print results to console
	print "guess params = [%.2f, %.2f, %.2f, %.2f]" % \
	   (guess_params['mat_00'], guess_params['mat_01'],
	    guess_params['mat_10'], guess_params['mat_11'])
	print "true  params = [%.2f, %.2f, %.2f, %.2f]" % \
	   (tru_proj_params['mat_00'], tru_proj_params['mat_01'],
	    tru_proj_params['mat_10'], tru_proj_params['mat_11'])
	print "est   params = [%.2f, %.2f, %.2f, %.2f]" % \
	   (est_params['mat_00'], est_params['mat_01'],
	    est_params['mat_10'], est_params['mat_11'])
	g = np.array([guess_params[k] for k in guess_params.keys()])
	t = np.array([tru_proj_params[k] for k in guess_params.keys()])
	e = np.array([est_params[k] for k in guess_params.keys()])
	g_err = np.linalg.norm(g - t)
	e_err =  np.linalg.norm(e - t)
	recov = 1 - (e_err / g_err)
	print "error recovery = %d%% (%.2f => %.2f)" % (round(recov*100), g_err, e_err)
	print "that took %f seconds" % (time() - time_start)