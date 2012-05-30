#!/usr/bin/python
import math, sys, getopt
from random import random
from simple_illustrations import illustrate_points, illustrate_forces, illustrate_assoc
import numpy as np
import copy
import pdb

"""
align_dup_param:

Matches points in the 2-d plane to translated points in the 2-d plane just like
align_no_rot.py.  I.e., the projection function is a simple translation.  The
twist in this version is that there are 4 translation parameters:
 {offset_x1, offset_x2, offset_y1, offset_y2}.  
Offset_x1 and offset_x2 do exactly the same thing.  How does the optimization
procedure deal with the dimensionality mismatch?  If the projection is off in x,
the jacobian will suggest modifying both parameters in the same way.
"""


def estimate_projection_params(
	img_pts,
	space_pts,
	projection_function = None,
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
		proj_pts = [projection_function(s_pt, proj_params) for s_pt in space_pts]
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
		# From each association, compute a suggested parameter update.
		param_update_list = assocs_to_updates(assocs = assocs, 
			img_pts = img_pts,
			proj_pts = proj_pts,
			space_pts = space_pts,
			move_portion = 0.5,
			n_params = n_params,
			jacobian_fcn = jacobian_fcn)
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
	jacobian_fcn = None):
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
			J = jacobian_fcn(space_pt) # derivative of projection WRT params
			Jinv = J.I # deriv of params WRT projection
			param_updates[assoc_idx, :] = move_portion * Jinv.dot(offset)
	return param_updates

# Associate space points to nearby image points.
def associate_points_all_to_all(img_pts = [], proj_pts = []):
	""" Returns a datastructure (dict) containing:
	pair_indicies : a list of pairs of associated points as (img_p, proj_p)
	offsets: a list of vectors img_p - proj_p 
	distances: the L2 magnitude of the offsets
	confidences: a weight to be placed on this association

	Within each list, the item at index i refers to the same association.
	"""
	pair_indicies = []
	offsets = []
	distances = []
	confidences = []
	for i_img in range(0, len(img_pts)):
		for i_proj in range(0, len(proj_pts)):
			pair_indicies.append((i_img, i_proj))
			offset = img_pt_subtract(img_pts[i_img], proj_pts[i_proj])
			offsets.append(offset)
			distances.append(np.linalg.norm((offset[0], offset[1])))
			confidences.append(1)
	a = {}
	a['pair_indicies'] = pair_indicies
	a['offsets'] = offsets
	a['distances'] = distances
	a['confidences'] = confidences
	return a

# Associate space points to nearby image points.
def associate_points_all_to_nearest(img_pts = [], proj_pts = []):
	""" Returns a datastructure (dict) containing:
	pair_indicies : a list of pairs of associated points as (img_p, proj_p)
	offsets: a list of vectors img_p - proj_p 
	distances: the L2 magnitude of the offsets
	confidences: a weight to be placed on this association

	Within each list, the item at index i refers to the same association.
	associates each projected space point to the nearest image point.
	"""
	pair_indicies = []
	offsets = []
	distances = []
	confidences = []
	for i_proj in range(0, len(proj_pts)):
		min_dist = np.inf
		for i_img in range(0, len(img_pts)):
			offset = img_pt_subtract(img_pts[i_img], proj_pts[i_proj])
			dist = np.linalg.norm((offset[0], offset[1]))
			if dist < min_dist:
				best_pair = (i_img, i_proj)
				best_offset = offset
				best_dist = dist
				best_conf = 1
				min_dist = dist
		pair_indicies.append(best_pair)
		offsets.append(best_offset)
		distances.append(best_dist)
		confidences.append(best_conf)
	a = {}
	a['pair_indicies'] = pair_indicies
	a['offsets'] = offsets
	a['distances'] = distances
	a['confidences'] = confidences
	return a

def gen_rand_space_pts(n_pts = 5):
	" picks n_pts randomly in 2d. unit box.  pts returned as list-of-tuples"
	space_pts = []
	for i in range(1, n_pts):
		x = random() - .5
		y = random() - .5
		space_pts.append((x, y))
	return space_pts

def img_pt_subtract(pt1, pt2):
	return (pt1[0] - pt2[0], pt1[1] - pt2[1] )

def proj_jacobian(space_pt = None):
	"""array element i,j is partial derivative projection dim j WRT
	param i"""
	(space_x, space_y) = space_pt
	if space_x > 0:
		y2_sign = -1
	else:
		y2_sign = 1
	if space_y > 0:
		x2_sign = -1
	else:
		x2_sign = 1
	return np.mat([[1.0, x2_sign * 1.0, 0.0, 0.0],
				   [0.0, 0.0, 1.0, y2_sign * 1.0]], np.double);

def project_translate(space_pt, proj_params = {}):
	"""forward projection transform, parameterized 
	offset_x1, offset_x2, offset_y1, offset_y2
	proj_params can be a dict or any sequence. """
	(space_x, space_y) = space_pt
	try:
		o_x1 = proj_params['offset_x1']
		o_x2 = proj_params['offset_x2']
		o_y1 = proj_params['offset_y1']
		o_y2 = proj_params['offset_y2']
	except TypeError:
		o_x1 = proj_params[0]
		o_x2 = proj_params[1]
		o_y1 = proj_params[2]
		o_y2 = proj_params[3]
	except ValueError:
		o_x1 = proj_params[0]
		o_x2 = proj_params[1]
		o_y1 = proj_params[2]
		o_y2 = proj_params[3]
	if (space_x > 0):
		o_y2 = -o_y2
	if (space_y > 0):
		o_x2 = -o_x2
	img_x = space_x + o_x1 + o_x2
	img_y = space_y + o_y1 + o_y2
	img_pt = (img_x, img_y)
	return img_pt

##########
## MAIN ##
##########

if __name__ == "__main__":
	from time import time

	# Simulation parameters
	time_start = time()
	noise_std = .025;
	n_pts = 25;
	tru_proj_params = {
	  'offset_x1': random() - .5,
	  'offset_x2': random() - .5,
	  'offset_y1': random() - .5,
	  'offset_y2': random() - .5}
	guess_params = {}
	for k in tru_proj_params.keys():
		guess_params[k] = tru_proj_params[k] + random() * .5
	# print to console the problem to solve.
	print "guess offset = [%.2f, %.2f, %.2f, %.2f]" % (
		guess_params['offset_x1'], guess_params['offset_x2'],
	    guess_params['offset_y1'], guess_params['offset_y2'])
	print "true offset = [%.2f, %.2f, %.2f, %.2f]" % (
		tru_proj_params['offset_x1'], tru_proj_params['offset_x2'],
	    tru_proj_params['offset_y1'], tru_proj_params['offset_y2'])

	# generate random space points
	space_pts = gen_rand_space_pts(n_pts)

	# img_pts = true projection of space_pts + some noise
	img_pts = []
	for space_pt in space_pts:
		(img_x, img_y) = project_translate(space_pt, tru_proj_params)
		noise_x = random() * noise_std
		noise_y = random() * noise_std
		img_pts.append((img_x + noise_x, img_y + noise_y))

	# perform optimization procedure to estimate true projection params
	est_params = estimate_projection_params(
			img_pts = img_pts, 
			space_pts = space_pts,
			projection_function = project_translate,
			jacobian_fcn = proj_jacobian,
			# associate_fcn = associate_points_all_to_all,
			associate_fcn = associate_points_all_to_nearest,
			guess_params = guess_params,
			iterations = 40,
			# valid illustrate includes 'projection', 'association'
			# illustrate = set(),
			illustrate = set(['projection', 'association']),
			verbose_on = False)

	# print results to console
	print "est  offset = [%.2f, %.2f, %.2f, %.2f]" % (
		est_params['offset_x1'], est_params['offset_x2'],
	    est_params['offset_y1'], est_params['offset_y2'])
	g = np.array([guess_params[k] for k in guess_params.keys()])
	t = np.array([tru_proj_params[k] for k in guess_params.keys()])
	e = np.array([est_params[k] for k in guess_params.keys()])
	g_err = np.linalg.norm(g - t)
	e_err =  np.linalg.norm(e - t)
	recov = 1 - (e_err / g_err)
	print "error recovery = %d%% (%.2f => %.2f)" % (round(recov*100), g_err, e_err)
	print "that took %f seconds" % (time() - time_start)


