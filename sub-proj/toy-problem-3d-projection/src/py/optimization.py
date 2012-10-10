import pdb
import numpy as np
import copy
from simple_illustrations import illustrate_points, illustrate_forces, illustrate_assoc, illustrate_jacobian, illustrate_toy_problem

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
		proj_pts = projection_fcn(space_pts, proj_params)
		# Draw current state of alignment
		if ('projection' in illustrate):
			illustrate_points(img_pts = img_pts, proj_pts = proj_pts)
		if (verbose_on):
			print "iteration = " + `i+1`
		# Associate projected space-points to nearby image-points
		# creating associated pairs.  Each pair infers an offset in image space.
		# This offset infers a confidence in the association and a parameter
		# manipulation to align the points.

		# these are lists of matricies
		print "img_pts = " + str(img_pts)
		print "proj_pts = " + str(proj_pts)
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
