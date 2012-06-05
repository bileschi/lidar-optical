#!/usr/bin/python
import math, sys, getopt
from random import random, uniform
#from simple_illustrations import illustrate_points, illustrate_forces, illustrate_assoc, illustrate_jacobian
import numpy as np
import copy
import pdb

"""
param_sweep.py:

based on align_rot_4_param.py

adds spurious points, and tests for robustness against these points.

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

#@CONSDIER http://projects.scipy.org/scikits/wiki/AnnWrapper
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

# Associate space points to image points at the same index.
def associate_points_cheating(img_pts = [], proj_pts = []):
	""" Returns a datastructure (dict) containing:
	pair_indicies : a list of pairs of associated points as (img_p, proj_p)
	offsets: a list of vectors img_p - proj_p 
	distances: the L2 magnitude of the offsets
	confidences: a weight to be placed on this association

	Within each list, the item at index i refers to the same association.
	associates each projected space point to the nearest image point.

	This degenerate association gives the exact correct association
	in the case where the img_pts are exactly the true projections of
	the incorrectly projected proj_pts
	"""
	pair_indicies = []
	offsets = []
	distances = []
	confidences = []
	for i_proj in range(0, len(proj_pts)):
		i_img = i_proj
		offset = img_pt_subtract(img_pts[i_img], proj_pts[i_proj])
		dist = np.linalg.norm((offset[0], offset[1]))
		pair_indicies.append((i_img, i_proj))
		offsets.append(offset)
		distances.append(dist)
		confidences.append(1)
	a = {}
	a['pair_indicies'] = pair_indicies
	a['offsets'] = offsets
	a['distances'] = distances
	a['confidences'] = confidences
	return a


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


def estimate_error_recovery(guess_params, est_params, tru_proj_params):
	"""
		@param guess_params 	are the intial placement in the space,
		@param est_params 		are the estimated parameters of alignment
		@param tru_proj_params 	are the actual offset that the projection 
							points from the img points
		@return [percent recovered, norm(guess_params), norm(est_params)]
	"""
	g = np.array([guess_params[k] for k in guess_params.keys()])
	t = np.array([tru_proj_params[k] for k in guess_params.keys()])
	e = np.array([est_params[k] for k in guess_params.keys()])
	g_err = np.linalg.norm(g - t)
	e_err =  np.linalg.norm(e - t)
	percnt = round((1 - (e_err / g_err))*100, 2)
	return [percnt, g_err, e_err]

def proceed_until(	img_pts, 
				space_pts, 
				tru_proj_params,
				projection_fcn = None,
				jacobian_fcn = None,
				associate_fcn = None,
				guess_params = None,
				illustrate = set(),
				verbose_on = True,
				max_itr = 20, 
				perc_recov_required = 90):

	"""
		Will run estimate_projection_params with 1 iteration over and over
		until either @param max_itr is hit, or @param perc_recov_required 
		is reached

		@return : upon success: ["finished", iterations completed]
		@return : upon failure: ["not finished", percent recovered]
	"""

	steps_done = 0
	perc_recov = 0
	max_perc_recov = -np.inf
	g_err = np.inf
	e_err = np.inf
	cur_est_params = guess_params
	while (steps_done < max_itr and perc_recov < perc_recov_required):
		cur_est_params = estimate_projection_params(	img_pts = img_pts,
											space_pts = space_pts,
											projection_fcn = projection_fcn,
											jacobian_fcn = jacobian_fcn,
											associate_fcn = associate_fcn,
											guess_params = cur_est_params,
											iterations = 1,
											illustrate = illustrate,
											verbose_on = False
											)

		[perc_recov, g_err, e_err]  = estimate_error_recovery(guess_params, cur_est_params, 
								tru_proj_params)
		if perc_recov > max_perc_recov:
			max_perc_recov = perc_recov
		steps_done+=1
		if verbose_on:
			print `steps_done` + " of " + `max_itr`
		
	if steps_done == max_itr:
		return ["not finished", perc_recov]
	elif perc_recov >= perc_recov_required:
		return ["finished", steps_done]
	else:
		return ["an error occurred"]

def run_simulation(max_itr = 5, n_pts = 10, frac_spurious = 0, noise_std = 0, req_recovery = 80):
	"""
		basically is the old main ported to over here, so i can run 
		many simulations in one go

		will create a new problem, and then run the method "proceed_until" on that 
		problem

	"""
	# Simulation parameters
	time_start = time()
	real_rot = random()
	tru_proj_params = {
	  'mat_00': np.cos(real_rot),
	  'mat_01': np.sin(real_rot),
	  'mat_10': -np.sin(real_rot),
	  'mat_11': np.cos(real_rot)
	}
	guess_rot = real_rot + .5 #@SB why?
	guess_params = {}
	guess_params = {
	  'mat_00': np.cos(guess_rot),
	  'mat_01': np.sin(guess_rot),
	  'mat_10': -np.sin(guess_rot),
	  'mat_11': np.cos(guess_rot)
	}


	# generate random space points
	space_pts = gen_rand_space_pts(n_pts)

	# img_pts = true projection of space_pts + some noise
	img_pts = []
	for space_pt in space_pts:
		(img_x, img_y) = project_rotate(space_pt, tru_proj_params)
		noise_x = random() * noise_std
		noise_y = random() * noise_std
		img_pts.append((img_x + noise_x, img_y + noise_y))


	#add spurious points
	n_spurious = int(round(n_pts*frac_spurious))
	spurious_pts = []
	for i in range(n_spurious):
		img_pts.append([uniform(-1,1), uniform(-1,1)])


	trial = proceed_until(	img_pts = img_pts, 
						space_pts = space_pts,
						tru_proj_params = tru_proj_params,
						projection_fcn = project_rotate,
						jacobian_fcn = proj_rotate_jacobian,
						associate_fcn = associate_points_all_to_nearest,
						guess_params = guess_params,
						illustrate = set(),
						verbose_on = False,
						max_itr = max_itr, 
						perc_recov_required = req_recovery);

	trial.append(time() - time_start)
	return trial

##########
## MAIN ##
##########

#the weird file io placements are intentional, so that i can see the file as it updates.
if __name__ == "__main__":
	from time import time

	#erase file
	f = open('./param_sweep.txt', 'w')
	f.close();

	req_recovery = 90
	num_trials = 20
	max_itr = 50 #@CONSIDER larger?	noise_std  = 0;
	n_pts_list = [10,20,50,100,250]
	#n_pts_list = [3]
	frac_spurious_list = [0,.25,.5,1,1.5]

	f = open('./param_sweep.txt', 'a')
	print "required recovery : " + `req_recovery`
	f.write("required recovery : " + `req_recovery`+"\n")
	print "max_itr : " + `max_itr`
	f.write("max_itr : " + `max_itr`+"\n")
	print "noise : " + `noise_std` + "\n"
	f.write("noise : " + `noise_std` + "\n"+"\n")
	f.close()
	
	for n_pts in n_pts_list:
		for frac_spurious in frac_spurious_list:
			f = open('./param_sweep.txt', 'a') # so that i can see progress

			print "num points = "  + `n_pts`
			f.write("num points = "  + `n_pts`+"\n")
			print "frac spurious = " + `frac_spurious`
			f.write("frac spurious = " + `frac_spurious`+"\n")
			print "------------------------------"
			f.write("------------------------------"+"\n")

			failures = []
			successes = []
			for i in range(num_trials):
				current_run = run_simulation(		max_itr = max_itr, 
											n_pts = n_pts, 
											frac_spurious = frac_spurious,
											noise_std = noise_std,
											req_recovery = req_recovery)
				if current_run[0] == 'not finished':
					failures.append(current_run)
				else:
					successes.append(current_run)
			print "% successful : " +  `((len(successes)+ 0.0)/num_trials)`
			f.write("% successful : " +  `((len(successes)+ 0.0)/num_trials)`+"\n")

			if len(successes) != 0:
				total_time = 0
				total_itr = 0
				for sum_idx in range(len(successes)):
					total_time += successes[sum_idx][2]
					total_itr += successes[sum_idx][1]
				print "avg iterations for successful: " + `(total_itr/len(successes))`
				f.write("avg iterations for successful: " + `(total_itr/len(successes))`+"\n")
				print "avg time for successful: " + `(total_time/len(successes))` 
				f.write("avg time for successful: " + `(total_time/len(successes))`+"\n\n")
			else:
				print "no successes"
				f.write("no successes"+"\n")

			print

			if len(failures) != 0:
				total_recovered = 0
				total_time = 0
				for sum_idx in range(len(failures)):
					total_time += failures[sum_idx][2]
					total_recovered += failures[sum_idx][1]
				print "avg recovery for failures: " + `(total_recovered/len(failures))`
				f.write("avg recovery for failures: " + `(total_recovered/len(failures))`+"\n")
				print "avg time for failures: " + `(total_time/len(failures))` 
				f.write("avg time for failures: " + `(total_time/len(failures))`+"\n")
			else:
				print "no failures"	
				f.write("no failures"+"\n")

			print "\n"	
			f.write("\n\n")

			f.close()


	"""
	print "\n"
	print successes
	print "\n"
	print failures
	"""

