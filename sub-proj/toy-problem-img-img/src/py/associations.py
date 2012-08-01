#!/usr/bin/python
from scipy import spatial
import numpy as np
import pdb

def img_pt_subtract(pt1, pt2):
	return (pt1[0] - pt2[0], pt1[1] - pt2[1] )


# Associate space points to nearby image points.
def all_to_all(img_pts = [], proj_pts = []):
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
def all_to_nearest(img_pts = [], proj_pts = []):
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
def cheating(img_pts = [], proj_pts = []):
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



def knn(proj_pts = [], kdtree = None, k = 1, eps = 0):
	""" Returns a datastructure (dict) containing:
	pair_indicies : a list of pairs of associated points as (img_p, proj_p)
	offsets: a list of vectors img_p - proj_p 
	distances: the L2 magnitude of the offsets
	confidences: a weight to be placed on this association

	Within each list, the item at index i refers to the same association.
	associates each projected space point to the k nearest img points (approximately)

	img_pts information is stored in the kdtree
	"""
	(dists_mat, tree_idx_mat) = kdtree.query(np.array([proj_pts]), k=k, eps=eps)
	n_assocs = np.size(dists_mat)
	pair_indicies = []
	offsets = []
	dists = []
	img_pts = kdtree.data
	for i_proj, proj_pt in enumerate(proj_pts):
		for i_k in range(0, k):
			# pdb.set_trace()
			if (k == 1):
				i_img = tree_idx_mat[0, i_proj]
				dists.append(dists_mat[0, i_proj])
			else:
				i_img = tree_idx_mat[0, i_proj, i_k]
				dists.append(dists_mat[0, i_proj, i_k])
			img_pt = img_pts[i_img]
			print "img_pt = " + str(img_pt)
			print "proj_pt = " + str(img_pt)
			pdb.set_trace()
			offset = img_pt_subtract(img_pt, proj_pt)
			offsets.append(offset)
			pair_indicies.append((i_img, i_proj))
	a = {}
	a['pair_indicies'] = pair_indicies
	a['offsets'] = offsets
	a['distances'] = dists
	a['confidences'] = np.ones(len(dists))
	return a

