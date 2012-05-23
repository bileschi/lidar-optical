#!/usr/bin/python
import math, sys, getopt
from random import random
from simple_illustrations import illustrate_points, illustrate_forces

def align(
	img_pts,
	space_pts,
	guess_params = {'offset_x': 0, 'offset_y': 0},
	iterations = 20,
	illustrate = set(),
	verbose_on = True):
	""" Determines a set of params to the projection function
	which align the projection of the space_pts onto the img_points, 
	according to some metric.  Caller must provide a guess as to the projection_params.
	Method returns optimized params.
	"""
	proj_params = guess_params;
	n_pts2 = len(space_pts)
	if(verbose_on):
		print "img_pts = " + `img_pts` + " and space_pts = " + `space_pts` + ".\n"\
				"aligning with " + `iterations` + "\iterations."
	accumulated_offset = [0, 0]
	for i in range(iterations):
		# project the space points into the image domain
		proj_pts = [project(s, proj_params) for s in space_pts]
		# Draw current state of alignment
		if ('projection' in illustrate):
			illustrate_points(img_pts = img_pts, proj_pts = proj_pts)
		if (verbose_on):
			print "iteration = " + `i+1`
		alf = aggregateForce(img_pts, proj_pts, illustrate = illustrate)
		proj_params['offset_x'] += .1*alf[0]/n_pts2
		proj_params['offset_y'] += .1*alf[1]/n_pts2
	return proj_params

#distance between two points
def distance(x1, y1, x2, y2):
	return math.sqrt(pow(x1 - x2,2) + pow(y1 - y2,2))

#linear force exterted between two points
#current model is 1/sqrt
def force_mag(d):
	if d <= .10:
		return (10 * d) * 0.316 # linear drop off within .1
	return 1/pow(d, .5)

# vector magnitude
def length(x, y):
	return distance(x, y, 0, 0)

#takes a vector in and returns the unit vector in the same direction
def unitVector(v):
	d = distance(0,0,v[0],v[1])
	if d==0:
		return [0,0]
	return [v[0]/d, v[1]/d]

# force exerted on point [x1, y1] by points in list space_pts
def force_on_one_point(x1, y1, space_pts):
	force_x = 0
	force_y = 0
	for (x2,y2) in space_pts:
		d = distance(x1, y1, x2, y2)
		f = force_mag(d)
		directionVec = unitVector([x1-x2,y1-y2])
		force_x += directionVec[0]*f
		force_y += directionVec[1]*f
	return (force_x, force_y)

#takes all linear forces between all sets of points, and adds them
def aggregateForce(img_pts, space_pts, illustrate = set()):
	forces = []
	for (x1,y1) in img_pts:
		forces.append(force_on_one_point(x1, y1, space_pts))
	if 'forces' in illustrate:
		illustrate_forces(forces)
	aggregate=[0,0]
	for (fx,fy) in forces:
		aggregate[0]+=fx
		aggregate[1]+=fy
	return aggregate

#returns [COM(img_pts),COM(space_pts)]
#will use as pivot points
def centerOfMass(points):
	xsum = 0.0
	ysum = 0.0
	for [x,y] in points:
		xsum+=x
		ysum+=y
	com = [xsum/len(points), ysum/len(points)]
	return com

def gen_rand_space_pts(n_pts = 5):
	" picks n_pts randomly in unit box"
	space_pts = []
	for i in range(1, n_pts):
		x = random()
		y = random()
		space_pts.append([x, y])
	return space_pts

def project(space_pt, proj_params = {'offset_x': 1, 'offset_y': 2}):
	"""forward projection transform, parameterized offset_x, offset_y """
	(space_x, space_y) = space_pt
	img_x = space_x + proj_params['offset_x']
	img_y = space_y + proj_params['offset_y']
	img_pt = (img_x, img_y)
	return img_pt

if __name__ == "__main__":
	from time import time
	# Simulation parameters
	time_start = time()
	noise_std = .25;
	n_pts = 5;
	tru_proj_params = {'offset_x': 1, 'offset_y': 2}
	guess_params = {'offset_x': 0, 'offset_y': 0}
	# generate random space points
	space_pts = gen_rand_space_pts(n_pts)
	# calculate images of those points
	img_pts = []
	for space_pt in space_pts:
		(img_x, img_y) = project(space_pt, tru_proj_params)
		noise_x = random() * noise_std
		noise_y = random() * noise_std
		img_pts.append((img_x + noise_x, img_y + noise_y))
	# perform optimization procedure
	est_params = align(
			img_pts = img_pts, 
			space_pts = space_pts,
			guess_params = guess_params,
			iterations = 20,
			# valid illustrate includes 'projection', 'forces'
			# illustrate = set(),
			illustrate = set(['projection']),
			verbose_on = False)
	# print results
	print "that took %f seconds" % (time() - time_start)
	print "true offset = (%f, %f)" % (tru_proj_params['offset_x'], tru_proj_params['offset_y'])
	print "estimated offset = (%f, %f)" % (est_params['offset_x'], est_params['offset_y'])


