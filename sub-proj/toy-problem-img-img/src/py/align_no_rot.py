#!/usr/bin/python
import math, sys, getopt
from random import random
from simple_illustrations import illustrate_points, illustrate_forces

def align(points1, points2, iterations, illustrate = set(), verbose_on = True):
	"""translate points2 to lie atop points2.  Return the learned transform"""
	n_pts2 = len(points2)
	if(verbose_on):
		print "points1 = " + `points1` + " and points2 = " + `points2` + ".\n"\
				"aligning with " + `iterations` + "\iterations."
	accumulated_offset = [0, 0]
	for i in range(iterations):
		# Draw current state of alignment
		if ('projection' in illustrate):
			illustrate_points(points1 = points1, points2 = points2)
		if (verbose_on):
			print "iteration = " + `i+1`
		alf = aggregateForce(points1, points2, illustrate = illustrate)
		unitAlf = unitVector(alf)
		if(verbose_on):
			print "unitalf  = " + `unitAlf`
		#print "alf = " + `alf`
		for j in range(len(points2)):
			offset = [alf[0]*.1/n_pts2, alf[1]*.1/n_pts2]
			points2[j] = [points2[j][0]+offset[0], points2[j][1]+offset[1]]
		accumulated_offset[0] += offset[0]
		accumulated_offset[1] += offset[1]
		if(verbose_on):
			print "------------------------------------------------------\npoints2 = " + `points2`
			print
	return accumulated_offset

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

# force exerted on point [x1, y1] by points in list points2
def force_on_one_point(x1, y1, points2):
	force_x = 0
	force_y = 0
	for (x2,y2) in points2:
		d = distance(x1, y1, x2, y2)
		f = force_mag(d)
		directionVec = unitVector([x1-x2,y1-y2])
		force_x += directionVec[0]*f
		force_y += directionVec[1]*f
	return (force_x, force_y)

#takes all linear forces between all sets of points, and adds them
def aggregateForce(points1, points2, illustrate = set()):
	forces = []
	for (x1,y1) in points1:
		forces.append(force_on_one_point(x1, y1, points2))
	if 'forces' in illustrate:
		illustrate_forces(forces)
	aggregate=[0,0]
	for (fx,fy) in forces:
		aggregate[0]+=fx
		aggregate[1]+=fy
	return aggregate

#returns [COM(points1),COM(points2)]
#will use as pivot points
def centerOfMass(points):
	xsum = 0.0
	ysum = 0.0
	for [x,y] in points:
		xsum+=x
		ysum+=y
	com = [xsum/len(points), ysum/len(points)]
	return com


if __name__ == "__main__":
	from time import time
	time_start = time()
	err_x = random() + 1
	err_y = random() + 1
	noise_std = .25;
	points1 = []
	points2 = []
	n_pts = 5;
	for i in range(1, n_pts):
		x = random()
		y = random()
		noise_x = random() * noise_std
		noise_y = random() * noise_std
		points1.append([x, y])
		points2.append([x + err_x + noise_x, y + err_y + noise_y])

	est_offset = align(  points1 = points1, 
			points2 = points2,
			iterations = 20,
			# valid illustrate includes 'projection', 'forces'
			illustrate = set(['projection']),
			verbose_on = False)

	print "that took %f seconds" % (time() - time_start)
	print "true offset = (%f, %f)" % (-err_x, -err_y)
	print "estimated offset = (%f, %f)" % (est_offset[0], est_offset[1])


