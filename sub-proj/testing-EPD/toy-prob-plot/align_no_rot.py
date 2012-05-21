#!/usr/bin/python
import math, sys, getopt
from random import random
import matplotlib.pyplot as plt

def illustrate_points(points1, points2):
	"""scatter plot of two point sets"""
	plt.figure(0)
	plt.clf()
	plt.plot([x for [x,y] in points1], [y for [x,y] in points1], 'ro')
	plt.plot([x for [x,y] in points2], [y for [x,y] in points2], 'bx')
	plt.axis([-7, 7, -5, 5])
	plt.grid(True)
	plt.title('this is a title')
	plt.xlabel('this is x label')
	plt.ylabel('this is y label')
	plt.show()
	plt.pause(0.01)

def illustrate_forces(force_list):
	"""draws bar plot of force in x and y direction for a list of forces
	force_list is like [(f_x1, f_y1), (f_x2, f_y2), ... ]"""
	x_force = [x for (x,y) in force_list]
	y_force = [y for (x,y) in force_list]
	plt.figure(1)
	plt.clf()
	plt.subplot(2, 1, 1)
	plt.bar(range(len(x_force)), x_force)
	plt.ylim([-20, 20])
	plt.title('x direction')
	plt.subplot(2, 1, 2)
	plt.bar(range(len(y_force)), y_force)
	plt.ylim([-20, 20])
	plt.title('y direction')
	plt.show()
	plt.pause(0.01)

def align(points1, points2, iterations, illustrate_on = False, verbose_on = True):
	"""translate points2 to lie atop points2.  Return the learned transform"""
	n_pts2 = len(points2)
	if(verbose_on):
		print "points1 = " + `points1` + " and points2 = " + `points2` + ".\n"\
				"aligning with " + `iterations` + "\iterations."
	accumulated_offset = [0, 0]
	for i in range(iterations):
		# Draw current state of alignment
		if (illustrate_on):
			illustrate_points(points1 = points1, points2 = points2)
		if (verbose_on):
			print "iteration = " + `i+1`
		alf = aggregateForce(points1, points2, illustrate_on = illustrate_on)
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
def aggregateForce(points1, points2, illustrate_on = False):
	forces = []
	for (x1,y1) in points1:
		forces.append(force_on_one_point(x1, y1, points2))
	if illustrate_on:
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
	err_x = random() * 1
	err_y = random() * 1
	points1 = []
	points2 = []
	n_pts = 5;
	for i in range(1, n_pts):
		x = random()
		y = random()
		points1.append([x, y])
		points2.append([x + err_x, y + err_y])

	est_offset = align(  points1 = points1, 
			points2 = points2,
			iterations = 20,
			illustrate_on = False,
			verbose_on = False)

	print "that took %f seconds" % (time() - time_start)
	print "true offset = (%f, %f)" % (err_x, err_y)
	print "estimated offset = (%f, %f)" % (est_offset[0], est_offset[1])


