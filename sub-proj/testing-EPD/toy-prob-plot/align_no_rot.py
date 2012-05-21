#!/usr/bin/python
import math, sys, getopt
from random import random
import matplotlib.pyplot as plt

def illustrate_points(points1, points2):
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
	plt.pause(0.1)

def illustrate_forces(force_list):
	"force_list is like [(f_x1, f_y1), (f_x2, f_y2), ... ]"
	x_force = [x for (x,y) in force_list]
	y_force = [y for (x,y) in force_list]
	plt.figure(1)
	plt.clf()
	plt.subplot(2, 1, 1)
	plt.bar(range(len(x_force)), x_force)
	plt.subplot(2, 1, 2)
	plt.bar(range(len(y_force)), y_force)
	plt.show()
	plt.pause(0.05)

class PointAligner(object):


	def __init__(self):
		pass

	@staticmethod
	def align(points1, points2, iterations, illustrate_on = False, verbose_on = True):
		if(verbose_on):
			print "points1 = " + `points1` + " and points2 = " + `points2` + ".\n"\
					"aligning with " + `iterations` + "\iterations."

		for i in range(iterations):
			# Draw current state of alignment
			if (illustrate_on):
				illustrate_points(points1 = points1, points2 = points2)
			if (verbose_on):
				print "iteration = " + `i+1`

			alf = PointAligner.aggregateForce(points1, points2, illustrate_on = False)
			unitAlf = PointAligner.unitVector(alf)
			if(verbose_on):
				print "unitalf  = " + `unitAlf`
			#print "alf = " + `alf`
			for j in range(len(points2)):
				points2[j] = [points2[j][0]+unitAlf[0]*.1, points2[j][1]+unitAlf[1]*.1]
			
			if(verbose_on):
				print "------------------------------------------------------\npoints2 = " + `points2`
				print





	#distance between two points
	@staticmethod
	def distance(x1, y1, x2, y2):
		return math.sqrt(pow(x1 - x2,2) + pow(y1 - y2,2))

	#linear force exterted between two points
	#current model is 1/sqrt
	@staticmethod
	def force_mag(d):
		if d <= .10:
			return (10 * d) * 0.316 # linear drop off within .1
		return 1/pow(d, .5)

	@staticmethod
	def length(x, y):
		return PointAligner.distance(x, y, 0, 0)

	#takes a vector in and returns the unit vector in the same direction
	@staticmethod
	def unitVector(v):
		d = PointAligner.distance(0,0,v[0],v[1])
		if d==0:
			return [0,0]
		return [v[0]/d, v[1]/d]

	@staticmethod
	def force_on_one_point(x1, y1, points2):
		force_x = 0
		force_y = 0
		for (x2,y2) in points2:
			d = PointAligner.distance(x1, y1, x2, y2)
			f = PointAligner.force_mag(d)
			directionVec = PointAligner.unitVector([x1-x2,y1-y2])
			force_x += directionVec[0]*f
			force_y += directionVec[1]*f
		return (force_x, force_y)

	#takes all linear forces between all sets of points, and adds them
	@staticmethod
	def aggregateForce(points1, points2, illustrate_on = False):
		forces = []
		for (x1,y1) in points1:
			forces.append(PointAligner.force_on_one_point(x1, y1, points2))
		if illustrate_on:
			illustrate_forces(forces)
		aggregate=[0,0]
		for (fx,fy) in forces:
			aggregate[0]+=fx
			aggregate[1]+=fy
		return aggregate


	#returns [COM(points1),COM(points2)]
	#will use as pivot points
	@staticmethod
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
	n_pts = 15;
	for i in range(1, n_pts):
		x = random()
		y = random()
		points1.append([x, y])
		points2.append([x + err_x, y + err_y])

	PointAligner.align(
			points1 = points1, 
			points2 = points2,
			iterations = 20,
			illustrate_on = True,
			verbose_on = False)

	print "that took %f seconds" % (time() - time_start)

