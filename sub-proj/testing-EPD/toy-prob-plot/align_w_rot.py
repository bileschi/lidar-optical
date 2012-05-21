#!/usr/bin/python
import math, sys, getopt
from random import random
import matplotlib.pyplot as plt

#say points1 is fixed, points2 translates, and rotates

class PointAligner(object):


	def __init__(self):
		pass

	@staticmethod
	def align(points1, points2, iterations):
		print "points1 = " + `points1` + " and points2 = " + `points2` + ".\n"\
				"aligning with " + `iterations` + "\iterations."

		for i in range(iterations):
			print "iteration = " + `i+1`

			alf = PointAligner.aggregateForce(points1, points2)
			unitAlf = PointAligner.unitVector(alf)
			print "unitalf  = " + `unitAlf`
			#print "alf = " + `alf`
			for j in range(len(points2)):
				points2[j] = [points2[j][0]+unitAlf[0]*.1, points2[j][1]+unitAlf[1]*.1]
			
			com2 = PointAligner.centerOfMass(points2)
			arf = PointAligner.aggregateAngularForce(points1, points2)

			#print "aggregate rotational force = " + `arf`
			for j in range(len(points2)):
				points2[j] = PointAligner.rotate(points2[j][0], points2[j][1], com2[0], com2[1], arf*.1)
				#PointAligner.length(alf[0], alf[1])
			
			print "------------------------------------------------------\npoints2 = " + `points2`
			print

			# Draw current state of alignment
			plt.clf()
			plt.plot([x for [x,y] in points1], [y for [x,y] in points1], 'ro')
			plt.plot([x for [x,y] in points2], [y for [x,y] in points2], 'bx')
			plt.axis([-7, 7, -5, 5])
			plt.grid(True)
			plt.title('this is a title')
			plt.xlabel('this is x label')
			plt.ylabel('this is y label')
			plt.show()
			plt.pause(.1)




	#distance between two points
	@staticmethod
	def distance(x1, y1, x2, y2):
		return math.sqrt(pow(x1 - x2,2) + pow(y1 - y2,2))

	#linear force exterted between two points
	#current model is 1/sqrt
	@staticmethod
	def force(x1, y1, x2, y2):
		r = PointAligner.distance(x1, y1, x2, y2)
		if r <= .10:
			return 0
		return 1/pow(r,.5)

	@staticmethod
	def length(x, y):
		return PointAligner.distance(x, y, 0, 0)

	#takes a vector in and returns the unit vector in the same direction
	@staticmethod
	def unitVector(v):
		d = PointAligner.distance(0,0,v[0],v[1])
		if d==0:
			return [0,0]
		return [v[0]/d,v[1]/d]

	#takes all linear forces between all sets of points, and adds them
	@staticmethod
	def aggregateForce(points1, points2):
		forces = []
		for (x1,y1) in points1:
			for (x2,y2) in points2:
				r = PointAligner.distance(x1,y1,x2,y2)
				f = PointAligner.force(x1,y1,x2,y2)
				directionVec = PointAligner.unitVector([x1-x2,y1-y2])
				forces.append( [ directionVec[0]*f , directionVec[1]*f ] )

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

	#takes all rotational forces between all sets of points, and adds them
	@staticmethod
	def aggregateAngularForce(points1, points2):
		unitVecs = []
		com1 = PointAligner.centerOfMass(points1)
		com2 = PointAligner.centerOfMass(points2)

		#calculate unit vectors in each theta angle
		for (x1,y1) in points1:
			for (x2,y2) in points2:

				if (x1==x2 and y1==y2):
					continue

				r = PointAligner.distance(x1,y1,x2,y2)
				x3 = x1-x2
				y3 = y1-y2
				tanVec = [com2[1]-y2,-(com2[0]-x2)]

				if(tanVec[0]==0 and tanVec[1]==0):
					continue

				theta = math.acos( (x3*tanVec[0]-y3*tanVec[1]) / (PointAligner.length(x3,y3) * 
								PointAligner.length(tanVec[0],tanVec[1])))
				unitVecs.append( [[math.cos(theta), math.sin(theta)],r] )

		#average the unit vectors
		aggregate = [0,0]
		totalWeight = 0
		for [[x,y],r] in unitVecs:
			aggregate[0] += x*r
			aggregate[1] += y*r
			totalWeight += r
		print "agg = " + `aggregate`
		if aggregate==[0,0]:
			return 0;
		return math.acos( (1*aggregate[0]) / ( PointAligner.length(aggregate[0],aggregate[1])))



	#rotate x, y about center c1, c2, angle theta 
	# COUNTERCLOCKWISE
	#found at: http://content.gpwiki.org/index.php/VB:Tutorials:Rotating_A_Point_In_2D
	@staticmethod
	def rotate(x, y, c1, c2, theta):
		r = PointAligner.distance(x,y,c1,c2)
		return [c1 + ( math.cos(theta) * (x - c1) - math.sin(theta) * (y - c2) ), 
					c2 + ( math.sin(theta) * (x - c1) + math.cos(theta) * (y - c2) ) ]



if __name__ == "__main__":
	err_x = random() * 3
	err_y = random() * 3
	err_rot = -random()/2
	points1 = []
	points2 = []
	numPoints = 10
	for i in range(1, numPoints+1):
		x = random()
		y = random()
		points1.append([x, y])
		points2.append([x + err_x, y + err_y])

	com = PointAligner.centerOfMass(points2)
	for i in range(1,numPoints+1):
		points2[i-1] = PointAligner.rotate(points2[i-1][0],points2[i-1][1],com[0],com[1],err_rot)

	PointAligner.align(
			points1 = points1, 
			points2 = points2,
			iterations = 50)


