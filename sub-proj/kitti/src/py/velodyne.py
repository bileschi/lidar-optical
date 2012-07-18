#!/usr/bin/python

def load_file(filename):
	from array import array
	input_file = open(filename, 'r')
	float_array = array('f')
	float_array.fromstring(input_file.read())
	x = float_array[0::4]
	y = float_array[1::4]
	z = float_array[2::4]
	r = float_array[3::4]
	return (x,y,z,r)


if __name__ == "__main__":
	import os
	vel_dir = '/Users/stanleybileschi/ResearchData/KITTI/raw-data/2011_09_26_drive_0001/velodyne_points/data'
	vel_file = '0000000043.bin'
	(x,y,z,_) = load_file(os.path.join(vel_dir, vel_file))
	import matplotlib.pyplot as plt
	plt.figure(0)
	plt.clf()
	plt.plot(x, y, 'bx')
	plt.show()
	plt.figure(1)
	plt.clf()
	plt.plot(x, z, 'bx')
	plt.show()
	plt.figure(2)
	plt.clf()
	plt.plot(y, z, 'bx')
	plt.show()

	import numpy as np
	from mpl_toolkits.mplot3d import Axes3D
	

	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')
	n = 100
	ax.scatter(x[::100], y[::100], z[::100], c='b', marker='^')

	ax.set_xlabel('X Label')
	ax.set_ylabel('Y Label')
	ax.set_zlabel('Z Label')

	plt.show()
