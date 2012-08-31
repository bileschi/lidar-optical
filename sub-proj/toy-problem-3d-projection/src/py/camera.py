#!/usr/bin/python
import numpy as np
import pdb
from rotations_3d import rot_mat

class Camera(object):
	# Camera does not store the individual rotation parameters (r_x, r_y, r_z)
	# it only stores the computed rotation matrix R.
	def __init__(self, 
		t_x=0, 
		t_y=0, 
		t_z=0,
		R=rot_mat(0,0,0), 
		c_x=320, 
		c_y=240,
		k=1000):
		self.t_x = t_x
		self.t_y = t_y
		self.t_z = t_z
		self.c_x = c_x
		self.c_y = c_y
		self.k = k
		self.R = R
		self.P = self.proj_mat()

	def __str__(self):
		return self.pretty_str()

	def pretty_str(self):
		out_str = "camera object\n"
		out_str += "trans: %f %f %f\n" % (self.t_x, self.t_y, self.t_z)
		out_str += "Rot:\n %s\n" % str(self.R)
		out_str += "k: %f\n" % self.k
		out_str += "c: %f %f\n" % (self.c_x, self.c_y)
		return out_str

	def proj_mat(self):
		"build projection matrix from camera elements"
		RI = self.R.I
		t = np.mat([self.t_x, self.t_y, self.t_z])
		K = np.mat([[self.k, 0.0, self.c_x], [0.0, self.k, self.c_y], [0.0, 0.0, 1]])
		P = np.concatenate((RI, -RI.dot(t.T)), axis=1)
		return K.dot(P)

	def as_list(self):
		return self.P.reshape((1,12)).tolist()[0]

	def as_dict(self):
		"dictionary version of the projection for legacy needs"
		P = self.P
		return {
			'p00':P.item(0),
			'p01':P.item(1),
			'p02':P.item(2),
			'p03':P.item(3),

			'p10':P.item(4),
			'p11':P.item(5),
			'p12':P.item(6),
			'p13':P.item(7),

			'p20':P.item(8),
			'p21':P.item(9),
			'p22':P.item(10),
			'p23':P.item(11)}

	def mutate(self, 
		t_x=0.0,
		t_y=0.0,
		t_z=0.0,
		r_x=0.0,
		r_y=0.0,
		r_z=0.0,
		k_mult=1.0):
		""" Returns a new copoy of the camera, with modified parameters. 
			c_x and c_y are not modified """
		cam2 = Camera(
			t_x = self.t_x + t_x,
			t_y = self.t_y + t_y,
			t_z = self.t_z + t_z,
			R = rot_mat(r_x, r_y, r_z).dot(self.R),
			k = self.k * k_mult
			)
		return cam2;
