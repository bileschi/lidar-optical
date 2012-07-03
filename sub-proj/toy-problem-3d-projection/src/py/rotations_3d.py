#!/usr/bin/python
import numpy as np

def x_rot(rads):
  return np.matrix([
  [1,0,0], 
  [0, np.cos(rads), -np.sin(rads)],
  [0, np.sin(rads), np.cos(rads)]])

def y_rot(rads):
  return np.matrix([
  [np.cos(rads), 0, np.sin(rads)],
  [0, 1, 0], 
  [-np.sin(rads), 0, np.cos(rads)]])

def z_rot(rads):
  return np.matrix([
  [np.cos(rads), -np.sin(rads), 0],
  [np.sin(rads), np.cos(rads), 0],
  [0,0,1]]) 

def rot_mat(x_rads = 0 , y_rads = 0, z_rads = 0):
	x_mat = x_rot(x_rads)
	y_mat = y_rot(y_rads)
	z_mat = z_rot(z_rads)
	return x_mat.dot(y_mat.dot(z_mat))

