#stdlib
from py.camera import *
from numpy.testing import assert_allclose
import math

expected_default_P = 	np.mat(
	[[1000, 0, 320, 0], 
	[0, 1000, 240, 0], 
	[0, 0, 1, 0]])
expected_default_dict = {'p03': 0.0, 'p02': 320.0, 'p01': 0.0, 'p00': 1000.0, 
	'p21': 0.0, 'p20': 0.0, 'p23': 0.0, 'p22': 1.0, 
	'p10': 0.0, 'p11': 1000.0, 'p12': 240.0, 'p13': 0.0}
expected_default_list = [1000.0, 0.0, 320.0, 0.0, 0.0, 1000.0, 240.0, 0.0, 0.0, 0.0, 1.0, 0.0]

translated_cam_P = np.mat(
	[[ 1000,   0,   320, -330000],
	 [ 0, 1000, 240, -340000],
     [ 0, 0, 1, -1000]])

rotated_cam_P = np.mat(
	[[1000.0, 0, -320.0, 0.0], 
	 [0.0, -1000.0, -240.0, 0.0], 
	 [0.0, 0, -1.0, 0.0]])


def test_make_camera():
	cam = Camera()
	assert_allclose(cam.P, expected_default_P)

def test_pretty_str():
	"test is brittle.  maybe just check that the strings for two different cameras are different?"
	cam = Camera()
	str = cam.pretty_str()
	assert("trans: 0.000000 0.000000 0.000000" in cam.pretty_str())
	assert("""Rot:
 [[ 1.  0.  0.]
 [ 0.  1.  0.]
 [ 0.  0.  1.]]""" in cam.pretty_str())
	assert("k: 1000.000000" in cam.pretty_str())
	assert("c: 320.000000 240.000000" in cam.pretty_str())

def test_as_dict():
	d = Camera().as_dict()
	assert(d == expected_default_dict)

def test_as_list():
	l = Camera().as_list()
	assert(l == expected_default_list)

def test_mutate():
	"""Tests that projection matrices match expectations.
	also tests that cam_1 is unmodified by the mutate."""
	cam_1 = Camera()
	cam_t = cam_1.mutate(t_x = 10, t_y = 100, t_z = 1000)
	cam_r = cam_1.mutate(r_x = math.pi)
	assert_allclose(cam_t.P, translated_cam_P)	
	assert_allclose(cam_r.P, rotated_cam_P, atol=1e-6)	

def test_proj_mat():
	cam = Camera()
	assert_allclose(cam.P, cam.proj_mat())

def test_mutation_commutivity():
	"rotate then translate vs translate then rotate"
	cam1 = Camera()
	cam1 = cam1.mutate(t_x = 10)
	cam1 = cam1.mutate(r_x = .5, r_y = .5, r_z = .5)
	cam2 = Camera()
	cam2 = cam2.mutate(r_x = .5, r_y = .5, r_z = .5)
	cam2 = cam2.mutate(t_x = 10)
	assert_allclose(cam1.P, cam2.P)
