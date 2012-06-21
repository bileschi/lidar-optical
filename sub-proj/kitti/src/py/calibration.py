#!/usr/bin/python
import numpy as np

def load_projection_matrix(filename, want_matrix = 'P'):
	"""kitti calibrations are stored in files named like:
	KITTI/raw-data/2011_09_28_calib/calib_velo_to_cam.txt

	the contents of which are like:
	calib_time: 15-Mar-2012 11:37:52
	R: 6.927964e-03 -9.999722e-01 -2.757829e-03 -1.162982e-03 2.749836e-03 -9.999955e-01 9.999753e-01 6.931141e-03 -1.143899e-03
	T: -2.457729e-02 -6.127237e-02 -3.321029e-01
	delta_f: 0.000000e+00 0.000000e+00
	delta_c: 0.000000e+00 0.000000e+00

	R indicates a rotation matrix and T indicates a translation matrix.
	The other fields are unused.

	The default return value is the projection matrix P.
	P := K * [R | -T]
	"""
	# read in whole calibration file
	with open(filename) as f:
		content = f.readlines()
	# save desired matricies
	# K = np.array([[1.0, 0, 0], [0, 1.0, 0], [0, 0, 1.0]])
	for line in content:
		(key, val) = line.split(':',1)
		if key == 'R':
			R = np.fromstring(val, sep=' ')
			R = R.reshape(3, 3)
		if key == 'T':
			T = np.fromstring(val, sep=' ')
			T = T.reshape(3, 1)
	if want_matrix == 'R':
		return R
	if want_matrix == 'T':
		return T
	P = np.hstack((R,-T))
	P = np.vstack((P,np.array([0, 0, 0, 1])))
	if want_matrix == 'P':
		return P
	return None

if __name__ == "__main__":
	P = load_projection_matrix("/Users/stanleybileschi/ResearchData/KITTI/raw-data/2011_09_28_calib/calib_velo_to_cam.txt")