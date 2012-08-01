#!/usr/bin/python
import numpy as np
import pdb

"""
Tools to project 3d points into a 2d image.
Data is expected to be in np.mat format.  Some functions allow for array-like input.
"""

def proj_3d_jacobian(space_pt = None, current_params = None):
	"""
	Let i = f(s,p) where f is a function from space S to space I.
	f is the function calculating the projection of s, parameterized by p. 

	proj_3d_jacobian calculates the partial derivatives of i with respect to 
	p.  In general, the derivative may depend on the space_point s.

	
	Input:
	space_pt: the tuple representing the point in space to project
	current_params: p.  a dict indicating the current param settings p

	Output:
	Output array[m,n] is partial derivative of i[m] W.R.T. p[n]

	i in non-homogenous coordinates is [x, y]
	s in homogenous coordinates is [s0, s1, s2, s3]
	p is a 3x4 matrix. it has 12 elements of the form pjk

	d(x)/d(p0j)) = sj / sum_k(p2k * sk)
	d(x)/d(p1j)) = 0
	d(x)/d(p2j)) = - sj * sum_k(p0k * sk) / (sum_k(p2k * sk) ^ 2)

	d(y)/d(p0j)) = 0
	d(y)/d(p1j)) = sj / sum_k(p2k * sk)
	d(y)/d(p2j)) = - sj * sum_k(p1k * sk) / (sum_k(p2k * sk) ^ 2)



	"""
	(s0, s1, s2) = space_pt
	p00 = current_params[0]
	p01 = current_params[1]
	p02 = current_params[2]
	p03 = current_params[3]

	p10 = current_params[4]
	p11 = current_params[5]
	p12 = current_params[6]
	p13 = current_params[7]

	p20 = current_params[8]
	p21 = current_params[9]
	p22 = current_params[10]
	p23 = current_params[11]

	c0 = p00 * s0 + p01 * s1 + p02 * s2 + p03 
	c1 = p10 * s0 + p11 * s1 + p12 * s2 + p13 
	c2 = p20 * s0 + p21 * s1 + p22 * s2 + p23 
	c2c2 = c2 * c2

	return np.mat([
				[s0/c2, s1/c2, s2/c2, 1.0/c2, # dx wrt p0_
				0.0, 0.0, 0.0, 0.0, # dx wrt p1_
				-s0*c0/c2c2, -s1*c0/c2c2, -s2*c0/c2c2, -c0/c2c2], # dx wrt p2_

				[0.0, 0.0, 0.0, 0.0, # dy wrt p0_
				s0/c2, s1/c2, s2/c2, 1.0/c2, # dy wrt p1_
				-s0*c1/c2c2, -s1*c1/c2c2, -s2*c1/c2c2, -c1/c2c2]]) # dy wrt p2_

def project_3d_legacy(space_pts, cam_params):
	imgs = []
	RI = cam_params['R'].I
	k = cam_params['k']
	cx = cam_params['cx']
	cy = cam_params['cy']
	for pt in space_pts:
		tpt = RI.dot((pt - cam_params['t']).T)
		denom = tpt.item(2)
		if(denom <= 0):
			imgs.append(None) # point behind or to side of camera
		else:
			px = k * tpt.item(0) / denom + cx
			py = k * tpt.item(1) / denom + cy
			if (px < 0) or (py < 0) or (px > (cx * 2)) or (py > (cy * 2)):
				imgs.append(None) # point outside camera frame
			else:
				imgs.append(np.mat([px, py]))
	return imgs

def cam_as_list(cam_params):
	P = proj_mat_from_cam_params(cam_params=cam_params)
	return P.reshape((1,12)).tolist()[0]

def cam_P_to_dict(cam_params):
	P = proj_mat_from_cam_params(cam_params=cam_params)
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



def proj_mat_from_cam_params(RI=None, k=None, cx=None, cy=None, t=None, 
	cam_params=None):
	"build projection matrix either from camera elements or cam-as-dict"
	if cam_params is not None:
		RI = cam_params['R'].I
		k = cam_params['k']
		cx = cam_params['cx']
		cy = cam_params['cy']
		t = cam_params['t']
	K = np.mat([[k, 0.0, cx], [0.0, k, cy], [0.0, 0.0, 1]])
	P = np.concatenate((RI, -RI.dot(t.T)), axis=1)
	return K.dot(P)

def hspace_to_himg(hspace_pts, P):
	P = np.mat(P).reshape((3,4))
	h_img_pts = []
	try:
		for h_spc_pt in hspace_pts:
			h_img_pt = P.dot(np.mat(h_spc_pt).T)
			h_img_pts.append(h_img_pt)
	except ValueError:
		pdb.set_trace()
	return h_img_pts

def himg_to_pixel(himg_pts, cam_params):
	cx = cam_params['cx']
	cy = cam_params['cy']
	px_pts = []
	for himg_pt in himg_pts:
		denom = himg_pt.item(2)
		if(denom <= 0):
			px_pts.append(None) # point behind or to side of camera
		else:
			px = himg_pt.item(0) / denom
			py = himg_pt.item(1) / denom
			if (px < 0) or (py < 0) or (px > (cx * 2)) or (py > (cy * 2)):
				px_pts.append(None) # point outside camera frame
			else:
				px_pts.append(np.mat([px, py]))
	return px_pts

def project_3d(space_pts, cam_params):
	" using matrix multiplication"
	px_pts = []
	try: # cam_params as dictionary
		RI = cam_params['R'].I
		k = cam_params['k']
		cx = cam_params['cx']
		cy = cam_params['cy']
		t = cam_params['t']
		if cam_params.has_key('P'):
			P = cam_params['P']
		else:
			P = proj_mat_from_cam_params(RI, k, cx, cy, t)
			cam_params['P'] = P
	except TypeError: # cam params as list
		P = np.mat(cam_params).reshape((3,4))
	hspace_pts = [to_homg(pt) for pt in space_pts]
	himg_pts = hspace_to_himg(hspace_pts, P)
	try:
		px_pts = himg_to_pixel(himg_pts, cam_params)
		return px_pts
	except TypeError:
		#if the cam_params is just a projection matrix (no cx, cy) use def
		return  himg_to_pixel(himg_pts, {'cx':320, 'cy':240})

def to_unhomg(p3):
	return np.take(p3, range(0, p3.size-1), mode='wrap') / np.take(p3, [-1], mode='wrap')
	 
def to_homg(pt):
	return np.append(np.mat(pt), np.mat(1), 1)
