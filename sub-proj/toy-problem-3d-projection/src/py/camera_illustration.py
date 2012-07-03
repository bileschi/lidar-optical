import numpy as np

def get_default_camera(size=2, k=1):
	""" lines describing a camera oriented down the z axis"""
	cam = []
	k = k * size/2.0
	cam.append([(0,0,0), (+k,+k,size)]) # edges
	cam.append([(0,0,0), (+k,-k,size)])
	cam.append([(0,0,0), (-k,-k,size)])
	cam.append([(0,0,0), (-k,+k,size)])
	cam.append([(+k,+k,size), (+k,-k,size)])
	cam.append([(+k,-k,size), (-k,-k,size)])
	cam.append([(-k,-k,size), (-k,+k,size)])
	cam.append([(-k,+k,size), (+k,+k,size)])
	cam.append([(0,0,0), (0,+k,size)]) # top marker
	return cam

def move_camera(in_cam, t=np.mat([0,0,0]), R=np.mat(np.diag([1,1,1]))):
	out_cam = []
	for line in in_cam:
		transformed_line = map(lambda pt: R.dot(pt) + t, line) # rotate then transalte
		out_cam.append(transformed_line)
	return out_cam

