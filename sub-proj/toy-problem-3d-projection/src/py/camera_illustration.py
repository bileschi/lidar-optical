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

def draw_camera(ax, cam, color=[0,0,0]):
	for line in cam:
		ax.plot([line[0].item(0), line[1].item(0)], 
			[line[0].item(1), line[1].item(1)],
			[line[0].item(2), line[1].item(2)],
			c=color)

def draw_3d_axis(ax, x_max = 12, y_max = 12, z_max = 12):
	for direction, color in zip((1, -1), [[0.3, 0.3, 0.3], [0.7, 0.3, 0.3]]):
		for point in np.diag(direction * np.array([x_max, y_max, z_max])):
			ax.plot([point[0]], [point[1]], [point[2]], 'w')
			ax.plot([0, point[0]], [0, point[1]], [0, point[2]], c=color)
	ax.set_xlabel('X')
	ax.set_ylabel('Y')
	ax.set_zlabel('Z')

