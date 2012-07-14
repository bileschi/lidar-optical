#!/usr/bin/python
import matplotlib.pyplot as plt
from simple_illustrations import illustrate_toy_problem
from align import gen_3d_toy_problem, mutate_camera_params

[space_pts, img_pts, true_params, guess_params] = gen_3d_toy_problem()
# illustrate rotations
for i in range(0,4):
	guess_params = mutate_camera_params(guess_params, z_rot=0.2)
	illustrate_toy_problem(space_pts, true_params, guess_params
		, title1='z_rot %d' % i)
	plt.pause(.2)
guess_params = true_params
for i in range(0,4):
	guess_params = mutate_camera_params(guess_params, y_rot=0.2)
	illustrate_toy_problem(space_pts, true_params, guess_params
		, title1='y_rot %d' % i)
	plt.pause(.2)
guess_params = true_params
for i in range(0,4):
	guess_params = mutate_camera_params(guess_params, x_rot=0.2)
	illustrate_toy_problem(space_pts, true_params, guess_params
	, title1='x_rot %d' % i)
	plt.pause(.2)
guess_params = true_params

#illustrate translations
for i in range(0,4):
	guess_params = mutate_camera_params(guess_params, t_z=3)
	illustrate_toy_problem(space_pts, true_params, guess_params
	, title1='z translation %d' % i)
	plt.pause(.2)
guess_params = true_params
for i in range(0,4):
	guess_params = mutate_camera_params(guess_params, t_y=3)
	illustrate_toy_problem(space_pts, true_params, guess_params
	, title1='y translation %d' % i)
	plt.pause(.2)
guess_params = true_params
for i in range(0,4):
	guess_params = mutate_camera_params(guess_params, t_x=3)
	illustrate_toy_problem(space_pts, true_params, guess_params
	, title1='x translation %d' % i)
	plt.pause(.2)
guess_params = true_params

# illustrate scaling
guess_params = mutate_camera_params(guess_params, k_mult=0.5)
for i in range(0,6):
	guess_params = mutate_camera_params(guess_params, k_mult=1.4)
	illustrate_toy_problem(space_pts, true_params, guess_params
	, title1='k %d' % i)
	plt.pause(.2)









