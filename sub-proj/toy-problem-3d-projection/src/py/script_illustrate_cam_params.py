#!/usr/bin/python
import matplotlib.pyplot as plt
from simple_illustrations import illustrate_toy_problem
from align import gen_3d_toy_problem, mutate_camera_params

[space_pts, img_pts, true_params, guess_params] = gen_3d_toy_problem()
# illustrate rotations
for i in range(0,4):
	guess_params = mutate_camera_params(guess_params, z_rot=0.1)
	illustrate_toy_problem(space_pts, true_params, guess_params)
	plt.pause(.2)
guess_params = true_params
for i in range(0,4):
	guess_params = mutate_camera_params(guess_params, y_rot=0.1)
	illustrate_toy_problem(space_pts, true_params, guess_params)
	plt.pause(.2)
guess_params = true_params
for i in range(0,4):
	guess_params = mutate_camera_params(guess_params, x_rot=0.1)
	illustrate_toy_problem(space_pts, true_params, guess_params)
	plt.pause(.2)
guess_params = true_params
#illustrate translations
for i in range(0,4):
	guess_params = mutate_camera_params(guess_params, t_z=2)
	illustrate_toy_problem(space_pts, true_params, guess_params)
	plt.pause(.2)
guess_params = true_params
for i in range(0,4):
	guess_params = mutate_camera_params(guess_params, t_y=2)
	illustrate_toy_problem(space_pts, true_params, guess_params)
	plt.pause(.2)
guess_params = true_params
for i in range(0,4):
	guess_params = mutate_camera_params(guess_params, t_x=2)
	illustrate_toy_problem(space_pts, true_params, guess_params)
	plt.pause(.2)
guess_params = true_params
# illustrate scaling
guess_params = mutate_camera_params(guess_params, k_mult=2)
for i in range(0,6):
	guess_params = mutate_camera_params(guess_params, k_mult=0.8)
	illustrate_toy_problem(space_pts, true_params, guess_params)
	plt.pause(.2)









