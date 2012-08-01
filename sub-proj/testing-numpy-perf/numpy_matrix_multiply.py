#!/usr/bin/python
import sys, random, time, numpy as np
from datetime import datetime
from copy import copy

#create matrix of size 
def create_matrix_lists(rows, cols):
	to_ret = []
	for i in range(rows):
		to_push = []
		random.seed(i+time.time())
		for j in range(cols):
			to_push = to_push + [random.uniform(0, 100)] #+ [random.randint(1,10)]
		to_ret = to_ret + [to_push]
	return to_ret

def zero_matrix(m,n):
    # Create zero matrix
    new_matrix = [[0 for row in range(n)] for col in range(m)]
    return new_matrix

#no error checking for now
def matrix_multiply(mat1, mat2):
	to_ret = []
	for i in range(len(mat1)):
		to_push = []
		for j in range(len(mat1[i])):
			to_push = to_push + [mat1[i][j] * mat2[j][i]]
		to_ret = to_ret + [to_push]
	return to_ret

def matrix_multiply(matrix1,matrix2):
    # Matrix multiplicatizon
    if len(matrix1[0]) != len(matrix2):
        # Check matrix dimensions
        print 'Matrices must be m*n and n*p to multiply!'
    else:
        # Multiply if correct dimensions
        new_matrix = zero_matrix(len(matrix1),len(matrix2[0]))
        for i in range(len(matrix1)):
            for j in range(len(matrix2[0])):
                for k in range(len(matrix2)):
                    new_matrix[i][j] += matrix1[i][k]*matrix2[k][j]
        return new_matrix

def timeThis(func, params):
	t1 = datetime.now()
	func(params[0], params[1])
	t2 = datetime.now()
	elapsed = (t2 - t1)
	return elapsed

if __name__ == "__main__":
	sizes = [2**n for n in range(2,10)]
	for n in sizes:
		# first, select the numbers to fill the matricies
		x1 = create_matrix_lists(n, n)
		x2 = create_matrix_lists(n, n)
		# calculate the time to copy the list of lists object.  
		t1 = datetime.now()
		m1 = copy(x1)
		m2 = copy(x2)
		t2 = datetime.now()
		t_make_list = t2 - t1
		# time numpy's matrix generation.
		t1 = datetime.now()
		M1 = np.mat(x1)
		M2 = np.mat(x2)
		t2 = datetime.now()
		t_make_np = t2 - t1
		# time the two multiplication routines
		t_mult_list = timeThis(matrix_multiply, [m1, m2]);
		t_mult_np = timeThis(np.dot, [M1, M2])
		# pretty print output
		if (t_make_np < t_make_list):
			victor = "numpy"
		else:
			victor = "list"
		print "  size %d make: %s wins by %fs" % (n, victor, (abs(t_make_list - t_make_np)).total_seconds())
		if (t_mult_np < t_mult_list):
			victor = "numpy"
		else:
			victor = "list"
		print "  size %d mult: %s wins by %fs\n" % (n, victor, (abs(t_mult_list - t_mult_np)).total_seconds())



