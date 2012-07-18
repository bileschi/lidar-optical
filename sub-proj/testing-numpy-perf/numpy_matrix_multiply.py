#!/usr/bin/python
import sys, random, time, numpy as np
from datetime import datetime


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
	
	n=2
	while( n <= 1024):
		m1 = create_matrix_lists(n,n)
		m2 = create_matrix_lists(n,n)

		print "np was ", timeThis(matrix_multiply, [m1, m2]) - timeThis(np.dot, [np.mat(m1), np.mat(m2)]), " slower on matrices of size ", n, "\n"
		n = n * 2



