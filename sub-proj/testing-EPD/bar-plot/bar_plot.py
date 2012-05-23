# just an illustration how to draw a bar plot.

import matplotlib.pyplot as plt

def illustrate_forces(x_val):
	plt.figure(1)
	plt.clf()
	plt.bar(range(len(x_val)), x_val)
	plt.show()

def main():
	illustrate_forces([x**2 for x in range(1,9) ])

if __name__ == '__main__':
	main()
