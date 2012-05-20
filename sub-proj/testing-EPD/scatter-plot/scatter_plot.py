import matplotlib.pyplot as plt
plt.plot([1,2,3,4], [1,4,9,16], 'ro')
plt.plot([1,2,3,4], [16,9,4,1], 'bx')
plt.axis([0, 6, 0, 20])
plt.grid(True)
plt.title('this is a title')
plt.xlabel('this is x label')
plt.ylabel('this is y label')
plt.show()