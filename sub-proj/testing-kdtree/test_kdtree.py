from scipy import spatial
import numpy as np

if __name__ == "__main__":
    x, y = np.mgrid[0:5, 2:8]
    tree = spatial.KDTree(zip(x.ravel(), y.ravel()))
    tree.data
    pts = np.array([[0, 0], [2.1, 2.9]])
    tree.query(pts)
