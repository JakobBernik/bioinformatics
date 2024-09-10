import unittest

from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage
import numpy as np

from helper_functions import UPGMA


class TestLocalAlignment(unittest.TestCase):
    def test_upgma1(self):
        distances = np.array([[0,4,5,6],
                                [4,0,3,4],
                                [5,3,0,3],
                                [6,4,3,0]])

        valid_alignments = [
            np.array([[1,2,3,2],[3,4,3.5,3],[0,5,5,4]]),
            np.array([[2,3,3,2],[1,4,3.5,3],[0,5,5,4]])
        ]

        res = UPGMA(distances)
        a = False
        for valid in valid_alignments:
            a |= np.alltrue(res == valid)
        self.assertTrue(a)
    
    def test_upgma2(self):
        distances = np.array([[0,17,21,31,23],
                                [17,0,30,34,21],
                                [21,30,0,28,39],
                                [31,34,28,0,43],
                                [23,21,39,43,0]])

        valid_alignments = [
            np.array([[0,1,17,2],[4,5,22,3],[2,3,28,2],[6,7,33,5]])
        ]

        res = UPGMA(distances)
        a = False
        for valid in valid_alignments:
            a |= np.alltrue(res == valid)
        self.assertTrue(a)
    
    def test_upgma3(self):
        # Create random distances matrix and solve using linkage and compare
        # with our implementation
        n = 100
        for _ in range(10):
            distances = np.random.rand(n,n)
            distances = distances * distances.T # Make the matrix symmetric
            for i in range(n):
                distances[i,i] = 0
            np.testing.assert_allclose(UPGMA(distances), linkage(squareform(distances), method='average'))



if __name__ == "__main__":
    unittest.main()