import unittest
from Bio import Phylo

from MinPD import GetMinPD
from main import getTreeFromFile
from preprocessing import preprocessTree
from MaxPD import GetMaxPD
from avg import GetnewAvgPD

class TestPDStats(unittest.TestCase):

    def test_10_min(self):
        tree = getTreeFromFile("t_11_1000")
        x = preprocessTree(tree, 10)
        y = GetMinPD(tree, 10, x[0], x[1], x[3])

        self.assertEqual(round(y[10]),14 )  # add assertion here

    def test_10_max(self):
        tree = getTreeFromFile("t_11_1000")
        x = preprocessTree(tree, 10)
        y = GetMaxPD(tree, 10, x[0], x[1], x[3])

        self.assertEqual(round(y[10]),102 ) # add assertion here

    def test_10_avg(self):
        tree = getTreeFromFile("t_11_1000")
        x = preprocessTree(tree, 10)
        y=GetnewAvgPD(tree,10,x[0],x[1],x[2],x[3])

        self.assertEqual(round(y[0][10]),85) # add assertion here


if __name__ == '__main__':
    unittest.main()
