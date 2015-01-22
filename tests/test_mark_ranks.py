#!/usr/bin/env python

__author__ = "Donovan Park"
__copyright__ = "Copyright 2015, The tax2tree project"
__credits__ = ["Donovan Park"]
__license__ = "BSD"
__version__ = "1.0"
__maintainer__ = "Donovan Park"
__email__ = "donovan.parks@gmail.com"
__status__ = "Development"

from unittest import TestCase, main
from io import StringIO

from t2t.mark_ranks import MarkRanks

from skbio import TreeNode


class MarkRanksTests(TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_average_dist_to_leaves(self):
        """Test calculation of average distance to leaf nodes"""

        mr = MarkRanks()

        # test ultrametric tree (i.e., equal branch lengths)
        tree = TreeNode.read(StringIO(u'((a:1,b:1):1,(d:1,e:1):1);'))
        avg_dist = mr.average_dist_to_leaves(tree)
        self.assertAlmostEqual(avg_dist, 2.0)

        # test calculation under varying branch lengths
        tree = TreeNode.read(StringIO(u'((a:2,b:2):2,(d:1,e:1):1);'))
        avg_dist = mr.average_dist_to_leaves(tree)
        self.assertAlmostEqual(avg_dist, 3.0)

        # test with branch length of zero
        tree = TreeNode.read(StringIO(u'((a:1,b:2):0,(d:2,e:1):1);'))
        avg_dist = mr.average_dist_to_leaves(tree)
        self.assertAlmostEqual(avg_dist, 2.0)

        # test calculation of leaf nodes
        for leaf_node in tree.tips():
            avg_dist = mr.average_dist_to_leaves(leaf_node)
            self.assertAlmostEqual(avg_dist, 0.0)

        # test calculation on internal nodes
        for child in tree.children:
            avg_dist = mr.average_dist_to_leaves(child)
            self.assertAlmostEqual(avg_dist, 6.0 / 4)




if __name__ == '__main__':
    main()
