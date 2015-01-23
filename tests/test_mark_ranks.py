"""Test marking of taxonomic ranks"""  #!/usr/bin/env python

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
        tree = TreeNode.read(StringIO(u'((a:1,b:1):1,(c:1,d:1):1);'))
        avg_dist = mr.average_dist_to_leaves(tree)
        self.assertAlmostEqual(avg_dist, 2.0)

        # test calculation under varying branch lengths
        tree = TreeNode.read(StringIO(u'((a:2,b:2):2,(c:1,d:1):1);'))
        avg_dist = mr.average_dist_to_leaves(tree)
        self.assertAlmostEqual(avg_dist, 3.0)

        # test with branch length of zero
        tree = TreeNode.read(StringIO(u'((a:1,b:2):0,(c:2,d:1):1);'))
        avg_dist = mr.average_dist_to_leaves(tree)
        self.assertAlmostEqual(avg_dist, 2.0)

        # test calculation of leaf nodes
        for leaf_node in tree.tips():
            avg_dist = mr.average_dist_to_leaves(leaf_node)
            self.assertAlmostEqual(avg_dist, 0.0)

        # test calculation on internal nodes
        for child in tree.children:
            avg_dist = mr.average_dist_to_leaves(child)
            self.assertAlmostEqual(avg_dist, 3.0 / 2)

    def test_mark_ultrametric(self):
        """Test marking of taxonomic ranks on ultrametric tree"""

        mr = MarkRanks()

        tree = StringIO(u'((a:5,b:5)t1:1,(c:1,d:1)t2:5)root;')
        root = mr.mark(tree, ['root'], 0)

        # verify rank assignments of nodes
        pass_test = True
        try:
            root.find('S__|a')
            root.find('S__|b')
            root.find('S__|c')
            root.find('S__|d')
            root.find('P__|t1')
            root.find('S__|t2')
            root.find('D__|root')
        except:
            pass_test = False

        self.assertTrue(pass_test, "Unexpected taxonomic rank assigned to an internal node.")

        self.assertEqual(root.count(), 7)
        self.assertEqual(root.count(tips=True), 4)

    def test_mark_skewed_tree(self):
        """Test marking of taxonomic ranks on tree with highly varying rates of evolution"""

        mr = MarkRanks()

        tree = StringIO(u'((a:5,b:5)t1:1,(c:10,d:10)t2:5)root;')
        root = mr.mark(tree, ['root'], 0)

        # verify rank assignments of nodes
        pass_test = True
        try:
            root.find('S__|a')
            root.find('S__|b')
            root.find('S__|c')
            root.find('S__|d')
            root.find('P__|t1')
            root.find('O__|t2')
            root.find('D__|root')
        except:
            pass_test = False

        self.assertTrue(pass_test, "Unexpected taxonomic rank assigned to an internal node.")

        self.assertEqual(root.count(), 7)
        self.assertEqual(root.count(tips=True), 4)

    def test_mark_skewed_leaves(self):
        """Test marking of taxonomic ranks on tree with highly varying rates at leaf nodes"""

        mr = MarkRanks()

        tree = StringIO(u'((a:5,b:5)t1:1,(c:10,d:1)t2:5)root;')
        root = mr.mark(tree, ['root'], 0)

        # verify rank assignments of nodes
        pass_test = True
        try:
            root.find('S__|a')
            root.find('S__|b')
            root.find('S__|c')
            root.find('S__|d')
            root.find('P__|t1')
            root.find('F__|t2')
            root.find('D__|root')
        except:
            pass_test = False

        self.assertTrue(pass_test, "Unexpected taxonomic rank assigned to an internal node.")

        self.assertEqual(root.count(), 7)
        self.assertEqual(root.count(tips=True), 4)

    def test_mark_missing_label(self):
        """Test proper handling of missing start label"""

        mr = MarkRanks()

        # test missing label
        tree = StringIO(u'((a:1,b:1):1,(c:1,d:1):1)root;')
        root = mr.mark(tree, ['bad_label'], 0)
        self.assertIsNone(root)

        # test valid labels
        tree = StringIO(u'((a:1,b:1)r1:1,(c:1,d:1)r2:1)root;')
        root = mr.mark(tree, ['r1', 'r2'], 0)
        self.assertIsNotNone(root)

    def test_remove_dummy_nodes(self):
        """Test removal of dummy nodes (nodes with a single child)."""

        mr = MarkRanks()

        root = TreeNode.read(StringIO(u'(((a:1,b:1)c:1)d:1,(e:1,f:1)g:1)r;'))
        mr._remove_dummy_nodes(root)

        # verify node d has been removed
        try:
            node = root.find('d')
        except:
            node = None

        self.assertIsNone(node)
        self.assertEqual(root.count(), 7)
        self.assertEqual(root.count(tips=True), 4)


if __name__ == '__main__':
    main()
