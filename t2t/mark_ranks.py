#!/usr/bin/env python

import sys
from collections import defaultdict

from skbio import TreeNode

__author__ = "Donovan Park"
__copyright__ = "Copyright 2015, The tax2tree project"
__credits__ = ["Donovan Park"]
__license__ = "BSD"
__version__ = "1.0"
__maintainer__ = "Donovan Park"
__email__ = "donovan.parks@gmail.com"
__status__ = "Development"


class MarkRanks(object):
    """Marks internal nodes with predicted taxonomic rank.

    Predicted taxonomic ranks are based on the average distance to
    leaf nodes.
    """
    def __init__(self):
        """Initialize class."""
        self.rank_prefixes = ['D__', 'P__', 'C__', 'O__', 'F__', 'G__', 'S__']

    def average_dist_to_leaves(self, node):
        """Calculate average distance from an internal node to its leaves.

        Parameters
        ----------
        node : TreeNode
            Internal node of interest.

        Returns
        -------
        float
            Average distance from node to its leaves.
        """

        if node.is_tip():
            return 0

        sum_dist = 0.0
        num_tips = 0
        for leaf in node.tips():
            sum_dist += leaf.accumulate_to_ancestor(node)
            num_tips += 1

        return sum_dist / num_tips

    def mark(self, tree, start_labels, start_rank, output):
        """Determine taxonomic ranks for all internal nodes.

        Each internal node is assigned a predicted taxonomic
        rank based on the branch length from internal nodes
        to leaves.

        Parameters
        ----------
        tree : str or TreeNode
            A newick string or a TreeNode
        start_label: str
            Label of node to start marking internal nodes.
        start_rank: int
            Rank of node with start_label.
        ouput_tree : Output stream
            Name of tree with marked taxonomic ranks.
        """

        # make sure we have a TreeNode object
        root = tree
        if not isinstance(root, TreeNode):
            root = TreeNode.read(root, convert_underscores=False)

        # find specified starting node
        start_nodes = []
        try:
            for label in start_labels:
                start_nodes.append(root.find(label))
        except:
            print 'Unable to locate node with label: ' + label
            return

        cur_nodes = set(start_nodes)
        for cur_rank in xrange(start_rank + 1, len(self.rank_prefixes) - 1):
            # mark putative nodes at current rank
            nodes_at_next_rank = set()
            for node in cur_nodes:
                # get average distance to ancestral nodes
                # marked at parental rank
                parent_node = node.parent
                dist = node.length
                dist_to_phyla = []
                while parent_node != None and node not in start_nodes:
                    if parent_node.name.startswith(self.rank_prefixes[cur_rank - 1]):
                        dist_to_phyla.append(dist)

                    if parent_node.length:
                        dist += parent_node.length
                    parent_node = parent_node.parent

                if len(dist_to_phyla) != 0:
                    mean_dist_to_phyla = sum(dist_to_phyla) / len(dist_to_phyla)
                else:
                    mean_dist_to_phyla = 0

                # get average distance to leaf nodes
                average_dist = self.average_dist_to_leaves(node)

                # get distance from node to rank of class
                rank_dist = (mean_dist_to_phyla + average_dist) / (len(self.rank_prefixes) - cur_rank - 1)

                # mark putative nodes at the current rank and get children
                # nodes below this rank
                nodes_at_rank = set()
                children_nodes = set()
                for n in node.non_tips():
                    if n.accumulate_to_ancestor(node) + mean_dist_to_phyla < rank_dist:
                        n.name = self.rank_prefixes[cur_rank] + '|' + n.name
                        print self.rank_prefixes[cur_rank] + ': ' + n.name

                        nodes_at_rank.add(n)
                        for child in n.children:
                            children_nodes.add(child)

                # the children of a node marked at the current rank
                # may also be marked at the current rank so should be
                # removed from the set of putative nodes to be marked
                # at the next rank
                nodes_at_next_rank |= (children_nodes - nodes_at_rank)

            cur_nodes = nodes_at_next_rank

        # save marked tree
        root.write(output)
