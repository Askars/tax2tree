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

    def _remove_dummy_nodes(self, node):
        """Remove all dummy nodes from a tree.

        A dummy node is define as a node with only one child.

        This function modifies the tree specified by node.

        Parameters
        ----------
        node : TreeNode
            Internal node of interest.
        """

        # get list of nodes in preorder as we will be
        # modifying the tree as we go
        preorder = list(node.preorder())

        for n in preorder:
            if len(n.children) == 1:
                # remove dummy node and connect its single
                # child to its parent; adjusting the branch
                # length as necessary
                length = n.length + n.children[0].length
                n.children[0].length = length
                n.parent.append(n.children[0])
                n.parent.remove(n)

    def mark(self, tree, start_labels, start_rank, verbose=False):
        """Determine taxonomic ranks for all internal nodes.

        Each internal node is assigned a predicted taxonomic
        rank based on the branch length from internal nodes
        to leaves.

        Parameters
        ----------
        tree : str or TreeNode
            A newick string or a TreeNode
        start_labels: str
            Labels of nodes to start marking taxonomic ranks.
        start_rank: int
            Rank of node with start_label.
        """

        # make sure we have a TreeNode object
        root = tree
        if not isinstance(root, TreeNode):
            root = TreeNode.read(root, convert_underscores=False)

        # find specified starting node
        start_nodes = []
        try:
            for label in start_labels:
                start_node = root.find(label)
                start_node.name = self.rank_prefixes[start_rank] + '|' + start_node.name
                start_nodes.append(start_node)
        except:
            print 'Unable to locate node with label: ' + label
            return None

        nodes_at_rank = set(start_nodes)

        # assign internal nodes with ranks from [phylum, genus],
        # as it is assumed the terminal taxa are dereplicated at
        # the species-level
        for cur_rank in xrange(start_rank, len(self.rank_prefixes) - 1):
            print '******************* RANK: ' + str(cur_rank)

            # get immediate children of nodes at current rank
            child_lineages = set()
            for n in nodes_at_rank:
                for child in n.children:
                    child_lineages.add(child)

            child_lineages = child_lineages - nodes_at_rank

            # mark putative nodes at child rank
            next_nodes = set()
            for node in child_lineages:
                print 'processing node: ' + node.name

                # get average distance to ancestral nodes marked at current rank
                parent_node = node.parent
                dist = node.length
                dist_to_parent_rank = []
                while parent_node != None and node not in start_nodes:
                    if parent_node.name.startswith(self.rank_prefixes[cur_rank]):
                        dist_to_parent_rank.append(dist)

                    if parent_node.length:
                        dist += parent_node.length
                    parent_node = parent_node.parent

                if len(dist_to_parent_rank) != 0:
                    mean_dist_to_parent_rank = sum(dist_to_parent_rank) / len(dist_to_parent_rank)
                else:
                    mean_dist_to_parent_rank = 0

                print '  mean_dist_to_parent_rank: ' + str(mean_dist_to_parent_rank)

                # get average distance to leaf nodes
                average_dist = self.average_dist_to_leaves(node)
                print '  average_dist: ' + str(average_dist)

                # get distance from current node for descendant nodes
                # to be marked at the current rank
                rank_dist = (mean_dist_to_parent_rank + average_dist) / (len(self.rank_prefixes) - cur_rank - 1)
                print '  rank_dist: ' + str(rank_dist)

                # mark putative nodes at the child rank
                nodes_to_process = [node]
                nodes_at_child_rank = set()
                while nodes_to_process:
                    n = nodes_to_process.pop()

                    if n == node:
                        dist = mean_dist_to_parent_rank
                    else:
                        dist = n.accumulate_to_ancestor(node) + mean_dist_to_parent_rank

                    if dist <= rank_dist:
                        n.name = self.rank_prefixes[cur_rank + 1] + '|' + n.name

                        if verbose:
                            print '  ' + self.rank_prefixes[cur_rank + 1] + ': ' + n.name

                        nodes_at_child_rank.add(n)
                        for child in n.children:
                            nodes_to_process.append(child)

                # add in dummy nodes so all children lineages contain a
                # node at the child rank
                if not nodes_at_child_rank:
                    name = self.rank_prefixes[cur_rank + 1] + '|' + node.name
                    parent = node.parent
                    dummy = TreeNode(name, 0.5 * rank_dist, parent, [node])
                    parent.append(dummy)
                    node.length = node.length - 0.5 * rank_dist

                    nodes_at_child_rank.add(dummy)

                    print '  name of dummy node: ' + dummy.name
                    print '  children of dummy node: ' + ','.join([x.name for x in dummy.children])
                    print '  parent of dummy node: ' + dummy.parent.name
                    print '  children of parent: ' + ','.join([x.name for x in dummy.parent.children])
                    print '  parent of children: ' + ','.join([x.name for x in dummy.children])

                print '  nodes_at_child_rank: ' + ','.join([x.name for x in nodes_at_child_rank])

                next_nodes |= nodes_at_child_rank
                print '  next_nodes: ' + ','.join([x.name for x in next_nodes])
                print ''

            nodes_at_rank = next_nodes

        # remove dummy nodes
        self._remove_dummy_nodes(root)

        # save marked tree
        return root
