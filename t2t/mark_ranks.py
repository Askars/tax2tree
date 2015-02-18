#!/usr/bin/env python

import sys

from skbio import TreeNode

from numpy import mean as np_mean

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

    Predicted taxonomic ranks are based on the average distance
    to leaf nodes.
    """
    def __init__(self):
        """Initialize class."""
        self.rank_prefixes = ['D__', 'P__', 'C__', 'O__', 'F__', 'G__', 'S__', 'ST__']

    def _mean_dist_to_leaves(self, node):
        """Calculate mean distance from an internal node to its leaves.

        Parameters
        ----------
        node : TreeNode
            Internal node of interest.

        Returns
        -------
        float
            Mean distance to leaves.
        """

        if node.is_tip():
            return 0

        dist = []
        for leaf in node.tips():
            dist.append(leaf.accumulate_to_ancestor(node))

        return np_mean(dist)

    def write_consistency(self, consistency_table, root):
        """Write table of labeled ranks versus predicted ranks

        Parameters
        ----------
        consistency_table : str
            Name of table to create.
        root: TreeNode
            Tree of interest.
        """

        if consistency_table:
            fout = open(consistency_table, 'w')

            for n in root.preorder():
                if ':' in n.name:
                    label = n.name[n.name.find(':'):]
                else:
                    continue

                ranks = '<none>'
                if '|' in n.name:
                    ranks, _name = n.name.split('|')

                fout.write('%s\t%s\n' % label, ranks)

            fout.close()

    def __create_new_node(self, node, node_name, dist_to_parent):
        """Create a new node between specified node and its parent.

        Inserts a new node into the tree between the specified
        node and its parent. This node is placed at the specified
        distance from the parent. Nodes are updated to reflect
        the insertion of this node.

        This function modifies the tree specified by node.

        Parameters
        ----------
        node : TreeNode
            Node specifying branch to place dummy node.
        dist_to_parent: float
            Desired distance from new node to its parent.

        Returns
        -------
        TreeNode
            New node placed in tree.
        """
        parent = node.parent
        dummy_length = dist_to_parent
        dummy_node = TreeNode(node_name, dummy_length, parent, [node])
        parent.append(dummy_node)
        node.length = node.length - dummy_length

        return dummy_node

    def __remove_dummy_nodes(self, node):
        """Remove all dummy nodes from a tree.

        A dummy node is define as a node with only one child.

        This function modifies the tree specified by node.

        Parameters
        ----------
        node : TreeNode
            Root of tree to modify.
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

    def __assess_consistency(self, root):
        """Sanity check on placed taxonomic labels.

        The placed taxonomic labels should be taxonomically
        consistent from root to tip.

        Parameters
        ----------
        root : TreeNode
            Root of tree to evaluate.
        """

        print 'Consistency test started.'

        for n in root.tips():
            # walk from tip to root
            cur_rank = len(self.rank_prefixes) - 1
            for node in n.ancestors():
                if not '|' in node.name:
                    continue

                _name, ranks = node.name.split('|')
                ranks = ranks.split(';')
                highest_rank = self.rank_prefixes.index(ranks[0])

                if highest_rank > cur_rank:
                    print 'ranks: ', ranks
                    print 'cur_rank: ', cur_rank
                    print '[Error] A taxonomic inconsistency was identified.'
                    return

                cur_rank = highest_rank

        print 'Consistency test passed.'

    def _marked_parent(self, node):
        """Get first parent with a predicted rank.

        Parameters
        ----------
        node : TreeNode
            Node to search from.

        Returns
        -------
        str
            Predicted ranks or None.
        """

        parent = node.parent
        while True:
            if len(parent.children) == 1:
                # this is a dummy node inserted
                # for book keeping purposes
                pass
            else:
                if '|' in parent.name:
                    _name, ranks = parent.name.split('|')
                    return ranks

            parent = parent.parent
            if parent == None:
                break

        return None

    def _mark_rank(self, node, rank_prefix, min_support):
        """Mark node with a specific rank.

        Parameters
        ----------
        node : TreeNode
            Node to mark.
        rank_prefix : str
            Prefix of rank to mark node with.
        min_support : float
            Minimum support required to mark a node.
        """

        # do not modify labels of leaf nodes, and
        # do not mark species
        if node.is_tip():
            return

        # don't append isolates species ranks
        # (this is simply to aid in visualizing the tree in ARB)
        if rank_prefix == 'S__':
            return

        # extract name and bootstrap support
        # for node
        if '|' in node.name:
            print '[Error] Node has an invalid name: %s' % node.name
            sys.exit()

        taxa_str = None
        if ':' in node.name:
            bootstrap, taxa_str = node.name.split(':')
        else:
            try:
                bootstrap = float(node.name)
            except:
                bootstrap = 0
                taxa_str = node.name

        # only mark nodes with sufficient support
        if int(bootstrap) < min_support:
            return

        # get last predicted rank in this lineage and
        # back-fill any missing taxonomic ranks
        parent_ranks = self._marked_parent(node)
        cur_rank_index = self.rank_prefixes.index(rank_prefix)
        if parent_ranks:
            deepest_rank = parent_ranks.split(';')[-1]
            deepest_rank_index = self.rank_prefixes.index(deepest_rank)
            if deepest_rank_index < cur_rank_index:
                rank_prefix = ';'.join(self.rank_prefixes[deepest_rank_index + 1:cur_rank_index + 1])
        else:
            rank_prefix = ';'.join(self.rank_prefixes[0:cur_rank_index + 1])

        # append rank to node
        if bootstrap:
            node.name = str(bootstrap) + ':'

        if taxa_str:
            node.name += taxa_str + '|' + rank_prefix
        else:
            node.name += '|' + rank_prefix

    def mark(self, tree, start_labels, start_rank, min_support, verbose=False):
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
        min_bootstrap: int
            Minimum required support to place label on internal node.
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
                start_node.name = start_node.name + '|' + self.rank_prefixes[start_rank]
                start_node.dist_to_parent_rank = 0
                start_nodes.append(start_node)
        except:
            print 'Unable to locate node with label: ' + label
            return None

        nodes_at_cur_rank = set(start_nodes)

        # set name of unnamed nodes to an empty strings
        # to facilitate prefixing each node
        for n in root.preorder():
            if not n.name:
                n.name = ''

        # assign internal nodes with ranks from
        num_intervals = 5  # [domain, genus]
        for cur_rank in xrange(start_rank, num_intervals):
            if verbose:
                print 'Processing rank %s.' % self.rank_prefixes[cur_rank]

            # get immediate lineages below current rank
            child_lineages = set()
            for n in nodes_at_cur_rank:
                if n.dist_to_parent_rank <= 0:
                    # node is above or at the mean point for the current rank,
                    # so its children may be potential children lineages
                    for child in n.children:
                        if child not in nodes_at_cur_rank:
                            child.dist_to_parent_rank = n.dist_to_parent_rank + child.length
                            child_lineages.add(child)
                else:
                    # this node is below the current rank which should never happen!
                    print '[Error] Node is at an unexpected position relative to its parent rank.'
                    print n.name
                    return None

            # determine most basal node for each lineage
            # below the current rank
            immediate_child_lineages = set()
            for n in child_lineages:
                if n.parent not in child_lineages:
                    immediate_child_lineages.add(n)

            # mark nodes predicted to be at the current rank and determine mean
            # distance to child rank for each child lineage
            next_nodes = set()
            for child_lineage in immediate_child_lineages:
                dist_to_current_rank = child_lineage.dist_to_parent_rank

                # get mean distance from current rank to child rank
                average_dist = self._mean_dist_to_leaves(child_lineage)
                dist_to_child_rank = (dist_to_current_rank + average_dist) / (num_intervals - cur_rank)
                if verbose:
                    print child_lineage.name
                    print '  dist_to_parent_rank: ' + str(dist_to_current_rank)
                    print '  dist_to_child_rank: ' + str(dist_to_child_rank)
                    print '---'

                # mark putative nodes at the current or child rank
                nodes_to_process = [child_lineage]
                nodes_above_child_rank = set()

                while nodes_to_process:
                    n = nodes_to_process.pop()

                    if n == child_lineage:
                        dist_from_current_rank = dist_to_current_rank
                    else:
                        dist_from_current_rank = n.accumulate_to_ancestor(child_lineage) + dist_to_current_rank

                    if dist_from_current_rank <= dist_to_child_rank:
                        if dist_from_current_rank <= 0.5 * dist_to_child_rank:
                            # node belongs to the current rank
                            rank_prefix = self.rank_prefixes[cur_rank]
                        else:
                            # node belongs to the child rank
                            rank_prefix = self.rank_prefixes[cur_rank + 1]

                        self._mark_rank(n, rank_prefix, min_support)
                        n.dist_to_parent_rank = dist_from_current_rank - dist_to_child_rank

                        # process children as they may also be above the mean point of
                        # the child rank
                        nodes_above_child_rank.add(n)
                        for child in n.children:
                            nodes_to_process.append(child)

                # if required, add in a dummy node so all children lineages
                # contain a node at the child rank
                if not nodes_above_child_rank:
                    if not child_lineage.is_tip() and self.rank_prefixes[cur_rank + 1] != 'S__':
                        name = child_lineage.name + '|' + self.rank_prefixes[cur_rank + 1]
                    else:
                        name = child_lineage.name

                    dist_to_parent = dist_to_child_rank - child_lineage.parent.dist_to_parent_rank
                    dummy_node = self.__create_new_node(child_lineage, name, dist_to_parent)
                    nodes_above_child_rank.add(dummy_node)
                    dummy_node.dist_to_parent_rank = 0

                next_nodes |= nodes_above_child_rank

            nodes_at_cur_rank = next_nodes

        # remove dummy nodes
        self.__remove_dummy_nodes(root)

        # verify taxonomic consistency of placed ranks
        if True:
            self.__assess_consistency(root)

        # save marked tree
        return root
