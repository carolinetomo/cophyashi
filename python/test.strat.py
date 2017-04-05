import sys
import tree_utils
import stratoML
import node_opt

if len(sys.argv) != 3:
    print "usage: "+sys.argv[0]+ " <newick file> <range data file>"
    sys.exit(0)

tree = tree_utils.read_tree(sys.argv[1])
ranges = tree_utils.read_strat(sys.argv[2])
tree_utils.match_strat(tree,ranges)
ranges = None

tree_utils.init_heights_strat(tree)
print tree.get_newick_repr(True)

print -stratoML.hr97_loglike(tree,lam=1.0)
#print tree.get_newick_repr(True)

opt = stratoML.optim_lambda_heights(tree,ranges)

#print opt
t1= tree.get_newick_repr(True)
t2= opt[0].get_newick_repr(True)
print opt[1][1]
print t1
