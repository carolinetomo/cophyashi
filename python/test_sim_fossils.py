import calc_bm_likelihood as brownian
import phylo3
from scipy import optimize
import sys
import node_opt
import tree_utils

opth=True
z = "MEDUSA"
s = 1
tree = tree_utils.read_tree("../test_data/2tips.tre") 
traits = tree_utils.read_traits("../test_data/2tips_traits.csv")
tree_utils.assign_sigsq(tree)
tree_utils.tip_dates(tree,"../test_data/2tips_heights.csv",float(sys.argv[1]))
print brownian.bm_prune(tree,traits)

tree_utils.init_heights(tree)
rng = 100.
for i in range(1,int(rng)+1):
    brownian.init_heights(tree,[i/10.])
    print i/10.,brownian.bm_prune(tree,traits)
print tree.get_newick_repr(showbl=True)

print brownian.bm_prune(tree,traits)
#a = node_opt.bm_height_optim(tree,traits)
#a = brownian.find_shifts(tree, traits,opt_nodes=opth,search=z,stop=s)
#print a



