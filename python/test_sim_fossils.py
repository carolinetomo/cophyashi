import calc_bm_likelihood as brownian
import phylo3
from scipy import optimize
import sys
import node_opt
import tree_utils

opth=True
z = "MEDUSA"
s = 1
tree = tree_utils.read_tree("../test_data/paleo_sim.tre") 
traits = tree_utils.read_traits("../test_data/paleo_BM.csv")
tree_utils.assign_sigsq(tree,[0.7])
tree_utils.tip_dates(tree,"../test_data/paleo_heights.csv",float(sys.argv[1]))
print brownian.bm_prune(tree,traits)
print tree.get_newick_repr(True)
tree_utils.init_heights(tree)
"""
rng = 100.
for i in range(1,int(rng)+1):
    if i/10. <0.51:
        continue
    tree_utils.tip_dates(tree,"../test_data/3tips_heights.csv",i/10.)
    brownian.init_heights(tree)
    print i/10.,brownian.bm_prune(tree,traits)
"""
#print tree.get_newick_repr(showbl=True)

print brownian.bm_prune(tree,traits)
print tree.get_newick_repr(True)
print [i.sigsq for i in tree.iternodes()]
a = node_opt.bm_height_optim(tree,traits,1) #a = brownian.find_shifts(tree, traits,opt_nodes=opth,search=z,stop=s)
print a



