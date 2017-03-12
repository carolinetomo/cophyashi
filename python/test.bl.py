import calc_bm_likelihood as brownian
import phylo3
from scipy import optimize
import node_opt
import tree_utils

tree = brownian.read_tree("../test_data/hap.tre")
traits = brownian.read_traits("../test_data/hap.30.multi.csv")#traits.0.nex")#multi_traits_t6_t7_t9_t10.csv")

tree_utils.assign_sigsq(tree)
print brownian.bm_prune(tree,traits)
#brownian.tip_dates(tree,"../test_data/sim_tip_dates.csv",1)
print tree.get_newick_repr(True)
"""
for j in range(10):
    for i in tree.iternodes():
        i.length = 0.01
    a =node_opt.bm_brlen_optim(tree,traits,j/float(10)+0.1) 
    print j/float(10)+0.1,a
"""
for i in tree.iternodes():
    i.length =0.1

a =node_opt.bm_brlen_optim(tree,traits)
print a
