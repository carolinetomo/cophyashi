import calc_bm_likelihood as brownian
import phylo3
from scipy import optimize
import sys

opth=True
z = "MEDUSA"
s = 1
tree = brownian.read_tree("./test_data/2tips.tre")
traits = brownian.read_traits("./test_data/2tips_traits.csv")

brownian.tip_dates(tree,"./test_data/2tips_heights.csv",float(sys.argv[1]))

brownian.init_heights(tree,[0.9410645])
"""
rng = 100.
for i in range(1,int(rng)+1):
    brownian.init_heights(tree,[i/10.])
    print i/10.,brownian.bm_prune(tree,traits)
print tree.get_newick_repr(showbl=True)
"""
print brownian.bm_prune(tree,traits)
#a = brownian.find_shifts(tree, traits,opt_nodes=opth,search=z,stop=s)
#print a



