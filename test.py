import calc_bm_likelihood as brownian
from scipy import optimize

tree = brownian.read_tree("test_data/0.mcc.tre")
traits = brownian.read_traits("./test_data/traits.0.nex") #traits.0.nex")
newtree = brownian.assign_sigsq(tree)

print brownian.bm_prune(newtree,traits)
#print brownian.bm_prune(newtree,traits)
#a = [0.1]
#optimize.fmin_powell(brownian.calc_like_single,a,args=(tree,traits))


