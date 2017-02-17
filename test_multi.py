import calc_bm_likelihood as brownian
import phylo3
from scipy import optimize

tree = brownian.read_tree("0.mcc.tre")

for i in tree.iternodes(order=1):
    shifts ={}
    shifts[i] = 1
    brownian.paint_branches(tree,shifts)
    traits = brownian.read_traits("multi_traits.csv")
    a = [0.1,0.1]
    optimize.fmin_powell(brownian.calc_like_multi,a,args=(tree,traits))
    print [j.label for j in i.leaves()]

#mrca = phylo3.getMRCA(["t6","t7"],tree)
#shifts = {}
#shifts[mrca] = 1
#mrca = phylo3.getMRCA(["t5","t4"],tree)
#shifts[mrca] = 2

#brownian.paint_branches(tree,shifts)

traits = brownian.read_traits("multi_traits.csv")
#newtree = brownian.assign_sigsq(tree)

#print brownian.bm_prune(newtree,traits)
#print brownian.bm_prune(newtree,traits)
#a = [0.1,0.1]
#optimize.fmin_powell(brownian.calc_like_multi,a,args=(tree,traits))
a = [0.1]
optimize.fmin_powell(brownian.calc_like_single,a,args=(tree,traits))




