import calc_bm_likelihood as brownian
import phylo3
from scipy import optimize

tree = brownian.read_tree("test_data/hap.tre")
traits = brownian.read_traits("./test_data/multi_traits_primates.csv")#multi_traits_t6_t7_t9_t10.csv")

#for i in tree.iternodes(order=1):
#    shifts ={}
#    shifts[i] = 1
#    brownian.paint_branches(tree,shifts)
#    traits = brownian.read_traits("./test_data/multi_traits_t6_t7_t9_t10.csv")
#    a = [0.1,0.1]
#    opt = optimize.fmin_powell(brownian.calc_like_multi,a,args=(tree,traits),full_output = True)
#    print opt#[0]
#    print [j.label for j in i.leaves()]

a = brownian.find_shifts(tree, traits)
print a
#print [[i.label,i.rate_class] for i in a[1].iternodes()] 

#    print [j.label for j in i.leaves()]

#mrca = phylo3.getMRCA(["t6","t7"],tree)
#shifts = {}
#shifts[mrca] = 1
#mrca = phylo3.getMRCA(["t5","t4"],tree)
#shifts[mrca] = 2

#brownian.paint_branches(tree,shifts)

#traits = brownian.read_traits("multi_traits.csv")
#newtree = brownian.assign_sigsq(tree)

#print brownian.bm_prune(newtree,traits)
#print brownian.bm_prune(newtree,traits)
#a = [0.1,0.1]
#optimize.fmin_powell(brownian.calc_like_multi,a,args=(tree,traits))
#a = [0.1]
#optimize.fmin_powell(brownian.calc_like_single,a,args=(tree,traits))




