import calc_bm_likelihood as brownian
from scipy import optimize

def calc_like(params,tree,traits):
    brownian.assign_sigsq_p(params[0],tree)
    try:
        val = -brownian.bm_prune(tree,traits)
    except:
        return 1000000
    print val,params
    return val

tree = brownian.read_tree("test_data/0.mcc.tre")
traits = brownian.read_traits("./test_data/multi_traits_t6_t7_t9_t10.csv") #traits.0.nex")
newtree = brownian.assign_sigsq(tree)

print brownian.bm_prune(newtree,traits)
print brownian.bm_prune(newtree,traits)
a = [0.1]
optimize.fmin_powell(brownian.calc_like_single,a,args=(tree,traits))


