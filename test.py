import calc_bm_likelihood as brownian
from scipy import optimize

tree = brownian.read_tree("./test_data/hap.tre")
traits = brownian.read_traits("./test_data/multi_traits_primates.csv") #traits.0.nex")
tree = brownian.read_tree("./test_data/0.mcc.tre")
traits = brownian.read_traits("./test_data/traits.0.nex")

newtree = brownian.assign_sigsq(tree)

print brownian.bm_prune(newtree,traits)
#print brownian.bm_prune(newtree,traits)
brownian.assign_node_nums(tree)
nnodes=0
tree
for i in tree.iternodes():
    if i.istip or i.parent == None:
        continue
    else:
        nnodes+=1

#brownian.tip_dates(tree,"./hap.tip_dates.csv",34.0)
brownian.tip_dates(tree,"./test_data/sim_tip_dates.csv",1.0)

brownian.init_heights(tree)

nstart = [i.height for i in tree.iternodes(order=0) if i.istip==False and i.parent!=None]

a = [0.1]
a = a +nstart
#brownian.assign_node_heights(a[1:],tree)
#print tree.get_newick_repr(True)
print a
#optimize.fmin_powell(brownian.calc_like_nodes,(nstart),args=(tree,traits))
res = optimize.fmin_powell(brownian.calc_like_single,(a),args=(tree,traits))
print res
brownian.assign_node_heights(res[1:],tree)
#print [i.height for i in tree.iternodes() if i.istip==False]
print tree.get_newick_repr(True)
