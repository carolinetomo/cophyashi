import tree_utils
import stratoML
import calc_bm_likelihood as brownian
import node_opt

tree = tree_utils.read_tree("../test_data/fossils_strat.tre")
ranges = tree_utils.read_strat("../test_data/fossils_strat.txt")
tree_utils.match_strat(tree,ranges)

traits = tree_utils.read_traits("../test_data/paleo_BM.csv")

tree_utils.tip_dates(tree,"../test_data/fossils_strat.tips.txt",27)
tree_utils.assign_sigsq(tree,[0.7])

print brownian.bm_prune(tree,traits)
#print stratoML.hr97_loglike(tree,lam=1.0)
print tree.get_newick_repr(True)

tree_utils.init_heights(tree,strat = True)

print -brownian.bm_prune(tree,traits)
print -stratoML.hr97_loglike(tree,lam=1.0)
print tree.get_newick_repr(True)

#stratoML.optim_lambda_heights(tree,ranges)


opt = node_opt.optim_strat_bm(tree,ranges,traits)
print opt
t1= tree.get_newick_repr(True)
t2= opt[0].get_newick_repr(True)
print t1
