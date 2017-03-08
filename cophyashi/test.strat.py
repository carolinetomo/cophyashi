import tree_utils
import stratoML

tree = tree_utils.read_tree("../test_data/fossils_strat.tre")
ranges = tree_utils.read_strat("../test_data/fossils_strat.txt")
tree_utils.match_strat(tree,ranges)
tree_utils.tip_dates(tree,"../test_data/fossils_strat.tips.txt",27)
tree_utils.init_heights(tree)

print stratoML.hr97_loglike(tree,lam=0.5)

stratoML.optim_lambda_heights(tree,ranges)

