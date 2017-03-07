import tree_utils
import stratoML

tree = tree_utils.read_tree("test_data/fossils_strat.tre")
ranges = tree_utils.read_strat("test_data/fossils_strat.txt")
tree_utils.match_strat(tree,ranges)

print stratoML.strat_like(tree,ranges,lam=0.5)
