import sys,os
import tree_utils
import stratoML
import node_opt

if len(sys.argv) != 4:
    print "usage: " +sys.argv[0]+ "<file containing trees> <stratigraphic range file> <morphological alignment>"
    sys.exit(0)

trees = sys.argv[1]
strt = sys.argv[2]
aln = sys.argv[3]

#calculate likelihoods
cmd = "iqtree-omp -s "+aln+" -st MORPH -nt 2 -m MK+ASC -z " +trees
os.system(cmd)
likefl = aln+".trees"
aln_liks = {}

for i in open(likefl,"r"):
    spls = i.strip().split("]")
    nwk = spls[1]
    likespls = spls[0].strip().split("lh=")
    lik = float(likespls[1])
    aln_liks[nwk] = lik

best = None
for i in aln_liks.keys():
    tree = tree_utils.read_tree(i,nwk=True)
    ranges = tree_utils.read_strat(strt)
    tree_utils.match_strat(tree,ranges)
    tree_utils.init_heights_strat(tree)
    ranges = None 
    strat_opt = stratoML.optim_lambda_heights(tree,ranges)
    strat_heights_tree = strat_opt[0].get_newick_repr(True)
    strat_like = -float(strat_opt[1][1])
    #strat_liks[strat_heights_tree] = float(strat_like)
    combined_like = strat_like + aln_liks[i]
    if best == None:
        best = (combined_like,strat_heights_tree)
    else:
        if best[0] < combined_like:
            best = (combined_like,strat_heights_tree)
            print best

treefl = open(aln+".tre","w")
treefl.write(best[1])
treefl.close()
print best


