import math
from scipy import optimize
import tree_reader
import tree_utils
import calc_bm_likelihood

LARGE = 100000000000

def bm_height_optim(tree,traits,nrates=False):
    if nrates != False:
        start = [random.uniform(0.000001,2.0)]*nrates
        tree_utils.assign_node_nums(tree)
        tree_utils.init_heights(tree,start)
        tree_utils.assign_sigsq(tree,start)
        nstart = [i.height for i in tree.iternodes() if i.istip == False and i.parent!=None]
        start = start + nstart
        #single = optimize.fmin_bfgs(calc_like_nodes,start,args=(tree,traits,nrates),full_output=True,disp=True)
        bounds = [(0.0,100000.0)]*len(start)
        opt = optimize.fmin_l_bfgs_b(calc_like_sigsq_nodes,start,approx_grad = True,bounds =bounds,args=(tree,traits,nrates))
        aic1 = 2. * (1+single[1])
        aic[aic1]= tree.get_newick_repr(True)
        assign_node_heights(single[0][1:],tree)
        return opt 
    elif nrates == False: ##estimate node heights with fixed rate
        tree_utils.assign_node_nums(tree)
        tree_utils.init_heights(tree)
        tree_utils.assign_sigsq(tree)
        start = [i.height for i in tree.iternodes() if i.istip == False and i.parent!=None]
        print start
        opt = optimize.fmin_bfgs(calc_like_nodes,start,args=(tree,traits),full_output=True,disp=True)
        #print tree.get_newick_repr(True)
        return opt

def calc_like_sigsq_nodes(ht,tree,traits,nrates):
    for i in ht:
        if i < 0:
            return LARGE
    z = 0
    while z < nrates:
        if ht[z] > 500:
            return LARGE
        z += 1
    if nrates == 1:
        tree_utils.assign_sigsq(tree,ht[0])
    elif nrates > 1:
        tree_utils.assign_sigsq(tree,ht[0:nrates])
    bad = tree_utils.assign_node_heights(ht[nrates:],tree)
    if bad:
        return LARGE
    try:
        val = -calc_bm_likelihood.bm_prune(tree,traits)
    except:
        return LARGE
    #print ht[0] 
    #print (val,ht)
    return val


def calc_like_nodes(ht,tree,traits):
    for i in ht:
        if i < 0:
            return LARGE
    bad = tree_utils.assign_node_heights(ht,tree)
    #print [i.height for i in tree.iternodes()] 
    if bad:
        return LARGE
    try:
        val = -calc_bm_likelihood.bm_prune(tree,traits)
    except:
        return LARGE
    #print ht[0] 
    #print (val,ht)
    return val


