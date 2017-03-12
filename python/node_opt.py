import math
from scipy import optimize
import tree_reader
import tree_utils
import calc_bm_likelihood
import random
import stratoML

LARGE = 100000000000

def bm_height_optim(tree,traits,nrates=False):
    if nrates != False:

        start = [random.uniform(0.000001,2.0)]*nrates
        tree_utils.assign_node_nums(tree)
        tree_utils.init_heights(tree,start)
        tree_utils.assign_sigsq(tree,start)
        nstart = [i.height for i in tree.iternodes() if i.istip == False and i.parent!=None]
        start = start + nstart
        opt = optimize.fmin_bfgs(calc_like_sigsq_nodes,start,args=(tree,traits,nrates),full_output=True,disp=True)
        #bounds = [(0.0,100000.0)]*len(start)
        #opt = optimize.fmin_l_bfgs_b(calc_like_sigsq_nodes,start,approx_grad = True,bounds =bounds,args=(tree,traits,nrates))
        tree_utils.assign_node_heights(opt[0][1:],tree)
        return [tree.get_newick_repr(True),opt]
    elif nrates == False: ##estimate node heights with fixed rate
        tree_utils.assign_node_nums(tree)
        tree_utils.init_heights(tree)
        sigsq_i = calc_bm_likelihood.sigsqML(tree)
        tree_utils.assign_sigsq(tree,[sigsq_i])
        start = [i.height for i in tree.iternodes() if i.istip == False]# and i.parent!=None]
        #print len(start)
        opt = optimize.fmin_bfgs(calc_like_nodes,start,args=(tree,traits),full_output=True,disp=True)
        return [tree.get_newick_repr(True),opt]

def bm_brlen_optim(tree,traits,rate=1):
    tree_utils.assign_branch_nums(tree)
    tree_utils.assign_sigsq(tree,[rate])
    start = [i.length for i in tree.iternodes() if i != tree]
    opt = optimize.fmin_bfgs(calc_like_brlens, start, args = (tree,traits),full_output=False,disp=True)
    return [tree.get_newick_repr(True),opt]

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
        tree_utils.assign_sigsq(tree,[ht[0]])
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

def calc_like_brlens(l,tree,traits):
    for i in l:
        if i<0:
            return LARGE
        bad = tree_utils.assign_brlens(l,tree)
        if bad:
            return LARGE
        try:
            ll = -calc_bm_likelihood.bm_prune(tree,traits)
        except:
            return LARGE
        return ll


def calc_like_nodes(ht,tree,traits):
    for i in ht:
        if i < 0:
            return LARGE
    bad = tree_utils.assign_node_heights(ht,tree)
    #print [i.height for i in tree.iternodes()] 
    if bad:
        return LARGE
    try:
        ll = -calc_bm_likelihood.bm_prune(tree,traits)
    except:
        return LARGE
    #print ht[0] 
    #print (val,ht)
    return ll

def calc_like_strat_bm(p,tree,strat,traits,nrates=1): #p[0] should be lambda, p[1:nrates] = sig2
    for i in p:
        if i < 0:
            return LARGE
    hstart = nrates+1
    if nrates == 1:
        tree_utils.assign_sigsq(tree,[p[1]])
    elif nrates > 1:
        tree_utils.assign_sigsq(tree,p[1:nrates])
    bad = tree_utils.assign_node_heights(p[hstart:],tree)
    if bad:
        return LARGE
    try:
        sll = -stratoML.hr97_loglike(tree,p[0])
        bmll = -calc_bm_likelihood.bm_prune(tree,traits)
        ll = sll+bmll
    except:
        return LARGE
    return ll

def optim_strat_bm(tree,strat,traits,nrates = 1):
    hstart = nrates+1
    tree_utils.assign_node_nums(tree)
    lamstart = [random.uniform(0.1,1)]
    sigstart = [random.uniform(0.01,2) for i in range(nrates)]
    nstart = [i.height for i in tree.iternodes(order = 0) if i.istip == False and i != tree]
    start = lamstart+sigstart+nstart
    #opt = optimize.fmin_bfgs(calc_like_strat_bm,start,args =(tree,strat,traits),full_output=True,disp=True)
    bounds = [(0.0,100000.0)]*len(start)
    opt = optimize.fmin_l_bfgs_b(calc_like_strat_bm,start,approx_grad = True,bounds =bounds,args=(tree,strat,traits))
    tree_utils.assign_node_heights(opt[0][hstart:],tree)
    return [tree,opt]







