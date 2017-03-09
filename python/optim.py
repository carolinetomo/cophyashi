import math
from scipy import optimize
import tree_reader
import tree_utils
import calc_bm_likelihood as bm_like

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
        assign_node_heights(single[0][1:],tree)
        return opt 
    elif nrates == False: ##estimate node heights with fixed rate
        tree.utils. assign_node_nums(tree)
        tree_utils.init_heights(tree)
        tree_utils.assign_sigsq(tree)
        start = [i.height for i in tree.iternodes() if i.istip == False and i.parent!=None]
        opt = optimize.fmin_bfgs(calc_like_nodes,start,args=(tree,traits),full_output=True,disp=True)
        assign_node_heights(single[0][1:],tree)
        return opt

def tree_aic(num_p,loglike):
    aic = 2. * (num_p+loglike)
    return aic

def calc_like_nodes(ht,tree,traits):
    for i in ht:
        if i < 0:
            return LARGE
    z = 0
    while z < nrates:
        if ht[z] > 500:
            return LARGE
        z += 1
    bad = assign_node_heights(ht,tree)
    if bad:
        return LARGE
    try:
        val = -bm_like.bm_prune(tree,traits)
    except:
        return LARGE
    #print ht[0] 
    #print (val,ht)
    return val



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
    bad = assign_node_heights(ht[nrates:],tree)
    if bad:
        return LARGE
    try:
        val = -bm_like.bm_prune(tree,traits)
    except:
        return LARGE
    #print ht[0] 
    #print (val,ht)
    return val


