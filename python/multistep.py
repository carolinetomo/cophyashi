from calc_bm_likelihood import *
import node_opt
import tree_utils
import math

def paint_branches(tree,shift_nodes): #shift_nodes should be a dictionary with nodes as keys and rate regime as value
    nodes = {}
    for i in shift_nodes.keys():
        nodes[i]=i.descendants("PREORDER")
    curshift = None
    for i in tree.iternodes(order = 0):
        if i in shift_nodes.keys():
            curshift = i
            i.rate_class = shift_nodes[i]
        elif curshift == None:
            i.rate_class = 0
            continue
        elif i in nodes[curshift]:
            i.rate_class = shift_nodes[curshift]
        else:
            i.rate_class = 0
    return tree

def pop_dict(tree):
    d = {}
    for i in tree.iternodes():
        d[i] = (i.height,i.length,i.sigsq)
    return d

def calc_like_single(params,tree,traits):
    assign_sigsq_p(params[0],tree)
    try:
        val = -bm_prune(tree,traits)
    except:
        return LARGE
    #print val,params
    return val

def find_shifts(tree,traits,stop=2,aic_cutoff=4,opt_nodes=True,search="MEDUSA",est_sigsq=True):
    aic = {}
    nrates = 1
    if opt_nodes == False:
        single = optimize.fmin_bfgs(calc_like_single,start,args=(tree,traits),full_output=True,disp=False)
    elif opt_nodes == True:
        single = node_opt.bm_height_optim(tree,traits,nrates)
        aic1 = node_opt.tree_aic(nrates,single[1]) 
        aic[aic1]= tree.get_newick_repr(True)
    if nrates == stop:
        return aic
    curlike = 0.0
    curbest = LARGE
    best_node = None
    best_tree = None
    best_tree_obj=None
    best = {}
    for i in tree.iternodes(order=1):
        if best_node == None:
            best_node = i
        shifts = {}
        shifts[i] = 1
        paint_branches(tree,shifts)
        nrates = 2
        if opt_nodes == False:
            opt = optimize.fmin_bfgs(calc_like_multi,start,args=(tree,traits),full_output =True,disp=False)
        elif opt_nodes == True:
            opt = node_opt.bm_height_optim(tree,traits,nrates)
        curlike = opt[1]
        opt2= None
        if curlike < curbest:
            curbest = curlike
            best_node2 = i
            best_tree2 = tree.get_newick_repr(showbl=True,show_rate=False)
            best_nh = pop_dict(tree)
            best = pop_dict(tree)
            opt2 = opt[0]
        else:
            for j in tree.iternodes():
                tup = best[j]
                j.height = tup[0]
                j.length = tup[1]
                j.sigsq = tup[2]
    aic2 = 2.*(3+curbest)
    aic[aic2] = best_tree2
    likes=[single[1],curbest]
    if nrates == stop:
        print aic.keys()
        return aic
    curbest= LARGE
    nrates = 3
    if search == "FULL":
        for i in tree.iternodes(order=1):
            for j in tree.iternodes(order=1):
                if i == j:
                    continue
                shifts={}
                shifts[i]=1
                shifts[j]=2
                paint_branches(tree,shifts)
                if opt_nodes == False:
                    opt = optimize.fmin_bfgs(calc_like_multi,start,args=(tree,traits),full_output= True,disp=False)
                elif opt_nodes == True:
                    nrates = 3
                    opt = node_opt.bm_height_optim(tree,traits,nrates)
                curlike = opt[1]
                if curlike < curbest:
                    curbest = curlike
                    best_nodes = shifts
                    best_tree3 = tree.get_newick_repr(showbl=False,show_rate=True)
                    best_tree_obj3 = tree
                    opt3 = opt[0]
                else:
                    tree = best_tree_obj3
    elif search == "MEDUSA": #take the best single shift and tries to add another 
        for i in tree.iternodes(order = 1):
            if i == best_node2:
                continue
            shifts = {}
            shifts[best_node2]=1
            shifts[i] = 2
            paint_branches(tree,shifts)
            start = [random.uniform(0.0,2.0),random.uniform(0.0,2.0),random.uniform(0.0,2.0)]
            if opt_nodes == False:
                opt = optimize.fmin_bfgs(calc_like_multi,start,args=(tree,traits),full_output= True,disp=False)
            elif opt_nodes == True:
                nrates = 3
                opt = node_opt.bm_height_optim(tree,traits,nrates)
            curlike = opt[1]
            if curlike < curbest:
                curbest = curlike
                best_node3 = i
                best_tree3 = tree.get_newick_repr(showbl=False,show_rate=True)
                best_tree_obj3 = tree
                opt3=opt[0]
                best = pop_dict(tree)
            else:
                tup = best[j]
                j.height = tup[0]
                j.length = tup[1]
                j.sigsq = tup[2]
    likes.append(curbest)
    aic3 = 2.*(5+curbest)
    aic[aic3] = best_tree3
    sm = 1000000000000.
    for i in aic:
        if i < sm and abs(sm-i) >= aic_cutoff:
            sm = i
    print aic.keys()
    return aic[sm]


