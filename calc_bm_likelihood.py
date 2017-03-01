from scipy import optimize
import math
import tree_reader
from numpy import random
import sys

LARGE = 100000000

def bm_prune(tree,traits):
    trait_likes = []
    for i in range(len(traits.values()[0])):
        node_likes = []
        match_traits_tips(tree,traits,i)
        for j in tree.iternodes():
            j.old_length = j.length
        for j in tree.iternodes(order=1):
            if len(j.children) > 0: # do internal nodes only  
                child_charst = [k.charst for k in j.children]
                brlens = [k.length for k in j.children]
                contrast = child_charst[0]-child_charst[1]
                cur_var = brlens[0]+brlens[1]
                #print j.sigsq
                cur_like =((-0.5)* ((math.log(2*math.pi*j.sigsq))+(math.log(cur_var))+(math.pow(contrast,2)/(j.sigsq*cur_var))))
                node_likes.append(cur_like)
                temp_charst = (((1/brlens[0])*child_charst[0])+((1/brlens[1])*child_charst[1]))/((1/brlens[0])+(1/brlens[1]))
                #temp_charst = (((brlens[1]/(sum(brlens)))*child_charst[0]))+(((brlens[0]/(sum(brlens))))*child_charst[1])
                temp_brlen = j.length+((brlens[0]*brlens[1])/(brlens[0]+brlens[1]))
                [k.remove_child for k in j.children]
                j.charst = temp_charst
                j.length = temp_brlen
        for j in tree.iternodes():
            j.length = j.old_length
        trait_likes.append(sum(node_likes))
        #return(sum(node_likes))
    return sum(trait_likes)

def tip_dates(tree,dates,root_height):
    d = {}
    for i in open(dates,"r"):
        spls = i.strip().split()
        d[spls[0]] = float(spls[1])
    for i in tree.iternodes():
        if i.istip == True:
            i.height = d[i.label]
        elif i.parent==None:
            i.height = root_height

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

def assign_sigsq_p(p,tree):
    for i in tree.iternodes():
        i.sigsq = p
    return tree

def assign_node_nums(tree):
    num = 0
    for i in tree.iternodes(order=0):
        if i.istip or i == tree:
            continue
        else:
            i.number = num
            num += 1
    return tree


def init_heights(tree,sigsq):
    for i in tree.iternodes(order=1):
        if i.istip or i.parent == None:
            if len(sigsq)==1:
                i.sigsq = sigsq[0]
            elif len(sigsq) > 1:
                i.sigsq = sigsq[i.rate_class]
            continue
        i.height = max([j.height for j in i.children])+0.01
        if len(sigsq)==1:
            i.sigsq = sigsq[0]
        elif len(sigsq) > 1:
            i.sigsq = sigsq[i.rate_class]
    for i in tree.iternodes():
        if i == tree:
            continue
        i.length = i.parent.height-i.height

def assign_node_heights(h,tree):
    for i in tree.iternodes(order=0):
        if i.istip:
            i.length = i.parent.height-i.height
            if i.length < 0:
                return True
        elif i.parent == None:
            i.length = 2.
        else:
            i.height = h[i.number]
            i.length = i.parent.height-i.height
            if i.length < 0:
                #print i.get_newick_repr(False),i.number,h[i.number],i.parent.height,i.length
                return True
    return False

def assign_sigsq_multi(p,tree): #params should be vector ordered like the rate classes
    for i in tree.iternodes():
        i.sigsq = p[i.rate_class]

def calc_like_nodes(ht,tree,traits,nrates):
    for i in ht:
        if i < 0:
            return LARGE
    z = 0
    while z < nrates:
        if ht[z] > 500:
            return LARGE
        z += 1
    if nrates == 1:
        assign_sigsq_p(ht[0],tree)
    elif nrates > 1:
        assign_sigsq_multi(ht[0:nrates],tree)
    bad = assign_node_heights(ht[nrates:],tree)
    if bad:
        return LARGE
    try:
        val = -bm_prune(tree,traits)
    except:
        return LARGE
    #print (val,ht)
    return val

def calc_like_multi(params,tree,traits):
    assign_sigsq_multi(params,tree)
    try:
        val = -bm_prune(tree,traits)
    except:
        return LARGE
    #print val,params
    return val

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

def bm_like(sigsq,cur_var,contrast):
    return ((-0.5)* ((math.log(2*math.pi*sigsq))+(math.log(cur_var))+(math.pow(contrast,2)/(sigsq*cur_var))))

def find_shifts(tree,traits,stop=2,aic_cutoff=4,opt_nodes=True,search="MEDUSA"):
    start = [random.uniform(0.0,2.0)]
    aic = {}
    nrates = 1
    if opt_nodes == False:
        single = optimize.fmin_bfgs(calc_like_single,start,args=(tree,traits),full_output=True,disp=False)
    elif opt_nodes == True:
        assign_node_nums(tree)
        init_heights(tree,start)
        nstart = [i.height for i in tree.iternodes() if i.istip == False and i.parent!=None]
        start = start + nstart
        single = optimize.fmin_bfgs(calc_like_nodes,start,args=(tree,traits,nrates),full_output=True,disp=True)
        #bounds = [(0.0,100000.0)]*len(start)
        #single = optimize.fmin_l_bfgs_b(calc_like_nodes,start,approx_grad = True,bounds =bounds,args=(tree,traits,nrates))
        aic1 = 2. * (1+single[1])
        aic[aic1]= tree.get_newick_repr(True)
        assign_node_heights(single[0][1:],tree)
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
        rand = random.uniform(0.005,2.0)
        rand2 = random.uniform(0.005,2.0)
        start = [rand,rand2]
        nrates = 2
        if opt_nodes == False:
            opt = optimize.fmin_bfgs(calc_like_multi,start,args=(tree,traits),full_output =True,disp=False)
        elif opt_nodes == True:
            nstart = [j.height for j in tree.iternodes() if j.istip == False and j.parent!=None]
            start = start + nstart
            opt = optimize.fmin_bfgs(calc_like_nodes,start,args=(tree,traits,nrates),full_output =True,disp=False)
            #bounds = [(0.0,100000.0)]*len(start)
            #opt = optimize.fmin_l_bfgs_b(calc_like_nodes,start,approx_grad = True,bounds = bounds,args=(tree,traits,nrates))
            assign_node_heights(opt[0][nrates:],tree)
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
                start = [random.uniform(0.0,2.0),random.uniform(0.0,2.0),random.uniform(0.0,2.0)]
                if opt_nodes == False:
                    opt = optimize.fmin_bfgs(calc_like_multi,start,args=(tree,traits),full_output= True,disp=False)
                elif opt_nodes == True:
                    nrates = 3
                    nstart = [i.height for i in tree.iternodes() if i.istip == False and i.parent!=None]
                    start = start+nstart
                    opt = optimize.fmin_bfgs(calc_like_nodes,start,args=(tree,traits,nrates),full_output=True,disp=False)
                curlike = opt[1]
                if curlike < curbest:
                    curbest = curlike
                    best_nodes = shifts
                    best_tree3 = tree.get_newick_repr(showbl=False,show_rate=True)
                    best_tree_obj3 = tree
                    opt3 = opt[0]
                else:
                    tree = best_tree_obj3
    elif search == "MEDUSA": #this takes the best single shift and tries to add another 
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
                nstart = [i.height for i in tree.iternodes() if i.istip == False and i.parent!=None]
                start = start+nstart
                opt = optimize.fmin_bfgs(calc_like_nodes,start,args=(tree,traits,nrates),full_output=True,disp=False)
            curlike = opt[1]
            if curlike < curbest:
                curbest = curlike
                best_node3 = i
                best_tree3 = tree.get_newick_repr(showbl=False,show_rate=True)
                best_tree_obj3 = tree
                opt3=opt[0]
    likes.append(curbest)
    aic3 = 2.*(5+curbest)
    aic[aic3] = best_tree3
    sm = 1000000000000.
    for i in aic:
        if i < sm and abs(sm-i) >= aic_cutoff:
            sm = i
    print aic.keys()
    return aic[sm]

def match_traits_tips(tree,traits,number):
    for i in tree.leaves():
        i.charst = traits[i.label][number]
    return tree

def read_tree(treefl):
    nwk = open(treefl,"r").readlines()[0].strip()
    tree = tree_reader.read_tree_string(nwk)
    for i in tree.iternodes():
        i.old_length = i.length
    return tree

def read_traits(traitfl): #should be tab separated
    traits = {}
    for i in open(traitfl,"r"):
        if "\t" not in i:
            continue
        spls = i.strip().split("\t")
        nm = spls[0]
        traitls = [float(j) for j in spls[1:]]
        traits[nm]= traitls
    return traits

def assign_sigsq(tree):
    for i in tree.iternodes():
        i.sigsq = 0.36878333
    return tree
