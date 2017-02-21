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
        newtree = match_traits_tips(tree,traits,i)
        for j in newtree.iternodes(order=1):
            if len(j.children) > 0: # do internal nodes only  
                child_charst = [k.charst for k in j.children]
                brlens = [k.length for k in j.children]
                contrast = child_charst[0]-child_charst[1]
                cur_var = brlens[0]+brlens[1]
                temp_charst = (((1/brlens[0])*child_charst[0])+((1/brlens[1])*child_charst[1]))/((1/brlens[0])+(1/brlens[1]))
                #temp_charst = (((brlens[1]/(sum(brlens)))*child_charst[0]))+(((brlens[0]/(sum(brlens))))*child_charst[1])
                temp_brlen = j.length+((brlens[0]*brlens[1])/(brlens[0]+brlens[1]))
                cur_like =((-0.5)* ((math.log(2*math.pi*j.sigsq))+(math.log(cur_var))+(math.pow(contrast,2)/(j.sigsq*cur_var))))
                node_likes.append(cur_like)
                [k.remove_child for k in j.children]
                j.charst = temp_charst
                j.length = temp_brlen
        for j in newtree.iternodes():
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
    return tree

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


def init_heights(tree):
    for i in tree.iternodes(order=1):
        if i.istip or i.parent == None:
            continue
        else:
            temp = 0.0
            n=i.descendants("POSTORDER")[0]
            temp += n.height
            while n!=i:
                temp += n.length
                n = n.parent
            i.height = temp

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
                return True
    return False

def assign_sigsq_multi(p,tree): #params should be vector ordered like the rate classes
    for i in tree.iternodes():
        i.sigsq = p[i.rate_class]
    return tree

def calc_like_nodes(ht,tree,traits,nrates):
    for i in ht:
        if i < 0:
            return LARGE
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
    #sys.exit(0)
    #print val,ht 
    return val

def calc_like_multi(params,tree,traits):
    assign_sigsq_multi(params,tree)
    try:
        val = -bm_prune(tree,traits)
    except:
        return LARGE
    #print val,params
    return val

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

def find_shifts(tree,traits,stop=2,aic_cutoff=2,opt_nodes=True):
    start = [random.uniform(0.0,2.0)]
    if opt_nodes == False:
        single = optimize.fmin_powell(calc_like_single,start,args=(tree,traits),full_output=True,disp=False)
    elif opt_nodes == True:
        init_heights(tree)
        nrates = 1
        nstart = [i.height for i in tree.iternodes() if i.istip == False and i.parent!=None]
        start = start + nstart
        single = optimize.fmin_powell(calc_like_nodes,start,args=(tree,traits,nrates),full_output=True,disp=False)
    print single;sys.exit(0)
    curlike = 0.0
    curbest = 100000000
    best_node = None
    best_tree = None
    print [i.height for i in tree.iternodes() if i.istip == False or i == tree]
    for i in tree.iternodes(order=1):
        if best_node == None:
            best_node = i
        shifts = {}
        shifts[i] = 1
        paint_branches(tree,shifts)
        start = [random.uniform(0.0,2.0),random.uniform(0.0,2.0)]
        if opt_nodes == False:
            opt = optimize.fmin_powell(calc_like_multi,start,args=(tree,traits),full_output =True,disp=False)
        elif opt_nodes == True:
            nrates = 2
            nstart = [i.height for i in tree.iternodes() if i.istip == False and i.parent!=None]
            start = start + nstart
            opt = optimize.fmin_powell(calc_like_nodes,start,args=(tree,traits,nrates),full_output =True,disp=False)
        curlike = opt[1]
        if curlike < curbest:
            curbest = curlike
            best_node = i
            best_tree2 = tree.get_newick_repr(showbl=False,show_rate=True)
            opt2 = opt[0]
    single_aic = 2. * (1+single[1])
    two_aic = 2.*(3+curbest)
    likes=[single[1],curbest]
    aic = [single_aic,two_aic]
    """
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
                opt = optimize.fmin_powell(calc_like_multi,start,args=(tree,traits),full_output= True,disp=False)
            elif opt_nodes == True:
                nrates = 3
                nstart = [i.height for i in tree.iternodes() if i.istip == False and i.parent!=None]
                start = start+nstart
                opt = optimize.fmin_powell(calc_like_nodes,start,args=(tree,traits,nrates),full_output=True,disp=False)
            curlike = opt[1]
            if curlike < curbest:
                curbest = curlike
                best_nodes = shifts
                best_tree3 = tree.get_newick_repr(showbl=False,show_rate=True)
                opt3 = opt[0]
    likes.append(curbest)
    three_aic = 2.*(5+curbest)
    aic.append(three_aic)
    """
    sm =1000000000000.
    for i in aic:
        if i < sm and abs(sm-i) >= aic_cutoff:
            sm = i
    return [aic,opt2]#,opt3]

def match_traits_tips(tree,traits,number):
    for i in tree.leaves():
        i.charst = traits[i.label][number]
    return tree

def read_tree(treefl):
    nwk = open(treefl,"r").readlines()[0].strip()
    return tree_reader.read_tree_string(nwk)

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
