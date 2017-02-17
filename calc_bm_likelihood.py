from scipy import optimize
import math
import tree_reader
from numpy import random

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

def assign_sigsq_multi(p,tree): #params should be vector ordered like the rate classes
    for i in tree.iternodes():
        i.sigsq = p[i.rate_class]
    return tree

def calc_like_multi(params,tree,traits):
    assign_sigsq_multi(params,tree)
    try:
        val = -bm_prune(tree,traits)
    except:
        return 1000000000
    #print val,params
    return val

def calc_like_single(params,tree,traits):
    assign_sigsq_p(params,tree)
    try:
        val = -bm_prune(tree,traits)
    except:
        return 1000000000
    #print val,params
    return val

def bm_like(sigsq,cur_var,contrast):
    return ((-0.5)* ((math.log(2*math.pi*sigsq))+(math.log(cur_var))+(math.pow(contrast,2)/(sigsq*cur_var))))

def find_shifts(tree,traits,stop=2,cutoff=2):
    start = [random.uniform(0.0,2.0)]
    single = optimize.fmin_powell(calc_like_single,start,args=(tree,traits),full_output=True,disp=False)
    curlike = 0.0 
    curbest = 100000000
    best_node = None
    best_tree = None
    for i in tree.iternodes(order=1):
        if best_node == None:
            best_node = i
        shifts = {}
        shifts[i] = 1
        paint_branches(tree,shifts)
        start = [0.5,0.5]
        opt = optimize.fmin_powell(calc_like_multi,start,args=(tree,traits),full_output =True,disp=False)
        curlike = opt[1]
        if curlike < curbest:
            curbest = curlike
            best_node = i
            best_tree = tree.get_newick_repr(showbl=False,show_rate=True)
    single_aic = 2. * (1-single[1])
    two_aic = 2.*(3-curbest)
    likes=[single[1],curbest]
    aic = [single_aic,two_aic]
    for i in tree.iternodes(order=1):
        for j in tree.iternodes(order=1):
            if i == j:
                continue
            shifts={}
            shifts[i]=1
            shifts[j]=2
            paint_branches(tree,shifts)
            start = [0.5,0.5,0.5]
            opt = optimize.fmin_powell(calc_like_multi,start,args=(tree,traits),full_output= True,disp=False)
            curlike = opt[1]
            if curlike < curbest:
                curbest = curlike
                best_nodes = shifts
    likes.append(curbest)
    three_aic = 2.*(5-curbest)
    aic.append(three_aic)
    [i.label for i in best_nodes.keys()[0].leaves()]
    return [i.label for i in best_nodes.keys()[1].leaves()]#,curbest,best_tree]



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
