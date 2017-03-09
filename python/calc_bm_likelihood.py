from scipy import optimize
import math
import tree_reader
from tree_utils import *
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
            if j.istip == False:# and j.parent != None: # do internal nodes only  
                child_charst = [k.charst for k in j.children]
                brlens = [k.length for k in j.children]
                contrast = child_charst[0]-child_charst[1]
                cur_var = brlens[0]+brlens[1]
                #x = ((brlens[1]*child_charst[0])+(brlens[0]*child_charst[1]))/(sum(brlens))
                curlike =((-0.5)* ((math.log(2*math.pi*j.sigsq))+(math.log(cur_var))+(math.pow(contrast,2)/(j.sigsq*cur_var))))
                node_likes.append(curlike)
                temp_charst = (((1/brlens[0])*child_charst[0])+((1/brlens[1])*child_charst[1]))/((1/brlens[0])+(1/brlens[1]))
                #temp_charst = ((brlens[1]*child_charst[0])+(brlens[0]*child_charst[1]))/(sum(brlens))
                temp_brlen = j.length+((brlens[0]*brlens[1])/(brlens[0]+brlens[1]))
                #[k.remove_child for k in j.children]
                j.charst = temp_charst
                j.length = temp_brlen
        for j in tree.iternodes():
            j.length = j.old_length
        #first = (tree.nnodes("tips")*math.log(2*math.pi))
        #Ly = (first+sum(node_likes))*(-0.5)
        trait_likes.append(sum(node_likes))
    #print tree.get_newick_repr(True)
    #print -sum(trait_likes)
    return sum(trait_likes)

def sigsqML(tree): #tree must already have characters mapped to tips using match_traits_tips()
    n = len(tree.lvsnms())
    vals = [None]*(n-2)
    p = 0
    for i in tree.iternodes(order=1):
        if i.istip == False and i != tree:
            x = [j.charst for j in i.children]
            t = [j.length for j in i.children]
            ui = abs(x[0]-x[1])
            Vi = sum(t)
            vals[p] = (ui,Vi)
            add = (t[0]*t[1])/(t[0]+t[1])
            i.length = i.length + add
            p += 1
        if i == tree:
            t = [j.length for j in i.children]
            Vi = sum(t)
            V0 = (t[0]*t[1])/(t[0]+t[1])
    div = sum([math.pow(i[0],2)/i[1] for i in vals])+(0/V0)
    sig2 = (1./n) * div
    return sig2
