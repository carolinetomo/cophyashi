import tree_utils
import tree_reader
import math
from scipy import optimize
import random

LARGE = 10000000000

def hr97_loglike(tree,lam):
    lik = []
    for i in tree.leaves():
        if i.occurrences == []:
            continue
        if i.occurrences[0] == "NA" or i.occurrences == None: #TODO figure out what to do with ghost lineages
            continue
        #elif len(i.occurrences) == 1 and i.occurrences[0]: 
        #    continue
        tf = i.parent.height
        tl = i.height
        if len(i.occurrences) > 1:
            f = max(i.occurrences)
            l = min(i.occurrences)
            num = len(i.occurrences)
            a = math.log(math.pow((abs(l-f)),(num-2)))
            b = math.log(math.pow(lam,num))
            c = -lam*(abs(tl-tf))
            top = a+b+c
            bot = math.factorial(num-2)
            #print str(top)+" "+str(bot)
            loglik = top-bot
        elif len(i.occurrences) == 1:
            loglik = math.log(1/(abs(tl-tf)))
        lik.append(loglik)
    return sum(lik)

def calc_like_strat(p,tree,strat): #p[0] should be lambda, p[1:] is node heights
    for i in p[1:]:
        if i < 0:
            return LARGE
    bad = tree_utils.assign_node_heights(p[1:],tree)
    if bad:
        return LARGE
    try:
        val = -hr97_loglike(tree,p[0])
    except:
        return LARGE
    return val

def optim_lambda_heights(tree,strat):
    tree_utils.assign_node_nums(tree)
    start = [random.uniform(0.1,1)]
    start += [i.height for i in tree.iternodes() if i.istip == False and i != tree]
    opt = optimize.fmin_bfgs(calc_like_strat,start,args=(tree,strat),full_output = True, disp = True)
    tree_utils.assign_node_heights(opt[0][1:],tree)
    return [tree,opt]
