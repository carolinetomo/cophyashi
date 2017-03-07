import tree_utils
import tree_reader
import math
from scipy import optimize

def hr97_like(tree,strat,lam):
    lik = []
    for i in tree.leaves():
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
            a = math.pow((l-f),(num-2))
            b = math.pow(lam,num)
            c = math.exp(-lam*(tl-tf))
            top = a*b*c
            bot = math.factorial(num-2)
            print str(top)+" "+str(bot)
            loglik = math.log(top/bot)
        elif len(i.occurrences) == 1:
            loglik = math.log(1/(tl-tf))
        lik.append(loglik)
    return sum(lik)
