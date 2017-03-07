import tree_utils
import tree_reader
import math
from scipy import optimize

def strat_like(tree,strat,lam):
    lik = []
    for i in tree.leaves():
        if i.occurrences[0] == "NA" or i.occurrences == None:
            continue
        elif len(i.occurrences) == 1 and i.occurrences[0]: #TODO: figure out what to do with extant tips
            continue
        f = max(i.occurrences)
        l = min(i.occurrences)
        num = len(i.occurrences)
        tf = i.parent.height
        tl = i.height
        a = math.pow((l-f),(num-2))
        b = math.pow(lam,num)
        c = math.exp(-lam*(tl-tf))
        top = a*b*c
        bot = math.factorial(num-2)
        loglik = math.log(top/bot)
        lik.append(loglik)
    return sum(lik)
