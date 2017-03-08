import tree_reader

def init_heights(tree,sigsq=None): ##TODO: FIX BUG IN THIS
    for i in tree.iternodes(order=1):
        if i.istip or i.parent == None:
            if sigsq != None:
                if len(sigsq)==1:
                    i.sigsq = sigsq[0]
                elif len(sigsq) > 1:
                    i.sigsq = sigsq[i.rate_class]
            continue
        o = []
        start = False
        for j in i.children:
            if j.occurrences == []:
                continue
            if j.occurrences != None and j.occurrences[0]!="NA":
                if start == False:
                    o = j.occurrences
                    start = True
                else:
                    o = o+j.occurrences
        if len(o) != 0:
            i.height = max(o) + 0.1
            print i.height,[(j.label,j.height) for j in i.children]
            #print i.height,[(m.label,m.height) for m in i.children]
        else:
            i.height = max([j.height for j in i.children])+0.1
        if sigsq!=None:
            if len(sigsq)==1:
                i.sigsq = sigsq[0]
            elif len(sigsq) > 1:
                i.sigsq = sigsq[i.rate_class]
    for i in tree.iternodes():
        if i == tree:
            continue
        i.length = i.parent.height-i.height

def assign_node_nums(tree):
    num = 0
    for i in tree.iternodes(order=0):
        if i.istip or i == tree:
            continue
        else:
            i.number = num
            num += 1
    return tree

def assign_node_heights(h,tree):
    for i in tree.iternodes(order=0):
        if i.istip:
            i.length = i.parent.height-i.height
            if i.length < 0:
                return True
        elif i.parent == None:
            i.length = 0.01
        else:
            i.height = h[i.number]
            i.length = i.parent.height-i.height
            if i.length < 0:
                #print i.get_newick_repr(False),i.number,h[i.number],i.parent.height,i.length
                return True
    return False

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

def match_traits_tips(tree,traits,number):
    for i in tree.leaves():
        i.charst = traits[i.label][number]
    return tree

def read_strat(stratfl):
    ran = {}
    for i in open(stratfl,"r"):
        spls = i.strip().split("\t")
        t = spls[0]
        s = [float(i) for i in spls[1:] if i != "NA"]
        #print s
        ran[t]=s
    return ran

def match_strat(tree,strat):
    for i in tree.iternodes(order=1):
        if i.label in strat.keys():
            i.occurrences = strat[i.label]

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
        i.sigsq = 0.5
    return tree
