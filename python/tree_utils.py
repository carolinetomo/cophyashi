import tree_reader

def init_heights(tree,strat=False,fixed_root=False): 
    for i in tree.iternodes(order=1):
        if i.istip:# or i.parent == None:
            continue
        elif i == tree and fixed_root == True:
            continue
        if strat == True:
            o = []
            start = False
            for j in i.children:
                if j.occurrences != None and "NA" not in j.occurrences:
                    if start == False:
                        o = j.occurrences
                        #print o
                        start = True
                    else:
                        o = o+j.occurrences
            if len(o) > 0:
                i.height = max(o+[j.height for j in i.children]) + 0.1
                #print [j.height for j in i.children],i.height,[(j.label,j.height) for j in i.children]
            else:
                i.height = max([j.height for j in i.children])+0.1
        else:
            i.height = max([j.height for j in i.children])+0.1
    for i in tree.iternodes():
        if i == tree:
            continue
        i.length = i.parent.height-i.height

def init_heights_strat(tree,fixed_root=False): #TODO: fix negative branchlength problem
    for i in tree.iternodes(order=1):
        if i == tree and fixed_root == True:
            continue
        elif i.istip and i.occurrences[0] != "NA": 
            if min(i.occurrences) != 0.0:
                i.height = min(i.occurrences)- 0.01
            else:
                i.height = min(i.occurrences)
        elif i.istip == False:
            o = []
            start = False
            for j in i.children:
                if j.occurrences != None and "NA" not in j.occurrences:
                    if start == False:
                        o = j.occurrences
                        #print o
                        start = True
                    else:
                        o = o+j.occurrences
            if len(o) > 0:
                i.height = max(o+[j.height for j in i.children]) + 0.1
                #print [j.height for j in i.children],i.height,[(j.label,j.height) for j in i.children]
            else:
                i.height = max([j.height for j in i.children])+0.1
    for i in tree.iternodes():
        if i == tree:
            continue
        elif i.occurrences[0] == "NA" and i.istip:
            i.height = i.parent.height - 0.1
        i.length = i.parent.height-i.height

def assign_node_nums(tree,tips = True,fixed_root = False):
    num = 0
    for i in tree.iternodes(order=0):
        if i.istip and tips == False:
            continue
        elif i == tree and fixed_root == True:
            continue
        else:
            i.number = num
            num += 1

def assign_branch_nums(tree):
    num = 0
    for i in tree.iternodes(order=0):
        if i == tree:
            continue
        i.number = num
        num +=1

def assign_brlens(l,tree):
    for i in tree.iternodes(order = 0):
        if i.parent == None:
            continue
        i.length = l[i.number]
        if i.length <0 :
            return True
    return False

def assign_node_heights(h,tree,tips=True,fixed_root=False):
    for i in tree.iternodes(order=0):
        if i.istip and tips == False:
            i.length = i.parent.height-i.height
            if i.length < 0:
                return True
        elif fixed_root == True and i == tree:
            continue
        elif fixed_root == False and i == tree:
            i.height = h[i.number]
            continue
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

def read_strat(stratfl):
    ran = {}
    for i in open(stratfl,"r"):
        spls = i.strip().split("\t")
        t = spls[0]
        s = [float(i) for i in spls[1:] if i != "NA"]
        if s == []:
            ran[t] = ["NA"]
        else:
            ran[t]=s
    return ran

def match_strat(tree,strat):
    for i in tree.iternodes(order=1):
        if i.label in strat.keys():
            i.occurrences = strat[i.label]
        elif i.occurrences == None:
            i.occurrences = ["NA"]

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

def assign_sigsq(tree,sigsq=[1]):
    for i in tree.iternodes():
        if len(sigsq) == 1:
            i.sigsq =sigsq[0]
        else:
            i.sigsq = sigsq[i.number] 
    return tree
