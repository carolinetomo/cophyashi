from phylo3 import *
import random


def tip_dates(tree,dates):
    d = {}
    for i in open(dates,"r"):
        spls = i.strip().split("\t")
        d[spls[0]] = float(spls[1])
    for i in tree.leaves():
        i.height = d[i]


def nni(tree,nnodes):
	#pick a node
	#while a node is a tip continue picking
	cnode = None
	count = 0
	stop = random.randint(0,nnodes)
	for i in tree.iternodes():
		if count == stop:
			cnode = i
		count += 1
	while cnode.istip:
		count = 0
		stop = random.randint(0,nnodes)
		for i in tree.iternodes():
			if count == stop:
				cnode = i
			count += 1
	ndli = []
	# if the parent is not the root get the parent of the parent
	# if the parent is the root, there should be three children, two other than the child
	#	then just change which goes to which , list = [nodes],use random.shuffle(list)
	children = cnode.children
	swap1 = None
	swap2 = None
	count = 0
	ndlip = []
	for i in children:
		if count >= 2:
			break
		if len(i.children) > 0:
			ndli.append(random.sample(i.children,1)[0])
			ndlip.append(i)
		else:
			ndli.append(i)
			ndlip.append(cnode)
	for i in range(len(ndlip)):
		ndlip[i].remove_child(ndli[i])
	for i in range(len(ndlip)):
		ndlip[-len(ndlip)+1-i].add_child(ndli[i])
	return tree

