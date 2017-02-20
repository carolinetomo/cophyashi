import sys,os
import node
from alignments import pairwise_alignment

def upgma_tree_builder(seqs,matrix=None,aligned=True,verbose=False):
	cur_nodes = []
	distances = {}
	s_dist = None
	names = []
	if matrix == None:
		if aligned:
			s_dist = calc_aligned_seq_dist_mat(seqs)
		else:
			print "calculating pairwise alignments"
			s_dist = calc_unaligned_seq_dist_mat(seqs)
	else:
		s_dist = matrix
	if verbose:
		pretty_print_dist_mat_seqs(seqs,s_dist)
	for i in s_dist:
		t = node.Node()
		t.label = i
		cur_nodes.append(t)
		t.istip = True
	orig_distances = {}
	for i in cur_nodes:
		ndist = {}
		ndist2 = {}
		for j in cur_nodes:
			if i!=j:
				ndist[j] = s_dist[i.label][j.label]
				ndist2[j] = s_dist[i.label][j.label]
		distances[i] = ndist
		orig_distances[i] = ndist2
	while len(cur_nodes) > 2:
		mini = 9999999999999
		smallnode1 = None
		smallnode2 = None
		for i in distances:
			for j in distances[i]:
				if distances[i][j] < mini:
					mini = distances[i][j]
					smallnode1 = i
					smallnode2 = j
		t = node.Node()
		#set node 1 length
		if smallnode1.istip == True:
			smallnode1.length = mini/2.
		else:
			tlen = 0
			cn = smallnode1.leaves()[0]
			while cn != smallnode1:
				tlen += cn.length
				cn = cn.parent
			smallnode1.length = (mini/2.)-tlen
		#set node 2 length
		if smallnode2.istip == True:
			smallnode2.length = mini/2.
		else:
			tlen = 0
			cn = smallnode2.leaves()[0]
			while cn != smallnode2:
				tlen += cn.length
				cn = cn.parent
			smallnode2.length = (mini/2.)-tlen
		#add node 1 and node 2 to the new node
		t.add_child(smallnode1)
		t.add_child(smallnode2)
		cur_nodes.remove(smallnode1)
		cur_nodes.remove(smallnode2)
		cur_nodes.append(t)
		#calculate new distances
		del distances[smallnode1]
		del distances[smallnode2]
		tdists = {}
		tlvs = t.leaves()
		for i in distances:
			del distances[i][smallnode1]
			del distances[i][smallnode2]
			lvs = i.leaves()
			tsum = 0
			count = 0
			for j in lvs:
				for k in tlvs:
					tsum += orig_distances[j][k]
					count += 1
			tsum /= float(count)
			tdists[i] = tsum
			distances[i][t] = tdists[i]
		distances[t] = tdists
	root = node.Node()
	#get the last distance
	tdist = 0
	for j in distances:
		for k in distances[j]:
			tdist = distances[j][k]
	for i in cur_nodes:
		if i.istip == True:
			i.length = tdist/2.
		else:
			tlen = 0
			cn = i.leaves()[0]
			while cn != i:
				tlen += cn.length
				cn = cn.parent
			i.length = (tdist/2.)-tlen
		root.add_child(i)
	return root

"""
this calculates a neighbor joining tree
from raw distances
can send either sequences (default) or a matrix dictionary
"""
def nj_tree_builder(seqs,matrix = None,aligned = True,verbose=False):
	cur_nodes = [] #list of the nodes that can be found in the distance map
	distances = {} #key is node value is map of distances
	s_dist = None
	names = []
	if matrix == None:
		if aligned:
			s_dist = calc_aligned_seq_dist_mat(seqs)
		else:
			print "calculating pairwise alignments"
			s_dist = calc_unaligned_seq_dist_mat(seqs)
	else:
		s_dist = matrix
	if verbose:
		pretty_print_dist_mat_seqs(seqs,s_dist)
	root = node.Node()
	for i in s_dist:
		t = node.Node()
		t.label = i
		cur_nodes.append(t)
		t.istip = True
		root.add_child(t)
	for i in cur_nodes:
		ndist = {}
		for j in cur_nodes:
			if i != j:
				ndist[j] = s_dist[i.label][j.label]
		distances[i] = ndist
	while len(cur_nodes) > 2:
		#get_min_dist
		mini = 999999999 #arbitrary large
		smallnode1 = None
		smallnode2 = None
		qmat = calc_q_matrix(cur_nodes,distances)
		if verbose:
			pretty_print_dist_mat_nodes(cur_nodes,qmat)
		for i in distances:
			for j in distances[i]:
				if qmat[i][j] < mini:
					mini = qmat[i][j]
					smallnode1 = i
					smallnode2 = j
		#create a new node
		mini = distances[smallnode1][smallnode2]
		t = node.Node()
		tp = smallnode1.parent
		tp.add_child(t)
		tp.remove_child(smallnode1)
		tp.remove_child(smallnode2)
		#calculate branch length
		r = len(distances)
		tsum1 = 0
		for i in distances[smallnode1]:
			if i != smallnode2:
				tsum1 += distances[smallnode1][i]
		tsum2 = 0
		for i in distances[smallnode2]:
			if i != smallnode1:
				tsum2 += distances[smallnode2][i]
		newlength1 = (mini/2.) + (1./(2*(r-2))*(tsum1-tsum2))
		smallnode1.length = newlength1
		smallnode2.length = mini-newlength1
		print "joining "+smallnode1.label+":"+str(smallnode1.length),smallnode2.label+":"+str(smallnode2.length)
		t.add_child(smallnode1)
		t.add_child(smallnode2)
		del distances[smallnode1]
		del distances[smallnode2]
		cur_nodes.remove(smallnode1)
		cur_nodes.remove(smallnode2)
		tdists = {}
		for i in distances:
			d1 = distances[i][smallnode1]
			d2 = distances[i][smallnode2]
			del distances[i][smallnode1]
			del distances[i][smallnode2]
			tdists[i] = (d1+d2-mini)/2.
			distances[i][t] = tdists[i]
		distances[t] = tdists
		cur_nodes.append(t)
	for i in distances:
		for j in distances[i]:
			i.length = distances[i][j]/2.
	return root

