# Functions for Stepic assignments

# 40-9
DNA_BASES = ['A', 'C', 'G', 'T']

# 39-3
# 39-5
# 40-9
# 41-4
# 43-4
# Find the probability of a kmer matching a profile matrix
def profile_probability(kmer, profile):
	total = 1.0
	for idx, base in enumerate(kmer):
		if base in profile[idx]:
			total = total * profile[idx][base]
		else:
			return 0.0
	return total
 	
# 39-3
# 39-5 	
# 40-9
# 41-4
# 51-3
# 53-6
# Enumerate all kmers in a sequence
def all_kmers(sequence, k):
	for i in range(len(sequence)-k+1):
		yield sequence[i:i+k]

# 39-3
# 39-5
# 40-9
# 41-4
# Find the kmer in a sequence that best matches a profile matrix
def profile_most_probable_kmer(sequence, k, profile):
	best_kmer = sequence[:k]
	best_kmer_probability = 0
	for kmer in all_kmers(sequence, k):
		prob = profile_probability(kmer, profile)
		if prob > best_kmer_probability:
 			best_kmer = kmer
 			best_kmer_probability = prob
 	return best_kmer

# 39-5
# 40-9
# 41-4
# 43-4
# Count the number of each base for corresponding positions of a set of kmers
def count_matrix(motifs):
	from collections import Counter
	
	k = len(motifs[0])
	matrix = []
	for base_position in range(k):
		column = []
		for motif in motifs:
			column.append(motif[base_position])
		matrix.append(Counter(column).most_common())
	return matrix

# 39-5
# 40-9
# 41-4
# 43-4
# Count the number of differences among a set of kmers
def score_motifs(motifs):
	score = 0
	cmatrix = count_matrix(motifs)
	for base_position in cmatrix:
		for base_count_pair in base_position[1:]:
			score = score + base_count_pair[1]
	return score
	
# 39-5
# Convert a set of motifs to a profile matrix
def profile_matrix(motifs):
	kmer_count = len(motifs)
	cmatrix = count_matrix(motifs)
	matrix = []
	for base_position in cmatrix:
		matrix.append({base_count_pair[0]: float(base_count_pair[1]) / kmer_count for base_count_pair in base_position})
	return matrix

# 39-5
def greedy_motif_search(sequences, k, t):
	best_motifs = [x[:k] for x in sequences]
	best_motif_score = 100000000
	
	for kmer in all_kmers(sequences[0], k):
		motifs = [kmer]
		for i in range(1, t):
			profile = profile_matrix(motifs)
			motifs.append(profile_most_probable_kmer(sequences[i], k, profile))
		s = score_motifs(motifs)
		if s < best_motif_score:
			best_motifs = motifs
			best_motif_score = s
	return best_motifs
	
# 40-9
# 41-4
# 43-4
# Convert a set of motifs to a profile matrix
def profile_matrix_with_pseudocounts(motifs):
	kmer_count = len(motifs)
	cmatrix = count_matrix(motifs)
	matrix = []
	for base_position in cmatrix:
		count_dict = dict(base_position)
		profile_dict = {}
		for base in DNA_BASES:
			if base in count_dict:
				profile_dict[base] = (count_dict[base] + 1.0) / (2 * kmer_count)
			else:
				profile_dict[base] = 1.0 /  (2 * kmer_count)
		matrix.append(profile_dict)
	return matrix

# 40-9
def greedy_motif_search_with_pseudocounts(sequences, k, t):
	best_motifs = [x[:k] for x in sequences]
	best_motif_score = 100000000
	
	for kmer in all_kmers(sequences[0], k):
		motifs = [kmer]
		for i in range(1, t):
			profile = profile_matrix_with_pseudocounts(motifs)
			motifs.append(profile_most_probable_kmer(sequences[i], k, profile))
		s = score_motifs(motifs)
		if s < best_motif_score:
			best_motifs = motifs
			best_motif_score = s
	return best_motifs

# 41-4
def randomized_motif_search(sequences, k, t):
	from random import choice
	
	best_motif_score = 100000000
	best_motifs = []
	for sequence in sequences:
		best_motifs.append(choice(list(all_kmers(sequence, k))))
	
	motifs = best_motifs
	while True:
		profile = profile_matrix_with_pseudocounts(motifs)
		motifs = [profile_most_probable_kmer(sequence, k, profile) for sequence in sequences]
		score = score_motifs(motifs)
		if score < best_motif_score:
			best_motifs = motifs
			best_motif_score = score
		else:
			break
	return best_motif_score, best_motifs

# 41-4
def repeated_randomized_motif_search(sequences, k, t, n):
	best_motif_score = 100000000
	best_motifs = []
	for i in range(n):
		score, motifs = randomized_motif_search(sequences, k, t)
		if score < best_motif_score:
			best_motifs = motifs
			best_motif_score = score
	return best_motifs			

# 43-4
# Find a kmer in a sequence that probabilistically matches a profile matrix
def profile_random_kmer(sequence, k, profile):
	from random import uniform
	
	total_probability = 0
	probability_pairs = []
	for kmer in all_kmers(sequence, k):
		prob = profile_probability(kmer, profile)
		total_probability += prob
		probability_pairs.append([prob, kmer])
	
	randval = uniform(0, total_probability)
	probability_acc = 0 
	chosen_kmer = ''
	while probability_acc < randval:
		prob, chosen_kmer = probability_pairs.pop()
		probability_acc += prob
	
 	return chosen_kmer

# 43-4
def gibbs_sampler(sequences, k, t, n):
	from random import choice, randrange
	
	best_motif_score = 100000000
	best_motifs = []
	for sequence in sequences:
		best_motifs.append(choice(list(all_kmers(sequence, k))))
	
	motifs = best_motifs
	for i in range(n):
		j = randrange(len(sequences))
		motifs.pop(j)
		profile = profile_matrix_with_pseudocounts(motifs)
		motifs.insert(j, profile_random_kmer(sequences[j], k, profile))
		score = score_motifs(motifs)
		if score < best_motif_score:
			best_motifs = motifs
			best_motif_score = score
	return best_motif_score, best_motifs

# 43-4
def repeated_gibbs_sampler(sequences, k, t, n, runcount):

	best_motif_score = 100000000
	best_motifs = []
	for i in range(runcount):
		score, motifs = gibbs_sampler(sequences, k, t, n)
		if score < best_motif_score:
			best_motifs = motifs
			best_motif_score = score
			print (score, best_motifs)
	return best_motifs			

# 52-7
def overlap(a, b):
	return a[1:] == b[:-1]
	
# 52-7
def overlap_to_str(o):
	return "{} -> {}".format(*o)

# 52-7
# construct an overlap graph from arbitrary fragments
def overlap_graph(sequences):
	graph = []
	for s1 in sequences:
		for s2 in sequences:
			if overlap(s1, s2):
				graph.append([s1, s2])
	return graph

# 53-6
# 54-7
def debruijn_to_str(source, sinks):
	return "{} -> {}".format(source, ','.join(sinks))

# 53-6
# 54-7
def debruijn_graph(sequences):
	graph = {}
	for sequence in sequences:
		source, sink = sequence[:-1], sequence[1:]
		if source in graph:
			graph[source].add(sink)
		else:
			graph[source] = set([sink])
	return graph

# 57-2
# 57-5
def parse_graph_edges(edge_strs):
	graph = {}
	for edge_str in edge_strs:
		source, dummy, sink_str = edge_str.split(' ')
		sinks = sink_str.split(',')
		graph[source] = set(sinks)
	return graph
	
# 57-2
# 57-5
# 57-6		
# 57-10
def find_cycle_starting_at(graph, startnode):
	from random import sample
	
	cycle = [startnode]
	node = startnode
	complete = False
	while not complete:
		sinks = graph[node]
		next = sample(sinks, 1)[0]
		if len(sinks) > 1:
			graph[node] = sinks - set([next])
		else:
			del(graph[node])
		cycle.append(next)
		node = next
		complete = (node == startnode)
	return cycle, graph

# 57-2
# 57-10
def find_cycle(graph):
	from random import choice	
	startnode = choice(graph.keys())
	return find_cycle_starting_at(graph, startnode)

# 57-2
# 57-5
# 57-6		
# 57-10
def combine_cycles(cycle, index, new_cycle):
	cycle = cycle[:index] + new_cycle + cycle[index+1:]
	return cycle

# 57-2		
# 57-10
def find_eulerian_cycle(graph):
	cycle, remaining_graph = find_cycle(graph)
	while remaining_graph:
		for index, new_start in enumerate(cycle):
			if new_start in remaining_graph:
				new_cycle, remaining_graph = find_cycle_starting_at(remaining_graph, new_start)
				cycle = combine_cycles(cycle, index, new_cycle)
				break
		else:
			raise Exception("Cannot find any nodes from {} in remaining graph {}".format(cycle, remaining_graph))
	return cycle

# 57-5
# 57-6		
def path_degrees(graph):
	indegree = {}
	outdegree = {}
	
	for source in graph:
		if source not in indegree:
			indegree[source] = 0
		outdegree[source] = len(graph[source])
		for sink in graph[source]:
			if sink in indegree:
				indegree[sink] += 1
			else:
				indegree[sink] = 1
			if sink not in outdegree:
				outdegree[sink] = 0
	
	return indegree, outdegree

# 57-5		
# 57-6		
def find_eulerian_endpoints(graph):
	indegree, outdegree = path_degrees(graph)
	startnode, endnode = None, None
	
	for node in indegree:
		ins, outs = indegree[node], outdegree[node]
		if outs == ins + 1:
			if startnode == None:
				startnode = node
			else:
				raise Exception("Eulerian Path would have two start nodes: {} and {}".format(startnode, node))
		elif ins == outs + 1:
			if endnode == None:
				endnode = node
			else:
				raise Exception("Eulerian Path would have two end nodes: {} and {}".format(endnode, node))
				
	if startnode == None or endnode == None:
		raise Exception("Could not find Eulerian Path endpoints")
		
	return startnode, endnode

# 57-5		
# 57-6		
def find_eulerian_path(graph):
	startnode, endnode = find_eulerian_endpoints(graph)
	if endnode in graph:
		graph[endnode].add(startnode)
	else:
		graph[endnode] = set([startnode])
	
	cycle, remaining_graph = find_cycle_starting_at(graph, startnode)
	while remaining_graph:
		for index, new_start in enumerate(cycle):
			if new_start in remaining_graph:
				new_cycle, remaining_graph = find_cycle_starting_at(remaining_graph, new_start)
				cycle = combine_cycles(cycle, index, new_cycle)
				break
		else:
			raise Exception("Cannot find any nodes from {} in remaining graph {}".format(cycle, remaining_graph))
	return cycle[:-1]

# 57-6
# 57-10
def overlap_n(a, b, n):
	return a[-n:] == b[:n]

# 57-6		
# 57-10		
def assemble_path(path):
	sequence = path[0]
	for kmer in path[1:]:
		if overlap_n(sequence, kmer, len(kmer)-1):
			sequence += kmer[-1]
		else:
			raise Exception('kmer {} does not extend existing sequence {}'.format(kmer, sequence))
	return sequence
	
# 57-10
def universal_extend(inlist):
	outlist = []
	for elem in inlist:
		outlist.append(elem + '0')
		outlist.append(elem + '1')
	return outlist		

# 57-10
def universal_extend_n(n):
	counter = n
	outlist = ['']
	while counter > 0:
		outlist = universal_extend(outlist)
		counter -= 1
	return outlist  

# 57-10
def k_universal_graph(k):
	graph = {}
	for node in universal_extend_n(k-1):
		graph[node] = set(universal_extend([node[1:]]))
	return graph
		
	