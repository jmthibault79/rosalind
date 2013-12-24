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
    for i in range(len(sequence) - k + 1):
        yield sequence[i:i + k]

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
                profile_dict[base] = 1.0 / (2 * kmer_count)
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
# 58-14	
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
# 58-14	
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
# 58-14	
# 59-5
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
# 58-14	
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
# 58-14	
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
# 59-5
def overlap_n(a, b, n):
    return a[-n:] == b[:n]

# 57-6		
# 57-10		
# 59-5
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

# 58-14
def parse_graph_from_pairs(pairs):
    graph = {}
    for pair in pairs:
        read1, read2 = pair.split('|')
        source = '{}|{}'.format(read1[:-1], read2[:-1])
        sink = '{}|{}'.format(read1[1:], read2[1:])

        if source in graph:
            graph[source].add(sink)
        else:
            graph[source] = set([sink])
    return graph

# 58-14	
def assemble_path_from_pairs(path, d):
    first_reads = [x.split('|')[0] for x in path]
    second_reads = [x.split('|')[1] for x in path]
    first_reads_path = assemble_path(first_reads)
    second_reads_path = assemble_path(second_reads)
    gap = d + len(first_reads[0]) + 1
    if not overlap_n(first_reads_path, second_reads_path, len(first_reads_path)-gap):
        raise Exception('paths {} and {} do not overlap'.format(first_reads_path, second_reads_path))

    return first_reads_path + second_reads_path[-gap:]

# 59-5
def debruijn_graph_with_duplicates(sequences):
    graph = {}
    for sequence in sequences:
        source, sink = sequence[:-1], sequence[1:]
        if source in graph:
            graph[source].append(sink)
        else:
            graph[source] = [sink]
    return graph

# 59-5
def find_contigs(graph, node, indegree, outdegree):
    contigs = []
    for next in graph[node]:
        new_path = [node, next]
        ins, outs = indegree[next], outdegree[next]
        while ins == 1 and outs == 1:
            node = next
            next = graph[node][0]
            new_path.append(next)
            ins, outs = indegree[next], outdegree[next]
        contigs.append(new_path)
    return contigs

# 59-5
def debruijn_to_contigs(graph):
    outpaths = []
    indegree, outdegree = path_degrees(graph)
    for node in outdegree:
        ins, outs = indegree[node], outdegree[node]
        if outs > 0 and not (outs == 1 and ins == 1):
            outpaths.extend(find_contigs(graph, node, indegree, outdegree))
    return outpaths

# 71-8
def make_change(change, coins):
    change_map = { 0: 0, 1: 1 }
    for money in range(2, change + 1):
        min = 10000000
        for coin in coins:
            if money - coin in change_map:
                num_coins = change_map[money - coin] + 1
                if num_coins < min:
                    min = num_coins
        change_map[money] = min
    return change_map[change]

# 72-9
def parse_matrix(instrings, n, m):
    mat = []
    if len(instrings) != n:
        raise Exception('Expected n={} rows, saw {}'.format(n, len(instrings)))
    for instring in instrings:
        row = map(int, instring.split(' '))
        if len(row) != m:
            raise Exception('Expected m={} columns, saw {}'.format(m, len(instring)))
        mat.append(row)
    return mat

# 72-9
def longest_path(n, m, downmatrix, rightmatrix):
    pathmatrix = []
    for rows in range(n + 1):
        row = []
        for cols in range(m + 1):
            row.append(0)
        pathmatrix.append(row)

    for row in range(1, n + 1):
        pathmatrix[row][0] = pathmatrix[row - 1][0] + downmatrix[row - 1][0]
    for col in range(1, m + 1):
        pathmatrix[0][col] = pathmatrix[0][col - 1] + rightmatrix[0][col - 1]
    for row in range(1, n + 1):
        for col in range(1, m + 1):
            down = pathmatrix[row - 1][col] + downmatrix[row - 1][col]
            right = pathmatrix[row][col - 1] + rightmatrix[row][col - 1]
            pathmatrix[row][col] = max(down, right)

    return pathmatrix[n][m]