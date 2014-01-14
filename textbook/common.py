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
# 248-3
# 248-5
# 248-7
# 249-8
# 250-12
# 250-14
def init_matrix(rows, cols):
    matrix = []
    for row in range(rows):
        m_row = []
        for col in range(cols):
            m_row.append(0)
        matrix.append(m_row)
    return matrix

# 72-9
def longest_path(n, m, downmatrix, rightmatrix):
    pathmatrix = init_matrix(n + 1, m + 1)

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

# 74-5
# 76-3
# 76-9
# 248-3
# 248-5
# 248-7
def max_and_direction(down, right, diag):
    dir = 'down'
    max = down
    if right > max:
        dir = 'right'
        max = right
    if diag > max:
        dir = 'diag'
        max = diag
    return max, dir

# 74-5
def print_matrix(matrix):
    for row in matrix:
        print(row)

# 74-5
def longest_common_subsequence(seq1, seq2):
    v = len(seq1)
    w = len(seq2)
    pathmatrix = init_matrix(v + 1, w + 1)
    backtrack_matrix = init_matrix(v + 1, w + 1)

    for row in range(1, v + 1):
        for col in range(1, w + 1):
            down = pathmatrix[row - 1][col]
            right = pathmatrix[row][col - 1]
            diag = pathmatrix[row - 1][col - 1]
            if seq1[row - 1] == seq2[col - 1]:
                diag += 1

            max, dir = max_and_direction(down, right, diag)
            pathmatrix[row][col] = max
            backtrack_matrix[row][col] = dir

    return pathmatrix[v][w], backtrack_matrix

# 74-5
def output_longest_common_subsequence(backtrack_matrix, v, i, j):
    if i == 0 or j == 0:
        return ''

    dir = backtrack_matrix[i][j]
    if dir == 'down':
        return output_longest_common_subsequence(backtrack_matrix, v, i - 1, j)
    elif dir == 'right':
        return output_longest_common_subsequence(backtrack_matrix, v, i, j - 1)
    else:
        retstr = output_longest_common_subsequence(backtrack_matrix, v, i - 1, j - 1)
        return retstr + v[i - 1]

# 74-7
def parse_dag_edges(edge_strs):
    graph = {}
    for edge_str in edge_strs:
        source, rest = edge_str.split('->')
        sink, weightstr = rest.split(':')
        weight = int(weightstr)
        if source in graph:
            graph[source].append([sink, weight])
        else:
            graph[source] = [[sink, weight]]
    return graph

# 74-7
def wikipedia_depth_first_topological_sort_visit(dag,node,unmarked,temp_marked,ordered):
    if node in temp_marked:
        raise Exception("Not a DAG!")
    if node in unmarked:
        temp_marked.add(node)
        for sinkweight in dag[node]:
            sink, weight = sinkweight
            unmarked,temp_marked,ordered = wikipedia_depth_first_topological_sort_visit(dag,sink,unmarked,temp_marked,ordered)
        unmarked.remove(node)
        temp_marked.remove(node)
        ordered = [node] + ordered
    return unmarked,temp_marked,ordered

# 74-7
def wikipedia_depth_first_topological_sort(dag, sink):
    from random import sample

    # set of nodes only
    unmarked = {x for x in dag}
    temp_marked = set()
    ordered = []
    while unmarked:
        node = sample(unmarked, 1)[0]
        unmarked,temp_marked,ordered = wikipedia_depth_first_topological_sort_visit(dag,node,unmarked,temp_marked,ordered)

    return ordered

# 74-7
def longest_dag_weight(dag, ordering, source, final_sink):
    vals = {}
    backtrack = {}
    for node in ordering:
        if node in vals:
            nodeval = vals[node]
        else:
            if node != source:
                continue
            nodeval = 0
        for sinkweight in dag[node]:
            sink, weight = sinkweight
            sinkval = nodeval + weight
            if sink in vals:
                if sinkval > vals[sink]:
                    vals[sink] = sinkval
                    backtrack[sink] = node
            else:
                vals[sink] = sinkval
                backtrack[sink] = node

    return vals[final_sink], backtrack

# 74-7
def output_longest_dag_path(backtrack, source, sink):
    pred = backtrack[sink]
    path = [pred, sink]
    while pred != source:
        pred = backtrack[pred]
        path = [pred] + path
    return '->'.join(path)

# 76-3
# 76-9
def parse_scoring_matrix(lines):
    # non-empty only
    header = [x for x in lines[0].strip().split(' ') if x]

    matrix = {}
    for line in lines[1:]:
        elements = [x for x in line.strip().split(' ') if x]
        matrix_row = dict(zip(header, map(int, elements[1:])))
        matrix[elements[0]] = matrix_row
    return matrix

# 76-3
# 248-3
def scored_longest_common_subsequence(scoring_matrix, indel_penalty, seq1, seq2):
    v = len(seq1)
    w = len(seq2)
    pathmatrix = init_matrix(v + 1, w + 1)
    backtrack_matrix = init_matrix(v + 1, w + 1)

    for row in range(1, v + 1):
        pathmatrix[row][0] = pathmatrix[row - 1][0] + indel_penalty
    for col in range(1, w + 1):
        pathmatrix[0][col] = pathmatrix[0][col - 1] + indel_penalty
    for row in range(1, v + 1):
        for col in range(1, w + 1):
            down = pathmatrix[row - 1][col] + indel_penalty
            right = pathmatrix[row][col - 1] + indel_penalty
            diag = pathmatrix[row - 1][col - 1] + scoring_matrix[seq1[row - 1]][seq2[col - 1]]

            max, dir = max_and_direction(down, right, diag)
            pathmatrix[row][col] = max
            backtrack_matrix[row][col] = dir

    return pathmatrix[v][w], backtrack_matrix

# 76-3
def output_longest_common_subsequence_aligned(backtrack_matrix, v, w, i, j):
    if i == 0:
        return '-'*j, w[:j]
    if j == 0:
        return v[:i], '-'*i

    dir = backtrack_matrix[i][j]
    if dir == 'down':
        retstr1, retstr2 = output_longest_common_subsequence_aligned(backtrack_matrix, v, w, i - 1, j)
        return retstr1 + v[i - 1], retstr2 + '-'
    elif dir == 'right':
        retstr1, retstr2 = output_longest_common_subsequence_aligned(backtrack_matrix, v, w, i, j - 1)
        return retstr1 + '-', retstr2 + w[j - 1]
    else:
        retstr1, retstr2 = output_longest_common_subsequence_aligned(backtrack_matrix, v, w, i - 1, j - 1)
        return retstr1 + v[i - 1], retstr2 + w[j - 1]

# 76-9
def scored_longest_common_subsequence_local(scoring_matrix, indel_penalty, seq1, seq2):
    v = len(seq1)
    w = len(seq2)
    pathmatrix = init_matrix(v + 1, w + 1)
    backtrack_matrix = init_matrix(v + 1, w + 1)
    best_result, best_row, best_col = -1, -1, -1
    for row in range(1, v + 1):
        for col in range(1, w + 1):
            down = pathmatrix[row - 1][col] + indel_penalty
            right = pathmatrix[row][col - 1] + indel_penalty
            diag = pathmatrix[row - 1][col - 1] + scoring_matrix[seq1[row - 1]][seq2[col - 1]]

            max, dir = max_and_direction(down, right, diag)
            if max < 0:
                max, dir = 0, 'zero'    # "free ride"
            elif max > best_result:
                best_result = max
                best_row = row
                best_col = col
            pathmatrix[row][col] = max
            backtrack_matrix[row][col] = dir
    return pathmatrix[best_row][best_col], backtrack_matrix, best_row, best_col

# 76-9
# 248-5
# 248-7
def output_longest_common_subsequence_local(backtrack_matrix, v, w, i, j):
    if i == 0 or j == 0:
        return '', ''

    dir = backtrack_matrix[i][j]
    if dir == 'down':
        retstr1, retstr2 = output_longest_common_subsequence_local(backtrack_matrix, v, w, i - 1, j)
        return retstr1 + v[i - 1], retstr2 + '-'
    elif dir == 'right':
        retstr1, retstr2 = output_longest_common_subsequence_local(backtrack_matrix, v, w, i, j - 1)
        return retstr1 + '-', retstr2 + w[j - 1]
    elif dir == 'diag':
        retstr1, retstr2 = output_longest_common_subsequence_local(backtrack_matrix, v, w, i - 1, j - 1)
        return retstr1 + v[i - 1], retstr2 + w[j - 1]
    else:
        return '', ''

# 248-3
def mismatch_scoring_matrix(alphabet):
    matrix = {}
    for letter in alphabet:
        matrix_row = {}
        for other_letter in alphabet:
            if letter == other_letter:
                matrix_row[other_letter] = 0
            else:
                matrix_row[other_letter] = -1
        matrix[letter] = matrix_row
    return matrix

# 248-5
def mismatch_scoring_matrix_fitted(alphabet):
    matrix = {}
    for letter in alphabet:
        matrix_row = {}
        for other_letter in alphabet:
            if letter == other_letter:
                matrix_row[other_letter] = 1
            else:
                matrix_row[other_letter] = -1
        matrix[letter] = matrix_row
    return matrix

# 248-5
def scored_longest_common_subsequence_fitted(scoring_matrix, indel_penalty, seq1, seq2):
    v = len(seq1)
    w = len(seq2)
    pathmatrix = init_matrix(v + 1, w + 1)
    backtrack_matrix = init_matrix(v + 1, w + 1)
    best_row, best_col = -1000000, -1000000
    for col in range(1, w + 1):
        best_result_for_col = -1000000
        pathmatrix[0][col] = pathmatrix[0][col - 1] + indel_penalty
        for row in range(1, v + 1):
            down = pathmatrix[row - 1][col] + indel_penalty
            right = pathmatrix[row][col - 1] + indel_penalty
            diag = pathmatrix[row - 1][col - 1] + scoring_matrix[seq1[row - 1]][seq2[col - 1]]

            max, dir = max_and_direction(down, right, diag)
            if max > best_result_for_col:
                best_result_for_col = max
                best_row = row
                best_col = col
            pathmatrix[row][col] = max
            backtrack_matrix[row][col] = dir
    return pathmatrix[best_row][best_col], backtrack_matrix, best_row, best_col

# 248-7
def mismatch_scoring_matrix_overlap(alphabet):
    matrix = {}
    for letter in alphabet:
        matrix_row = {}
        for other_letter in alphabet:
            if letter == other_letter:
                matrix_row[other_letter] = 1
            else:
                matrix_row[other_letter] = -2
        matrix[letter] = matrix_row
    return matrix

# 248-7
def scored_longest_common_subsequence_overlap(scoring_matrix, indel_penalty, seq1, seq2):
    v = len(seq1)
    w = len(seq2)
    pathmatrix = init_matrix(v + 1, w + 1)
    backtrack_matrix = init_matrix(v + 1, w + 1)
    best_result, best_row, best_col = None, None, None
    for col in range(1, w + 1):
        for row in range(1, v + 1):
            down = pathmatrix[row - 1][col] + indel_penalty
            right = pathmatrix[row][col - 1] + indel_penalty
            diag = pathmatrix[row - 1][col - 1] + scoring_matrix[seq1[row - 1]][seq2[col - 1]]

            best, dir = max_and_direction(down, right, diag)
            pathmatrix[row][col] = best
            backtrack_matrix[row][col] = dir

    best_result = None
    best_row = v
    for col in range(1, w + 1):
        val = pathmatrix[best_row][col]
        if val > best_result:
            best_result = val
            best_col = col

    return pathmatrix[best_row][best_col], backtrack_matrix, best_row, best_col

# 249-8
def lower_max_and_direction(lower_down, lower_from_middle):
    dir = 'down'
    max = lower_down
    if lower_from_middle > max:
        dir = 'lower_from_middle'
        max = lower_from_middle
    return max, dir

# 249-8
def upper_max_and_direction(upper_right, upper_from_middle):
    dir = 'right'
    max = upper_right
    if upper_from_middle > max:
        dir = 'upper_from_middle'
        max = upper_from_middle
    return max, dir

# 249-8
def middle_max_and_direction(diag, lower_best, upper_best):
    dir = 'diag'
    max = diag
    if lower_best > max:
        dir = 'middle_from_lower'
        max = lower_best
    if upper_best > max:
        dir = 'middle_from_upper'
        max = upper_best
    return max, dir

# 249-8
def scored_longest_common_subsequence_affine(scoring_matrix, gap_open, gap_extend, seq1, seq2):
    v = len(seq1)
    w = len(seq2)
    lower = init_matrix(v + 1, w + 1)
    middle = init_matrix(v + 1, w + 1)
    upper = init_matrix(v + 1, w + 1)
    backtrack_matrix_lower = init_matrix(v + 1, w + 1)
    backtrack_matrix_middle = init_matrix(v + 1, w + 1)
    backtrack_matrix_upper = init_matrix(v + 1, w + 1)

    for col in range(1, w + 1):
        for row in range(1, v + 1):
            lower_down = lower[row - 1][col] + gap_extend
            lower_from_middle = middle[row - 1][col] + gap_open
            lower_best, lower_dir = lower_max_and_direction(lower_down, lower_from_middle)
            lower[row][col] = lower_best
            backtrack_matrix_lower[row][col] = lower_dir

            upper_right = upper[row][col - 1] + gap_extend
            upper_from_middle = middle[row][col - 1] + gap_open
            upper_best, upper_dir = upper_max_and_direction(upper_right, upper_from_middle)
            upper[row][col] = upper_best
            backtrack_matrix_upper[row][col] = upper_dir

            diag = middle[row - 1][col - 1] + scoring_matrix[seq1[row - 1]][seq2[col - 1]]

            middle_best, middle_dir = middle_max_and_direction(diag, lower_best, upper_best)

            middle[row][col] = middle_best
            backtrack_matrix_middle[row][col] = middle_dir

    return middle[v][w], backtrack_matrix_lower, backtrack_matrix_middle, backtrack_matrix_upper, v, w

# 249-8
def output_longest_common_subsequence_affine_lower(backtrack_matrix_lower, backtrack_matrix_middle, backtrack_matrix_upper, v, w, i, j):
    if i == 0:
        return '-'*j, w[:j]
    if j == 0:
        return v[:i], '-'*i

    dir = backtrack_matrix_lower[i][j]
    if dir == 'lower_from_middle':
        retstr1, retstr2 = output_longest_common_subsequence_affine_middle(backtrack_matrix_lower, backtrack_matrix_middle, backtrack_matrix_upper, v, w, i - 1, j)
    else:
        retstr1, retstr2 = output_longest_common_subsequence_affine_lower(backtrack_matrix_lower, backtrack_matrix_middle, backtrack_matrix_upper, v, w, i - 1, j)

    return retstr1 + v[i - 1], retstr2 + "-"

# 249-8
def output_longest_common_subsequence_affine_upper(backtrack_matrix_lower, backtrack_matrix_middle, backtrack_matrix_upper, v, w, i, j):
    if i == 0:
        return '-'*j, w[:j]
    if j == 0:
        return v[:i], '-'*i

    dir = backtrack_matrix_upper[i][j]
    if dir == 'upper_from_middle':
        retstr1, retstr2 = output_longest_common_subsequence_affine_middle(backtrack_matrix_lower, backtrack_matrix_middle, backtrack_matrix_upper, v, w, i, j - 1)
    else:
        retstr1, retstr2 = output_longest_common_subsequence_affine_upper(backtrack_matrix_lower, backtrack_matrix_middle, backtrack_matrix_upper, v, w, i, j - 1)

    return retstr1 + "-", retstr2 + w[j - 1]

# 249-8
def output_longest_common_subsequence_affine_middle(backtrack_matrix_lower, backtrack_matrix_middle, backtrack_matrix_upper, v, w, i, j):
    if i == 0:
        return '-'*j, w[:j]
    if j == 0:
        return v[:i], '-'*i

    dir = backtrack_matrix_middle[i][j]
    if dir == 'middle_from_lower':
        return output_longest_common_subsequence_affine_lower(backtrack_matrix_lower, backtrack_matrix_middle, backtrack_matrix_upper, v, w, i, j)
    elif dir == 'middle_from_upper':
        return output_longest_common_subsequence_affine_upper(backtrack_matrix_lower, backtrack_matrix_middle, backtrack_matrix_upper, v, w, i, j)
    else:
        retstr1, retstr2 = output_longest_common_subsequence_affine_middle(backtrack_matrix_lower, backtrack_matrix_middle, backtrack_matrix_upper, v, w, i - 1, j - 1)
        return retstr1 + v[i - 1], retstr2 + w[j - 1]

# 250-12
# 250-14
def middle_edge_half_matrix(scoring_matrix, indel_penalty, halfseq1, halfseq2):
    v = len(halfseq1)
    w = len(halfseq2)
    m = init_matrix(v + 1, w + 1)
    backtrack_col = { 0: 'right' }
    for row in range(1, v + 1):
        m[row][0] = m[row - 1][0] + indel_penalty
    for col in range(1, w + 1):
        m[0][col] = m[0][col - 1] + indel_penalty
        for row in range(1, v + 1):
            down = m[row - 1][col] + indel_penalty
            right = m[row][col - 1] + indel_penalty
            diag = m[row - 1][col - 1] + scoring_matrix[halfseq1[row - 1]][halfseq2[col - 1]]

            best, direction = max_and_direction(down, right, diag)
            m[row][col] = best

            if col == w:
                backtrack_col[row] = direction
    return m, backtrack_col

# 250-12
# 250-14
def revstr(fwd_str):
    return fwd_str[::-1]

# 250-12
# 250-14
def alignment_middle_edge(scoring_matrix, indel_penalty, seq1, seq2):
    v = len(seq1)
    w = len(seq2)
    middle_col = w/2

    forward_half, forward_backtrack = middle_edge_half_matrix(scoring_matrix, indel_penalty, seq1, seq2[:middle_col])
    backward_half, backward_backtrack = middle_edge_half_matrix(scoring_matrix, indel_penalty, revstr(seq1), revstr(seq2[middle_col:]))

    backward_col_offset = middle_col + 1
    if middle_col * 2 == w:
         backward_col_offset = middle_col

    best_row, best_total = None, None
    for row in range(0, v + 1):
        total = forward_half[row][middle_col] + backward_half[v - row][backward_col_offset]
        if total > best_total:
            best_total = total
            best_row = row

    # which direction to go?  which was best in backward_half
    direction = backward_backtrack[v - best_row]
    if direction == 'diag':
        to_row, to_col = best_row + 1, middle_col + 1
    elif direction == 'down':
        to_row, to_col = best_row + 1, middle_col
    elif direction == 'right':
        to_row, to_col = best_row, middle_col + 1

    return best_row, middle_col, to_row, to_col

# 250-14
def linear_space_alignment_part(scoring_matrix, indel_penalty, seq1, seq2, row_start, col_start, row_end, col_end):
    if row_start == row_end:
        right_path = {}
        for col in range(col_start, col_end):
            right_path[(row_start, col)] = (row_start, col + 1)
        return right_path
    elif col_start == col_end:
        down_path = {}
        for row in range(row_start, row_end):
            down_path[(row, col_start)] = (row + 1, col_start)
        return down_path

    subseq1 = seq1[row_start:row_end]
    subseq2 = seq2[col_start:col_end]

    from_row, from_col, to_row, to_col = alignment_middle_edge(scoring_matrix, indel_penalty, subseq1, subseq2)

    edge_start_row, edge_start_col = row_start + from_row, col_start + from_col
    edge_end_row, edge_end_col = row_start + to_row, col_start + to_col
    path = { (edge_start_row, edge_start_col): (edge_end_row, edge_end_col) }

    first_half_path = linear_space_alignment_part(scoring_matrix, indel_penalty, seq1, seq2, row_start, col_start, edge_start_row, edge_start_col)
    second_half_path = linear_space_alignment_part(scoring_matrix, indel_penalty, seq1, seq2, edge_end_row, edge_end_col, row_end, col_end)

    path.update(first_half_path)
    path.update(second_half_path)

    return path

# 250-14
def score_and_align_linear_space(scoring_matrix, indel_penalty, seq1, seq2, v, w, path):
    score, align1, align2 = 0, '', ''
    row, col = 0, 0
    while row != v or col != w:
        to_row, to_col = path[(row, col)]
        if to_row == row + 1 and to_col == col + 1:
            let1, let2 = seq1[row], seq2[col]
            score = score + scoring_matrix[let1][let2]
            align1, align2 = align1 + let1, align2 + let2
        elif to_row == row + 1 and to_col == col:
            score = score + indel_penalty
            let1 = seq1[row]
            align1, align2 = align1 + let1, align2 + '-'
        elif to_row == row and to_col == col + 1:
            score = score + indel_penalty
            let2 = seq2[col]
            align1, align2 = align1 + '-', align2 + let2
        row, col = to_row, to_col

    return score, align1, align2

# 250-14
def linear_space_alignment(scoring_matrix, indel_penalty, seq1, seq2):
    v = len(seq1)
    w = len(seq2)
    path = linear_space_alignment_part(scoring_matrix, indel_penalty, seq1, seq2, 0, 0, v, w)
    return score_and_align_linear_space(scoring_matrix, indel_penalty, seq1, seq2, v, w, path)

# 251-5
def init_matrix_3(x_dim, y_dim, z_dim):
    matrix = []
    for x in range(x_dim):
        m_x = []
        for y in range(y_dim):
            m_y = []
            for z in range(z_dim):
                m_y.append(0)
            m_x.append(m_y)
        matrix.append(m_x)
    return matrix

# 251-5
def align_3_maxdir(x_only, y_only, z_only, xy_only, yz_only, zx_only, xyz):
    dir = 'x'
    max = x_only
    if y_only > max:
        dir = 'y'
        max = y_only
    if z_only > max:
        dir = 'z'
        max = z_only
    if xy_only > max:
        dir = 'xy'
        max = xy_only
    if yz_only > max:
        dir = 'yz'
        max = yz_only
    if zx_only > max:
        dir = 'zx'
        max = zx_only
    if xyz > max:
        dir = 'xyz'
        max = xyz

    return max, dir

# 251-5
def align_3(seq1, seq2 ,seq3):
    x_len = len(seq1)
    y_len = len(seq2)
    z_len = len(seq3)
    pathmatrix = init_matrix_3(x_len + 1, y_len + 1, z_len + 1)
    backtrack_matrix = init_matrix_3(x_len + 1, y_len + 1, z_len + 1)
    best_result = None

    for x in range(1, x_len + 1):
        for y in range(1, y_len + 1):
            for z in range(1, z_len + 1):
                x_only = pathmatrix[x-1][y][z]
                y_only = pathmatrix[x][y-1][z]
                z_only = pathmatrix[x][y][z-1]
                xy_only = pathmatrix[x-1][y-1][z]
                yz_only = pathmatrix[x][y-1][z-1]
                zx_only = pathmatrix[x-1][y][z-1]
                xyz = pathmatrix[x-1][y-1][z-1]

                # score 1 for a perfect match, 0 for all others
                if seq1[x-1] == seq2[y-1] and seq1[x-1] == seq3[z-1]:
                    xyz = xyz + 1

                best, dir = align_3_maxdir(x_only, y_only, z_only, xy_only, yz_only, zx_only, xyz)
                if best > best_result:
                    best_result = best

                pathmatrix[x][y][z] = best
                backtrack_matrix[x][y][z] = dir
    return best, backtrack_matrix, x_len, y_len, z_len

# 251-5
def output_align_3(backtrack_matrix, seq1, seq2, seq3, x, y, z):
    if x == 0 and y == 0 and z == 0:
        return '', '', ''
    if x == 0 and y == 0:
        return '-'*z, '-'*z, seq3[:z]
    if y == 0 and z == 0:
        return seq1[:x], '-'*x, '-'*x
    if z == 0 and x == 0:
        return '-'*y, seq2[:y], '-'*y

    if x == 0:
        dir = 'yz'
    elif y == 0:
        dir = 'zx'
    elif z == 0:
        dir = 'xy'
    else:
        dir = backtrack_matrix[x][y][z]

    if dir == 'x':
        retstr1, retstr2, retstr3 = output_align_3(backtrack_matrix, seq1, seq2, seq3, x-1, y, z)
        return retstr1 + seq1[x - 1], retstr2 + '-', retstr3 + '-'
    elif dir == 'y':
        retstr1, retstr2, retstr3 = output_align_3(backtrack_matrix, seq1, seq2, seq3, x, y-1, z)
        return retstr1 + '-', retstr2 + seq2[y - 1], retstr3 + '-'
    elif dir == 'z':
        retstr1, retstr2, retstr3 = output_align_3(backtrack_matrix, seq1, seq2, seq3, x, y, z-1)
        return retstr1 + '-', retstr2 + '-', retstr3 + seq3[z - 1]
    elif dir == 'xy':
        retstr1, retstr2, retstr3 = output_align_3(backtrack_matrix, seq1, seq2, seq3, x-1, y-1, z)
        return retstr1 + seq1[x - 1], retstr2 + seq2[y - 1], retstr3 + '-'
    elif dir == 'yz':
        retstr1, retstr2, retstr3 = output_align_3(backtrack_matrix, seq1, seq2, seq3, x, y-1, z-1)
        return retstr1 + '-', retstr2 + seq2[y - 1], retstr3 + seq3[z - 1]
    elif dir == 'zx':
        retstr1, retstr2, retstr3 = output_align_3(backtrack_matrix, seq1, seq2, seq3, x-1, y, z-1)
        return retstr1 + seq1[x - 1], retstr2 + '-', retstr3 + seq3[z - 1]
    elif dir == 'xyz':
        retstr1, retstr2, retstr3 = output_align_3(backtrack_matrix, seq1, seq2, seq3, x-1, y-1, z-1)
        return retstr1 + seq1[x - 1], retstr2 + seq2[y - 1], retstr3 + seq3[z - 1]

# 286-2
def greedysorting_parse(perm_str):
    return map(int, perm_str.replace('(','').replace(')','').split())

# 286-2
def greedysorting_out(sequence):
    retlist = []
    for line in sequence:
        linestr = ''
        for val in line:
            if val < 0:
                linestr += str(val) + ' '
            else:
                linestr += '+' + str(val) + ' '
        retlist.append('({})'.format(linestr.strip()))
    return '\n'.join(retlist)

# 286-2
def greedysorted(permutation):
    for index, element in enumerate(permutation):
        if element != index + 1:
            return False
    return True

# 286-2
def greedysorting_find(permutation, goal):
    for index, element in enumerate(permutation):
        if abs(element) == abs(goal):
            return index
    return None

# 286-2
def greedysorting_reverse(permutation, from_index, to_idx):
    pre_middle = permutation[from_index:to_idx + 1]
    to_middle = []
    for element in pre_middle:
        to_middle.insert(0, -element)
    return permutation[:from_index] + to_middle + permutation[to_idx + 1:]

# 286-2
def greedysorting(permutation):
    retlist = []
    currperm = list(permutation)
    while (not greedysorted(currperm)):
        for index, element in enumerate(currperm):
            if element != index + 1:
                if element == -(index + 1):
                    currperm[index] = -element
                else:
                    dest_idx = greedysorting_find(currperm, index + 1)
                    currperm = greedysorting_reverse(currperm, index, dest_idx)
                retlist.append(list(currperm))  # copy to new list so it won't be overwritten
                break
    return retlist
