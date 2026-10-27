



    
def BackTrack(begin, end, K, minloop, seq, seps, partial = False):

    queue = {(begin, end), } # start from the last cell
    basepairs = []

    while queue:

        newq = set()

        # for each cell in the queue
        for i,j in queue:
            # if we have (k, j) base pair
            if (i,j) in K:
                k = K[(i,j)]
                # if we may have more base pairs before k  
                if ((k - 1) - i > minloop or (k-1)-i > 0 and any(ch in seps
                                                                 for ch in seq[i+1:k-1])) \
                   and not partial:
                    newq.add((i, k - 1))
                # if we may have more base pairs inside [k, j]
                if (j - 1) - (k + 1) > minloop or (j-1)-(k+1) > 0 and any(ch in seps
                                                                          for ch in seq[k+2:j-1]):
                    newq.add((k + 1, j - 1))
                # add (k, j) bp
                basepairs.append((k,j))
            # if j is unpaired
            else:
                # if we may have more base pairs inside [i, j - 1]
                if ((j - 1) - i > minloop or (j-1)-i > 0 and any(ch in seps
                                                                 for ch in seq[i+1:j-1])) \
                   and not partial:
                    newq.add((i, j - 1))
        # update queue
        queue = newq

    return sorted(basepairs)


def Nussinov(seq, stems, N, seps, minloop = 3, matrix = None):

    import numpy as np

    if matrix is None:
        SCORES = {(v,w):-stem[2] for stem in stems for v,w in stem[0]}
    else:
        SCORES = {(v,w):-matrix[v,w]
                  for v in range(N - 1)
                  for w in range(v + 1,N)
                  if matrix[v,w] > 0}

    # score matrix
    D = np.zeros((N,N))
    # base pair dictionary for traceback
    K = {}

    # going diagonale by diagonale,
    # skipping first [minloop] diagonals
    # to avoid too short hairpins
    #for h in range(minloop + 1, N):

    for h in range(1, N):
        # row index
        for i in range(N - h):
            # column index
            j = i + h

            bestk, bestscorek = -1, 10**9
            # searching for the optimal base pair
            #for k in range(i, j - minloop):
            for k in range(i, j - 1):

                if (k,j) in SCORES:
                    scorek = (D[i, k - 1] if k > i else 0) + D[k + 1, j - 1] + SCORES[(k,j)]
                    if scorek < bestscorek:
                        bestk, bestscorek = k, scorek
            # see if no base pair is not better
            if bestscorek <= D[i, j - 1]:

                K[(i, j)] = bestk
                D[i, j]   = bestscorek
            # if it is
            else:
                D[i, j] = D[i, j - 1]

    # traceback
    predictedbps = BackTrack(0, N - 1, K, minloop, seq, seps)
    
    return predictedbps


def PredictedStems(pairs, stems):
    """Map the predicted base pairs back onto the source stems:
    returns [predicted part of the stem, its length, stem score, stem index],
    each predicted base pair assigned to the best-scored stem containing it"""

    pairs    = set(pairs)
    assigned = set()
    result   = []

    for ind in sorted(range(len(stems)),
                      key = lambda x: (stems[x][2], stems[x][1]),
                      reverse = True):

        part = [bp for bp in stems[ind][0]
                if bp in pairs and bp not in assigned]
        if part:
            assigned |= set(part)
            result.append([part, len(part), stems[ind][2], ind])

    return result


def ILM2(seq, stems, N, seps, minloop = 3):
    """Two-pass Nussinov: the second pass is run without the stems
    predicted in the first pass and without the stems that cannot
    outscore the first-pass stems they overlap with, then the second-pass
    stems are added to the first-pass prediction one by one if they either
    do not overlap with it or outscore all the stems they overlap with in total"""
 
    def Weight(stem):
        return stem[1] * stem[2]
         
    def Nts(stem):
        return {pos for bp in stem[0] for pos in bp}
         
    # first iteration
    first = PredictedStems(Nussinov(seq, stems, N, seps, minloop), stems)
    for stem in first:
        stem.append(stems[stem[3]]) # keep the source stem for re-extension
          
    # second iteration without the first-iteration stems
    # and without the stems that can never replace what they overlap with
    excluded  = {stem[3] for stem in first}
    remaining = [stem for ind, stem in enumerate(stems)
                 if ind not in excluded and
                 Weight(stem) > sum(Weight(fst) for fst in first
                                    if Nts(fst) & Nts(stem))]
    second = PredictedStems(Nussinov(seq, remaining, N, seps, minloop),
                            remaining)
    for stem in second:     
        stem.append(remaining[stem[3]])
          
    # merging, heaviest second-iteration stems first
    current = first[:]
  
    for stem in sorted(second, key = Weight, reverse = True):
      
        nts = Nts(stem)

        overlapping = [ind for ind, curstem in enumerate(current)
                       if Nts(curstem) & nts]

        if not overlapping:
            current.append(stem)
        elif Weight(stem) > sum(Weight(current[ind]) for ind in overlapping):
            current = [curstem for ind, curstem in enumerate(current)
                       if ind not in overlapping] + [stem]
              
    # re-extending the stem fragments to their source stems
    # wherever the nucleotides are still free
    used = set().union(*(Nts(stem) for stem in current))

    for stem in sorted(current, key = Weight, reverse = True):
        extra = [(v,w) for v,w in stem[4][0]
                 if (v,w) not in stem[0] and v not in used and w not in used]
        stem[0] = sorted(stem[0] + extra)
        stem[1] = len(stem[0])
        used |= {pos for bp in extra for pos in bp}

    pairs = sorted(bp for stem in current for bp in stem[0])
         
    return pairs



def Edmonds(stems, power = 1.7, matrix = None):
    
    import networkx as nx

    if matrix is None:
        edges = [(v,w,stem[2]**power) for stem in stems for v,w in stem[0]]
    else:
        edges = [(v,w,matrix[v,w]**power) for v in range(matrix.shape[0] - 1)
                                   for w in range(v + 1, matrix.shape[0])
                                   if matrix[v,w] > 0]

    G = nx.Graph()
    G.add_weighted_edges_from(edges)
    pairs = sorted(nx.max_weight_matching(G))
    return pairs


def Hungarian(seq, stems, N, seps, minloop = 3, power = 1.7, matrix = None):

    import numpy as np
    from scipy.optimize import linear_sum_assignment

    if matrix is None:
        mat = np.zeros((N,N))
        for stem in stems:
            for v,w in stem[0]:
                mat[v,w] = -(stem[2]**power)
                mat[w,v] = -(stem[2]**power)
    else:
        mat = -(matrix**power)

    row_ind, col_ind = linear_sum_assignment(mat)
    sol = {i:j for i,j in zip(row_ind,col_ind)}

    pairs = [(k,sol[k]) for k in sol
             if (k < sol[k]-minloop or k < sol[k] and any(ch in seps
                                                          for ch in seq[k+1:sol[k]])) \
                and sol[k] in sol and sol[sol[k]] == k and mat[k,sol[k]] != 0]

    return pairs
    
