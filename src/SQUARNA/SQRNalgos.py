



    
def BackTrack(begin, end, K, minloop, partial = False):

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
                if (k - 1) - i > minloop and not partial:
                    newq.add((i, k - 1))
                # if we may have more base pairs inside [k, j]
                if (j - 1) - (k + 1) > minloop:
                    newq.add((k + 1, j - 1))
                # add (k, j) bp
                basepairs.append((k,j))
            # if j is unpaired
            else:
                # if we may have more base pairs inside [i, j - 1]
                if (j - 1) - i > minloop and not partial:
                    newq.add((i, j - 1))
        # update queue
        queue = newq

    return sorted(basepairs)


def Nussinov(stems, N, minloop = 3, matrix = None):

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
    for h in range(minloop + 1, N):
        # row index
        for i in range(N - h):
            # column index
            j = i + h

            bestk, bestscorek = -1, 10**9
            # searching for the optimal base pair
            for k in range(i, j - minloop):

                if (k,j) in SCORES:
                    scorek = D[i, k - 1] + D[k + 1, j - 1] + SCORES[(k,j)]
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
    predictedbps = BackTrack(0, N - 1, K, minloop)
    
    return predictedbps


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


def Hungarian(stems, N, minloop = 3, power = 1.7, matrix = None):

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
             if k < sol[k]-minloop and sol[k] in sol and sol[sol[k]] == k and mat[k,sol[k]] != 0]

    return pairs
    
