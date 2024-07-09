
from multiprocessing import Pool
import numpy as np
import io
import sys

try:
    from SQRNdbnseq import DBNToPairs, PairsToDBN, UnAlign, ParseRestraints,\
                           BPMatrix, AnnotateStems, RunSQRNdbnseq, GAPS,\
                           EncodedReactivities, SEPS
except:
    from .SQRNdbnseq import DBNToPairs, PairsToDBN, UnAlign, ParseRestraints,\
                            BPMatrix, AnnotateStems, RunSQRNdbnseq, GAPS,\
                            EncodedReactivities, SEPS


def ReAlignDict(shortseq, longseq):
    """ Return the dictionary between unaligned and realigned indices"""

    dct = {}

    i1, i2 = 0, 0

    while i1 < len(shortseq):

        # Shift the longseq index if gap
        if longseq[i2] in GAPS:
            i2 += 1
        else: # otherwise save the pair of matching indices and increment
            dct[i1] = i2
            i1 += 1
            i2 += 1

    return dct


def PrintMatrix(seq, matrix, dbn1='', dbn2=''):
    """Print matrix and sequence into tsv format"""
    print('',*list(seq),sep='\t')

    Pairs1 = DBNToPairs(dbn1) # bps from dbn1 will be framed with []
    Pairs2 = DBNToPairs(dbn2) # bps from dbn2 will be framed with <>

    for i in range(len(seq)):
        print(seq[i], end='\t')
        line = []
        for j in range(len(seq)):
            x = str(matrix[i][j])
            if (i,j) in Pairs1:
                x = '['+x+']'
            if (i,j) in Pairs2:
                x = '<'+x+'>'
            line.append(x)
        print(*line, sep='\t')


def YieldStems(seq, reactivities = None, restraints = None, 
               bpweights = {}, interchainonly = False,
               minlen = 2, minbpscore = 0):
    """Returns a list of stems based on a (possibly aligned) sequence"""

    # turn seq into UPPERCASE & replace T with U
    seq = seq.upper().replace("T", "U")

    # if not restraints - generate a string of dots
    if not restraints:
        restraints = '.'*len(seq)

    # assert the restraints are of matching length
    assert len(seq) == len(restraints)
    
    # Unalign seq, reactivities and restraints strings
    shortseq, shortrest = UnAlign(seq, restraints)

    if reactivities:
        shortreacts = [reactivities[i]
                       for i in range(len(seq))
                       if seq[i] not in GAPS]
    else:
        shortreacts = None
    
    # Parse restraints into unpaired bases (rxs) and base pairs (rbps)
    rbps, rxs, rlefts, rrights = ParseRestraints(shortrest)

    # Build the BP matrices based on the unaligned seq and restraints
    bpboolmatrix, bpscorematrix = BPMatrix(shortseq, bpweights, rxs,
                                           rlefts, rrights,
                                           interchainonly,
                                           reacts = shortreacts)

    # Annotate stems from the matrices and restraint bps
    # using the maxlen definition for stems (diff = 0)
    stems = AnnotateStems(bpboolmatrix, bpscorematrix,
                          rbps, rstems = [],
                          minlen = minlen,
                          minscore = minbpscore,
                          diff = 0)

    # aligned-unaligned indices dictionary
    radict = ReAlignDict(shortseq, seq)

    # Return the stems with re-aligned indices of their base pairs
    return [[[(radict[v], radict[w]) for v, w in stem[0]], stem[-1]]
            for stem in stems]


def mpYieldStems(args):
    """multiprocessing (single-parameter) version of YieldStems"""
    #Unpack args
    seq, reacts, rests, bpweights, interchainonly, minlen, minbpscore = args
    #Call YieldStems
    return YieldStems(seq, reacts, rests,
                      bpweights, interchainonly,
                      minlen, minbpscore)


def MatrixToDBNs(mat, score, depth, verbose = False, sink = sys.stdout):
    """Build dbn structures from the stem-scored matrix
       in a greedy manner"""
    N = mat.shape[0]

    # the bp score threshold
    thr = score*depth
    
    dct   = dict(enumerate(mat.flatten())) # 1-index key: score value
    # Cells sorted in the decreasing order of their scores
    cells = sorted(dct.items(), key = lambda x: x[1], reverse = True)

    # List of [list of base pairs, set of paired positions] lists
    res = [[[],set()],]

    if verbose:
        print(">Conserved base pairs (one by one)", file = sink)

    for cell in cells:
        # Obtain the 2-index base pair from the 1-index
        # E.g., for a square 2x2 matrix the scheme looks like this:
        # 0 -> (0, 0); 1 -> (0, 1); 2 -> (1, 0); 3 -> (1, 1)
        bp = np.unravel_index(cell[0], mat.shape)

        # Ignore all the bps below the threshold
        if cell[1] < thr:
            break

        # Ignore (v, w) bps with v > w as well as
        # the ones forming too short hairpins of <3nt
        if not bp[1]-bp[0] >= 4:
            continue

        # Whether the base pair was added to the existing struct
        # or is in conflict with everything and requires a new struct
        added = False

        for struct in res:
            if bp[0] not in struct[1] and bp[1] not in struct[1]:
                struct[0].append(bp) # Adding the base pair
                struct[1].add(bp[0]) # Adding its paired positions
                struct[1].add(bp[1]) #
                added = True
                break
        if not added:
            res.append([[bp,], set(bp)])

        if verbose:
            print(PairsToDBN([bp,], N), round(cell[1], 3), sep='\t', file = sink)

    # Convert bp-lists to dbns
    dbns = [PairsToDBN(struct[0], N) for struct in res]

    if verbose:
        print(">Conserved base pairs (assembled)", file = sink)
        for dbn in dbns:
            print(dbn, file = sink)
    # Return dbns
    return dbns


def Metrics(ref, pred):
    """Return the calculated metrics values for the two dbns"""
    if not ref:
        return [np.nan]*6

    rb = set(DBNToPairs(ref))   # reference base pairs
    pb  = set(DBNToPairs(pred)) # predicted base pairs
    TP  = len(pb & rb)
    FP  = len(pb - rb)
    FN  = len(rb - pb)
    PRC = (round(TP / (TP + FP), 3)) if (TP + FP) else 1
    RCL = (round(TP / (TP + FN), 3)) if (TP + FN) else 1
    FSC = (round(2*TP / (2*TP + FP + FN), 3)) if (2*TP + FP + FN) else 1
    return [TP, FP, FN, FSC, PRC, RCL]


def SQRNdbnali(objs, defrests = None, defreacts = None, defref = None,
               bpweights = {}, interchainonly = False,
               minlen = 2, minbpscore = 0,
               threads = 1, verbose = False,
               sink = sys.stdout):
    """Returns a predicted dbn for a set of aligned sequences"""

    # Placeholder LENxLEN matrix of zeros 
    stemmatrix = np.zeros((len(objs[0][1]),len(objs[0][1])))

    # multiprocessing stuff
    with Pool(threads) as pool:

        inputs = [(obj[1],                           # Sequence
                   obj[2],                           # Reactivities
                   defrests if defrests else obj[3], # Restraints
                   bpweights,
                   interchainonly,
                   minlen,
                   minbpscore) for obj in objs]

        # Sum up the stem scores over all sequences
        for stems in pool.imap(mpYieldStems, inputs):
            for stem in stems:
                for v,w in stem[0]:
                    stemmatrix[v, w] += stem[-1]
                    stemmatrix[w, v] += stem[-1]

    # Build the dbns from the assembled stem-scored matrix
    pred_dbns = MatrixToDBNs(stemmatrix, minbpscore, len(objs), verbose, sink = sink)
    # Return the first dbn and the stemmatrix
    return pred_dbns[0], stemmatrix


def mpRunSQRNdbnseq(args):
    """multiprocessing (single-parameter) version of RunSQRNdbnseq;
       here we parallelize over the sequences instead of parallelizing
       over alternative structures within each sequence prediction"""
    obj, paramsetnames, paramsets, threads,\
    rankbydiff, rankby, hardrest, interchainonly,\
    toplim, outplim, conslim, reactformat, verbose, smat, poollim = args

    evalonly = False

    name, seq, reacts, rests, ref = obj

    # We use a printing buffer so that the output is ordered
    # instead of being mixed up due to parallelization
    with io.StringIO() as buffer:
    
        cons, dbns, nums1, nums2 = RunSQRNdbnseq(name, seq, reacts, rests, ref, paramsetnames,
                                                 paramsets, threads, rankbydiff, rankby,
                                                 hardrest, interchainonly, toplim, outplim,
                                                 conslim, reactformat, evalonly, poollim,
                                                 mp = False, sink = buffer, stemmatrix = smat)
        return cons, buffer.getvalue()


def Consensus(structs, freqlimit = 0.0, verbose = False, sink = sys.stdout):
    """Builds a Consensus dbn from a list of dbns"""
    bps = {}

    # In how many sequences a base pair
    # has to be present to be added
    # to the consensus
    freqlimit *= len(structs)

    # Count numbers of sequences for each bp
    for struct in structs:
        for bp in DBNToPairs(struct):
            if bp not in bps:
                bps[bp] = 0
            bps[bp] += 1

    resbps = []

    seen = set()

    if verbose:
        print(">Step 2, Populated base pairs", file = sink)

    # Greedily assemble the set of the most frequent non-conflicting bps
    for bp in sorted(bps.keys(), key = lambda x: bps[x], reverse = True):
        if verbose:
            print(PairsToDBN([bp,], len(structs[0])), bps[bp], file = sink)
        if bps[bp] >= freqlimit and bp[0] not in seen and bp[1] not in seen:
            seen.add(bp[0])
            seen.add(bp[1])
            resbps.append(bp)

    # Convert the bps list to dbn and return
    return PairsToDBN(list(set(resbps)), len(structs[0]))


def ReactScore(reacts, seq, dbn):
    """Returns the reactivity score for a predicted structure"""
    if not reacts:
        return 0.5

    # set of paired positions
    paired = set()
    for v,w in DBNToPairs(dbn):
        paired.add(v)
        paired.add(w)

    # how many chain separators we have in the sequence
    sepnum = sum(1 for _ in seq if _ in SEPS)

    # calculate the score as (1 - average-error)
    # where the error for each position is calculated
    # as reactivity for paired positions
    # and (1 - reactivity) for unpaired positions
    # the score ranges from 0.0 (worst) to 1.0 (best)
    reactscore = 1 - sum(reacts[i] if i in paired else 1 - reacts[i]
                         for i in range(len(seq))
                         if seq[i] not in SEPS) / (len(seq) - sepnum)
    return reactscore


def RunSQRNdbnali(objs, defreacts, defrests, defref,
                  levellimit, freqlimit, verbose, step3,
                  paramsetnames, paramsets, threads, rankbydiff, rankby,
                  hardrest, interchainonly, toplim, outplim,
                  conslim, reactformat, poollim, sink = sys.stdout):

    # Alignment length is derived as the length
    # of the first aligned sequence
    N = len(objs[0][1])

    bpweights  = paramsets[0]['bpweights']
    minlen     = paramsets[0]['minlen']
    minbpscore = paramsets[0]['minbpscore']

    if verbose:
        print(">Step 1, Iteration 1", file = sink)

    # Iteration 1    
    pred_dbn, smat = SQRNdbnali(objs, defrests, defreacts, defref, 
                                bpweights, interchainonly,
                                minlen, minbpscore,
                                threads, verbose, sink = sink)

    if verbose:
        print(">Step 1, Iteration 2", file = sink)

    # Iteration 2 - here we use the pred_dbn from the
    # previous iteration instead of defrests
    pred_dbn = SQRNdbnali(objs, pred_dbn, defreacts, defref, 
                          bpweights, interchainonly,
                          minlen, minbpscore,
                          threads, verbose, sink = sink)[0]

    # Truncate pseudoknotted bps of order higher than levellimit for the step1 structure
    step1dbn = PairsToDBN(DBNToPairs(pred_dbn), N,
                          levellimit = levellimit)

    # Normalize the stemmatrix obtained from the first iteration
    smat = smat / np.max(smat) * 5

    if verbose:
        print(">Step 1, Result", file = sink)
        print(step1dbn, file = sink)

    structs = []
    if step3 != '1':
        if verbose:
            print(">Step 2, Individuals", file = sink)
        # Run stemmatrix-weighted single-sequence predictions
        with Pool(threads) as pool:
            inputs = [(obj, paramsetnames, paramsets, threads,
                       rankbydiff, rankby, hardrest, interchainonly,
                       toplim, outplim, conslim, reactformat,
                       verbose, smat, poollim) for obj in objs]
            for cons, output in pool.imap(mpRunSQRNdbnseq, inputs):
                if verbose:
                    print(output, end = '', file = sink)
                structs.append(cons) # Collect the derived consensus predictions
        # Build the overall consensus from all the single-seq consensus dbns
        step2dbn = Consensus(structs, freqlimit, verbose, sink = sink)
        if verbose:
            print(">Step 2, Consensus", file = sink)
            for lim in range(0, 101, 5):
                print(Consensus(structs, lim / 100), str(lim)+'%', sep='\t', file = sink)
    else:
        step2dbn = '.' * N

    # Truncate pseudoknotted bps of order higher than levellimit for the step2 structure
    step2dbn = PairsToDBN(DBNToPairs(step2dbn), N,
                          levellimit = levellimit)

    if verbose:
        print("=" * N, file = sink) # Separator line 1 if verbose 

    # Print default input lines if any
    if defreacts:
        print(EncodedReactivities(objs[0][1],
                                  defreacts,
                                  reactformat),
              "reactivities", sep='\t', file = sink)
    if defrests:
        print(''.join([defrests[i]
                       if   objs[0][1][i] not in SEPS
                       else objs[0][1][i]
                       for i in range(N)]), "restraints", sep='\t', file = sink)
    if defref:
        print(''.join([defref[i]
                       if   objs[0][1][i] not in SEPS
                       else objs[0][1][i]
                       for i in range(N)]), "reference", sep='\t', file = sink)

    if defreacts or defref or defrests:
        print("_" * N, file = sink) # Separator line 2 if default input lines

        
    print(step1dbn,
          "Step-1" + ('\t'+str(round(ReactScore(defreacts, objs[0][1], step1dbn), 2))) * bool(defreacts),
          "TP={},FP={},FN={},FS={},PR={},RC={}".format(*Metrics(defref, step1dbn)) * bool(defref),
          sep = '\t', file = sink)

    print(step2dbn,
          "Step-2" + "(skipped)" * (step3 == '1') +\
           ('\t'+str(round(ReactScore(defreacts, objs[0][1], step2dbn), 2))) * bool(defreacts) * (step3 != '1'),
          "TP={},FP={},FN={},FS={},PR={},RC={}".format(*Metrics(defref, step2dbn))*bool(defref) * (step3 != '1'),
          sep = '\t', file = sink)

    # Derive the step3dbn based on the step3 param
    if step3 == '1':
        step3dbn = step1dbn
    elif step3 == '2':
        step3dbn = step2dbn
    elif step3 == 'i':
        step3dbn = PairsToDBN(sorted(set(DBNToPairs(step1dbn)) & set(DBNToPairs(step2dbn))), N)
    else: # elif step3 == 'u'
        step1pairs = DBNToPairs(step1dbn)
        seen_pos = set(pos for bp in step1pairs for pos in bp)
        for v, w in DBNToPairs(step2dbn):
            if v not in seen_pos and w not in seen_pos:
                step1pairs.append((v, w))
        step3dbn = PairsToDBN(sorted(step1pairs), N)

    print(step3dbn,
          "Step-3({})".format(step3) +\
           ('\t'+str(round(ReactScore(defreacts, objs[0][1], step3dbn), 2))) * bool(defreacts),
          "TP={},FP={},FN={},FS={},PR={},RC={}".format(*Metrics(defref, step3dbn))*bool(defref),
          sep = '\t', file = sink)


    
