
from multiprocessing import Pool
import numpy as np
import io

from SQRNdbnseq import DBNToPairs, PairsToDBN, UnAlign, ParseRestraints,\
                       BPMatrix, AnnotateStems, RunSQRNdbnseq, GAPS,\
                       EncodedReactivities, SEPS


def ReAlignDict(shortseq, longseq):
    """ Return the dictionary between unaligned and realigned indices"""

    dct = {}

    i1, i2 = 0, 0

    while i1 < len(shortseq):

        if longseq[i2] in ('.','-','~'):
            i2 += 1
        else:
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

    # turn seq into UPPERCASE & replace T with U
    seq = seq.upper().replace("T", "U")

    # if not restraints - generate a string of dots
    if not restraints:
        restraints = '.'*len(seq)

    assert len(seq) == len(restraints)
    
    # Unalign seq, reacts and restraints strings
    shortseq, shortrest = UnAlign(seq, restraints)

    if reactivities:
        shortreacts = [reactivities[i] for i in range(len(seq)) if seq[i] not in GAPS]
    else:
        shortreacts = None
    
    # Parse restraints into unpaired bases (rxs) and base pairs (rbps)
    rbps, rxs, rlefts, rrights = ParseRestraints(shortrest)

    bpboolmatrix, bpscorematrix = BPMatrix(shortseq, bpweights, rxs,
                                           rlefts, rrights,
                                           interchainonly)

    stems = AnnotateStems(bpboolmatrix, bpscorematrix, rbps,
                          rstems=[], minlen=minlen, minscore=minbpscore,
                          diff=0)

    radict = ReAlignDict(shortseq, seq)

    # Weight the score with reactivities if needed
    if shortreacts:
        for stem in stems:
            reactfactor = sum(1 - shortreacts[pos]
                              for bp in stem[0]
                              for pos in bp) / len(stem[0])
            stem[-1] *= reactfactor

    return [[[(radict[v],radict[w]) for v,w in stem[0]], stem[-1]]
            for stem in stems if stem[-1] >= minbpscore]


def mpYieldStems(args):

    seq, reacts, rests, bpweights, interchainonly, minlen, minbpscore = args
    return YieldStems(seq, reacts, rests,
                      bpweights, interchainonly,
                      minlen, minbpscore)


def MatrixToDBNs(mat, score, depth, verbose = False):

    N = mat.shape[0]

    thr = score*depth
    
    dct   = dict(enumerate(mat.flatten()))
    cells = sorted(dct.items(), key = lambda x: x[1], reverse = True)

    res = [[[],set()],]

    if verbose:
        print(">Conserved base pairs (one by one)")

    for cell in cells:

        bp = np.unravel_index(cell[0],mat.shape)

        if cell[1] < thr:
            break

        if not bp[1]-bp[0] >= 4:
            continue

        added = False

        for struct in res:
            if bp[0] not in struct[1] and bp[1] not in struct[1]:
                struct[0].append(bp)
                struct[1].add(bp[0])
                struct[1].add(bp[1])
                added = True
                break
        if not added:
            res.append([[bp,],set(bp)])

        if verbose:
            print(PairsToDBN([bp,],N),round(cell[1], 3),sep='\t')

    dbns = [PairsToDBN(struct[0], N) for struct in res]

    if verbose:
        print(">Conserved base pairs (assembled)")
        for dbn in dbns:
            print(dbn)

    return dbns


def Metrics(ref, pred):

    if not ref:
        return [np.nan]*6

    rb = set(DBNToPairs(ref))
    pb  = set(DBNToPairs(pred))
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
               threads = 1, verbose = False):

    stemmatrix = np.zeros((len(objs[0][1]),len(objs[0][1])))

    with Pool(threads) as pool:

        inputs = [(obj[1],                           # Sequence
                   obj[2],                           # Reactivities
                   defrests if defrests else obj[3], # Restraints
                   bpweights,
                   interchainonly,
                   minlen,
                   minbpscore) for obj in objs]
                    
        for stems in pool.imap(mpYieldStems, inputs):
            for stem in stems:
                for v,w in stem[0]:
                    stemmatrix[v, w] += stem[-1]
                    stemmatrix[w, v] += stem[-1]

    pred_dbns = MatrixToDBNs(stemmatrix, minbpscore, len(objs), verbose)
    return pred_dbns[0]


def mpRunSQRNdbnseq(args):

    obj, step1, paramsetnames, paramsets, threads,\
    rankbydiff, rankby, hardrest, interchainonly,\
    toplim, outplim, conslim, reactformat, verbose = args

    name, seq, reacts, restrs, ref = obj

    with io.StringIO() as buffer:
    
        cons, dbns, nums1, nums2 = RunSQRNdbnseq(name, seq, reacts, step1, ref, paramsetnames,
                                                 paramsets, threads, rankbydiff, rankby,
                                                 hardrest, interchainonly, toplim, outplim,
                                                 conslim, reactformat,
                                                 mp = False, sink = buffer)
        return cons, buffer.getvalue()


def Consensus(structs, freqlimit = 0.0, verbose = False):

    bps = {}

    freqlimit *= len(structs)

    for struct in structs:
        for bp in DBNToPairs(struct):
            if bp not in bps:
                bps[bp] = 0
            bps[bp] += 1

    resbps = []

    seen = set()

    if verbose:
        print(">Step 2, Populated base pairs")

    for bp in sorted(bps.keys(), key = lambda x: bps[x], reverse = True):
        if verbose:
            print(PairsToDBN([bp,],len(structs[0])), bps[bp])
        if bps[bp] >= freqlimit and bp[0] not in seen and bp[1] not in seen:
            seen.add(bp[0])
            seen.add(bp[1])
            resbps.append(bp)

    return PairsToDBN(list(set(resbps)), len(structs[0]))


def ReactScore(reacts, seq, dbn):

    if not reacts:
        return 0

    paired = set()
    for v,w in DBNToPairs(dbn):
        paired.add(v)
        paired.add(w)
        
    sepnum = sum(1 for _ in seq if _ in SEPS)

    reactscore = 1 - sum(reacts[i] if i in paired else 1 - reacts[i]
                         for i in range(len(seq))
                         if seq[i] not in SEPS) / (len(seq) - sepnum)
    return reactscore


def RunSQRNdbnali(objs, defreacts, defrests, defref,
                  levellimit, freqlimit, verbose, step3,
                  paramsetnames, paramsets, threads, rankbydiff, rankby,
                  hardrest, interchainonly, toplim, outplim,
                  conslim, reactformat):

    N = len(objs[0][1])

    bpweights  = paramsets[0]['bpweights']
    minlen     = paramsets[0]['minlen']
    minbpscore = paramsets[0]['minbpscore']

    if verbose:
        print(">Step 1, Iteration 1")
        
    pred_dbn = SQRNdbnali(objs, defrests, defreacts, defref, 
                          bpweights, interchainonly,
                          minlen, minbpscore,
                          threads, verbose)

    if verbose:
        print(">Step 1, Iteration 2")
        
    pred_dbn = SQRNdbnali(objs, pred_dbn, defreacts, defref, 
                          bpweights, interchainonly,
                          minlen, minbpscore,
                          threads, verbose)

    # Truncate pseudoknotted bps of order higher than levellimit
    step1dbn = PairsToDBN(DBNToPairs(pred_dbn), N,
                          levellimit = levellimit)

    if verbose:
        print(">Step 1, Result")
        print(step1dbn)

    structs = []
    if step3 != '1':
        if verbose:
            print(">Step 2, Individuals")

        with Pool(threads) as pool:
            inputs = [(obj, step1dbn, paramsetnames, paramsets, threads,
                       rankbydiff, rankby, hardrest, interchainonly,
                       toplim, outplim, conslim, reactformat,
                       verbose) for obj in objs]
            for cons, output in pool.imap(mpRunSQRNdbnseq, inputs):
                if verbose:
                    print(output, end = '')
                structs.append(cons)

        step2dbn = Consensus(structs, freqlimit, verbose)
        if verbose:
            print(">Step 2, Consensus")
            for lim in range(0, 101, 5):
                print(Consensus(structs, lim / 100), str(lim)+'%', sep='\t')
    else:
        step2dbn = '.' * N

    if verbose:
        print("=" * N)

    if defreacts:
        print(EncodedReactivities(objs[0][1],
                                  defreacts,
                                  reactformat),
              "reactivities", sep='\t')
    if defrests:
        print(''.join([defrests[i]
                       if   objs[0][1][i] not in SEPS
                       else objs[0][1][i]
                       for i in range(N)]), "restraints", sep='\t')
    if defref:
        print(''.join([defref[i]
                       if   objs[0][1][i] not in SEPS
                       else objs[0][1][i]
                       for i in range(N)]), "reference", sep='\t')

    if defreacts or defref or defrests:
        print("_" * N)

        
    print(step1dbn,
          "Step-1" + ('\t'+str(round(ReactScore(defreacts, objs[0][1], step1dbn), 2))) * bool(defreacts),
          "TP={},FP={},FN={},FS={},PR={},RC={}".format(*Metrics(defref, step1dbn)) * bool(defref),
          sep = '\t')

    print(step2dbn,
          "Step-2" + "(skipped)" * (step3 == '1') +\
           ('\t'+str(round(ReactScore(defreacts, objs[0][1], step2dbn), 2))) * bool(defreacts) * (step3 != '1'),
          "TP={},FP={},FN={},FS={},PR={},RC={}".format(*Metrics(defref, step2dbn))*bool(defref) * (step3 != '1'),
          sep = '\t')

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
          sep = '\t')
