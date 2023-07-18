
from multiprocessing import Pool
import numpy as np

from SQRNdbnseq import DBNToPairs, PairsToDBN, UnAlign, ParseRestraints, BPMatrix, AnnotateStems, SQRNdbnseq, GAPS


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

        if defref:
            rb = set(DBNToPairs(defref))
            pb  = set(DBNToPairs(pred_dbns[0]))
            TP  = len(pb & rb)
            FP  = len(pb - rb)
            FN  = len(rb - pb)
            PRC = (round(TP / (TP + FP), 3)) if (TP + FP) else 0
            RCL = (round(TP / (TP + FN), 3)) if (TP + FN) else 0
            FSC = (round(2*TP / (2*TP + FP + FN), 3)) if (2*TP + FP + FN) else 0
            return pred_dbns[0], [TP, FP, FN, FSC, PRC, RCL]
        
        return pred_dbns[0], [np.nan]*6


def RunSQRNdbnali(objs, defreacts, defrests, defref,
                  levellimit, freqlimit, verbose, step3,
                  paramsetnames, paramsets, threads, rankbydiff, rankby,
                  hardrest, interchainonly, toplim, outplim,
                  conslim, reactformat):

    bpweights  = paramsets[0]['bpweights']
    minlen     = paramsets[0]['minlen']
    minbpscore = paramsets[0]['minbpscore']

    if verbose:
        print(">Step 1, Iteration 1")
        
    pred_dbn, metrics = SQRNdbnali(objs, defrests, defreacts, defref, 
                                   bpweights, interchainonly,
                                   minlen, minbpscore,
                                   threads, verbose)

    if verbose:
        print(">Step 1, Iteration 2")
        
    pred_dbn, metrics = SQRNdbnali(objs, pred_dbn, defreacts, defref, 
                                   bpweights, interchainonly,
                                   minlen, minbpscore,
                                   threads, verbose)


    if verbose:
        print("="*len(objs[0][1]))
        
    if defref:
        print(pred_dbn,
              "Step-1",
              "TP={},FP={},FN={},FS={},PR={},RC={}".format(*metrics),
              sep = '\t')
    else:
        print(pred_dbn,
              "Step-1",
              sep = '\t')


