
import numpy as np

from SQRNdbnseq import DBNToPairs, PairsToDBN, UnAlign, ParseRestraints, BPMatrix, AnnotateStems, nonmpSQRNdbnseq
from multiprocessing import Pool


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



   
def YieldStems(seq, restraints = None, dbn = None, 
               bpweights = {}, interchainonly = False,
               minlen=2, minscore=0):

    # turn seq into UPPERCASE & replace T with U
    seq = seq.upper().replace("T", "U")

    # if not restraints - generate a string of dots
    if not restraints:
        restraints = '.'*len(seq)

    assert len(seq) == len(restraints)
    
    # Unalign seq, dbn, and restraints strings
    shortseq, shortrest = UnAlign(seq, restraints)

    if dbn:
        assert len(seq) == len(dbn)
        shortseq, shortdbn  = UnAlign(seq, dbn)

    predbps = DBNToPairs(shortdbn)
    
    # Parse restraints into unpaired bases (rxs) and base pairs (rbps)
    rbps, rxs, rlefts, rrights = ParseRestraints(shortrest)

    bpboolmatrix, bpscorematrix = BPMatrix(shortseq, bpweights, rxs,
                                           rlefts, rrights,
                                           interchainonly)

    stems = AnnotateStems(bpboolmatrix, bpscorematrix, rbps + predbps,
                          rstems=[], minlen=minlen, minscore=minscore, diff=0)

    radict = ReAlignDict(shortseq, seq)

    return [[[(radict[v],radict[w]) for v,w in stem[0]], stem[-1]] for stem in stems]


def Alt(mat, score, depth):

    N = mat.shape[0]

    thr = score*depth
    
    dct   = dict(enumerate(mat.flatten()))
    cells = sorted(dct.items(), key = lambda x: x[1], reverse = True)

    res = [[[],set()],]

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

    #for struct in res:
    #    print(PairsToDBN(struct[0],N))

    return [PairsToDBN(struct[0],N) for struct in res]


def mpYieldStems(args):

    name, seq, rst, consrest, dbn, bpweights,interchainonly,minlen,minscore = args

    return YieldStems(seq, restraints = rst if rst else consrest,
                      dbn = dbn,
                      bpweights=bpweights,
                      interchainonly=interchainonly,
                      minlen=minlen,
                      minscore=minscore)


def SQRNdbnali(objs, consrest = None, ref = None,
               bpweights = {}, interchainonly = False,
               minlen=2, minscore=0, limitscore = 0,
               iterative = False, threads = 1):

    dbn = '.'*len(objs[0][1])

    best_fscore = -1
    best_ls = -1

    rb = set(DBNToPairs(ref))

    seen_preds = set()

    first = True

    while True:

        stemmatrix = np.zeros((len(queue[0][1]),len(queue[0][1])))

        with Pool(24) as pool:

            inputs = [(obj[0],obj[1],obj[2],consrest,dbn,bpweights,interchainonly,minlen,minscore) for obj in queue]
                    
            for stems in pool.imap(mpYieldStems,inputs):
                for stem in stems:
                    for v,w in stem[0]:
                        stemmatrix[v, w] += stem[-1]
                        stemmatrix[w, v] += stem[-1]

        predbps = DBNToPairs(dbn)

        return Alt(stemmatrix, limitscore, len(objs))

        #########################
        if iterative:

            preds = Alt(stemmatrix, limitscore, len(objs))

            if preds[0] == consrest or preds[0] in seen_preds:
                print('last')
                return preds
            else:
                pb = set(DBNToPairs(preds[0]))
                TP = len(pb & rb)
                FP = len(pb - rb)
                FN = len(rb - pb)
                FS = 2*TP / (2*TP + FP + FN)
                print('next', round(FS,3), preds[0])
                seen_preds.add(preds[0])
                consrest = preds[0]
                if first:
                    first = False
                    continue
                else:
                    return preds
        #######################

        for v,w in predbps:
            stemmatrix[v, w] = 0
            stemmatrix[w, v] = 0

        mx = stemmatrix.max()

        

        if mx < limitscore * len(objs):
            break

        newbp = tuple(sorted(np.unravel_index(stemmatrix.argmax(),
                                              stemmatrix.shape)))

        dbn = PairsToDBN(predbps + [newbp,], len(dbn))
        print(dbn)

        pb = set(DBNToPairs(dbn))

        TP = len(pb & rb)
        FP = len(pb - rb)
        FN = len(rb - pb)

        FS = 2*TP / (2*TP + FP + FN)

        if FS > best_fscore:

            best_ls = mx / len(objs)
            best_fscore = FS
        
        print(mx, round(FS,3))
        
        #if len(predbps)==1:
        #    PrintMatrix('N'*len(queue[0][1]), stemmatrix, dbn1 = ref)
                       
    return dbn, round(best_ls,3), round(best_fscore,3), round(FS,3)


def Consensus(rst, structs, thr = 0.0, truncate = True):

    bps = {}

    thr *= len(structs)

    for struct in structs:
        for bp in DBNToPairs(struct):
            if bp not in bps:
                bps[bp] = 0
            bps[bp] += 1

    resbps = []

    seen = set()

    for bp in sorted(bps.keys(),key = lambda x: bps[x], reverse = True):

        if bps[bp] >= thr and bp[0] not in seen and bp[1] not in seen:
            seen.add(bp[0])
            seen.add(bp[1])
            resbps.append(bp)

    if truncate:
        return PairsToDBN(list(set(resbps) & set(DBNToPairs(rst))), len(structs[0]))
    return PairsToDBN(list(set(resbps)), len(structs[0]))

def mpSQRNdbnseq(args):

    name, seq, rst,pred,paramset,threads,interchainonly = args
    cons, dbns, nums1, nums2 = nonmpSQRNdbnseq(seq, reacts = None, restraints = pred, dbn = None,
                                              paramsets = [paramset,], conslim = 1,
                                              rankbydiff = False, interchainonly = interchainonly,
                                              threads = threads)
    return cons
    
if __name__ == "__main__":

    bpweights = {'GU' :  -1.25,
                 'AU' :  1.25,
                 'GC' :  3.25,}
    
    threads = 1

    threads1 = 32

    interchainonly = False

    queue = []

    with open("rfam10/afa/RF00023.afa") as file:
        lines = file.readlines()

        ref = lines[0].strip()

        for ii in range(1,len(lines)-1,2):

            nm = lines[ii].strip()[1:]
            sq = lines[ii+1].strip()
            queue.append([nm, sq, None])

    iters = True

    for bpweights in ({'GU' :  -1.25,'AU' :  1.25,'GC' :  3.25,},):

        minlen     = 2
        minscore   = bpweights['GC']+bpweights['AU']
        limitscore = bpweights['GC']+bpweights['AU']

        preds = SQRNdbnali(queue, consrest = None, ref = ref,
                           bpweights = bpweights, interchainonly = interchainonly,
                           minlen=minlen, minscore=minscore,
                           limitscore = limitscore, iterative = iters, threads = threads1)
        if iters:
            preds = SQRNdbnali(queue, consrest = preds[0], ref = ref,
                               bpweights = bpweights, interchainonly = interchainonly,
                               minlen=minlen, minscore=minscore,
                               limitscore = limitscore, iterative = iters, threads = threads1)
        
        preds[0] = PairsToDBN(DBNToPairs(preds[0]),len(ref),limitlevel=3 - int(len(ref)>500))####

        rb = set(DBNToPairs(ref))
        pb = set(DBNToPairs(preds[0]))

        TP = len(pb & rb)
        FP = len(pb - rb)
        FN = len(rb - pb)

        FS = 2*TP / (2*TP + FP + FN)
        print(preds[0],round(FS,3))

        prev_fs = FS

        pred = preds[0]

        for orderpen in (0.75,):
        
            paramset  = {"bpweights" :  bpweights,
                        "suboptmax": 0.99,
                        "suboptmin": 0.99,
                        "suboptsteps": 1,
                        "minlen": minlen,
                        "minbpscore": minscore,
                        "minfinscorefactor": 1.0,
                        "distcoef" : 0.09,
                        "bracketweight": -2,
                        "orderpenalty": orderpen,
                        "loopbonus": 0.125,
                        "maxstemnum": 10**6}

            structs = []

            with Pool(threads1) as pool:

                inputs = [(obj[0],obj[1],obj[2],pred,paramset,threads,interchainonly) for obj in queue]
                    
                for cons in pool.imap(mpSQRNdbnseq,inputs):

                    structs.append(cons)

            #print("Consensus")

                for thr in (0.35,):
                    for truncate in (False,):

                        consensus = Consensus(pred, structs, thr, truncate)

                        rb = set(DBNToPairs(ref))
                        pb = set(DBNToPairs(consensus))

                        TP = len(pb & rb)
                        FP = len(pb - rb)
                        FN = len(rb - pb)

                        FS = 2*TP / (2*TP + FP + FN)

                        predpairs = set(DBNToPairs(pred))
                        seenpos = set(x for bp in predpairs for x in bp)
                        for bp in DBNToPairs(consensus):
                            if bp[0] not in seenpos and bp[1] not in seenpos:
                                predpairs.add(bp)
                        TP2 = len(predpairs & rb)
                        FP2 = len(predpairs - rb)
                        FN2 = len(rb - predpairs)
                        FS2 = 2*TP2 / (2*TP2 + FP2 + FN2)

                        res = PairsToDBN(sorted(predpairs),len(ref))

                        #print(pred, ls, bfscore, fscore)
                        print(consensus, bpweights['GU'], orderpen, round(thr,3), truncate, round(prev_fs,3),
                              round(FS,3), round(FS2,3), sep='\t')# )
                        #print(altbps)
                        print(res)
    

    






