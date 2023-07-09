
import numpy as np

from SQRNdbnseq import DBNToPairs, PairsToDBN, UnAlign, ParseRestraints, BPMatrix, AnnotateStems, SQRNdbnseq


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

    for stem in stems:
        yield [[(radict[v],radict[w]) for v,w in stem[0]], stem[-1]]


def SQRNdbnali(objs, consrest = None, ref = None,
               bpweights = {}, interchainonly = False,
               minlen=2, minscore=0, limitscore = 0):

    dbn = '.'*len(objs[0][1])

    best_fscore = -1
    best_ls = -1

    rb = set(DBNToPairs(ref))

    while True:

        stemmatrix = np.zeros((len(queue[0][1]),len(queue[0][1])))

        for obj in objs:

            name, seq, rst = obj

            for stem in YieldStems(seq, restraints = rst if rst else consrest,
                                   dbn = dbn,
                                   bpweights=bpweights,
                                   interchainonly=interchainonly,
                                   minlen=minlen,
                                   minscore=minscore):

                for v,w in stem[0]:
                    stemmatrix[v, w] += stem[-1]
                    stemmatrix[w, v] += stem[-1]

        predbps = DBNToPairs(dbn)

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


def Consensus(rst, structs):

    bps = {}

    for struct in structs:
        for bp in DBNToPairs(struct):
            if bp not in bps:
                bps[bp] = 0
            bps[bp] += 1

    resbps = []

    seen = set()

    for bp in sorted(bps.keys(),key = lambda x: bps[x], reverse = True):

        if bp[0] not in seen and bp[1] not in seen:
            seen.add(bp[0])
            seen.add(bp[1])
            resbps.append(bp)
        
    return PairsToDBN(list(set(resbps) & set(DBNToPairs(rst))), len(structs[0]))


if __name__ == "__main__":

    bpweights = {'GU' :  -1.25,
                 'AU' :  1.25,
                 'GC' :  3.25,}
    minlen   = 2
    minscore = 0
    limitscore = bpweights['GC']+bpweights['AU']
    threads = 1

    interchainonly = False


    queue = []

    with open("rfam10/RF00001.afa") as file:
        lines = file.readlines()

        ref = lines[0].strip()

        for ii in range(1,len(lines)-1,2):

            nm = lines[ii].strip()[1:]
            sq = lines[ii+1].strip()
            queue.append([nm, sq, None])

    pred, ls, bfscore, fscore = SQRNdbnali(queue, consrest = None, ref = ref,
                                           bpweights = bpweights, interchainonly = interchainonly,
                                           minlen=minlen, minscore=minscore,
                                           limitscore = limitscore)
    print(ls, bfscore, fscore)
    
    paramset  = {"bpweights" :  {'GU' :  -1.25,
                 'AU' :  1.25,
                 'GC' :  3.25,},
              "suboptmax": 0.9,
              "suboptmin": 0.65,
              "suboptsteps": 1,
              "minlen": 2,
              "minbpscore": 4.5,
              "minfinscorefactor": 1.0,
              "distcoef" : 0.09,
              "bracketweight": -2,
              "orderpenalty": 1.0,
              "loopbonus": 0.125,
              "maxstemnum": 10**6}

    structs = []

    for obj in queue:

        name, seq, rst = obj
        print(seq)
        cons, dbns, nums1, nums2 = SQRNdbnseq(seq, reacts = None, restraints = pred, dbn = None,
                                              paramsets = [paramset,], conslim = 1,
                                              rankbydiff = False, interchainonly = interchainonly,
                                              threads = threads)
        print(cons)
        structs.append(cons)

    print("Consensus")

    consensus = Consensus(pred, structs)

    rb = set(DBNToPairs(ref))
    pb = set(DBNToPairs(consensus))

    TP = len(pb & rb)
    FP = len(pb - rb)
    FN = len(rb - pb)

    FS = 2*TP / (2*TP + FP + FN)

    print(pred, ls, bfscore, fscore)
    print(consensus, round(FS,3))

    






