
from SQRNdbnseq import DBNToPairs, PairsToDBN, UnAlign, ParseRestraints, BPMatrix, AnnotateStems, SQRNdbnseq
from multiprocessing import Pool

def SQRNdbnali(objs, defrests = None, defreacts = None, defref = None,
               bpweights = {}, interchainonly = False,
               minlen = 2, minbpscore = 0,
               threads = 1, verbose = False):

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

    
    return dbns, [np.nan]*6


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
        
    pred_dbns = SQRNdbnali(objs, defrests, defreacts, defref, 
                           bpweights, interchainonly,
                           minlen, minbpscore,
                           threads, verbose)



