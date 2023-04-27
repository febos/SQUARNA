
import numpy as np

def PairsToDBN(newpairs, length = 0, returnlevels = False):
    """Convert a list of base pairs into a dbn string of the given length"""

    # Initialize the dbn string
    dbn = ['.']*length

    # Define "brackets" for 30 pseudoknot levels
    levels = ['()','[]','{}','<>','Aa','Bb','Cc','Dd','Ee','Ff','Gg',
              'Hh','Ii','Jj','Kk','Ll','Mm','Nn','Oo','Pp','Qq','Rr',
              'Ss','Tt','Uu','Vv','Ww','Xx','Yy','Zz']

    # groups of non-conflicting base pairs
    groups = [[],]

    # normalize the pairs (i.e. ensure v < w)
    pairs = set((min(v, w), max(v, w)) for v, w in newpairs)
    
    for pair in sorted(pairs):

        level = 0

        # find the minimum level where the pair is not in conflict
        # with any base pair of that level
        while any(v[0]<=pair[0]<=v[1]<=pair[1] or
                  pair[0]<=v[0]<=pair[1]<=v[1] for v in groups[level]):
            level += 1
            if level == len(groups):
                groups.append([])

        # add the pair to the determined level
        groups[level].append(pair)

    # kind of a bubble sort of the base pairs among the levels
    # to maximize the number of base pairs of the lowest levels
    # e.g. to turn (..[[[...)...]]] into [..(((...]...)))
    for times in range(len(groups)-1):
        for i in range(len(groups)-1):

            # take the base pairs of the current level that are in conflict with
            # the next level
            conflicting = [v for v in groups[i] if any(v[0]<=w[0]<=v[1]<=w[1] or
                                                       w[0]<=v[0]<=w[1]<=v[1]
                                                       for w in groups[i+1])]
            # if the next level group is smaller than 
            # the conflicting group -> we switch them 
            if len(conflicting) < len(groups[i+1]):
                groups[i]   = [p for p in groups[i]
                               if p not in conflicting] + groups[i+1]
                groups[i+1] = conflicting

    if returnlevels:
        levels = {}
        for lev, group in enumerate(groups):
            for bp in group:
                levels[bp] = lev + 1
        return levels

    # add all the pairs to the dbn string
    # according to their levels
    for i, group in enumerate(groups):
        for pair in group:
            dbn[pair[0]] = levels[i][0]
            dbn[pair[1]] = levels[i][1]
            
    return ''.join(dbn)


def StemsToDBN(stems, seq):
    """Convert a list of stems (lists of bps) into a dbn string"""
    return PairsToDBN([bp for stem in stems for bp in stem[0]],len(seq))


def DBNToPairs(dbn):
    """Convert the dbn string into a sorted list of base pairs"""
    pairs = set()

    # keys == closing brackets, values == matching opening brackets
    closing = {'>':'<',']':'[',')':'(','}':'{','a':'A','b':'B','c':'C','d':'D',
               'e':'E','f':'F','g':'G','h':'H','i':'I','j':'J','k':'K','l':'L',
               'm':'M','n':'N','o':'O','p':'P','q':'Q','r':'R','s':'S','t':'T',
               'u':'U','v':'V','w':'W','x':'X','y':'Y','z':'Z',}
    # 30 bp stacks for 30 allowed pseudoknot levels
    stack = {'<':[],'(':[],'{':[],'[':[],'A':[],'B':[],'C':[],'D':[],'E':[],
             'F':[],'G':[],'H':[],'I':[],'J':[],'K':[],'L':[],'M':[],'N':[],
             'O':[],'P':[],'Q':[],'R':[],'S':[],'T':[],'U':[],'V':[],'W':[],
             'X':[],'Y':[],'Z':[],}

    for i,v in enumerate(dbn):
        # if we observe an opening bracket
        # then add its index into the matching stack
        if v in stack: 
            stack[v].append(i)
        # else if we observe the closing bracket
        # take the opening index from the matching stack
        # and add the base pair to the pairs set
        elif v in closing:
            # this is to handle closing brackets with no
            # opening partner - they will be ignored
            if stack[closing[v]]:
                pairs.add((stack[closing[v]].pop(), i))

    return sorted(pairs)


def ReAlign(shortdbn, longseq, seqmode=False):
    """ Realign the shortdbn string according to the longseq string
    if seqmode==True the gaps will be hyphens, otherwise - dots"""

    # convert shortdbn into a list
    dbn = [x for x in shortdbn]

    assert len(shortdbn) + longseq.count('.') + longseq.count('-') == len(longseq),\
    "Cannot ReAlign dbn string - wrong number of gaps:\n{}\n{}".format(longseq,shortdbn)

    # initialize the returning string
    newdbn = []

    for x in longseq:
        if x in ('.','-','~'): # if gap (hyphen/dot/tilda) -> add gap
            if seqmode:
                newdbn.append('-')
            else:
                newdbn.append('.')
        else: # otherwise -> add the next symbol
            newdbn.append(dbn.pop(0))

    # return the aligned string
    return ''.join(newdbn)


def UnAlign2(seq, dbn):
    """Removes gaps from a pair of seq & dbn strings
    hyphens, dots, and tildas are treated as gaps"""

    # Unalign the dbn based on the seq string first
    newdbn = ''.join([dbn[i] for i in range(len(seq))
                      if seq[i] not in ('-','.','~')])
    # Then unalign the seq string itself
    newseq = ''.join([x for x in seq if x not in ('-','.','~')])
    # return the unaligned seq & dbn strings
    return newseq, newdbn


def BPMatrix(seq, weights):
    """Return a matrix with bp-weights in cells"""

    # Set of different symbols in seq
    bases = set(seq)

    # Dictionary of all possible bp weights
    bps = {}

    # for (x,y) from weights add bps[(x,y)] & bps[(y,x)]
    for bp in weights:
        bps[bp] = weights[bp]
        bps[bp[1]+bp[0]] = weights[bp]

    # add 0 values for all possible base pairs
    # not yet in bps dictionary
    for b1 in bases:
        for b2 in bases:
            if b1+b2 not in bps:
                bps[b1+b2] = 0

    # Initialize the matrix
    mat = np.zeros((len(seq), len(seq)), dtype = float)

    # Fill the upper triangle of the matrix with bp weights
    # Ignore the bps that would form a hairpin of len < 3
    for i in range(len(seq)-1):
        for j in range(i+4,len(seq)):
            mat[i,j] = bps[seq[i]+seq[j]]

    return mat


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


def ParseRestraints(restraints):
    """ Convert a dbn string into base pairs and unpaired bases"""
    rbps = DBNToPairs(restraints) # Base pairs
    rxs = [i for i in range(len(restraints)) if restraints[i]=='x'] # Unpaired bases

    return rbps, rxs


def AnnotateStems(bpmatrix, rbps, rxs, rstems, minlen, minscore):

    # copy the bpmatrix
    matrix = bpmatrix.copy()
    # obtain the matrix length
    N = bpmatrix.shape[0]
    # initialize result
    stems = []

    # apply bp-restraints - zero out the entire row and column except the cell for each bp
    for v, w in rbps:
        matrix[v, :] *= 0
        matrix[:, v] *= 0
        matrix[w, :] *= 0
        matrix[:, w] *= 0
        matrix[v, w]  = bpmatrix[v, w]

    # apply stem-restraints - zero out the entire rows and columns for each stem
    for stem in rstems:
        for v, w in stem[0]:
            matrix[v, :] *= 0
            matrix[:, v] *= 0
            matrix[w, :] *= 0
            matrix[:, w] *= 0

    # apply unpaired-restraints - zero out the entire row and column for each x
    for i in rxs:
        matrix[i,:] *= 0
        matrix[:,i] *= 0

    # start cells of all the diagonals, except the ones producing hairpins of len < 3
    diagstarts  = [(0, x)   for x in range(4, N)]
    diagstarts += [(y, N-1) for y in range(1, N-4)]

    # for each diagonal select the stems
    # as maximal score sub-arrays
    for x, y in diagstarts:

        curstem  = []
        curscore = 0
        bestscore = 0

        i, j = x, y
        
        while i <= j - 4: # this is to forbid hairpins of len < 3
            
            val = matrix[i, j]

            if val > 0: # if bp has positive score - add it to the current stem                  
                curstem.append((i, j))
                lastpositive = len(curstem)
                curscore = curscore + val
                bestscore = curscore

            elif val < 0: # if bp has negative score
                if curstem:
                    curscore = curscore + val
                    if curscore <= 0: # if it is better to stop the stem before the last negatives
                        if lastpositive >= minlen and bestscore >= minscore:
                            stems.append([curstem[:lastpositive], lastpositive, bestscore])
                        curstem = []
                        curscore = 0
                    else: # otherwise - growing the stem
                        curstem.append((i,j))
            elif curstem: # if no bp but we had the stem before
                if lastpositive >= minlen and bestscore >= minscore:
                    stems.append([curstem[:lastpositive], lastpositive, bestscore])
                curstem = []
                curscore = 0
            
            i += 1
            j -= 1

        # in case we stopped at a positive value
        if curstem:
            if lastpositive >= minlen and bestscore >= minscore:
                stems.append([curstem[:lastpositive], lastpositive, bestscore])

    return stems


def ScoreStems(seq, stems, rstems, minscore,
               bracketweight, distcoef,
               orderpenalty, fiveprime):
    """Adjust the scores based on the stem distance, distance from 5'-end, and pseudoknot level"""

    # values == indices of the bp partners or -1 for unpaired bases
    bppartners = [-1 for _ in range(len(seq))]

    # set of restraint-bps
    rbps = set()

    for stem in rstems:
        for bp in stem[0]:
            rbps.add(bp)
            bppartners[bp[0]] = bp[1]
            bppartners[bp[1]] = bp[0]

    # keys = bps, values = their pseudoknot level (starting with 1)
    bplevels = PairsToDBN(rbps, returnlevels = True)

    # calculating adjusted scores for each stem
    for stem in stems:

        bps = stem[0]
        descr = "len={},bps={}".format(stem[1], stem[2]) # for debugging only

        levelset = set() # with base pairs of how many pseudoknot levels the stem is in conflict
        dots = 0 # number of confined unpaired bases
        brackets = 0 # number of bases belonging to the pseudoknotted wings (see URS defs)

        # confined region
        stemstart, stemend = bps[-1]

        # if we are inside the sub-ECR - where it will end
        inblockend = -1

        for pos in range(stemstart+1, stemend):

            partner = bppartners[pos]

            # if dot and not inside sub-ECR - increment dots
            if partner == -1:

                if pos > inblockend:
                    dots += 1

            # if the pseudoknotted wing and not inside sub-ECR - increment brackets
            elif partner < stemstart or partner > stemend:

                if pos > inblockend:
                    brackets += 1
                    levelset.add(bplevels[(min(pos, partner),
                                           max(pos, partner))])

            # if sub-ECR face - increase inblockend
            elif pos < partner:
                inblockend = partner
        

        stemdist = dots + bracketweight*brackets # stem distance
        stemdistfactor = (1/(1 + abs(stemdist - 4)))**distcoef

        order = len(levelset) # number of conflicting pseudoknot levels
        orderfactor = (1/(1 + order))**orderpenalty

        fiveprimedist   = (bps[0][0] - 0)/len(seq) #how far from 5'-end 
        fiveprimefactor = (1 - fiveprimedist)**fiveprime

        initscore  = stem[2] # initial bp score
        finalscore = initscore * stemdistfactor * orderfactor * fiveprimefactor
        
        descr += ",dt={},br={},sd={},sdf={}".format(dots, brackets,
                                                    round(stemdist,2),
                                                    round(stemdistfactor,2))
        descr += ",or={},orf={}".format(order,
                                        round(orderfactor,2))
        descr += ",fpd={},fpf={}".format(round(fiveprimedist,2),
                                         round(fiveprimefactor,2))

        stem.append(finalscore)
        stem.append(descr)

    # remove all stems with finscore < minfinscore
    return [stem for stem in stems if stem[3] >= minscore]


def ChooseStems(allstems, subopt = 1.0):
    """Choose the optimal stems from the stem set"""

    # sort in decreasing order of the adjusted score
    sortedstems = sorted(allstems, key = lambda x: x[3], reverse = True)

    resultstems = sortedstems[:1] # start with the top ranked stem

    # if no stems
    if not resultstems:
        return []

    # the adjusted score of the top ranked stem
    bestscore = resultstems[0][3]

    # add all stems with the subopt range that are
    # in conflict with all the better stems
    for stem in sortedstems[1:]:

        bps = stem[0]
        finscore = stem[3]

        # in conflict == the intersection of their sets
        # of the paired bases is not empty
        if all({p for bp1 in bps
                for p in bp1} & {p for bp2 in betterstem[0]
                                 for p in bp2}
               for betterstem in resultstems) and finscore >= subopt * bestscore:
            resultstems.append(stem)

    return resultstems

 
def OptimalStems(seq, bpmatrix, rbps = set(), rxs = set(), rstems = [],
                 subopt = 1.0, minlen = 2,
                 minbpscore = 6, minfinscore = 0,
                 bracketweight = 1.0, distcoef = 0.1,
                 orderpenalty = 0.0, fiveprime = 0.0):
    """Return the top stems"""

    # Remove already predicted bps from bp-restraints
    rbps = set(rbps) - {bp for stem in rstems for bp in stem[0]}

    allstems = AnnotateStems(bpmatrix, rbps, rxs, rstems, minlen, minbpscore)

    allstems = ScoreStems(seq, allstems, rstems, minfinscore,
                          bracketweight, distcoef,
                          orderpenalty, fiveprime)
    """
    # TEMPORARY PRINTING
    print('##################################################')
    for stem in sorted(allstems, key = lambda x: x[3], reverse = True):
        print(seq)
        print(StemsToDBN(rstems, seq))
        print(StemsToDBN([stem,],seq))
        print(stem[1:])
    """
    
    resultstems = ChooseStems(allstems, subopt)

    """
    # TEMPORARY PRINTING
    print('##################################################')
    for stem in resultstems:
        print(seq)
        print(StemsToDBN(rstems, seq))
        print(StemsToDBN([stem,],seq))
        print(stem[1:])
    """

    return resultstems


def ConsensusStemSet(stemsets):
    """Returns the set of bps present in all the stemsets"""
    
    if not stemsets:
        return set()

    # start with bps of the first stemset
    bps = {bp for stem in stemsets[0] for bp in stem[0]}

    # and intersect with all the other stemsets
    for stemset in stemsets[1:]:
        bps &= {bp for stem in stemset for bp in stem[0]}

    return bps

    
def SQRNdbnseq(seq, bpweights, restraints = None, dbn = None,
               subopt = 1.0, minlen = 2, minbpscore = 6,
               minfinscorefactor = 0.0, bracketweight = 1.0,
               distcoef = 0.1, orderpenalty = 0.0,
               fiveprime = 0.0, maxstemnum = 10**6):
    """seq == sequence (possibly with gaps or any other non-ACGU symbols
    bpweights == dictionary with canonical bps as keys and their weights as values
    restraints == string in dbn format; x symbols are forced to be unpaired
    dbn == known secondary structure in dbn format
    subopt == what share of the top stem score will be searched for alternative stems
    minlen == minimum stem length to predict
    minbpscore == minimum stem bp-score to predict
    minfinscorefactor == the percentage of minbpscore that will serve as the adjusted score threshold
    bracketweight == mult factor for pseudoknotted brackets in the stem distance calculation
    distcoef == how much the stem distance affects the stem score (0 == no effect)
    orderpenalty == how much we prioritize lower pseudoknot levels
    fiveprime == how much we prioritize 5'-close stems
    maxstemnum == how many stems we allow in a single predicted structure
    
    SQRNdbnseq returns a list of alternative predicted secondary structures in dbn format"""

    # turn seq into UPPERCASE & replace T with U
    seq = seq.upper().replace("T", "U")

    minfinscore = minbpscore * minfinscorefactor

    # if not restraints - generate a string of dots
    if not restraints:
        restraints = '.'*len(seq)

    assert len(seq) == len(restraints)
    
    # Unalign seq, dbn, and restraints strings
    shortseq, shortrest = UnAlign2(seq, restraints)

    if dbn:
        assert len(seq) == len(dbn)
        shortseq, shortdbn  = UnAlign2(seq, dbn)

    # Generate initial bp-matrix (with bp weights in cells)
    bpmatrix = BPMatrix(shortseq, bpweights)
    
    # Parse restraints into unpaired bases (rxs) and base pairs (rbps)
    rbps, rxs = ParseRestraints(shortrest.lower())

    # List of lists of stems (each stem list is a currently predicted secondary structure
    curstemsets = [[],]

    # List of finalized stem lists
    finstemsets = []

    # while list of intermediate stem lists is not empty
    while curstemsets:

        # new iteration
        newcurstemsets = []

        for stems in curstemsets:

            if len(stems) == maxstemnum:
                finstemsets.append(stems)
                continue

            # new optimal stems based on the current stem list    
            newstems = OptimalStems(shortseq, bpmatrix.copy(),
                                    rbps.copy(), rxs.copy(), stems,
                                    subopt, minlen, minbpscore,
                                    minfinscore, bracketweight,
                                    distcoef, orderpenalty,
                                    fiveprime)
            # append new intermediate stem lists
            if newstems:
                for newstem in newstems:
                    newcurstemsets.append(stems + [newstem,])
            # if no newstems returned - the stem list is considered final
            else:
                finstemsets.append(stems)

        # update the current stem lists
        curstemsets = newcurstemsets

    # list of dbn strings
    dbns = []

    # sort the final stem lists in descresing order of their total bp-score
    # and convert all final stem lists into dbn strings
    # and not forget about non-predicted bps from restraints
    finstemsets = sorted(finstemsets, key = lambda x: sum(y[2] for y in x), reverse = True)
    for stems in finstemsets:
        dbns.append(PairsToDBN({bp for stem in stems for bp in stem[0]} | set(rbps),
                    len(shortseq)))

    consbps = ConsensusStemSet(finstemsets)

    # ReAlign the dbn strings accoring to seq
    dbns = [ReAlign(x, seq) for x in dbns]
    cons = ReAlign(PairsToDBN(consbps, len(shortseq)), seq)

    # if input dbn is known - calculate the quality metrics
    if dbn:
        knownbps = set(DBNToPairs(shortdbn))

        constp = len(consbps & knownbps)
        consfp = len(consbps - knownbps)
        consfn = len(knownbps - consbps)

        consprc = round(constp / (constp + consfp), 3)
        consrcl = round(constp / (constp + consfn), 3)
        consfsc = round(2*constp / (2*constp + consfp + consfn), 3)

        consresult = [constp, consfp, consfn, consfsc, consprc, consrcl]

        bestfsc = 0.0
        result = []

        for rank, stemset in enumerate(finstemsets):

            setbps = {bp for stem in stemset for bp in stem[0]}

            tp = len(setbps & knownbps)
            fp = len(setbps - knownbps)
            fn = len(knownbps - setbps)

            prc = round(tp / (tp + fp), 3)
            rcl = round(tp / (tp + fn), 3)
            fsc = round(2*tp / (2*tp + fp + fn), 3)

            if fsc > bestfsc:
                bestfsc = fsc 
                result  = [tp, fp, fn, fsc, prc, rcl, rank + 1]

        return cons, dbns, consresult, result
    return cons, dbns, [], []


if __name__ == "__main__":

    from collections import Counter

    rst = None

    #SAM riboswitch
    seq = "GUUCUUAUCAAGAGAAGCAGAGGGACUGGCCCGACGAAGCUUCAGCAACCGGUGUAAUGGCGAAAGCCAUGACCAAGGUGCUAAAUCCAGCAAGCUCGAACAGCUUGGAAGAUAAGAACA"
    dbn = "(((((((((....(((((...(((.[[[[)))......)))))(((..(((((...(((((....))))).)))..)).)))...(]]]](((((.......)))))..))))))))))." 
    #TRNA
    #seq = "GGGCGGCUAGCUCAGCGGAAGAGCGCUCGCCUCACACGCGAGAGGUCGUAGGUUCAAGUCCUACGCCGCCCACCA"
    #dbn = "(((((((..((((....[..)))).(((((.......))))).....(((((..]....))))))))))))...."
    #Twister ribozyme
    #seq = "GAAAUAAUGCUACCAGACUGCACAAAAAGUCUGCCCGUUCCAAGUCGGGAUGAAAGCAGAGGAUUUC"
    #dbn = "((((...(((({{((((((........))))))(((..[[[..}}.))).....))))..]]]))))"
    #rst = "(x((xxx....xx..............,,,,,,..............................)).)"
    #rst = "..................................................................." 
    
    queue = []

    with open("CoRToise.fas") as file:
        lines = file.readlines()

        for ii in range(0,len(lines)-2,3):

            nm = lines[ii].strip()[1:]
            sq = lines[ii+1].strip()
            db = lines[ii+2].strip()

            queue.append([nm, sq, db, rst])

    #queue  = [["default", seq, dbn, rst],]
    #queue += [["default", seq, dbn, rst],]

    bpweights = {
                 'GU' : -1,
                 'AU' :  2,
                 'GC' :  4,
                 }

    subopt = 0.95
    minlen = 2
    minbpscore = 6
    minfinscorefactor = 0.5
    bracketweight = 1.1
    distcoef = 0.05
    orderpenalty = 0.05
    fiveprime = 0.01
    maxstemnum = 1

    resultsB = []
    resultsC = []

    for obj in queue:

        name, seq, dbn, rst = obj

        result = SQRNdbnseq(seq, bpweights, rst, dbn,
                            subopt, minlen, minbpscore,
                            minfinscorefactor, bracketweight,
                            distcoef, orderpenalty, fiveprime,
                            maxstemnum)

        print(name)
        print(seq)
        print(dbn)
        print('_'*len(seq))
        print(result[0],result[2])
        print('_'*len(seq))
        for rank, pred in enumerate(result[1]):
            if rank == result[3][-1]-1:
                print(pred, result[3])
            else:
                print(pred)
        print("#"*len(seq))

        resultsC.append(result[2])
        resultsB.append(result[3])

    print(resultsC)
    print(resultsB)

    tpC = sum(x[0] for x in resultsC)
    fpC = sum(x[1] for x in resultsC)
    fnC = sum(x[2] for x in resultsC)
    fsC = [x[3] for x in resultsC]
    prC = [x[4] for x in resultsC]
    rcC = [x[5] for x in resultsC]

    print(round(2*tpC / (2*tpC + fpC + fnC), 3), round(np.mean(fsC), 3),
          round(np.mean(prC), 3), round(np.mean(rcC), 3))

    tpB = sum(x[0] for x in resultsB)
    fpB = sum(x[1] for x in resultsB)
    fnB = sum(x[2] for x in resultsB)
    fsB = [x[3] for x in resultsB]
    prB = [x[4] for x in resultsB]
    rcB = [x[5] for x in resultsB]
    rkB = [x[6] for x in resultsB]
    
    print(round(2*tpB / (2*tpB + fpB + fnB), 3), round(np.mean(fsB), 3),
          round(np.mean(prB), 3), round(np.mean(rcB), 3))
    print(Counter(rkB))







    
