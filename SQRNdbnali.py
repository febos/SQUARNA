
import numpy as np
from multiprocessing import Pool

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

    boolmat = np.zeros((len(seq), len(seq)), dtype = float)

    # Fill the upper triangle of the boolmatrix with ones for bps
    # Ignore the bps that would form a hairpin of len < 3
    for i in range(len(seq)-1):
        for j in range(i+4,len(seq)):
            boolmat[i,j] = int(seq[i]+seq[j] in bps)

    # add 0 values for all possible base pairs
    # not yet in bps dictionary
    for b1 in bases:
        for b2 in bases:
            if b1+b2 not in bps:
                bps[b1+b2] = 0

    # Initialize the bp-score matrix
    scoremat = np.zeros((len(seq), len(seq)), dtype = float)

    # Fill the upper triangle of the matrix with bp weights
    # Ignore the bps that would form a hairpin of len < 3
    for i in range(len(seq)-1):
        for j in range(i+4,len(seq)):
            scoremat[i,j] = bps[seq[i]+seq[j]]

    return boolmat, scoremat


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
    rxs = [i for i in range(len(restraints)) if restraints[i]=='_'] # Unpaired bases
    return rbps, rxs


def StemsFromDiagTheFirstVersion(diag):
    """going from the outer pairs to the inner ones"""
    stems = []
    
    curstem  = []
    curscore = 0
    bestscore = 0

    for isbp, val, bp in diag:

        if isbp:
            if val >= 0: # if bp has non-negative score - add it to the current stem                  
                curstem.append(bp)
                lastpositive = len(curstem)
                curscore = curscore + val
                bestscore = curscore

            elif val < 0: # if bp has negative score
                if curstem:
                    curscore = curscore + val
                    if curscore <= 0: # if it is better to stop the stem before the last negatives
                        stems.append([curstem[:lastpositive], lastpositive, bestscore])
                        curstem = []
                        curscore = 0
                    else: # otherwise - growing the stem
                        curstem.append(bp)
        elif curstem: # if no bp but we had the stem before
            stems.append([curstem[:lastpositive], lastpositive, bestscore])
            curstem = []
            curscore = 0

    # in case we stopped at a positive value
    if curstem:
        stems.append([curstem[:lastpositive], lastpositive, bestscore])

    return stems


def StemsFromDiagTheFirstVersionReversed(diag):
    """going from the inner pairs to the outer pairs"""
    stems = []
    
    curstem  = []
    curscore = 0
    bestscore = 0

    for isbp, val, bp in diag[::-1]:

        if isbp:
            if val >= 0: # if bp has non-negative score - add it to the current stem                  
                curstem.insert(0, bp)
                lastpositive = len(curstem)
                curscore = curscore + val
                bestscore = curscore

            elif val < 0: # if bp has negative score
                if curstem:
                    curscore = curscore + val
                    if curscore <= 0: # if it is better to stop the stem before the last negatives
                        stems.append([curstem[-lastpositive:], lastpositive, bestscore])
                        curstem = []
                        curscore = 0
                    else: # otherwise - growing the stem
                        curstem.insert(0, bp)
        elif curstem: # if no bp but we had the stem before
            stems.append([curstem[-lastpositive:], lastpositive, bestscore])
            curstem = []
            curscore = 0

    # in case we stopped at a positive value
    if curstem:
        stems.append([curstem[-lastpositive:], lastpositive, bestscore])

    return stems


def PreStemsFromDiag(diag):
    "Return longest sequences of consecutive bps"""
    prestems = []
    
    firstbp = -1

    for i, pos in enumerate(diag):

        isbp, val, bp = pos

        # if first bp of the stem
        if isbp and firstbp < 0:
            firstbp = i

        # if first non-bp after the stem
        if not isbp and firstbp >= 0:
            prestems.append(diag[firstbp:i])
            firstbp = -1

    #if the diagonal ends with a stem
    if firstbp >= 0:
        prestems.append(diag[firstbp:])

    return prestems


def StemsFromPreStemsMaxLen(prestems):
    """Define stems as the longest sequence of consecutive bps"""
    stems = []
    for prestem in prestems:
        bps = [x[2] for x in prestem]
        score = sum(x[1] for x in prestem)
        stems.append([bps, len(bps), score])
    return stems


def StemsFromPreStemsTopScore(prestems):
    """Define stems as the top-scored sub-sequence of consecutive bps"""

    def MaxSubArrays(prestem):
        """split into max-subarrays"""
        if not prestem:
            return []

        cursum  = -10**6
        bestsum = -10**6

        begin    = -1
        curbegin = -1 
        end      = -1
        
        for i, pos in enumerate(prestem):
            val = pos[1]
            if val >= cursum + val:
                curbegin = i
            cursum = max(val, cursum + val)
            if cursum > bestsum:
                begin = curbegin
                end = i
                bestsum = cursum

        return MaxSubArrays(prestem[:begin]) + [prestem[begin : end + 1],] + MaxSubArrays(prestem[end + 1:])

    stems = []
    for prestem in prestems:
        for substem in MaxSubArrays(prestem):
            bps = [x[2] for x in substem]
            score = sum(x[1] for x in substem)
            stems.append([bps, len(bps), score])
    return stems


def StemsFromPreStemsAll(prestems):
    """Define stems as any sub-sequence of consecutive bps"""
    stems = []
    for prestem in prestems:
        for i in range(len(prestem)):
            for j in range(i + 1, len(prestem) + 1):
                substem = prestem[i:j]
                bps = [x[2] for x in substem]
                score = sum(x[1] for x in substem)
                stems.append([bps, len(bps), score])
    return stems


def StemsFromDiag(diag, mode):
    """Annotate stems at the given diagonal"""
    if mode == "ver1": # My first bugged version of the max-subarray stems
        return StemsFromDiagTheFirstVersion(diag)
    if mode == "ver1rev": # The reversed first bugged version of the max-subarray stems
        return StemsFromDiagTheFirstVersionReversed(diag)
    if mode in {"all", "maxlen", "topscore"}:
        prestems = PreStemsFromDiag(diag)
        if mode == "maxlen": # Stems as longest sequences of consecutive bps
            return StemsFromPreStemsMaxLen(prestems)
        if mode == "topscore":
            return StemsFromPreStemsTopScore(prestems)
        if mode == "all":
            return StemsFromPreStemsAll(prestems)


    assert True == False, "check your mode value" ### SHOULD NEVER HAPPEN


def AnnotateStems(bpboolmatrix, bpscorematrix, rbps, rxs,
                  rstems, minlen, minscore, mode):

    # copy the bpboolmatrix
    matrix = bpboolmatrix.copy()
    # obtain the matrix length
    N = bpboolmatrix.shape[0]
    # initialize result
    stems = []

    # apply bp-restraints - zero out the entire row and column except the cell for each bp
    for v, w in rbps:
        matrix[v, :] *= 0
        matrix[:, v] *= 0
        matrix[w, :] *= 0
        matrix[:, w] *= 0
        matrix[v, w]  = bpboolmatrix[v, w]

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

        diag = []

        i, j = x, y
        
        while i <= j - 4: # this is to forbid hairpins of len < 3
            diag.append([matrix[i, j], bpscorematrix[i, j], (i, j)])
            i += 1
            j -= 1

        for stem in StemsFromDiag(diag, mode):
            if stem[1] >= minlen and stem[2] >= minscore:
                stems.append(stem)

    return stems


def ScoreStems(seq, stems, rstems, minscore,
               bracketweight, distcoef,
               orderpenalty, fiveprimeleft,
               fiveprimeright,
               idealdist1, idealdist2,
               loopbonus):
    """Adjust the scores based on the stem distance, distance from 5'-end, and pseudoknot level"""

    # short near-symmetric internal loops
    goodloops = {(0, 0), (0, 1), (1, 0),
                 (1, 1), (1, 2), (2, 1),
                 (2, 2), (2, 3), (3, 2),
                 (3, 3), (3, 4), (4, 3),
                 (4, 4), (4, 2), (2, 4),
                 (3, 1), (1, 3),
                 (0, 2), (2, 0),
                 (5, 5)}

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

        descr = ''
        bps = stem[0]

        descr = "len={},bps={}".format(stem[1], stem[2]) # for debugging only

        levelset = set() # with base pairs of how many pseudoknot levels the stem is in conflict
        dots = 0 # number of confined unpaired bases
        brackets = 0 # number of bases belonging to the pseudoknotted wings (see URS defs)

        # confined region
        stemstart, stemend = bps[-1]

        block_edges = []

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
            elif pos < partner and partner > inblockend:
                inblockend = partner
                block_edges.append((pos, partner))

        # to prioritize short near-symmetric internal loops
        goodloop = False
        if len(block_edges) == 1:
            if (block_edges[0][0] - stemstart - 1,
                stemend - block_edges[0][1] - 1) in goodloops:
                goodloop = True
        loopfactor = 1 + loopbonus*goodloop

        # idealdist1 for hairpins, idealdist2 for everything else
        idealdist = idealdist1 if inblockend == -1 else idealdist2

        stemdist = dots + bracketweight*brackets # stem distance
        stemdistfactor = (1 / (1 + abs(stemdist - idealdist)))**distcoef

        order = len(levelset) # number of conflicting pseudoknot levels
        orderfactor = (1 / (1 + order))**orderpenalty

        fiveprimeleftdist  = (bps[0][0]   - 0)/len(seq) #how far the left wing is from the 5'-end
        fiveprimerightdist = (bps[-1][-1] - 0)/len(seq) #how far the right wing is from the 5'-end

        fiveprimeleftfactor  = (1 / (1 + fiveprimeleftdist))**fiveprimeleft
        fiveprimerightfactor = (1 / (1 + fiveprimerightdist))**fiveprimeright

        initscore  = stem[2] # initial bp score
        finalscore = initscore * stemdistfactor * orderfactor * fiveprimeleftfactor * fiveprimerightfactor * loopfactor
        
        descr += ",dt={},br={},sd={},sdf={}".format(dots, brackets,
                                                    round(stemdist,2),
                                                    round(stemdistfactor,2))
        descr += ",or={},orf={}".format(order,
                                        round(orderfactor,2))
        descr += ",fpld={},fplf={}".format(round(fiveprimeleftdist,2),
                                           round(fiveprimeleftfactor,2))
        descr += ",fprd={},fprf={},lf={}".format(round(fiveprimerightdist,2),
                                                 round(fiveprimerightfactor,2),
                                                 round(loopfactor,2))

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

    suboptrange = subopt * bestscore

    # add all stems within the subopt range that are
    # in conflict with all the better stems
    for stem in sortedstems[1:]:

        bps = stem[0]
        finscore = stem[3]

        if finscore < suboptrange:
            return resultstems

        # in conflict == the intersection of their sets
        # of the paired bases is not empty
        if all({p for bp1 in bps
                for p in bp1} & {p for bp2 in betterstem[0]
                                 for p in bp2}
               for betterstem in resultstems):
            resultstems.append(stem)

    return resultstems

 
def OptimalStems(seq, rstems, bpboolmatrix, bpscorematrix,
                 rbps = set(), rxs = set(),
                 subopt = 1.0, minlen = 2,
                 minbpscore = 6, minfinscore = 0,
                 bracketweight = 1.0, distcoef = 0.1,
                 orderpenalty = 0.0, fiveprimeleft = 0.0,
                 fiveprimeright = 0.0, mode = "ver1",
                 idealdist1 = 4, idealdist2 = 4,
                 loopbonus = 0.0,):
    """Return the top stems"""

    # Remove already predicted bps from bp-restraints
    rbps = set(rbps) - {bp for stem in rstems for bp in stem[0]}

    allstems = AnnotateStems(bpboolmatrix, bpscorematrix, rbps, rxs,
                             rstems, minlen, minbpscore, mode)

    allstems = ScoreStems(seq, allstems, rstems, minfinscore,
                          bracketweight, distcoef,
                          orderpenalty, fiveprimeleft,
                          fiveprimeright,
                          idealdist1, idealdist2,
                          loopbonus,
                          )
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


def mpOptimalStems(args):
    """OptimalStems version with a single input parameter for mp.pool.imap"""
    rstems   = args[1] # current stem set
    newstems = OptimalStems(*args)
    return newstems, rstems


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


def ScoreStruct(seq, stemset):
    """ Return the overall structure score based on the stemset"""
    bpscores = {"GU": -0.25,
                "AU": 1.5,
                "GC": 4.0,
               }
    for bp in sorted(bpscores.keys()): # CG for GC, etc.
        bpscores[bp[::-1]] = bpscores[bp]

    power = 1.75 # we will sum up these powers of the stem scores
    thescore = 0

    for stem in stemset:
        bpsum = sum(bpscores[seq[bp[0]]+seq[bp[1]]] for bp in stem[0])
        if bpsum >= 0:
            thescore += bpsum**power if bpsum >= 0 else -(abs(bpsum)**power)
    return round(thescore, 3)

   
def SQRNdbnali(seq, restraints = None, dbn = None,
               bpweights = {}):

    # turn seq into UPPERCASE & replace T with U
    seq = seq.upper().replace("T", "U")

    # if not restraints - generate a string of dots
    if not restraints:
        restraints = '.'*len(seq)

    assert len(seq) == len(restraints)
    
    # Unalign seq, dbn, and restraints strings
    shortseq, shortrest = UnAlign2(seq, restraints)

    if dbn:
        assert len(seq) == len(dbn)
        shortseq, shortdbn  = UnAlign2(seq, dbn)
    
    # Parse restraints into unpaired bases (rxs) and base pairs (rbps)
    rbps, rxs = ParseRestraints(shortrest)

    bpboolmatrix, bpscorematrix = BPMatrix(shortseq, bpweights)

    stems = AnnotateStems(bpboolmatrix, bpscorematrix, rbps, rxs,
                          rstems=[], minlen=2, minscore=0, mode="maxlen")

    
    radict = ReAlignDict(shortseq, seq)

    for stem in stems:

        yield [[(radict[v],radict[w]) for v,w in stem[0]], stem[-1]]


if __name__ == "__main__":

    from collections import Counter

    cons = "ACYURRACUAAGYCCARRARYAURRURYYN-----YARARYUR-AYRRGYCAURURCUUA-NCYAARAYRA-GAUYYYNYYYNYYYUNNN-NNNNNNNNNNNNNRNNNNYNNNYCUCAURACYYNUGYUURGARANA--UCUGGRYU-NRRARUYNGRRR-NNNNNNAYYYUGCCARYYRNNRYUAR"
    rst  = ".............[[[......................(......................................(....((......................(((...)))..................))......]]]......)............)......................."
    rst  = ".............[[[[............(.....(((((...................................(((((..((...................((((((...))))))<<<.....{.{{{{.)).....]]]]....)))))....>>>..)))))).....}}}}.}........" #rst2


    #rst = ".............[[[[............(.....(((((......................[[[[..{{{....(((((..((...................((((((...)))))){{{.......]]]].)).....]]]]....)))))....}}}..))))))..}}}.............." #rst1
    

    threads = 1

    queue = []

    with open("RF04188_selected_SS.fas") as file:
        lines = file.readlines()

        for ii in range(0,len(lines)-1,2):

            nm = lines[ii].strip()[1:]
            sq = lines[ii+1].strip()
            db = None
            queue.append([nm, sq, db, rst])

    queue = queue[:-6]

    bpweights = {'GU' : -1,
                 'AU' :  2,
                 'GC' :  4,}

    stemmatrix = np.zeros((len(queue[0][1]),len(queue[0][1])),dtype=int)



    for obj in queue:

        name, seq, dbn, rst = obj

        for stem in SQRNdbnali(seq, rst, dbn,
                               bpweights):
            for v,w in stem[0]:
                stemmatrix[v, w] += stem[-1]
                stemmatrix[w, v] += stem[-1]

    PrintMatrix(cons, stemmatrix)





    
