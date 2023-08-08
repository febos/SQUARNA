
import numpy as np
from multiprocessing import Pool
import sys

# Gapped values
GAPS = {'-', '.', '~'}
# Separator values (~chain breaks)
SEPS = {';', '&'}

# Reactivities encoding
ReactDict = {"_" : 0.00,  "+" : 0.50,  "#" : 1.00,
             "0" : 0.05,  "1" : 0.15,  "2" : 0.25,
             "3" : 0.35,  "4" : 0.45,  "5" : 0.55,
             "6" : 0.65,  "7" : 0.75,  "8" : 0.85,
             "9" : 0.95,  "a" : 0.00,  "b" : 0.04,
             "c" : 0.08,  "d" : 0.12,  "e" : 0.16,
             "f" : 0.20,  "g" : 0.24,  "h" : 0.28,
             "i" : 0.32,  "j" : 0.36,  "k" : 0.40,
             "l" : 0.44,  "m" : 0.48,  "n" : 0.52,
             "o" : 0.56,  "p" : 0.60,  "q" : 0.64,
             "r" : 0.68,  "s" : 0.72,  "t" : 0.76,
             "u" : 0.80,  "v" : 0.84,  "w" : 0.88,
             "x" : 0.92,  "y" : 0.96,  "z" : 1.00,
             }


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
        

def EncodedReactivities(seq, reacts, reactformat):
    """Takes a list of floats, returns a string"""

    # Put the reactivities into (0, 1) range
    reacts = [x if 0 <= x <= 1 else 0 if x < 0 else 1 for x in reacts]
    
    if reactformat == 3:
        reactline = ''.join(["_+##"[int(x * 3)]
                             for x in reacts])
    elif reactformat == 10:
        reactline = ''.join(['01234567899'[int(x * 10)]
                             for x in reacts])
    else:
        reactline = ''.join(['abcdefghijklmnopqrstuvwxyz'[int(x * 25) + 0.5]
                             for x in reacts])

    # Introduce the chain separators
    reactline = ''.join([reactline[i] if seq[i] not in SEPS else seq[i]
                         for i in range(len(seq))])
    return reactline


def PairsToDBN(newpairs, length = 0, returnlevels = False, levellimit = -1):
    """Convert a list of base pairs into a dbn string of the given length"""

    # Initialize the dbn string
    dbn = ['.']*length

    # Define "brackets" for 30 pseudoknot levels (and 19 more encoded with cyrillic letters)
    # Higher levels will be simply ignored
    levels = ['()','[]','{}','<>','Aa','Bb','Cc','Dd','Ee','Ff','Gg',
              'Hh','Ii','Jj','Kk','Ll','Mm','Nn','Oo','Pp','Qq','Rr',
              'Ss','Tt','Uu','Vv','Ww','Xx','Yy','Zz',
              'Бб','Гг','Дд','Ёё','Жж','Йй','Лл','Пп',
              'Фф','Цц','Чч','Шш','Щщ','Ьь','Ыы','Ъъ','Ээ','Юю','Яя']

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
            if level == len(levels):
                levels.append('..')

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

    # remove all levels higher than levellimit (if specified)
    if levellimit >= 0:
        groups = groups[:levellimit]

    # add all the pairs to the dbn string
    # according to their levels  
    for i, group in enumerate(groups):
        for pair in group:
            dbn[pair[0]] = levels[i][0]
            dbn[pair[1]] = levels[i][1]
            
    return ''.join(dbn)


def StemsToDBN(stems, seq):
    """Convert a list of stems (having lists
    of bps as the first elements) into a dbn string"""
    return PairsToDBN([bp for stem in stems for bp in stem[0]], len(seq))


def DBNToPairs(dbn):
    """Convert the dbn string into a sorted list of base pairs"""
    pairs = set()

    # keys == closing brackets, values == matching opening brackets
    closing = {'>':'<',']':'[',')':'(','}':'{','a':'A','b':'B','c':'C','d':'D',
               'e':'E','f':'F','g':'G','h':'H','i':'I','j':'J','k':'K','l':'L',
               'm':'M','n':'N','o':'O','p':'P','q':'Q','r':'R','s':'S','t':'T',
               'u':'U','v':'V','w':'W','x':'X','y':'Y','z':'Z',
               'б':'Б','г':'Г','д':'Д','ё':'Ё','ж':'Ж','й':'Й','л':'Л','п':'П',
               'ф':'Ф','ц':'Ц','ч':'Ч','ш':'Ш','щ':'Щ','ь':'Ь','ы':'Ы','ъ':'Ъ',
               'э':'Э','ю':'Ю','я':'Я',}
    # 30+19 bp stacks for 30+19 allowed pseudoknot levels
    stack = {'<':[],'(':[],'{':[],'[':[],'A':[],'B':[],'C':[],'D':[],'E':[],
             'F':[],'G':[],'H':[],'I':[],'J':[],'K':[],'L':[],'M':[],'N':[],
             'O':[],'P':[],'Q':[],'R':[],'S':[],'T':[],'U':[],'V':[],'W':[],
             'X':[],'Y':[],'Z':[],
             'Б':[],'Г':[],'Д':[],'Ё':[],'Ж':[],'Й':[],'Л':[],'П':[],'Ф':[],
             'Ц':[],'Ч':[],'Ш':[],'Щ':[],'Ь':[],'Ы':[],'Ъ':[],'Э':[],'Ю':[],
             'Я':[],}
              
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

    assert len(shortdbn) + sum(longseq.count(GAP) for GAP in GAPS) == len(longseq),\
    "Cannot ReAlign dbn string - wrong number of gaps:\n{}\n{}".format(longseq,shortdbn)

    # initialize the returning string
    newdbn = []

    for x in longseq:
        if x in GAPS: # if gap (hyphen/dot/tilda) -> add gap
            if seqmode:
                newdbn.append('-')
            else:
                newdbn.append('.')
        else: # otherwise -> add the next symbol
            newdbn.append(dbn.pop(0))

    # return the aligned string
    return ''.join(newdbn)


def UnAlign(seq, dbn):
    """Removes gaps from a pair of seq & dbn strings
    hyphens, dots, and tildas are treated as gaps"""

    # First, remove base pairs of gap columns
    # to avoid any shifts
    cleandbn = list(dbn)
    pairs = DBNToPairs(dbn)
    for v,w in pairs:
        if seq[v] in GAPS or seq[w] in GAPS:
            cleandbn[v] = '.'
            cleandbn[w] = '.'

    # Unalign the dbn based on the seq string
    newdbn = ''.join([cleandbn[i] for i in range(len(seq))
                      if seq[i] not in GAPS])
    # Then unalign the seq string itself
    newseq = ''.join([x for x in seq if x not in GAPS])
    # return the unaligned seq & dbn strings
    return newseq, newdbn


def BPMatrix(seq, weights, rxs, rlefts, rrights, interchainonly = False, reacts = None):
    """Return a matrix with bp-weights in cells"""

    # list of ordinal numbers of chains
    chains = [0 for _ in range(len(seq))]
    if interchainonly:
        curr   =  0
        for i in range(len(seq)):
            if seq[i] in SEPS:
                curr += 1
            else:
                chains[i] = curr

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
    # And satisfy the restraints
    for i in range(len(seq) - 1):

        # Allow short "hairpins" between different RNA chains
        inc4 = 4
        for chk in (1, 2):
            if i + chk < len(seq) and seq[i + chk] in SEPS:
                inc4 = chk + 1

        for j in range(i + inc4, len(seq)):
            boolmat[i, j] = int(seq[i] + seq[j] in bps) * \
                            (not interchainonly or chains[i] != chains[j]) * \
                            (i not in rxs and j not in rxs) * \
                            (j not in rlefts) * \
                            (i not in rrights)

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
    for i in range(len(seq) - 1):

        # Allow short "hairpins" between different RNA chains
        inc4 = 4
        for chk in (1, 2):
            if i + chk < len(seq) and seq[i + chk] in SEPS:
                inc4 = chk + 1
        
        for j in range(i + inc4, len(seq)):

            # Weight with reactivities if any
            reactfactor = 1 if reacts is None else (1 - (reacts[i] + reacts[j]) / 2 ) * 2
            if bps[seq[i] + seq[j]] <= 0:
                reactfactor = 1 / max(reactfactor, 0.01)

            scoremat[i, j] = bps[seq[i] + seq[j]] * boolmat[i, j] * reactfactor

    return boolmat, scoremat


def ParseRestraints(restraints):
    """ Convert a dbn string into base pairs and unpaired bases"""
    rbps = DBNToPairs(restraints) # Base pairs
    rxs     = {i for i in range(len(restraints)) if restraints[i]=='_'} # Unpaired bases
    rlefts  = {i for i in range(len(restraints)) if restraints[i]=='/'}  # 5'-ends of bps
    rrights = {i for i in range(len(restraints)) if restraints[i]=='\\'} # 3'-ends of bps
    return rbps, rxs, rlefts, rrights


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


def StemsFromPreStemsDiffEdge(prestems, diff = 1):
    """Define stems as maximum sub-sequence of consecutive bps with up to
    <diff> trimmed base pairs from either of the edges"""
    stems = []
    for prestem in prestems:
        for i in range(diff + 1):
            for j in range(len(prestem) - diff, len(prestem) + 1):
                if j > i:
                    substem = prestem[i:j]
                    bps = [x[2] for x in substem]
                    score = sum(x[1] for x in substem)
                    stems.append([bps, len(bps), score])
    return stems


def StemsFromDiag(diag, diff = 1):
    """Annotate stems at the given diagonal"""
    prestems = PreStemsFromDiag(diag)
    return StemsFromPreStemsDiffEdge(prestems, diff)


def AnnotateStems(bpboolmatrix, bpscorematrix, rbps,
                  rstems, minlen, minscore, diff = 1):

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

    # start cells of all the diagonals, except the ones producing hairpins of len < 3
    diagstarts  = [(0, x)   for x in range(4, N)]
    diagstarts += [(y, N-1) for y in range(1, N-4)]

    # for each diagonal select the stems
    # as maximal score sub-arrays
    for x, y in diagstarts:

        diag = []
        i, j = x, y
        
        while i <= j - 1: 
            diag.append([matrix[i, j], bpscorematrix[i, j], (i, j)])
            i += 1
            j -= 1

        for stem in StemsFromDiag(diag, diff):
            if stem[1] >= minlen and stem[2] >= minscore:
                stems.append(stem)
    
    return stems


def IsGNRA(seq):
    """Whether the input sequence is of GNRA pattern"""
    if len(seq) != 4:
        return False
    if seq[0] == 'G' and seq[2] in ('G','A') and seq[3] == 'A':
        return True
    return False


def ScoreStems(seq, stems, rstems,
               reacts, minscore,
               bracketweight, distcoef,
               orderpenalty, loopbonus):
    """Adjust the scores based on the stem distance, pseudoknot level, and
       possibly loop bonus"""

    # short near-symmetric internal loops
    goodloops = {(0, 0), (0, 1), (1, 0),
                 (1, 1), (0, 2), (2, 0),
                 (2, 2), (1, 2), (2, 1),
                 (3, 1), (1, 3),
                 (2, 3), (3, 2),
                 (3, 3), (3, 4), (4, 3), 
                 (4, 4), (4, 2), (2, 4),
                }

    # values == indices of the bp partners or -1 for unpaired bases
    bppartners = [-1 for _ in range(len(seq))]

    # set of restraint-bps
    rbps = set()

    # Fill the bppartners list
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

        reactfactor =  1 #(sum((1 - (reacts[bp[0]] + reacts[bp[1]]) / 2 ) for bp in bps) / len(bps) * 2) ** 1.5

        descr = "len={},bps={}".format(stem[1], stem[2]) # for debugging only

        levelset = set() # with base pairs of how many pseudoknot levels the stem is in conflict
        dots = 0 # number of confined unpaired bases
        brackets = 0 # number of bases belonging to the pseudoknotted wings (see URS defs)

        # confined region
        stemstart, stemend = bps[-1]

        block_edges = []

        # if we are inside the sub-ECR - where it will end
        inblockend = -1

        # Whether the stem confines an external loop
        between_chains = False

        for pos in range(stemstart+1, stemend):

            partner = bppartners[pos]

            # if dot and not inside sub-ECR - increment dots
            if partner == -1:

                if pos > inblockend:
                    dots += 1
                # If a separator instead of a dot
                if seq[pos] in SEPS:
                    between_chains = True

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

        # to prioritize short near-symmetric internal loops inside
        goodloop = False
        diff1 = 0
        if len(block_edges) == 1:
            if (block_edges[0][0] - stemstart - 1,
                stemend - block_edges[0][1] - 1) in goodloops:
                goodloop = True
                diff1 = abs((block_edges[0][0] - stemstart - 1) - (stemend - block_edges[0][1] - 1))
        # to prioritize short near-symmetric internal loops outside
        goodloopout = False
        diff2 = 0
        outerstart, outerend = bps[0]
        vv,ww = outerstart - 1, outerend + 1
        while vv >= 0 and outerstart - vv - 1 < 5 and bppartners[vv] == -1:
            vv -= 1
        while ww < len(seq) and ww - outerend - 1 < 5 and bppartners[ww] == -1:
            ww += 1
        if bppartners[vv] == ww and bppartners[ww] == vv and\
           (outerstart - vv - 1, ww - outerend - 1) in goodloops:
            goodloopout = True
            diff2 = abs((outerstart - vv - 1) - (ww - outerend - 1))

        # Double bonus for symmetric loops, 1.5-bonus for (x-1,x) loops
        # and a single bonus for (x-2,x) loops
        loopfactor = 1 + loopbonus*goodloop*(2-diff1/2) + loopbonus*goodloopout*(2-diff2/2)

        # Bonus for tetraloops
        tetrafactor = 1+0.25*IsGNRA(seq[stemstart+1:stemend])

        # 4 for hairpins, 2 for everything else
        idealdist = 4 if inblockend == -1 else 2

        stemdist = dots + bracketweight*brackets # stem distance

        # Ignore stemdistfactor for external loops
        stemdistfactor = (1 / (1 + abs(stemdist - idealdist)))**distcoef if not between_chains else 1

        order = len(levelset) # number of conflicting pseudoknot levels
        orderfactor = (1 / (1 + order))**orderpenalty

        initscore  = stem[2] # initial bp score
        finalscore = initscore * stemdistfactor * orderfactor * loopfactor * reactfactor * tetrafactor
        
        descr += ",dt={},br={},sd={},sdf={}".format(dots, brackets,
                                                    round(stemdist,2),
                                                    round(stemdistfactor,2))
        descr += ",or={},orf={}".format(order,
                                        round(orderfactor,2))
        descr += ",lf={},rf={}".format(round(loopfactor,2),
                                       round(reactfactor,2))

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
                 reacts, rbps = set(), 
                 subopt = 1.0, minlen = 2,
                 minbpscore = 6, minfinscore = 0,
                 bracketweight = 1.0, distcoef = 0.1,
                 orderpenalty = 0.0, loopbonus = 0.0,):
    """Return the top stems"""

    # Remove already predicted bps from bp-restraints
    rbps = set(rbps) - {bp for stem in rstems for bp in stem[0]}

    allstems = AnnotateStems(bpboolmatrix, bpscorematrix, rbps,
                             rstems, minlen, minbpscore)

    allstems = ScoreStems(seq, allstems, rstems, reacts,
                          minfinscore,
                          bracketweight, distcoef,
                          orderpenalty, loopbonus,
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


def ScoreStruct(seq, stemset, reacts):
    """ Return the overall structure score based on the stemset"""
    bpscores = {"GU": -1.0,
                "AU":  1.5,
                "GC":  4.0,
               }
    for bp in sorted(bpscores.keys()): # CG for GC, etc.
        bpscores[bp[::-1]] = bpscores[bp]

    power = 1.7 # we will sum up these powers of the stem scores
    thescore = 0
    
    paired = set()

    for stem in stemset:
        
        bpsum = 0
        for v,w in stem[0]:
            bpsum += bpscores[seq[v]+seq[w]]
            paired.add(v)
            paired.add(w)
            
        if bpsum > 0:
            thescore += bpsum**power 

    # seq length without the separator positions
    sepnum = sum(1 for _ in seq if _ in SEPS)

    # reactscore - sum of "errors" for each position
    # treating paired positions as predicted to be 0
    # and unpaired positions as predicted to be 1
    # then for each paired position the error is the
    # reactivity, and for unpaired = 1-reactivity
    reactscore = 1 - sum(reacts[i] if i in paired else 1 - reacts[i]
                         for i in range(len(seq))
                         if seq[i] not in SEPS) / (len(seq) - sepnum)

    # Return three scores - totalscore, structscore, reactscore
    return round(thescore * reactscore, 3), round(thescore, 3), round(reactscore, 3)


def RankStructs(stemsets, rankbydiff = False, rankby = (0, 2, 1)):
    """Rank the predicted structures"""

    # First sort according to the rankby order of scores
    finstemsets = sorted(stemsets,
                         key = lambda x: [x[1][rb] for rb in rankby],
                         reverse = True)

    # if not rankbydiff - that is all
    if not rankbydiff or len(finstemsets) < 3:
        return finstemsets

    # Otherwise - iteratively select the furthest from the
    # already selected structures breaking ties by the scores
    allbps = set()

    # Compile the entire set of observed bps among all
    # the structures
    for stemset in finstemsets:
        bps = {bp for stem in stemset[0] for bp in stem[0]}
        stemset.append(bps)
        allbps |= bps

    seenbps = {bp for bp in finstemsets[0][-1]}
    curind  = 1

    # While the pre-selected subset of structures
    # does not include all the observed bps
    # or until the last but one structure
    while seenbps != allbps and curind < len(finstemsets) - 1:

        finstemsets = finstemsets[:curind] +\
                      sorted(finstemsets[curind:],
                             key = lambda x: (len(x[-1] - seenbps),
                                              [x[1][rb] for rb in rankby]),
                             reverse = True)
        seenbps |= finstemsets[curind][-1]
        curind += 1

    # Sort the remaining tail according to the rankby
    # order of the scores
    finstemsets = finstemsets[:curind] +\
                  sorted(finstemsets[curind:],
                         key = lambda x: [x[1][rb] for rb in rankby],
                         reverse = True)
    # Remove the appended bps sets and return
    return [x[:-1] for x in finstemsets]  

   
def SQRNdbnseq(seq, reacts = None, restraints = None, dbn = None,
               paramsets = [], conslim = 1, toplim = 5,
               hardrest = False, rankbydiff = False,
               rankby = (0, 2, 1), interchainonly = False,
               threads = 1, mp = True, stemmatrix = None):
    """seq == sequence (possibly with gaps or any other non-ACGU symbols
    reacts == list of reactivities (values from 0 to 1)
    restraints == string in dbn format; x symbols are forced to be unpaired
    dbn == known secondary structure in dbn format
    conslim == how many alternative structures are used to derive the consensus
    toplim == how many top ranked structures are used to measure the performance
    hardrest == force restraint base pairs to be present in the predicted dbns or not
    rankbydiff == Rank the most different structures on top
    rankby == indices of the scores used for rankings
    interchainonly = allow only bps between different chains
    threads == number of CPUs to use
    mp == whether to use multiprocessing at all
    stemmatrix == if specified the bpscorematrix will be weighted with the stemmatrix
                  that was derived from the sequence alignment
    SQRNdbnseq returns a list of alternative predicted secondary structures in dbn format"""

    assert set(rankby) == {0, 1, 2} and len(rankby) == 3, "Invalid ranking indices"

    # turn seq into UPPERCASE & replace T with U
    seq = seq.upper().replace("T", "U")

    # if not restraints - generate a string of dots
    if not restraints:
        restraints = '.'*len(seq)

    assert len(seq) == len(restraints), "Invalid restraints given"

    # default reacts == list of 0.5 values
    if not reacts:
        reacts = [0.5 for i in range(len(seq))]

    assert len(reacts) == len(seq), "Invalid reactivities given"

    # if we need to decode the reactivities from a string
    if type(reacts) == str:
        reacts = [ReactDict[ch] for ch in reacts]
    
    # Unalign seq, dbn, reacts, and restraints strings
    shortseq, shortrest   = UnAlign(seq, restraints)
    shortreacts = [reacts[i] for i in range(len(seq)) if seq[i] not in GAPS]

    if dbn:
        assert len(seq) == len(dbn)
        shortseq, shortdbn  = UnAlign(seq, dbn)

    # Unalign stemmatrix
    if stemmatrix is not None:
        gapinds = [i for i in range(len(seq)) if seq[i] in GAPS]
        shortsmat = np.delete(stemmatrix, gapinds, 0)
        shortsmat = np.delete(shortsmat, gapinds, 1)
    
    # Parse restraints into unpaired bases (rxs) and base pairs (rbps)
    rbps, rxs, rlefts, rrights = ParseRestraints(shortrest)

    # final stem sets among all the parameter sets
    finfinstemsets = []
    seen_structures = set() # set of bp-sets (to avoid repeats in finfinstemsets)

    for psi, paramset in enumerate(paramsets):

        bpweights         = paramset["bpweights"]
        suboptmax         = paramset["suboptmax"]
        suboptmin         = paramset["suboptmin"]
        suboptsteps       = paramset["suboptsteps"]
        minlen            = paramset["minlen"]
        minbpscore        = paramset["minbpscore"]
        minfinscorefactor = paramset["minfinscorefactor"]
        bracketweight     = paramset["bracketweight"]
        distcoef          = paramset["distcoef"]
        orderpenalty      = paramset["orderpenalty"]
        maxstemnum        = paramset["maxstemnum"]
        loopbonus         = paramset["loopbonus"]

        # starting subopt value
        cursubopt = suboptmin
        # increment for the subopt value
        suboptinc = (suboptmax - suboptmin)/suboptsteps

        minfinscore = minbpscore * minfinscorefactor

        # Generate initial bp-matrix (with bp weights in cells)
        bpboolmatrix, bpscorematrix = BPMatrix(shortseq, bpweights, rxs,
                                               rlefts, rrights,
                                               interchainonly,
                                               reacts = shortreacts)

        #Weight the bpscorematrix with the alignment-derived scores
        if stemmatrix is not None:
            bpscorematrix = bpscorematrix * shortsmat

        # List of lists of stems (each stem list is a currently predicted secondary structure
        curstemsets = [[],]

        # List of finalized stem lists
        finstemsets = []

        # Starting with a single empty structure
        cursize = len(curstemsets)

        if mp: # multiprocessing version
            with Pool(threads) as pool:

                # while list of intermediate stem lists is not empty
                while curstemsets:

                    # Each time we diverge - cursubopt is increased by suboptinc
                    if len(curstemsets) > cursize and cursubopt < suboptmax:
                        cursize = len(curstemsets)
                        cursubopt += suboptinc

                    # filtering by len(stems) == maxstemnum
                    newcurstemsets = []
                    for stems in curstemsets:
                        if len(stems) == maxstemnum:
                            finstemsets.append(stems)
                        else:
                            newcurstemsets.append(stems)
                    curstemsets = newcurstemsets

                    # new iteration
                    newcurstemsets = []

                    inputs = ((shortseq, stems,
                               bpboolmatrix.copy(), bpscorematrix,
                               shortreacts, rbps.copy(), 
                               cursubopt, minlen, minbpscore, 
                               minfinscore, bracketweight,
                               distcoef, orderpenalty,
                               loopbonus,
                               ) for stems in curstemsets)

                    # new optimal stems based on the current stem list 
                    for newstems, stems in pool.imap(mpOptimalStems, inputs):

                        # append new intermediate stem lists
                        if newstems:
                            for newstem in newstems:
                                newcurstemsets.append(stems + [newstem,])
                        # if no newstems returned - the stem list is considered final
                        else:
                            finstemsets.append(stems)

                    # update the current stem lists
                    curstemsets = newcurstemsets

        else: # non-multiprocessing version
            # while list of intermediate stem lists is not empty
            while curstemsets:

                # Each time we diverge - cursubopt is increased by suboptinc
                if len(curstemsets) > cursize and cursubopt < suboptmax:
                    cursize = len(curstemsets)
                    cursubopt += suboptinc

                # filtering by len(stems) == maxstemnum
                newcurstemsets = []
                for stems in curstemsets:
                    if len(stems) == maxstemnum:
                        finstemsets.append(stems)
                    else:
                        newcurstemsets.append(stems)
                curstemsets = newcurstemsets

                # new iteration
                newcurstemsets = []

                for stems in curstemsets:

                    # new optimal stems based on the current stem list
                    newstems = OptimalStems(shortseq, stems, bpboolmatrix.copy(),
                                            bpscorematrix, shortreacts, rbps.copy(),
                                            cursubopt, minlen, minbpscore, 
                                            minfinscore, bracketweight,
                                            distcoef, orderpenalty,
                                            loopbonus)

                    # append new intermediate stem lists
                    if newstems:
                        for newstem in newstems:
                            newcurstemsets.append(stems + [newstem,])
                    # if no newstems returned - the stem list is considered final
                    else:
                        finstemsets.append(stems)

                # update the current stem lists
                curstemsets = newcurstemsets

        for finstemset in finstemsets:

            bpsset = tuple(sorted([bp for stem in finstemset for bp in stem[0]]))

            if bpsset not in seen_structures:
                # append [stemset, structscore, paramsetind]
                finfinstemsets.append([finstemset,
                                       ScoreStruct(shortseq, finstemset, shortreacts),
                                       psi])
                seen_structures.add(bpsset)

    # list of dbn strings
    dbns = []

    # sort the final stem lists in descreasing order of their total bp-score by default
    # and if rankbydiff - prioritize the most diverged structures
    finstemsets = RankStructs(finfinstemsets, rankbydiff, rankby)

    forcedbps = {(v,w) for v,w in rbps
                 if shortseq[v]+shortseq[w] in bpweights or
                    shortseq[w]+shortseq[v] in bpweights} if hardrest else set()

    # convert all final stem lists into dbn strings
    # and not forget about non-predicted bps from restraints
    for stems, structscore, paramsetind in finstemsets:
        dbns.append(PairsToDBN({bp for stem in stems for bp in stem[0]} | forcedbps,
                    len(shortseq)))

    consbps = ConsensusStemSet([xx[0] for xx in finstemsets[:conslim]]) | forcedbps # Consensus of the Top-ranked

    # ReAlign the dbn strings accoring to seq
    dbns = [ReAlign(x, seq) for x in dbns]
    cons = ReAlign(PairsToDBN(consbps, len(shortseq)), seq)

    # Introducing chain separators into the predicted structures
    dbns = [''.join([_[i] if seq[i] not in SEPS else seq[i]
                     for i in range(len(seq))]) for _ in dbns]
    cons = ''.join([cons[i] if seq[i] not in SEPS else seq[i]
                    for i in range(len(seq))])

    # if input dbn is known - calculate the quality metrics
    if dbn:
        knownbps = set(DBNToPairs(shortdbn))

        constp = len(consbps & knownbps)
        consfp = len(consbps - knownbps)
        consfn = len(knownbps - consbps)

        consprc = (round(constp / (constp + consfp), 3)) if (constp + consfp) else 1
        consrcl = (round(constp / (constp + consfn), 3)) if (constp + consfn) else 1
        consfsc = (round(2*constp / (2*constp + consfp + consfn), 3)) if (2*constp + consfp + consfn) else 1

        consresult = [constp, consfp, consfn, consfsc, consprc, consrcl]

        bestfsc = -1
        result = []

        for rank, stemset in enumerate(finstemsets):

            setbps = {bp for stem in stemset[0] for bp in stem[0]} | forcedbps

            tp = len(setbps & knownbps) # Correctly predicted bps
            fp = len(setbps - knownbps) # Wrongly predicted bps
            fn = len(knownbps - setbps) # Missed bps

            prc = (round(tp / (tp + fp), 3)) if (tp + fp) else 1
            rcl = (round(tp / (tp + fn), 3)) if (tp + fn) else 1
            fsc = (round(2*tp / (2*tp + fp + fn), 3)) if (2*tp + fp + fn) else 1

            if fsc > bestfsc:
                bestfsc = fsc 
                result  = [tp, fp, fn, fsc, prc, rcl, rank + 1]

            # TopN only (N == toplim)
            if rank + 1 >= toplim:
                break

        return cons, [(dbns[jj],*finstemsets[jj][1:]) for jj in range(len(dbns))], consresult, result
    return cons, [(dbns[jj],*finstemsets[jj][1:]) for jj in range(len(dbns))], [np.nan]*6, [np.nan]*7


def RunSQRNdbnseq(name, sequence, reactivities, restraints,
                  reference, paramsetnames,
                  paramsets, threads, rankbydiff, rankby,
                  hardrest, interchainonly, toplim, outplim,
                  conslim, reactformat, mp = True,
                  sink = sys.stdout, stemmatrix = None):
    """Main-like function;
       sink param is standard system output by default,
       we need it to use the buffer in alignment-mode parallelizations"""

    # Run prediction
    prediction = SQRNdbnseq(sequence, reactivities, restraints, reference,
                            paramsets, conslim, toplim, hardrest,
                            rankbydiff, rankby, interchainonly, threads, mp, stemmatrix)
    
    # Unpack the results
    consensus, predicted_structures, consensus_metrics, topN_metrics = prediction

    print(name, file = sink)
    print(sequence, file = sink)

    # Printing everything observed in the input
    if reactivities:
        print(EncodedReactivities(sequence,
                                  reactivities,
                                  reactformat),
              "reactivities", sep = '\t', file = sink)
    if restraints:
        print(''.join([restraints[i]
                       if sequence[i] not in SEPS
                       else sequence[i]
                       for i in range(len(sequence))]),
              "restraints", sep = '\t', file = sink)
    if reference:
        print(''.join([reference[i]
                       if sequence[i] not in SEPS
                       else sequence[i]
                       for i in range(len(sequence))]),
              "reference", sep = '\t', file = sink)

    # Separator line 1
    print('_'*len(sequence), file = sink)

    # Printing consensus
    # along with its metrics if reference is present
    if reference:
        print(consensus,
              "top-{}_consensus".format(conslim),
              "TP={},FP={},FN={},FS={},PR={},RC={}".format(*consensus_metrics),
              sep = '\t', file = sink)
    else:
        print(consensus,
              "top-{}_consensus".format(conslim), sep = '\t', file = sink)

    # Separator line 2
    print('='*len(sequence), file = sink)

    # Printing up to outplim predicted structures
    # along with their scores and
    # metrics of the best one if reference is present
    for i, pred in enumerate(predicted_structures[:outplim]):
        
        struct, scores, paramsetind = pred
        totalscore, structscore, reactscore = scores

        if reference and i + 1 == topN_metrics[-1]:
            print(struct, "#{}".format(i+1), totalscore,
                  structscore, reactscore,
                  paramsetnames[paramsetind],
                  "TP={},FP={},FN={},FS={},PR={},RC={},RK={}".format(*topN_metrics),
                  sep='\t', file = sink)
        else:
            print(struct, "#{}".format(i+1), totalscore,
                  structscore, reactscore,
                  paramsetnames[paramsetind], sep='\t', file = sink)

    return consensus, predicted_structures, consensus_metrics, topN_metrics



if __name__ == "__main__":

    from collections import Counter
  
    queue = []

    with open("datasets/S01_200.fas") as file:
        lines = file.readlines()

        for ii in range(0,len(lines)-2,4):

            nm  = lines[ii].strip()[1:]
            sq  = lines[ii+1].strip()
            db  = lines[ii+2].strip()
            rct = list(map(float,lines[ii+3].strip().split()))
            queue.append([nm, sq, db, None, rct])

    outp = open('temp.tsv','w')

    title = '\t'.join("suboptmin suboptmax suboptsteps minlen minbpscore minfinscorefactor distcoef orderpenalty def2 rankbydiff rankby tpc fpc fnc fstotc fsc prc rcc tp5 fp5 fn5 fstot5 fs5 pr5 rc5".split())
    print(title)
    outp.write(title+'\n')
    outp.close()

    ######################################

    threads = 32

    conslim = 1
    toplim  = 5
    maxstemnum = 10**6

    rbdic = {"rs":(0,2,1),"r":(2,0,1),"s":(1,2,0)}

    for suboptmin in (0.7, 0.75, 0.8, 0.85):
        for suboptmax in (0.9, 0.95, 0.99):
            for suboptsteps in (2, 3):
                for minlen in (2, ):
                    for minbpscore in (5.75, 7, 8.25, 9.5):
                        for minfinscorefactor in (1.0, 1.25, 1.5):
                            for distcoef in (0.09, 0.07, 0.11):
                                for orderpenalty in (1.0, 0.75, 0.5):
                                    for def2 in (True,):
                                        for rankbydiff in (False, True):
                                            for rankby in ("rs","s"):
                                                            
                                                print(suboptmin, suboptmax, suboptsteps, minlen, minbpscore, minfinscorefactor, distcoef, orderpenalty, def2,
                                                      rankbydiff, rankby, sep='\t', end='\t')

                                                paramsets = []
                                                paramsets.append({"bpweights" : {'GU' : -1.25,
                                                                                'AU' : 1.25,
                                                                                'GC' : 3.25,},
                                                                  "suboptmax" : suboptmax,
                                                                  "suboptmin" : suboptmin,
                                                                  "suboptsteps": suboptsteps,
                                                                  "minlen" : minlen,
                                                                  "minbpscore" : minbpscore,
                                                                  "minfinscorefactor" : minfinscorefactor,
                                                                  "distcoef" : distcoef,
                                                                  "bracketweight" :  -2,
                                                                  "orderpenalty"  : orderpenalty,
                                                                  "loopbonus": 0.125,
                                                                  "maxstemnum" : maxstemnum,
                                                                  })

                                                if def2:
                                                    paramsets.append({"bpweights" : {'GU' : 1,
                                                                                    'AU' :  1,
                                                                                    'GC' :  2,},
                                                                      "suboptmax" : suboptmax,
                                                                      "suboptmin" : suboptmin,
                                                                      "suboptsteps": suboptsteps,
                                                                      "minlen" : minlen,
                                                                      "minbpscore" : 3,
                                                                      "minfinscorefactor" : 0.99,
                                                                      "distcoef" : 0.1,
                                                                      "bracketweight" :  -2,
                                                                      "orderpenalty"  : 1.35,
                                                                      "loopbonus": 0.125,
                                                                      "maxstemnum" : maxstemnum,
                                                                      })

                                                resultsB = []
                                                resultsC = []

                                                for obj in queue:

                                                    name, seq, dbn, rst, react = obj

                                                    result = SQRNdbnseq(seq, react, rst, dbn,
                                                                        paramsets, conslim, toplim,
                                                                        rankbydiff = rankbydiff,
                                                                        rankby = rbdic[rankby],
                                                                        threads = threads)

                                                    resultsC.append(result[2])
                                                    resultsB.append(result[3])

                                                tpC = sum(x[0] for x in resultsC)
                                                fpC = sum(x[1] for x in resultsC)
                                                fnC = sum(x[2] for x in resultsC)
                                                fsC = [x[3] for x in resultsC]
                                                prC = [x[4] for x in resultsC]
                                                rcC = [x[5] for x in resultsC]

                                                print(tpC, fpC, fnC, round(2*tpC / (2*tpC + fpC + fnC), 3), round(np.mean(fsC), 3),
                                                      round(np.mean(prC), 3), round(np.mean(rcC), 3), sep='\t', end='\t')

                                                tpB = sum(x[0] for x in resultsB)
                                                fpB = sum(x[1] for x in resultsB)
                                                fnB = sum(x[2] for x in resultsB)
                                                fsB = [x[3] for x in resultsB]
                                                prB = [x[4] for x in resultsB]
                                                rcB = [x[5] for x in resultsB]
                                                rkB = [x[6] for x in resultsB]
                                                                        
                                                print(tpB, fpB, fnB, round(2*tpB / (2*tpB + fpB + fnB), 3), round(np.mean(fsB), 3),
                                                      round(np.mean(prB), 3), round(np.mean(rcB), 3), fsB, sep = '\t')

                                                outp = open('temp.tsv','a')
                                                toprint = '\t'.join([str(xx) for xx in [suboptmin, suboptmax, suboptsteps, minlen, minbpscore, minfinscorefactor, distcoef, orderpenalty,
                                                                                        def2,rankbydiff, rankby,
                                                                                        tpC, fpC, fnC,
                                                                                        round(2*tpC / (2*tpC + fpC + fnC), 3),
                                                                                        round(np.mean(fsC), 3),
                                                                                        round(np.mean(prC), 3),
                                                                                        round(np.mean(rcC), 3),
                                                                                        tpB, fpB, fnB,
                                                                                        round(2*tpB / (2*tpB + fpB + fnB), 3),
                                                                                        round(np.mean(fsB), 3),
                                                                                        round(np.mean(prB), 3),
                                                                                        round(np.mean(rcB), 3)]])
                                                outp.write(toprint+'\n')
                                                outp.close()


