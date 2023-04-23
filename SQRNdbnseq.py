
import numpy as np

def PairsToDBN(newpairs, length):
    """Convert a list of base pairs into a dbn string of the given length"""

    # Initialize the dbn string
    dbn = ['.']*length

    # Define "brackets" for 30 pseudoknot levels
    levels = ['()','[]','{}','<>','Aa','Bb','Cc','Dd','Ee','Ff','Gg',
              'Hh','Ii','Jj','Kk','Ll','Mm','Nn','Oo','Pp','Qq','Rr',
              'Ss','Tt','Uu','Vv','Ww','Xx','Yy','Zz']

    # groups of non-conflicting base pairs
    groups = [[],]

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

    # add all the pairs to the dbn string
    # according to their levels
    for i, group in enumerate(groups):
        for pair in group:
            dbn[pair[0]] = levels[i][0]
            dbn[pair[1]] = levels[i][1]
            
    return ''.join(dbn)


def StemsToDBN(stems, seq):
    """Convert a list of stems (lists of bps) into a dbn string"""
    return PairsToDBN([bp for stem in stems for bp in stem],len(seq))


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
    mat = np.zeros((len(seq),len(seq)),dtype=int)

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


def AnnotateStems(bpmatrix, rbps, rxs, rstembps):

    matrix = bpmatrix.copy()
    N = bpmatrix.shape[0]
    stems = []

    for v,w in rbps:

        matrix[v,:] *= 0
        matrix[:,v] *= 0
        matrix[w,:] *= 0
        matrix[:,w] *= 0
        matrix[v,w] = bpmatrix[v,w]

    for stem in rstembps:
        for v,w in stem:
            matrix[v,:] *= 0
            matrix[:,v] *= 0
            matrix[w,:] *= 0
            matrix[:,w] *= 0

    for i in rxs:
        matrix[i,:] *= 0
        matrix[:,i] *= 0

    diagstarts  = [(0, x)   for x in range(1, N)]
    diagstarts += [(y, N-1) for y in range(1,N-1)]

    for x, y in diagstarts:

        i, j = x, y
        
        diag = [[],[]]
        
        while i < j:
            diag[0].append(i)
            diag[1].append(j)
            i += 1
            j -= 1

        arr = matrix[diag]

        candarrs = []

        instem = False
        for k in range(arr.shape[0]):
            if arr[k] != 0:
                if not instem:
                    instem = True
                    candarrs.append([k,])
            elif instem:
                candarrs[-1].append(k)
                instem = False
        if instem:
            candarrs[-1].append(arr.shape[0])

        print(arr)
        print(candarrs,[arr[x[0]:x[1]] for x in candarrs])

        for s,t in candarrs:
            for v,w in MaxSubarrays(arr[s:t]):
                stems.append([(diag[0][z], diag[1][z]) for z in range(s + v, s + w)])

    return stems


def DistFactor(distance, coef):
    """Multiplication factor for a stem score based on the confining distance"""
    return (1/(1+abs(distance-4)))**coef

 
def OptimalStems(seq, bpmatrix, rbps = set(), rxs = set(), stems = [],
                 subopt = 1.0, minlen = 2, minscore = 6,
                 bracketweight = 1.0, distcoef = 0.1):

    ## remove from rbps all bps that are in stems

    stems = AnnotateStems(bpmatrix, rbps, rxs, rstembps)
    #stems = ScoreStems(seq, stems, rstembps)

    

    return stems
            

        

    


def SQRNdbnseq(seq, bpweights, restraints = None, dbn = None,
               subopt = 1.0, minlen = 2, minscore = 6,
               bracketweight = 1.0, distcoef = 0.1):
    """seq == sequence (possibly with gaps or any other non-ACGU symbols
    bpweights == dictionary with canonical bps as keys and their weights as values
    restraints == string in dbn format; x symbols are forced to be unpaired
    dbn == known secondary structure in dbn format
    subopt == what share of the top stem score will be searched for alternative stems
    minlen == minimum stem length to predict
    minscore == minimum stem bp-score to predict
    bracketweight == mult factor for pseudoknotted brackets in the stem distance calculation
    distcoef == how much the stem distance affects the stem score (0 == no effect)
    
    SQRNdbnseq returns a list of alternative predicted secondary structures in dbn format"""

    # turn seq into UPPERCASE & replace T with U
    seq = seq.upper().replace("T", "U")

    # if not restraints - generate a string of dots
    if not restraints:
        restraints = '.'*len(seq)
    # if not known dbn - generate a string of dots
    if not dbn:
        dbn = '.'*len(seq)

    # Unalign seq, dbn, and restraints strings
    shortseq, shortrest = UnAlign2(seq, restraints)
    shortseq, shortdbn  = UnAlign2(seq, dbn)

    # Generate initial bp-matrix (with bp weights in cells)
    bpmatrix = BPMatrix(shortseq, bpweights)
    
    # Parse restraints into unpaired bases (rxs) and base pairs (rbps)
    rbps, rxs = ParseRestraints(shortrest.lower())

    newstems = OptimalStems(shortseq, bpmatrix.copy(),
                            rbps.copy(), rxs.copy(), [],
                            subopt, minlen, minscore,
                            bracketweight, distcoef)

    # List of lists of stems (each stem list is a currently predicted secondary structure
    curstemsets = [[stem,] for stem in newstems]

    # List of finalized stem lists
    finstemsets = []

    # while list of intermediate stem lists is not empty
    while curstemsets:

        # new iteration
        newcurstemsets = []

        for stems in curstemsets:

            # new optimal stems based on the current stem list    
            newstems = OptimalStems(shortseq, bpmatrix.copy(),
                                    rbps.copy(), rxs.copy(), stems,
                                    subopt, minlen, minscore,
                                    bracketweight, distcoef):
            # append new intermediate stem lists
            if newstems:
                for newstem in newstems:
                    newcurstemsets.append(stems + [newstem,])
            # if no newstems returned - the stem list is considered final
            else:
                finstemsets.append(stems)

    # list of dbn strings
    dbns = []

    # convert all final stem lists into dbn strings
    # and not forget about non-predicted bps from restraints
    for stems in finstemsets:

        dbns.append({bp for stem in stems for bp in stem} | set(rbps),
                    len(shortseq))

    # ReAlign the dbn strings accoring to seq
    dbns = [ReAlign(x, seq) for x in dbns]

    sorteddbns = SortDBNs(dbns) #Sort in decreasing order of structure "quality"
    cons       = DBNsToConsensus(sorteddbns) #Consensus dbn

    return cons, sorteddbns

    ### Finalize the OptimalStems func
    ###MaxSubarrays func


if __name__ == "__main__":

    #SAM riboswitch
    #seq = "GUUCUUAUCAAGAGAAGCAGAGGGACUGGCCCGACGAAGCUUCAGCAACCGGUGUAAUGGCGAAAGCCAUGACCAAGGUGCUAAAUCCAGCAAGCUCGAACAGCUUGGAAGAUAAGAACA"
    #dbn = "(((((((((....(((((...(((.[[[[)))......)))))(((..(((((...(((((....))))).)))..)).)))...(]]]](((((.......)))))..))))))))))." 
    #TRNA
    #seq = "GGGCGGCUAGCUCAGCGGAAGAGCGCUCGCCUCACACGCGAGAGGUCGUAGGUUCAAGUCCUACGCCGCCCACCA"
    #dbn = "(((((((..((((....[..)))).(((((.......))))).....(((((..]....))))))))))))...."
    #Twister ribozyme
    seq = "GAAAUAAUGCUACCAGACUGCACAAAAAGUCUGCCCGUUCCAAGUCGGGAUGAAAGCAGAGGAUUUC"
    dbn = "((((...(((({{((((((........))))))(((..[[[..}}.))).....))))..]]]))))"
    rst = "(.((xxx....xx..............,,,,,,..............................)).)"

    bpweights = {'GU':-1,'AU':2,'GC':4}

    subopt = 1.0
    minlen = 2
    minscore = 6
    bracketweight = 1.0
    distcoef = 0.1

    dbns = SQRNdbnseq(seq, bpweights, rst, dbn,
                      subopt, minlen, minscore,
                      bracketweight, distcoef)














    
