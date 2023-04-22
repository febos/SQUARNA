
import numpy as np

def PairsToDBN(newpairs, length):
    
    dbn = ['.']*length

    levels = ['()','[]','{}','<>','Aa','Bb','Cc','Dd','Ee','Ff','Gg',
              'Hh','Ii','Jj','Kk','Ll','Mm','Nn','Oo','Pp','Qq','Rr',
              'Ss','Tt','Uu','Vv','Ww','Xx','Yy','Zz']

    groups = [[],]

    seen = set()
    pairs = set()

    for p in newpairs:

        if p[0] not in seen and p[1] not in seen:
            pairs.add(p)
            seen.add(p[0])
            seen.add(p[1])
    
    for pair in sorted(pairs):

        level = 0

        while any(v[0]<=pair[0]<=v[1]<=pair[1] or pair[0]<=v[0]<=pair[1]<=v[1] for v in groups[level]):
            level += 1
            if level == len(groups):
                groups.append([])

        groups[level].append(pair)

    for times in range(len(groups)-1):
    
        for i in range(len(groups)-1):

            test = [v for v in groups[i] if any(v[0]<=w[0]<=v[1]<=w[1] or w[0]<=v[0]<=w[1]<=v[1] for w in groups[i+1])]

            if len(test) < len(groups[i+1]):

                groups[i]   = [p for p in groups[i] if p not in test] + groups[i+1]
                groups[i+1] = test

    for i,g in enumerate(groups):

        for p in g:

            dbn[p[0]] = levels[i][0]
            dbn[p[1]] = levels[i][1]
            
    return ''.join(dbn)


def StemsToDBN(stems, seq):
    """Convert a list of stems (lists of bps) into a dbn string"""
    return PairsToDBN([bp for stem in stems for bp in stem],len(seq))


def DBNToPairs(dbn):

    pairs = set()

    closing = {'>':'<',']':'[',')':'(','}':'{','a':'A','b':'B','c':'C','d':'D',
               'e':'E','f':'F','g':'G','h':'H','i':'I','j':'J','k':'K','l':'L',
               'm':'M','n':'N','o':'O','p':'P','q':'Q','r':'R','s':'S','t':'T',
               'u':'U','v':'V','w':'W','x':'X','y':'Y','z':'Z',}
    stack = {'<':[],'(':[],'{':[],'[':[],'A':[],'B':[],'C':[],'D':[],'E':[],
             'F':[],'G':[],'H':[],'I':[],'J':[],'K':[],'L':[],'M':[],'N':[],
             'O':[],'P':[],'Q':[],'R':[],'S':[],'T':[],'U':[],'V':[],'W':[],
             'X':[],'Y':[],'Z':[],}

    for i,v in enumerate(dbn):

        if v in stack:
            stack[v].append(i)
        if v in closing:
            if stack[closing[v]]:
                pairs.add((stack[closing[v]].pop(),i))

    return sorted(pairs)


def ReAlign(shortdbn,longseq,seqmode=False):

    dbn = [x for x in shortdbn]

    assert len(shortdbn)+longseq.count('.')+longseq.count('-') == len(longseq),\
    "Cannot ReAlign dbn string - wrong number of gaps:\n{}\n{}".format(longseq,shortdbn)

    newdbn = []

    for x in longseq:
        if x in ('.','-','~'):
            if seqmode:
                newdbn.append('-')
            else:
                newdbn.append('.')
        else:
            newdbn.append(dbn.pop(0))

    return ''.join(newdbn)


def UnAlign2(seq, dbn):

    newdbn = ''.join([dbn[i] for i in range(len(seq)) if seq[i] not in ('-','.','~')])
    newseq = ''.join([x for x in seq if x not in ('-','.','~')])

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

                    
def OptimalStems(seq, bpmatrix, rbps, rxs, stems):

    ## remove from rbps all bps that are in stems

    stems = AnnotateStems(bpmatrix, rbps, rxs, rstembps)
    #stems = ScoreStems(seq, stems, rstembps)

    

    return stems
            

        

    


def SQUARNA(seq, bpweights, restraints = None, dbn = None):
    """seq == sequence (possibly with gaps or any other non-ACGU symbols
    bpweights == dictionary with canonical bps as keys and their weights as values
    restraints == string in dbn format; x symbols are forced to be unpaired
    dbn == known secondary structure in dbn format

    SQUARNA returns a list of alternative predicted secondary structures in dbn format"""

    # turn seq into UPPERCASE & replace T with U
    seq = seq.upper().replace("T","U")

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
                            rbps.copy(), rxs.copy(), [])

    # List of lists of stems (each stem list is a currently predicted secondary structure
    curstemsets = [[stem,] for stem in newstems]

    # List of finalized stem lists
    finstemsets = []

    while curstemsets:

        newcurstemsets = []
    
        for stems in curstemsets:
    
            newstems = OptimalStems(shortseq, bpmatrix.copy(),
                                    rbps.copy(), rxs.copy(), stems):
            if newstems:
                for newstem in newstems:
                    newcurstemsets.append(stems + [newstem,])
            else:
                finstemsets.append(stems)

    dbns = []

    for stems in finstemsets:

        dbns.append({bp for stem in stems for bp in stem} | set(rbps),
                    len(shortseq))

    dbns = [ReAlign(x, seq) for x in dbns]

    sorteddbns = SortDBNs(dbns) #Sort in decreasing order of structure "quality"
    cons       = DBNsToConsensus(sorteddbns) #Consensus dbn

    return cons, sorteddbns

    ### finalize the SQUARNA func + comments!!
    ###parameters: subopt minlen minscore bracketweight distcoef
    ### comment PairsToDBN & DBNToPairs
    ### ReAlign & UnAlign2
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

    dbns = SQUARNA(seq, bpweights, rst, dbn)














    
