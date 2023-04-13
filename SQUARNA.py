
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

    bases = set(seq)

    bps = {}

    for bp in weights:
        bps[bp] = weights[bp]
        bps[bp[1]+bp[0]] = weights[bp]

    for b1 in bases:
        for b2 in bases:
            if b1+b2 not in bps:
                bps[b1+b2] = 0

    mat = np.zeros((len(seq),len(seq)),dtype=int)

    for i in range(len(seq)-1):
        for j in range(i+1,len(seq)):
            mat[i,j] = bps[seq[i]+seq[j]]

    return mat


def PrintMatrix(seq, matrix, dbn1='', dbn2=''):

    print('',*list(seq),sep='\t')

    Pairs1 = DBNToPairs(dbn1)
    Pairs2 = DBNToPairs(dbn2)

    for i in range(len(seq)):
        print(seq[i], end='\t')
        line = []
        for j in range(len(seq)):
            x = str(matrix[i][j])
            if (i,j) in Pairs1:
                x = '['+x+']'
            if (i,j) in Pairs2:
                x = '('+x+')'
            line.append(x)
        print(*line, sep='\t')


def ParseRestraints(restraints):

    rbps = DBNToPairs(restraints)
    rxs = [i for i in range(len(restraints)) if restraints[i]=='x']

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

                    
def Stems(seq, bpmatrix, rbps, rxs, rstembps):

    stems = AnnotateStems(bpmatrix, rbps, rxs, rstembps)
    #stems = ScoreStems(seq, stems, rstembps)

    return stems
            

        

    


def SQUARNA(seq, restraints, dbn, bpweights):

    seq = seq.upper()
    if not restraints:
        restraints = '.'*len(seq)
    if not dbn:
        dbn = '.'*len(seq)

    shortseq, shortrest = UnAlign2(seq, restraints)
    shortseq, shortdbn  = UnAlign2(seq, dbn)

    dbns = []
    rstembps = []
    predstems = []

    bpmatrix = BPMatrix(shortseq, bpweights)
    
    rbps, rxs = ParseRestraints(shortrest.lower())

    stems = Stems(shortseq, bpmatrix, rbps, rxs, rstembps)

    #PrintMatrix(shortseq, bpmatrix, dbn)
    #PrintMatrix(shortseq, stems, dbn)

    #curstems = PredictStem(stems, stemmatrix)

    #while curstems:

    #    pass 

    dbns = []

    dbns = [ReAlign(x,seq) for x in dbns]

    return dbns


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
    rst = "((((xxx....................xxxxxx..............................))))"

    bpweights = {'GU':-1,'AU':2,'GC':4}

    dbns = SQUARNA(seq, rst, dbn, bpweights)














    
