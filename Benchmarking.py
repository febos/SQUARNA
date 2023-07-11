import os, time

def StemmedIsolated(sorted_pairs):

    sp = sorted_pairs

    stemmed, isolated = [], []

    if len(sp) < 2:
        if not sp:
            return [], []
        else:
            return [],[sp[0][0],sp[0][1]]

    for i in range(len(sp)):

        if i==0:
            if sp[i][0] + 1 == sp[i+1][0] and sp[i][1] == sp[i+1][1] + 1:
                stemmed.append(sp[i])
            else:
                isolated.append(sp[i][0])
                isolated.append(sp[i][1])
                
        elif i==len(sp)-1:
            if sp[i-1][0] + 1 == sp[i][0] and sp[i-1][1] == sp[i][1] + 1:
                stemmed.append(sp[i])
            else:
                isolated.append(sp[i][0])
                isolated.append(sp[i][1])
            
        else:
            if sp[i][0] + 1 == sp[i+1][0] and sp[i][1] == sp[i+1][1] + 1 or \
               sp[i-1][0] + 1 == sp[i][0] and sp[i-1][1] == sp[i][1] + 1:
                stemmed.append(sp[i])
            else:
                isolated.append(sp[i][0])
                isolated.append(sp[i][1])

    return stemmed, isolated  

def CombinePairsToDBN(newpairs, length, initpairs=()):
    """order of keeping pairs:
        1) stemmed mapped pairs
        2) all predicted pairs removing base triples
        3) isolated mapped pairs removing base triples"""

    #print(newpairs)
    #print(length)

    dbn = ['.']*length

    levels = ['()','[]','{}','<>','Aa','Bb','Cc','Dd','Ee','Ff','Gg','Hh','Ii','Jj','Kk','Ll','Mm','Nn']

    groups = [[],]

    stemmed  = set(StemmedIsolated(sorted(initpairs))[0])
    isolated = set(initpairs) - stemmed

    seen = set()
    pairs = set()

    for p in stemmed:

        pairs.add(p)
        seen.add(p[0])
        seen.add(p[1])

    for p in newpairs:

        if p[0] not in seen and p[1] not in seen:
            pairs.add(p)
            seen.add(p[0])
            seen.add(p[1])

    for p in isolated:

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


def GetPairs(dbn):

    pairs = set()

    closing = {'>':'<',']':'[',')':'(','}':'{','a':'A','b':'B','c':'C','d':'D','e':'E'}
    stack = {'<':[],'(':[],'{':[],'[':[],'A':[],'B':[],'C':[],'D':[],'E':[]}

    for i,v in enumerate(dbn):

        if v in stack:
            stack[v].append(i)
        if v in closing:
            if stack[closing[v]]:
                pairs.add((stack[closing[v]].pop(),i))

    return sorted(pairs)
          

def NoLone(dbn):

    pairs = GetPairs(dbn)
    lone_pos_set = set(StemmedIsolated(pairs)[1])

    return ''.join([dbn[i] if i not in lone_pos_set else '.' for i in range(len(dbn))])


def PredictRNAfold(seq):

    inpf = "inp.tmp"
    with open(inpf,'w') as inp:
        inp.write('>seq\n')
        inp.write(seq+'\n')

    os.system("RNAfold --noPS < inp.tmp > outp.tmp")
    with open("outp.tmp") as outp:
        dbn = outp.readlines()[2].split()[0]
    return [dbn,]        

def PredictIPknot(seq):

    inpf = "inp.tmp"
    with open(inpf,'w') as inp:
        inp.write('>seq\n')
        inp.write(seq+'\n')

    os.system("~/software/ipknot-master/build/ipknot inp.tmp > outp.tmp")
    with open("outp.tmp") as outp:
        dbn = outp.readlines()[2].strip()
    return [dbn,] 

def PredictMXfold2(seq):

    inpf = "inp.tmp"
    with open(inpf,'w') as inp:
        inp.write('>seq\n')
        inp.write(seq+'\n')

    os.system("mxfold2 predict inp.tmp > outp.tmp")
    with open("outp.tmp") as outp:
        dbn = outp.readlines()[2].split()[0]
    return [dbn,] 

def BPSEQtoDBN(bpseqfile):

    with open(bpseqfile) as outp:
        pairs = []
        lines = outp.readlines()[1:]
        for line in lines:
            x,y,z = line.strip().split()
            x = int(x)
            z = int(z)
            if x < z:
                pairs.append((x-1,z-1))
        return CombinePairsToDBN(pairs,x)

def PredictSPOTRNA(seq):

    seqs = [seq,]

    inpf = "inp.tmp"
    names = []
    with open(inpf,'w') as inp:
        for ii,seq in enumerate(seqs):
            names.append("seq%s"%ii)
            inp.write('>seq%s\n'%ii)
            inp.write(seq+'\n')

    os.system('conda run -n spotrna python ~/software/SPOT-RNA/SPOT-RNA.py --inputs inp.tmp --outputs "tmp" --cpu 32')

    res = []

    for name in names:
        res.append(BPSEQtoDBN("tmp/{}.bpseq".format(name)))
    return res


def PredictSQUARNA(seq, conf = "def.conf", top = 1):

    with open("inp.tmp","w") as inp:
        inp.write(">seq"+'\n')
        inp.write(seq+'\n')

    if len(seq) > 500:
        conf = "long.conf"

    os.system("python3 SQUARNA.py i=inp.tmp c={} toplim={} > outp.tmp".format(conf,top))

    cnt = 0
    flag = False
    res = []

    with open("outp.tmp") as outp:
        for line in outp:

            if line.startswith("="):
                flag = True
                continue

            if flag:
                res.append(line.strip().split()[0])

    return res[:top]


def PredictSQUARNA5(seq):
    return PredictSQUARNA(seq, top = 5)

def PredictSQUARNAN(seq):
    return PredictSQUARNA(seq, top = 10**6)

def PredictSQUARNAalt(seq):
    return PredictSQUARNA(seq, conf = "alt.conf")

def PredictSQUARNAalt5(seq):
    return PredictSQUARNA(seq, conf = "alt.conf", top = 5)

def PredictSQUARNAaltN(seq):
    return PredictSQUARNA(seq, conf = "alt.conf", top = 10**6)

def PredictSQUARNAlong(seq):
    return PredictSQUARNA(seq, conf = "long.conf")

def PredictSQUARNAlong5(seq):
    return PredictSQUARNA(seq, conf = "long.conf", top = 5)

def PredictSQUARNAlongN(seq):
    return PredictSQUARNA(seq, conf = "long.conf", top = 10**6)

     

if __name__ == "__main__":

    dataset = "SRtrain"
    tool    = "SPOT-RNA"
    NL      =  False

    if NL:
        dataset += "NL"

    with open('datasets/{}.fas'.format(dataset)) as file:

        outname = "{}_{}".format(dataset,tool)
        title = '\t'.join("NAME LEN TIME RANK TP FP FN PRC RCL FS SEQ DBN PRED".split())
        outp1 = open(outname+'.fas','w')
        outp2 = open(outname+'.tsv','w')
        outp2.write(title+'\n')
        lines = file.readlines()

        t0 = time.time()
        
        for i in range(0,len(lines)-2,3):

            name = lines[i].strip()
            print(name,end='')
            seq = lines[i+1].strip().upper()
            dbn = lines[i+2].strip()

            structs = {"RNAfold":PredictRNAfold,
                       "IPknot": PredictIPknot,
                       "MXfold2":PredictMXfold2,
                       "SPOT-RNA":PredictSPOTRNA,
                       "SQUARNA": PredictSQUARNA,
                       "SQUARNA5": PredictSQUARNA5,
                       "SQUARNAN": PredictSQUARNAN,
                       "SQUARNAalt": PredictSQUARNAalt,
                       "SQUARNAalt5": PredictSQUARNAalt5,
                       "SQUARNAaltN": PredictSQUARNAaltN,
                       "SQUARNAlong": PredictSQUARNAlong,
                       "SQUARNAlong5": PredictSQUARNAlong5,
                       "SQUARNAlongN": PredictSQUARNAlongN,
                       }[tool](seq)

            t1 = time.time()-t0

            print("...COMPLETE ({}sec)".format(round(t1,3)))

            

            # Clean <3nt hairpins and non-canonical pairs
            structs = [CombinePairsToDBN([(v,w) for v,w in GetPairs(_)
                                          if w-v >= 4 and seq[v]+seq[w] in {'GC','CG',
                                                                            'GU','UG',
                                                                            'AU','UA'}],
                                         len(seq))
                       for _ in structs]

            # Clean lone bps if NL
            if NL:
                structs = [NoLone(_) for _ in structs]

            best_ind    = -1
            best_fscore = -1
            BTP, BFP, BFN = -1, -1, -1

            pairsr = set(GetPairs(dbn))

            for i,pred in enumerate(structs):
            
                pairsq = set(GetPairs(pred))

                TP = len(pairsr & pairsq)
                FP = len(pairsq - pairsr)
                FN = len(pairsr - pairsq)
            
                FS = 2*TP / (2*TP + FN + FP) if (TP + FN + FP) else 1
                PRC = (TP / (TP + FP)) if (TP+FP) else 1
                RCL = (TP / (TP + FN)) if (TP+FN) else 1

                if FS > best_fscore:
                    best_ind = i
                    best_fscore = FS
                    BTP, BFP, BFN = TP, FP, FN

            best_prc = (BTP / (BTP + BFP)) if (BTP+BFP) else 1
            best_rcl = (BTP / (BTP + BFN)) if (BTP+BFN) else 1

            outp1.write(name+'\n')
            outp1.write(seq+'\n')
            outp1.write(dbn+'\n')
            for pred in structs:
                outp1.write(pred+'\n')
            outp1.write("LEN={} TIME={}sec RANK={} TP={} FP={} FN={} PRC={} RCL={} FS={}\n"\
                        .format(len(seq),round(t1,3),best_ind+1,
                                BTP,BFP,BFN,
                                round(best_prc,3),
                                round(best_rcl,3),
                                round(best_fscore,3)))
            res = [name[1:], len(seq), round(t1,3), best_ind+1, BTP,BFP,BFN,
                   round(best_prc,3), round(best_rcl,3), round(best_fscore,3),
                   seq,dbn,structs[best_ind]]
            outp2.write('\t'.join([str(g) for g in res])+'\n')

    outp1.close()
    outp2.close()


