import glob, os

def ReadStockholm(stkfile):

    seqnames = []
    seqdict  = {}
    gcnames  = []
    gcdict   = {}
    headers  = []

    try:
        file = open(stkfile)
    except:
        file = open(stkfile,encoding="iso8859-15")

    for line in file:
        if line.startswith('#=GC '):

            linesplit = line.strip().split()
            seq = linesplit[-1]
            name = ' '.join(linesplit[1:-1])

            if name not in gcdict:
                gcnames.append(name)
                gcdict[name] = seq
            else:
                gcdict[name] += seq

        elif line.startswith('#'):

            headers.append(line)

        elif line.startswith('//'):
            pass
        elif not line.strip():
            pass
        else:

            linesplit = line.strip().split()
            seq = linesplit[-1]
            name = ' '.join(linesplit[:-1])

            if name not in seqdict:
                seqnames.append(name)
                seqdict[name] = seq
            else:
                seqdict[name] += seq

    file.close()

    headers1 = [x for x in headers if not x.startswith("#=GF SQ")]
    headers2 = [x for x in headers if x.startswith("#=GF SQ")]

    headers = headers1 + headers2

    return headers, seqnames, seqdict, gcnames, gcdict


def WriteStockholm(stkfile,headers,seqnames,seqdict,gcnames,gcdict):

    famname = os.path.basename(stkfile).split('_')[0]

    newseqnames = []

    for seqname in seqnames:

        newseqname = famname+'__'+seqname if not seqname.startswith(famname) else seqname
        newseqnames.append(newseqname)
        seqdict[newseqname] = seqdict[seqname]

    for i in range(len(headers)):
        if headers[i].startswith("#=GR "):
            grname = headers[i].strip().split()[1]
            if not grname.startswith(famname):
                headers[i] = headers[i].replace(grname,famname+'__'+grname)

    seqnames = newseqnames

    grlens = [len(' '.join(h.strip().split()[:-1])) for h in headers if h.startswith("#=GR ")]

    longest = max([4,] + [len(x) for x in seqnames] + [len(x)+len('#=GC ') for x in gcnames] + grlens) + 1

    if not headers or 'STOCKHOLM' not in headers[0]:
        headers.insert(0,'# STOCKHOLM 1.0\n')

    if headers[0] != '# STOCKHOLM 1.0\n':
        headers[0] = '# STOCKHOLM 1.0\n'

    seqnum = len(seqnames)

    if headers[-1].strip() != '#=GF SQ {}'.format(seqnum):
        headers.append('#=GF SQ {}\n'.format(seqnum))

    with open(stkfile,'w') as file:

        for line in headers:
            if not line.startswith("#=GR "):
                file.write(line)

        file.write('\n')

        for name in seqnames:
            file.write(name + ' '*(longest-len(name)) + seqdict[name] + '\n')

        for line in headers:
            if line.startswith("#=GR "):
                file.write(line)
        
        for name in gcnames:
            file.write('#=GC ' + name + ' '*(longest-len(name)-len('#=GC ')) + gcdict[name] + '\n')

        file.write('//\n')

def Process(inseq, data):

    res1 = '-'*(int(data[0][0][0])-1)
    res3 = '-'*(len(inseq)-int(data[0][-1][-1]))
    res2 = ''.join([x[1] for x in data[1]])

    q = inseq[:len(res1)]+''.join([x[1] for x in data[0]])+inseq[-len(res3):]
    seq, gaps = res1+res2+res3, [len(q[:i].replace('-','')) for i in range(len(q)) if q[i]=='-']
    d = {}
    for g in gaps:
        if g not in d:
            d[g] = 0
        d[g] += 1

    return seq.upper().replace('T','U'), d, seq.replace('-','')


def Sim(seq1,seq2):

    return sum(seq1[h]==seq2[h] for h in range(len(seq1)))/len(seq1)


infiles = sorted(glob.glob("blastn/*"),key = lambda x: int(os.path.basename(x).split('_')[0]))

objs = []

with open("../S01raw.fas") as infile:
    lines = infile.readlines()
    for i in range(0,len(lines)-4,4):
        name = lines[i].strip()[1:]
        seq = lines[i+1].strip()
        ref = lines[i+2].strip()
        react = lines[i+3].strip().split()
        objs.append((name,seq,ref,react))

mincov = 0.9

for k, obj in enumerate(objs):

    infile = infiles[k]
    obj = objs[k]
    inseq = obj[1]
    ref = obj[2]
    react = obj[3]
    seqs = [(inseq,{}),]
    gaps = {}
    names = [os.path.basename(infile)[:-4],]

    seen = set()
    seen.add(inseq)
    
    with open(infile) as file:
        data = [[],[]]
        for line in file:
            if line.startswith("Sequence ID:"):
                if data[0]:
                    seq, gap, unaliseq = Process(inseq,data)
                    if unaliseq not in seen and len(unaliseq)/len(inseq)>=mincov and name not in names:
                        seen.add(unaliseq)
                        names.append(name)
                        seqs.append((seq,gap))
                        for g in gap:
                            if g not in gaps or gap[g] > gaps[g]:
                                gaps[g] = gap[g]
                data = [[],[]]
                name = line.strip().split()[2]
            elif line.startswith("Range"):
                linesplit = line.strip().split()
                if data[0]:
                    seq, gap, unaliseq = Process(inseq,data)
                    if unaliseq not in seen and len(unaliseq)/len(inseq)>=mincov and name not in names:
                        seen.add(unaliseq)
                        names.append(name)
                        seqs.append((seq,gap))
                        for g in gap:
                            if g not in gaps or gap[g] > gaps[g]:
                                gaps[g] = gap[g]
                name = name.split('_')[0]+'_'+linesplit[2]+'_'+linesplit[4]
                data = [[],[]]
            elif line.startswith("Query "):
                data[0].append(line.strip().split()[1:])
            elif line.startswith("Sbjct"):
                data[1].append(line.strip().split()[1:])

        seq, gap, unaliseq = Process(inseq,data)
        
        if unaliseq not in seen and len(unaliseq)/len(inseq)>=mincov and name not in names:
            seen.add(unaliseq)
            names.append(name)
            seqs.append((seq,gap))
            for g in gap:
                if g not in gaps or gap[g] > gaps[g]:
                    gaps[g] = gap[g]

    aliseqs = []

    for kk, tok in enumerate(seqs):
        seq, gap = tok
        newseq = ''
        cur = 0
        for g in sorted(gaps.keys()):
            newseq += seq[cur:g]
            k = gaps[g] if g not in gap else gaps[g]-gap[g]
            newseq += '-'*k
            cur = g
        newseq += seq[cur:]
        #print(">"+names[kk])
        #print(newseq)
        aliseqs.append(newseq)

    print(names[0])
    print(len(aliseqs))
    if len(aliseqs) > 100:

        newseqs  = []
        newnames = []

        i,jj = 0, -1

        mn = 100

        for j in range(i+1,len(aliseqs)):
            if Sim(aliseqs[i],aliseqs[j]) < mn:
                ii = i
                jj = j
                mn = Sim(aliseqs[i],aliseqs[j])

        curst = set([ii,jj])
        cur = [aliseqs[ii],aliseqs[jj]]


        while len(cur) < 100:
            print(len(cur))
            mn = 1000
            kk = -1
            for kkk, useq in enumerate(aliseqs):
                cmn = max([Sim(useq,cur[l]) for l in range(len(cur))])
                if cmn < mn and kkk not in curst:
                    kk = kkk
                    mn = cmn
            cur.append(aliseqs[kk])
            curst.add(kk)

        aliseqs = [aliseqs[i] for i in range(len(names)) if i in curst]
        names = [names[i] for i in range(len(names)) if i in curst]

    print(len(aliseqs))

    newref = ''
    newreact = []
    cur = 0

    ln = len(aliseqs[0])

    allgap = {i for i in range(ln) if all(sq[i]=='-' for sq in aliseqs)}

    aliseqs = [''.join([seq[i] for i in range(ln) if i not in allgap]) for seq in aliseqs] 

    for x in aliseqs[0]:
        if x == '-':
            newref += '.'
            newreact.append('0.5')
        else:
            newref += ref[cur]
            newreact.append(react[cur])
            cur += 1

    fasta   = 'afa_raw/{}.afa'.format(names[0])
    stofile = 'sto2/{}.sto'.format(names[0])
    with open(fasta,'w') as outp:
        outp.write(' '.join(newreact)+'\n')
        outp.write('\n')
        outp.write(newref+'\n')
        for kkk,seq in enumerate(aliseqs):
            outp.write('>'+names[kkk]+'\n')
            outp.write(seq+'\n')

    headers = []
    seqnames = names
    seqdict = {names[i]:aliseqs[i] for i in range(len(names))}
    gcnames = ['SS_cons']
    gcdict = {'SS_cons':newref}

    WriteStockholm(stofile,headers,seqnames,seqdict,gcnames,gcdict)
    
    
