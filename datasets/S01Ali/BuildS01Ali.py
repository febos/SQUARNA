

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
    return seq, d

infile = "1_PreQ1_riboswitch_B_subtilis.txt"

mincov = 90

seqs = []
gaps = {}
names = [infile[:-4],]

with open(infile) as file:
    data = [[],[]]
    for line in file:
        if line.startswith(">>"):
            if data[0]:
                if cov >= mincov:
                    seq, gap = Process(inseq,data)
                    names.append(name)
                    seqs.append((seq,gap))
                    for g in gap:
                        if g not in gaps or gap[g] > gaps[g]:
                            gaps[g] = gap[g]
            data = [[],[]]
            name = line.strip().split()[1]
        elif line.startswith("Query: "):
            inseq = line.strip().split()[1]
            seqs.append((inseq,{}))
        elif line.startswith("Query"):
            data[0].append(line.strip().split()[1:])
        elif line.startswith("Sbjct"):
            data[1].append(line.strip().split()[1:])
        elif line.startswith("E-value"):
            cov = float(line.strip().split()[6][:-1])
    if cov >= mincov:
        seq, gap = Process(inseq,data)
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
    print(">"+names[kk])
    print(newseq)
    aliseqs.append(newseq)

sims = []

seqs = aliseqs
print(len(seqs))
print(len(set(seqs)))

st = set()

useqs = []
unames = []


for kk,seq in enumerate(seqs):
    if seq not in st:
        useqs.append(seq)
        unames.append(names[kk])
        st.add(seq)
    

def Sim(seq1,seq2):

    return sum(seq1[k]==seq2[k] for k in range(len(seq1)))/len(seq1)

i,jj = 0, -1

mn = 100

for j in range(i+1,len(useqs)):
    if Sim(useqs[i],useqs[j]) < mn:
        ii = i
        jj = j
        mn = Sim(useqs[i],useqs[j])

curst = set([ii,jj])
cur = [useqs[ii],useqs[jj]]


while len(cur) < 20:

    mn = 1000
    kk = -1
    for k, useq in enumerate(useqs):
        cmn = max([Sim(useq,cur[l]) for l in range(len(cur))])
        if cmn < mn and k not in curst:
            kk = k
            mn = cmn
    cur.append(useqs[kk])
    curst.add(kk)

finseqs = []
finnames = []

for i in sorted(curst):
    finseqs.append(useqs[i])
    finnames.append(unames[i])

allgaps = set([k for k in range(len(finseqs[0]))
               if all(seq[k]=='-' for seq in finseqs)])

for k, seq in enumerate(finseqs):
    print(">"+finnames[k])
    print(''.join([seq[f] for f in range(len(seq)) if f not in allgaps]))

