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

infiles = sorted(glob.glob("../S01Ali_init/afa_clean/*"),key = lambda x: int(os.path.basename(x).split('_')[0]))

objs = []

with open("../S01.fas") as infile:
    lines = infile.readlines()
    for i in range(0,len(lines)-4,4):
        name = lines[i].strip()[1:]
        seq = lines[i+1].strip()
        ref = lines[i+2].strip()
        react = lines[i+3].strip().split()
        objs.append((name,seq,ref,react))



for k, obj in enumerate(objs):

    infile = infiles[k]
    obj = objs[k]
    inseq = obj[1]
    ref = obj[2]
    react = obj[3]
    gaps = {}

    outfile = infile.replace("S01Ali_init","S01Ali").replace("afa_clean","mafft")

    command = "mafft --inputorder --anysymbol --auto {} > {}".format(infile,outfile)
    os.system(command)

    seqs = []
    names = []

    with open(outfile) as inp:
        for line in inp:
            if line.startswith('>'):
                name = line.strip()[1:]
                names.append(name)
                seqs.append('')
            else:
                seqs[-1] += line.strip()

    newref = ''
    newreact = []
    cur = 0

    ln = len(seqs[0])

    for x in seqs[0]:
        if x == '-':
            newref += '.'
            newreact.append('0.5')
        else:
            newref += ref[cur]
            newreact.append(react[cur])
            cur += 1

    
    fasta   = 'afa/{}.afa'.format(names[0])
    stofile = 'sto/{}.sto'.format(names[0])
    with open(fasta,'w') as outp:
        outp.write(' '.join(newreact)+'\n')
        outp.write('\n')
        outp.write(newref+'\n')
        for kkk,seq in enumerate(seqs):
            outp.write('>'+names[kkk]+'\n')
            outp.write(seq+'\n')

    headers = []
    seqnames = names
    seqdict = {names[i]:seqs[i] for i in range(len(names))}
    gcnames = ['SS_cons']
    gcdict = {'SS_cons':newref}

    WriteStockholm(stofile,headers,seqnames,seqdict,gcnames,gcdict)
    
    
