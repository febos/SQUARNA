import os, glob
import random

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

def PrintFasta(fasta,seqnames,seqdict):
    with open(fasta,'w') as outp:
        for seqname in seqnames:
            outp.write('>'+seqname+'\n')
            outp.write(seqdict[seqname]+'\n')


if __name__ == "__main__":

    
    stofiles = glob.glob("RNAStralignExt/sto/*")

    for stofile in stofiles:

        headers, seqnames, seqdict, gcnames, gcdict = ReadStockholm(stofile)
        fam = os.path.basename(stofile).split('.')[0]

        for depth in (2,5,10,20,30,50,75,100,150,200,300):
            if depth > len(seqnames):
                break
            for tryN in range(1,11):
                sampled = random.sample(seqnames, depth)
                outstk = os.path.join("SubAli","sto","{}_{}_{}.sto".format(fam,depth,tryN))
                outafa = os.path.join("SubAli","afa","{}_{}_{}.afa".format(fam,depth,tryN))
                WriteStockholm(outstk,headers,sampled,seqdict,gcnames,gcdict)
                PrintFasta(outafa, sampled, seqdict)
