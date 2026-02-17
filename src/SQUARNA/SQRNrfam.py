
import os
import re
import urllib.request
import gzip
import shutil

try:
    from SQRNdbnseq import SEPS, GAPS, ReAlign
    from SQRNdbnseq import PairsToDBN, DBNToPairs, UnAlign
except:
    from .SQRNdbnseq import SEPS, GAPS, ReAlign
    from .SQRNdbnseq import PairsToDBN, DBNToPairs, UnAlign



def ProcessSeqLineDBNLine(start, end, origseqlen, seq, dbn):

    #handle insertions
    seqlst = seq.split('*')

    positions = [i for i in range(len(seq)) if seq[i]=='*']

    cur = -1
    dbnlst = []
    for pos in positions:
        dbnlst.append(dbn[cur+1:pos])
        cur = pos
    dbnlst.append(dbn[cur+1:])

    for k in range(len(seqlst)):

        if not seqlst[k]:
            continue

        if seqlst[k][0] == '[' or seqlst[k][-1] == ']':
            seqlst[k] = seqlst[k].strip('><')
            length = int(seqlst[k][1:-1])
            seqlst[k] = 'N'*length
            dbnlst[k] = 'N'*length

    seq = ''.join(seqlst)
    dbn = ''.join(dbnlst)

    # remove gaps
    seq, dbn = UnAlign(seq, dbn)

    # extend tails
    seq = '-'*(start-1) + seq + '-'*(origseqlen-end)
    dbn = '-'*(start-1) + dbn + '-'*(origseqlen-end)

    return DBNToPairs(dbn)


def CMScan(seq, homedir):

    shortseq = ''.join([x if x not in SEPS else "N" for x in seq if x not in GAPS])
    
    rfamcm = os.path.join(homedir,'Rfam.cm')
    infile = "squarna_cmscan.fas"
    outfile = "squarna_cmscan.out"

    illegal = {'e','f','j','l','o','p','q','z'}

    with open(infile,'w') as inp:
        inp.write(">seq\n")
        inp.write(''.join(['N' if (ch.lower() in illegal or ord(ch) > 127) else ch
                           for ch in shortseq]) + '\n') 
    
    os.system("cmscan --notextw -E 1e-4 -o {} --rfam --toponly {} {}".format(outfile, rfamcm, infile))

    with open(outfile) as file:
        flag = False
        lines = []
        for line in file:
            if line.startswith("Hit alignments:"):
                flag = True
            elif flag:
                if line.startswith("Internal HMM-only"):
                    break
                lines.append(line.strip('\n'))

    fams = []

    paired = set()
    pairs = []
    
    for k,line in enumerate(lines):
        if line.startswith(">>"):
            fam = line.split()[1]
            ls = lines[k+3].split()
            if ls[11] == '-': #ignore minus-strand hits
                continue
            desc = fam+'('+ls[9]+'-'+ls[10]+')'
            fams.append(desc)
            start = int(ls[9])
            end   = int(ls[10])

            dbnline = lines[k+6].split()[0]
            dbnline_ind = lines[k+6].find(dbnline)
            seqline = lines[k+9][dbnline_ind:dbnline_ind+len(dbnline)]            

            newpairs = ProcessSeqLineDBNLine(start, end, len(shortseq), seqline, dbnline)

            for v,w in newpairs:
                if v not in paired and w not in paired:
                    pairs.append((v,w))
                    paired.add(v)
                    paired.add(w)

    shortdbn = PairsToDBN(pairs,len(shortseq))

    dbn = ReAlign(shortdbn, seq)

    return dbn, ','.join(fams)


def G4Hscore(match):

    splt = []
    cur  = 0
    prev = 0

    N = len(match)

    while cur < N:

        if match[cur] not in {'G','C'}:
            if prev < cur:
                splt.append(match[prev:cur])
            cur += 1
            prev = cur
        elif match[cur] != match[prev]:
            splt.append(match[prev:cur])
            prev = cur
            cur += 1
        else:
            cur += 1
    if prev < cur:
        splt.append(match[prev:cur])

    tot = 0
    score = 0

    for chunk in splt:

        score += (1 - 2*(chunk[0]=='C')) * len(chunk) * min(len(chunk), 4)

    score /= N
            
    return score
    

def FindG4(seq, g4sym, scorelim = 1.2):

    patterns = (r'(?=((G{2,5})(\w{1,2}?)(G{2,5})(\w{1,2}?)(G{2,5})(\w{1,2}?)(G{2,5})))',
                r'(?=((G{3,5})(\w{1,12}?)(G{3,5})(\w{1,12}?)(G{3,5})(\w{1,12}?)(G{3,5})))',)
    patterns = (re.compile(p) for p in patterns)

    found = False

    g4 = ['.' for _ in seq]

    for k, pattern in enumerate(patterns):
        for match in pattern.finditer(seq):
            start = match.start()
            matched_string = match.group(1)
            if G4Hscore(matched_string) >= scorelim:
                found = True
                cur = start
                for i in range(2,9):
                    isG = not i % 2
                    for ch in match.group(i):
                        if isG:
                            g4[cur] = g4sym
                        cur += 1

    return ''.join(g4), found


def SearchG4(seq, rfamdbn, rfamfound, g4sym = '+'):

    shortseq = ''.join([x if x not in SEPS else "N" for x in seq if x not in GAPS]).upper()

    #####
    #shortg4 = ''.join([g4sym if ch=='G' else '.' for ch in shortseq])
    #g4found = g4sym in shortg4
    shortg4, g4found = FindG4(shortseq, g4sym)
    #####

    if not g4found:
        return rfamdbn, rfamfound

    g4 = ReAlign(shortg4, seq)

    if not rfamfound:
        return g4, "G4(+)"

    pairs = [(v,w) for v,w in DBNToPairs(rfamdbn) if g4[v] != g4sym and g4[w] != g4sym]
    res = PairsToDBN(pairs,len(seq))
    res = ''.join([ch if g4[i] != g4sym else g4sym for i,ch in enumerate(res)])
    return res, "G4(+),"+rfamfound  
    

def SearchRfamG4(seq, homedir, write_to, rfam, g4):
    """ return restraints,rfam-families """

    if not rfam:
        return SearchG4(seq, None, False)

    path = shutil.which("cmscan")

    if path is None:
        print("ERROR: could not find cmscan, rfam search disabled; to fix this, install Infernal: eddylab.org/infernal/",
              file = write_to)
        if not g4:
            return None, False
        else:
            return SearchG4(seq, None, False)

    else:
        if not os.path.exists(os.path.join(homedir,"Rfam.cm")) and\
           not os.path.exists(os.path.join(homedir,"Rfam.cm.i1f")):
            print("ERROR: could not find Rfam.cm, rfam search disabled; to fix this, run SQUARNA-build-rfam",
                  file = write_to)
            if not g4:
                return None, False
            else:
                return SearchG4(seq, None, False)
        else:
            print("Running Rfam search...", end = '', file = write_to)
            dbn, fams = CMScan(seq, homedir)
            if fams:
                print(': '+fams, file = write_to)
                if not g4:
                    return dbn, fams
                else:
                    return SearchG4(seq, dbn, fams)
            else:
                print(': no hits.', file = write_to)
                if not g4:
                    return None, False
                else:
                    return SearchG4(seq, None, False)

def BuildRfam():

    homedir = os.path.dirname(os.path.abspath(__file__))
    
    url = "https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz"

    gz_file_path     = os.path.join(homedir,"Rfam.cm.gz")
    output_file_path = os.path.join(homedir,"Rfam.cm")

    urllib.request.urlretrieve(url, gz_file_path)

    with gzip.open(gz_file_path, 'rb') as f_in:
        with open(output_file_path, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)

    os.system("cmpress {}".format(output_file_path))




'''
if __name__ == "__main__":

    seq = "GGGCCAUUGGGUGGGAUCUGGGGGGG"
    seq = "GGGCAAGGGAAAGGGCCCGGG"
    #seq = "GGCUGGUGAUUGGGACCGGGCAGGGCGGGCACGGGCCAGCC"
    g4, found = FindG4(seq, '+')

    print(seq)
    print(g4)
    print(found)
'''

