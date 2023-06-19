import os
import sys


from SQRNdbnseq import SQRNdbnseq, ReactDict, SEPS


def ParseConfig(configfile):

    params = {"bpweights",
              "suboptmax",
              "suboptmin",
              "suboptsteps",
              "minlen",
              "minbpscore",
              "minfinscorefactor",
              "distcoef",
              "bracketweight",
              "orderpenalty",
              "loopbonus",
              "maxstemnum"}

    paramsets = []
    names = []
    cnt = 0

    with open(configfile) as file:
        for line in file:
            cleanline = line.split('#', 1)[0].strip()
            if cleanline:
                if cleanline.startswith('>'):
                    names.append(cleanline[1:])
                    cnt += 1
                    if cnt == 1:
                        paramset = {}
                    else:  
                        paramsets.append(paramset)
                        paramset = {k:v for k, v in paramsets[0].items()}
                else:
                    key, val = cleanline.split(maxsplit = 1)
                    if key == "bpweights":
                        paramset[key] = {}
                        for kv in val.split(','):
                            k, v = kv.strip().split('=')
                            paramset[key][k] = float(v)
                    else:
                        paramset[key] = float(val)
    paramsets.append(paramset)

    if not all([_ in paramsets[0] for _ in params]):
        raise ValueError("Missing some of the parameters in"+\
                         " the first parameter set"+\
                         " of the config file: {}"\
                         .format(', '.join([_ for _ in params
                                            if _ not in paramset])))
    return names, paramsets


def EncodedReactivities(seq, reacts, reactformat):

    if reactformat == 3:
        reactline = ''.join(["_+##"[int(x * 3)]
                             for x in reacts])
    elif reactformat == 10:
        reactline = ''.join(['0123456789XX'[int(x * 11)]
                             for x in reacts])
    else:
        reactline = ''.join(['abcdefghijklmnopqrstuvwxyzz'[int(x * 26)]
                             for x in reacts])

    # Introduce the chain separators
    reactline = ''.join([reactline[i] if seq[i] not in SEPS else seq[i]
                         for i in range(len(seq))])
    return reactline


def RunSQRNdbnseq(name, data, inputformat, paramsetnames,
                  paramsets, threads, rankbydiff, rankby,
                  hardrest, interchainonly, toplim, outplim,
                  conslim, reactformat):

    while len(data) < len(inputformat):
        data.append(None)

    sequence     = data[inputformat.index('q')]
    reactivities = data[inputformat.index('t')] if 't' in inputformat else None
    restraints   = data[inputformat.index('r')] if 'r' in inputformat else None
    reference    = data[inputformat.index('f')] if 'f' in inputformat else None

    try:
        if reactivities:
            if len(reactivities) != len(sequence):
                reactivities = list(map(float, reactivities.split()))
            else:
                reactivities = [ReactDict[char] for char in reactivities]

        assert not reactivities or len(reactivities) == len(sequence)
    except:
        raise ValueError('Inappropriate reactivities line for entry "{}":\n {}'\
                         .format(name[1:], data[inputformat.index('t')]))

    assert not restraints or len(restraints) == len(sequence),\
           'Inappropriate restraints line for entry "{}":\n {}'\
           .format(name[1:], data[inputformat.index('r')])

    assert not reference or len(reference) == len(sequence),\
           'Inappropriate reference line for entry "{}":\n {}'\
           .format(name[1:], data[inputformat.index('f')])

    prediction = SQRNdbnseq(sequence, reactivities, restraints, reference,
                            paramsets, conslim, toplim, hardrest,
                            rankbydiff, rankby, interchainonly, threads)

    consensus, predicted_structures, consensus_metrics, topN_metrics = prediction

    print(name)
    print(sequence)
    #########################################################################
    if reactivities:
        print(EncodedReactivities(sequence,
                                  reactivities,
                                  reactformat),
              "reactivities", sep = '\t')
    if restraints:
        print(restraints,
              "restraints", sep = '\t')
    if reference:
        print(reference,
              "reference", sep = '\t')

    print('_'*len(sequence))

    if reference:
        print(consensus,
              "top-{}_consensus".format(conslim),
              "TP={},FP={},FN={},FS={},PR={},RC={}".format(*consensus_metrics),
              sep = '\t')
    else:
        print(consensus,
              "top-{}_consensus".format(conslim), sep = '\t')

    print('='*len(sequence))

    for i, pred in enumerate(predicted_structures[:outplim]):
        
        struct, scores, paramsetind = pred
        totalscore, structscore, reactscore = scores

        if reference and i + 1 == topN_metrics[-1]:
            print(struct, "#{}".format(i+1), totalscore,
                  structscore, reactscore,
                  paramsetnames[paramsetind],
                  "TP={},FP={},FN={},FS={},PR={},RC={},RK={}".format(*topN_metrics),
                  sep='\t')
        else:
            print(struct, "#{}".format(i+1), totalscore,
                  structscore, reactscore,
                  paramsetnames[paramsetind], sep='\t')


if __name__ == "__main__":

    def PrintUsage():
        print()
        print("Usage:")
        print()
        print('pathto/python3 pathto/SQUARNA.py i=inputfile [OPTIONS]')
        print()
        print("For further details read the help message:")
        print()
        print('pathto/python3 pathto/SQUARNA.py --help')
        print()
        exit(1)

    HOME_DIR = os.path.dirname(os.path.abspath(__file__))
    
    args = sys.argv[1:]

    if not args:
        PrintUsage()
    if "--help" in args or "-help" in args or "help" in args or\
       "--h" in args or "-h" in args or "h" in args or\
       "--H" in args or "-H" in args or "H" in args:
        with open(os.path.join(HOME_DIR,"README.md")) as helpfile:
            print(helpfile.read())
        exit(0)

    # DEFAULTS
    inputfile  = os.path.join(HOME_DIR, "examples", "seq_input.fas")
    configfile = os.path.join(HOME_DIR, "def.conf")

    inputformat = "qtrf"               # Input line order, q=seQuence,t=reacTivities,r=Restraints,f=reFerence

    maxstemnumset = False              # do we overwrite the maxstemnum from configfile
    maxstemnum    = 10**6              # maximum number of stems for each structure

    threads       = os.cpu_count()     # Number of cpus to use

    rankbydiff     = False             # Output diverse structures first
    rankby         = "rs"              # Rank by, r / s / rs / dr / ds / drs, r=reactscore,s=structscore,d=rankbydiff
    hardrest       = False             # Force bp-restraints into predicted structures 
    interchainonly = False             # Forbid intra-chain base pairs

    toplim        = 5                  # Top-N to print
    outplim       = toplim             # Top-N structs used for metrics calculations if reference
    conslim       = 1                  # Top-N structs used for consensus

    reactformat   = 3                  # 3 / 10 / 26

    # Parsing arguments
    for arg in args:
        # inputfile
        if arg.lower().startswith("i=") or\
           arg.lower().startswith("input="):
            inputfile = arg.split('=', 1)[1]
            assert os.path.exists(inputfile), "Input file does not exist."
        # configfile
        elif arg.lower().startswith("c=") or\
           arg.lower().startswith("config="):
            configfile = arg.split('=', 1)[1]
            assert os.path.exists(configfile), "Config file does not exist."
        # inputformat
        elif arg.lower().startswith("if=") or\
             arg.lower().startswith("inputformat="):
            inputformat = arg.split('=', 1)[1].lower()
            assert ''.join(sorted(inputformat)) in {"q","fq","qr","qt", "qrt",
                                                    "fqr", "fqt", "fqrt"}, \
                   'Inappropriate inputformat value (subset of "fqrt" with "q" being mandatory): {}'\
                   .format(arg.split('=', 1)[1])
        # maxstemnum
        elif arg.lower().startswith("msn=") or\
             arg.lower().startswith("maxstemnum="):
            try:
                maxstemnum = int(float(arg.split('=', 1)[1]))
                assert maxstemnum >= 0
                maxstemnumset = True
            except:
                raise ValueError("Inappropriate maxstemnum value (non-negative integer): {}"\
                                 .format(arg.split('=', 1)[1]))
        # threads
        elif arg.lower().startswith("t=") or\
             arg.lower().startswith("threads="):
            try:
                threads = int(float(arg.split('=', 1)[1]))
                threads = max(1, threads)
                threads = min(threads, os.cpu_count())
            except:
                raise ValueError("Inappropriate threads value (integer): {}"\
                                 .format(arg.split('=', 1)[1]))
        # rankby
        elif arg.lower().startswith("rb=") or\
             arg.lower().startswith("rankby="):
            rankby = ''.join(sorted(arg.split('=', 1)[1].lower()))
            assert rankby in {"r", "s", "rs", "dr", "ds", "drs"}, \
                   'Inappropriate rankby value (r/s/rs/dr/ds/drs): {}'\
                   .format(arg.split('=', 1)[1])
        # hardrest
        elif arg.lower() in {"hr", "hardrest"}:
            hardrest = True
        # interchainonly
        elif arg.lower() in {"ico", "interchainonly"}:
            interchainonly = True
        # toplim
        elif arg.lower().startswith("tl=") or\
             arg.lower().startswith("toplim="):
            try:
                toplim = int(float(arg.split('=', 1)[1]))
                assert toplim > 0
            except:
                raise ValueError("Inappropriate toplim value (positive integer): {}"\
                                 .format(arg.split('=', 1)[1]))
        # outplim
        elif arg.lower().startswith("ol=") or\
             arg.lower().startswith("outplim="):
            try:
                outplim = int(float(arg.split('=', 1)[1]))
                assert outplim > 0
            except:
                raise ValueError("Inappropriate outplim value (positive integer): {}"\
                                 .format(arg.split('=', 1)[1]))
        # conslim
        elif arg.lower().startswith("cl=") or\
             arg.lower().startswith("conslim="):
            try:
                conslim = int(float(arg.split('=', 1)[1]))
                assert conslim > 0
            except:
                raise ValueError("Inappropriate conslim value (positive integer): {}"\
                                 .format(arg.split('=', 1)[1]))
        # reactformat
        elif arg.lower().startswith("rf=") or\
             arg.lower().startswith("reactformat="):
            reactformat = arg.split('=', 1)[1]
            assert reactformat in {"3", "10", "26"},\
                   "Inappropriate reactformat value (3/10/26): {}"\
                   .format(arg.split('=', 1)[1])
            reactformat = int(float(reactformat))

    # Process rankby
    if "d" in rankby:
        rankbydiff = True
    if "r" in rankby and "s" in rankby:
        rankby = (0, 2, 1)
    elif "r" in rankby:
        rankby = (2, 0, 1)
    elif "s"  in rankby:
        rankby = (1, 2, 0)

    # Parse config
    paramsetnames, paramsets = ParseConfig(configfile)

    # Overwrite maxstemnum
    if maxstemnumset:
        for i in range(len(paramsets)):
            paramsets[i]['maxstemnum'] = maxstemnum

    # Reading + Running + Writing
    name = None
    data = []
    with open(inputfile) as file:
        for line in file:
            if line.startswith('>'):
                if name:
                    RunSQRNdbnseq(name, data, inputformat, paramsetnames,
                                  paramsets, threads, rankbydiff, rankby,
                                  hardrest, interchainonly, toplim, outplim,
                                  conslim, reactformat)
                name = line.strip()
                data = []
            else:
                data.append(line.strip())

    RunSQRNdbnseq(name, data, inputformat, paramsetnames,
                  paramsets, threads, rankbydiff, rankby,
                  hardrest, interchainonly, toplim, outplim,
                  conslim, reactformat)


    
