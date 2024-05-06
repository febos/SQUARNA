import os
import sys
import io
from multiprocessing import Pool

from SQRNdbnseq import RunSQRNdbnseq, ReactDict, SEPS, GAPS
from SQRNdbnali import RunSQRNdbnali


def ParseConfig(configfile):
    """Parses the config file"""
    # Set of mandatory parameters
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
            # Ignore everythin after # symbol
            cleanline = line.split('#', 1)[0].strip()
            # If non-empty line
            if cleanline:
                # If paramset name
                if cleanline.startswith('>'):
                    names.append(cleanline[1:])
                    cnt += 1
                    # If this is the first param set
                    if cnt == 1:
                        paramset = {}
                    # Otherwise 
                    else:  
                        paramsets.append(paramset)
                        # Init all the following sets with the first set values
                        paramset = {k:v for k, v in paramsets[0].items()}
                else:
                    key, val = cleanline.split(maxsplit = 1)
                    # bpweights require parsing values like GC=3,AU=2,GU=1
                    if key == "bpweights":
                        paramset[key] = {}
                        for kv in val.split(','):
                            k, v = kv.strip().split('=')
                            paramset[key][k] = float(v)
                    # all the other params are simply float values
                    else:
                        paramset[key] = float(val)
    # don't forget the last one
    paramsets.append(paramset)

    # Confirm the first param set contains all the params
    if not all([_ in paramsets[0] for _ in params]):
        raise ValueError("Missing some of the parameters in"+\
                         " the first parameter set"+\
                         " of the config file: {}"\
                         .format(', '.join([_ for _ in params
                                            if _ not in paramset])))
    return names, paramsets


def ParseDefaultInput(inputname, inputformat, returndefaults = False, ignore = False):
    """Returns object lists of format [name,sequence,reactivities,restraints,reference]
       or the list [default-reactivities, default-restraints, default reference] if
       returndefaults param is True"""

    warningsT = False
    warningsR = False
    warningsF = False

    def ProcessIndividual(data):
        """Returns a single [seq,reacts,rests,ref] list"""

        nonlocal warningsT, warningsR, warningsF

        while len(data) < len(inputformat):
            data.append(None)

        # split()[0] to allow space-separated comments
        # after the input within the line
        # for sequence, restraints, reference, but not for reactivities
        sequence     = data[q_ind].split()[0]
        reactivities = data[t_ind] if t_ind > 0 else None
        restraints   = data[r_ind].split()[0] if r_ind > 0 and data[r_ind] else None
        reference    = data[f_ind].split()[0] if f_ind > 0 and data[f_ind] else None

        N = len(sequence)

        # Fill features with default values if applicable
        if not reactivities and defT:
            if (len(defT) == N or len(defT.split()) == N):
                reactivities = defT
            elif not warningsT:
                warningsT = True
                if ignore:
                    print("WARNING: some sequences differ in length from the default reactivities line",
                          file=sys.stderr)
                else:
                    raise ValueError("WARNING: some sequences differ in length from the default reactivities line "+\
                                     "[Switch on the iw/ignore parameter to proceed anyway]")
        if not restraints and defR:
            if len(defR) == N:
                restraints   = defR
            elif not warningsR:
                warningsR = True
                if ignore:
                    print("WARNING: some sequences differ in length from the default restraints line",
                          file=sys.stderr)
                else:
                    raise ValueError("WARNING: some sequences differ in length from the default restraints line "+\
                                     "[Switch on the iw/ignore parameter to proceed anyway]")
        if not reference and defF:
            if len(defF) == N:
                reference    = defF
            elif not warningsF:
                warningsF = True
                if ignore:
                    print("WARNING: some sequences differ in length from the default reference line",
                          file=sys.stderr)
                else:
                    raise ValueError("WARNING: some sequences differ in length from the default reference line "+\
                                     "[Switch on the iw/ignore parameter to proceed anyway]")

        # Check reactivities for consistency and resolve them if needed
        try:
            if reactivities:
                if len(reactivities) != len(sequence):
                    reactivities = list(map(float, reactivities.split()))
                else:
                    reactivities = [ReactDict[char] for char in reactivities]

            assert not reactivities or len(reactivities) == len(sequence)
        except:
            raise ValueError('Inappropriate reactivities line for entry "{}":\n {}'\
                             .format(name[1:], reactivities))

        # Assert restraints and reference are of the consistent length
        # or empty line / None
        assert not restraints or len(restraints) == len(sequence),\
               'Inappropriate restraints line for entry "{}":\n {}'\
               .format(name[1:], restraints)
        assert not reference or len(reference) == len(sequence),\
               'Inappropriate reference line for entry "{}":\n {}'\
               .format(name[1:], reference)

        return sequence, reactivities, restraints, reference

    name = None
    defT = None
    defR = None
    defF = None
    data = []

    q_ind = inputformat.index('q')
    t_ind = inputformat.find('t')
    r_ind = inputformat.find('r')
    f_ind = inputformat.find('f')
    
    with open(inputname) as file:
        for line in file:
            if line.startswith('>'):
                # If not the first entry - process the previous one
                if name:
                    yield (name, *ProcessIndividual(data))
                else:
                    # Default TRF lines
                    defdata = data
                    while len(defdata) < len(inputformat) - 1:
                        defdata.append(None)
                    defdata.insert(q_ind, None)
                    defT = defdata[t_ind] if t_ind > 0 else None
                    defR = defdata[r_ind] if r_ind > 0 else None
                    defF = defdata[f_ind] if f_ind > 0 else None

                    if returndefaults:
                        yield (defT, defR, defF)
                        return None
                        
                name = line.strip()
                data = []
            else:
                data.append(line.strip())

    if name:
        yield (name, *ProcessIndividual(data))


def GuessFormat(inp):
    """Identify the input file format: default / fasta / stockholm / clustal"""
    with open(inp) as file:
        line1 = file.readline()

        entry_lines = 0
        seq_lines   = 0

        if line1.startswith('#') and "STOCKHOLM" in line1:
            return "stockholm"
        
        if line1.startswith("CLUSTAL"):
            return "clustal"

        if line1.startswith(">"):
            entry_lines += 1

        for line in file:
            if line.startswith(">"):
                entry_lines += 1
            else:
                acgut = sum(1 for ch in line.upper() if ch in {'A','C','G','U','T'})
                if acgut > len(line) / 2:
                    seq_lines += 1
                if seq_lines > 1000:
                    break

        if seq_lines > entry_lines and entry_lines > 0:
            return "fasta"

    return "default"


def ParseFasta(inp, returndefaults = False):

    if returndefaults:
        yield (None, None, None)
        return None

    name, seq = None, ''

    with open(inp) as file:
        for line in file:
            if line.startswith('>'):
                if name:
                    yield (name, seq, None, None, None)
                name = line.strip()
                seq  = ''
            elif line.strip():
                seq += line.strip()
    yield (name, seq, None, None, None)


def ReadStockholm(stkfile):
    """Parses Stockholm format into three lists and two dicts"""

    seqnames = [] # Sequence names
    seqdict  = {} # Sequence dict with name keys and sequence values
    gcnames  = [] # Structure names
    gcdict   = {} # Structure dict with name keys and structure values
    headers  = [] # Headers list

    try:
        file = open(stkfile)
    except:
        # Non-standard encoding found in some
        # of the Rfam families
        file = open(stkfile, encoding="iso8859-15")

    for line in file:
        if line.startswith('#=GC '): # Structure lines

            linesplit = line.strip().split()
            seq = linesplit[-1]
            name = ' '.join(linesplit[1:-1])

            if name not in gcdict:
                gcnames.append(name)
                gcdict[name] = seq
            else:
                gcdict[name] += seq

        elif line.startswith('#'):
            # Header lines
            headers.append(line)

        elif line.startswith('//'):
            pass
        elif not line.strip():
            pass
        else:
            # Sequence lines
            linesplit = line.strip().split()
            seq = linesplit[-1]
            name = ' '.join(linesplit[:-1])

            if name not in seqdict:
                seqnames.append(name)
                seqdict[name] = seq
            else:
                seqdict[name] += seq

    file.close()

    # Put #=GF lines to the end of the headers
    headers1 = [x for x in headers if not x.startswith("#=GF SQ")]
    headers2 = [x for x in headers if x.startswith("#=GF SQ")]
    headers = headers1 + headers2

    return headers, seqnames, seqdict, gcnames, gcdict


def ParseStockholm(inp, returndefaults = False):
    """Treats SS_cons dbn as default reference"""
    headers, seqnames, seqdict, gcnames, gcdict = ReadStockholm(inp)

    if returndefaults:
        return None, None, gcdict["SS_cons"] if "SS_cons" in gcnames else None

    return [('>'+seqname, seqdict[seqname], None, None,
             gcdict["SS_cons"] if "SS_cons" in gcnames else None)
            for seqname in seqnames]


def ParseClustal(inp, returndefaults = False):

    if returndefaults:
        return None, None, None

    objs = {}
    names = []

    with open(inp) as file:
        for line in file:
            if line.strip() and not line.startswith("CLUSTAL") and not line.startswith(' '):
                name, seq = line.strip().split()
                if name not in objs:
                    names.append(name)
                    objs[name] = ''
                objs[name] += seq

    return [('>'+name, objs[name], None, None, None) for name in names]


def ParseSeq(inputseq, returndefaults):

    if returndefaults:
        return None, None, None
    return [('>inputseq',inputseq,None,None,None),]


def ParseInput(inputseq, inputname, inputformat, returndefaults = False,
               fmt = "unknown", ignore = False):
    """Parser selector"""

    if inputseq:
        return ParseSeq(inputseq, returndefaults), fmt
    
    if fmt == "unknown":
        fmt = GuessFormat(inputname)
        if fmt != "default":
            print("Non-default input file format is recognized: {}".format(fmt.upper()))

    if fmt == "default":
        if returndefaults:
            return next(ParseDefaultInput(inputname, inputformat, returndefaults)), fmt
        return ParseDefaultInput(inputname, inputformat,
                                 returndefaults, ignore = ignore), fmt
    elif fmt == "fasta":
        if returndefaults:
            return next(ParseFasta(inputname, returndefaults)), fmt
        return ParseFasta(inputname, returndefaults), fmt
    elif fmt == "stockholm":
        return ParseStockholm(inputname, returndefaults), fmt
    elif fmt == "clustal":
        return ParseClustal(inputname, returndefaults), fmt

def byseqRunSQRNdbnseq(args):
    """multiprocessing (single-parameter) version of RunSQRNdbnseq;
       here we parallelize over the sequences instead of parallelizing
       over alternative structures within each sequence prediction"""
    name, seq, reacts, restrs, ref, theparamsetnames,\
    theparamsets, threads, rankbydiff, rankby,\
    hardrest, interchainonly, toplim, outplim,\
    conslim, reactformat, evalonly, poollim = args

    # We use a printing buffer so that the output is ordered
    # instead of being mixed up due to parallelization
    with io.StringIO() as buffer:
    
        RunSQRNdbnseq(name, seq, reacts, restrs, ref, theparamsetnames,
                      theparamsets, threads, rankbydiff, rankby,
                      hardrest, interchainonly, toplim, outplim,
                      conslim, reactformat, evalonly, poollim,
                      mp = False, sink = buffer)
        return buffer.getvalue()


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

    # If no arguments - print the short usage
    if not args:
        PrintUsage()

    # If asking for help message
    if "--help" in args or "-help" in args or "help" in args or\
       "--h" in args or "-h" in args or "h" in args or\
       "--H" in args or "-H" in args or "H" in args:
        with open(os.path.join(HOME_DIR,"README.md")) as helpfile:
            print(helpfile.read())
        exit(0)

    # DEFAULTS
    inputfile      = None
    fileformat     = "unknown"
    inputseq       = None
    configfile     = os.path.join(HOME_DIR, "def.conf")
    configfile500  = os.path.join(HOME_DIR, "500.conf")
    configfile1000 = os.path.join(HOME_DIR, "1000.conf")
    configfileset  = False             # Whether the user defined the config file

    inputformat = "qtrf"               # Input line order, q=seQuence,t=reacTivities,r=Restraints,f=reFerence

    maxstemnumset = False              # do we overwrite the maxstemnum from configfile
    maxstemnum    = 10**6              # maximum number of stems for each structure

    threads       = os.cpu_count()     # Number of cpus to use
    byseq         = False              # Parallelize by input sequences, not by structure pool

    rankbydiff     = False             # Output diverse structures first
    rankby         = "s"               # Rank by, r / s / rs / dr / ds / drs, r=reactscore,s=structscore,d=rankbydiff
    evalonly       = False             # Just evaluate the reference and do not predict anything
    hardrest       = False             # Force bp-restraints into predicted structures 
    interchainonly = False             # Forbid intra-chain base pairs

    toplim        = 5                  # Top-N to print
    outplimset    = False              # if the user specified the outplim value 
    outplim       = toplim             # Top-N structs used for metrics calculations if reference
    conslim       = 1                  # Top-N structs used for consensus
    poollim       = 1000               # Maximum number of structures allowed to populate the current
                                       # structure pool (if exceeded, no bifurcation will occur anymore)

    reactformat   = 3                  # 3 / 10 / 26

    alignment   = False                # Alignment mode
    levellimit  = None                 # Pseudoknot level threshold
    freqlimit   = 0.35                 # The percentage of sequences required to have the base pair
                                       # to include it into the step-2 result (alignment mode)
    verbose     = False                # Print the intermediate output or not
    step3       = "u"                  # i(intersection)/u(union)/1(step1)/2(step2) - what should be the
                                       # step-3 result dbn (alignment mode)

    ignorewarn  = False                # Ignore warnings


    # Allow standard parameter input
    formatted_args = []
    cnt = 0
    while cnt < len(args):
        if args[cnt].lower() in {"-i", "--i", "-input", "--input",
                                 "-c", "--c", "-config", "--config",
                                 "-if", "--if", "-inputformat", "--inputformat",
                                 "-rb", "--rb", "-rankby", "--rankby",
                                 "-ff", "--ff", "-fileformat", "--fileformat",
                                 "-fl", "--fl", "-freqlim", "--freqlim",
                                 "-ll", "--ll", "-levlim", "--levlim",
                                 "-tl", "--tl", "-toplim", "--toplim",
                                 "-ol", "--ol", "-outplim", "--outplim",
                                 "-cl", "--cl", "-conslim", "--conslim",
                                 "-pl", "--pl", "-poollim", "--poollim",
                                 "-s3", "--s3", "-step3", "--step3",
                                 "-msn", "--msn", "-maxstemnum", "--maxstemnum",
                                 "-rf", "--rf", "-reactformat", "--reactformat",
                                 "-s", "--s", "-seq", "--seq","-sequence", "--sequence",
                                 "-t", "--t", "-threads", "--threads",}:
            formatted_args.append(args[cnt].lstrip('-')+'='+args[cnt+1])
            cnt += 1
        elif args[cnt].lower() in {"-a", "--a", "-ali", "--ali", "-alignment", "--alignment",
                                   "-bs", "--bs", "-byseq", "--byseq",
                                   "-eo", "--eo", "-evalonly", "--evalonly",
                                   "-hr", "--hr", "-hardrest", "--hardrest",
                                   "-iw", "--iw", "-ignore", "--ignore",
                                   "-ico", "--ico", "-interchainonly", "--interchainonly",
                                   "-v", "--v", "-verbose", "--verbose",}:
            formatted_args.append(args[cnt].lstrip('-'))
        else:
            formatted_args.append(args[cnt])
        cnt += 1

    args = formatted_args
    # Parsing arguments
    for arg in args:
        # inputseq
        if arg.lower().startswith("s=") or\
           arg.lower().startswith("seq=") or\
           arg.lower().startswith("sequence="):
            inputseq = arg.split('=', 1)[1]
        # inputfile
        elif arg.lower().startswith("i=") or\
           arg.lower().startswith("input="):
            inputfile = arg.split('=', 1)[1]
            assert os.path.exists(inputfile), "Input file does not exist."
        # fileformat
        elif arg.lower().startswith("ff=") or\
           arg.lower().startswith("fileformat="):
            fileformat = arg.split('=', 1)[1].lower()
            assert fileformat in {'unknown','fasta','default','stockholm','clustal'},\
            "Wrong fileformat, choose one of these: default,fasta,stockholm,clustal"
        # configfile
        elif arg.lower().startswith("c=") or\
           arg.lower().startswith("config="):
            configfile = arg.split('=', 1)[1]
            configfileset = True
            assert os.path.exists(configfile), "Config file does not exist."
        # inputformat
        elif arg.lower().startswith("if=") or\
             arg.lower().startswith("inputformat="):
            inputformat = arg.split('=', 1)[1].lower()
            assert ''.join(sorted(inputformat.replace('x',''))) in {"q","fq","qr","qt", "qrt",
                                                                 "fqr", "fqt", "fqrt"}, \
                   'Inappropriate inputformat value (subset of "fqrtx" with "q" being mandatory): {}'\
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
        # byseq
        elif arg.lower() in {"bs", "byseq"}:
            byseq = True
        # rankby
        elif arg.lower().startswith("rb=") or\
             arg.lower().startswith("rankby="):
            rankby = ''.join(sorted(arg.split('=', 1)[1].lower()))
            assert rankby in {"r", "s", "rs", "dr", "ds", "drs"}, \
                   'Inappropriate rankby value (r/s/rs/dr/ds/drs): {}'\
                   .format(arg.split('=', 1)[1])
        # evalonly
        elif arg.lower() in {"eo", "evalonly"}:
            evalonly = True
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
                if not outplimset:
                    outplim = toplim
            except:
                raise ValueError("Inappropriate toplim value (positive integer): {}"\
                                 .format(arg.split('=', 1)[1]))
        # outplim
        elif arg.lower().startswith("ol=") or\
             arg.lower().startswith("outplim="):
            try:
                outplim = int(float(arg.split('=', 1)[1]))
                assert outplim > 0
                outplimset = True
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
        # poollim
        elif arg.lower().startswith("pl=") or\
             arg.lower().startswith("poollim="):
            try:
                poollim = int(float(arg.split('=', 1)[1]))
                assert poollim > 0
            except:
                raise ValueError("Inappropriate poollim value (positive integer): {}"\
                                 .format(arg.split('=', 1)[1]))
        # reactformat
        elif arg.lower().startswith("rf=") or\
             arg.lower().startswith("reactformat="):
            reactformat = arg.split('=', 1)[1]
            assert reactformat in {"3", "10", "26"},\
                   "Inappropriate reactformat value (3/10/26): {}"\
                   .format(arg.split('=', 1)[1])
            reactformat = int(float(reactformat))
        # alignment
        elif arg.lower() in {"a", "ali", "alignment"}:
            alignment = True
        # levellimit
        elif arg.lower().startswith("ll=") or\
             arg.lower().startswith("levlim=") or\
             arg.lower().startswith("levellim=") or\
             arg.lower().startswith("levlimit=") or\
             arg.lower().startswith("levellimit="):
            try:
                levellimit = int(float(arg.split('=', 1)[1]))
            except:
                raise ValueError("Inappropriate levellimit value (integer): {}"\
                                 .format(arg.split('=', 1)[1]))
        # freqlimit
        elif arg.lower().startswith("fl=") or\
             arg.lower().startswith("freqlim=") or\
             arg.lower().startswith("freqlimit=") or\
             arg.lower().startswith("frequencylim=") or\
             arg.lower().startswith("frequencylimit="):
            try:
                freqlimit = float(arg.split('=', 1)[1])
                assert 0 <= freqlimit <= 1
            except:
                raise ValueError("Inappropriate freqlimit value (float between 0.0 and 1.0): {}"\
                                 .format(arg.split('=', 1)[1]))
        # verbose
        elif arg.lower() in {"v", "verbose"}:
            verbose = True
        # ignore
        elif arg.lower() in {"iw", "ignore"}:
            ignorewarn = True
        # step3
        elif arg.lower().startswith("s3=") or\
             arg.lower().startswith("step3="):
            try:
                step3 = arg.split('=', 1)[1].lower()
                assert step3 in {'u', 'i', '1', '2'}
            except:
                raise ValueError("Inappropriate freqlimit value (float between 0.0 and 1.0): {}"\
                                 .format(arg.split('=', 1)[1]))
        else:
            if len(args) == 1:
                if os.path.exists(arg):
                    inputfile = arg
                elif sum(arg.lower().count(x) for x in (GAPS | set("acgut"))) > len(arg) / 2:
                    inputseq  = arg
                else:
                    inputfile = arg
            else:
                print("Unrecognized option: {}".format(arg))

    print(inputfile)
    assert os.path.exists(str(inputfile)) or inputseq, "Input file does not exist."

    # Process rankby
    if "d" in rankby:
        rankbydiff = True
    if "r" in rankby and "s" in rankby:
        rankby = (0, 2, 1)
    elif "r" in rankby:
        rankby = (2, 0, 1)
    elif "s"  in rankby:
        rankby = (1, 2, 0)

    # If alignment mode - use ali.conf by default
    if alignment and not configfileset:
        configfile = os.path.join(HOME_DIR, "ali.conf")

    # Parse config
    paramsetnames, paramsets = ParseConfig(configfile)

    # prepare 500.conf & 1000.conf for autoconfig
    if not configfileset:
        paramsetnames500,  paramsets500  = ParseConfig(configfile500)
        paramsetnames1000, paramsets1000 = ParseConfig(configfile1000)

    # Overwrite maxstemnum
    if maxstemnumset:
        for i in range(len(paramsets)):
            paramsets[i]['maxstemnum'] = maxstemnum
        if not configfileset:
            for i in range(len(paramsets500)):
                paramsets500[i]['maxstemnum'] = maxstemnum
            for i in range(len(paramsets1000)):
                paramsets1000[i]['maxstemnum'] = maxstemnum

    # Running single-sequence SQUARNA
    if not alignment:
        # Parallelizing over a structure pool for each sequence
        if not byseq:
            for name, seq, reacts, restrs, ref in ParseInput(inputseq, inputfile, inputformat,
                                                             fmt = fileformat, ignore = ignorewarn)[0]:
                # no autoconfig    
                if configfileset:
                    theparamsetnames, theparamsets = paramsetnames, paramsets
                # apply autoconfig
                else:
                    theparamsetnames, theparamsets = paramsetnames, paramsets
                    if len(seq) >= 500:
                        theparamsetnames, theparamsets = paramsetnames500, paramsets500
                    if len(seq) >= 1000:
                        theparamsetnames, theparamsets = paramsetnames1000, paramsets1000
                
                RunSQRNdbnseq(name, seq, reacts, restrs, ref, theparamsetnames,
                              theparamsets, threads, rankbydiff, rankby,
                              hardrest, interchainonly, toplim, outplim,
                              conslim, reactformat, evalonly, poollim)
        # Parallelizing over input sequences
        else:
            batchsize = threads*10
            with Pool(threads) as pool:
                inputs_batch = []
                for name, seq, reacts, restrs, ref in ParseInput(inputseq, inputfile, inputformat,
                                                                 fmt = fileformat, ignore = ignorewarn)[0]:
                    # no autoconfig    
                    if configfileset:
                        theparamsetnames, theparamsets = paramsetnames, paramsets
                    # apply autoconfig
                    else:
                        theparamsetnames, theparamsets = paramsetnames, paramsets
                        if len(seq) >= 500:
                            theparamsetnames, theparamsets = paramsetnames500, paramsets500
                        if len(seq) >= 1000:
                            theparamsetnames, theparamsets = paramsetnames1000, paramsets1000
                    # Collecting inputs
                    inputs_batch.append((name, seq, reacts, restrs, ref, theparamsetnames,
                                  theparamsets, threads, rankbydiff, rankby,
                                  hardrest, interchainonly, toplim, outplim,
                                  conslim, reactformat, evalonly, poollim))
                    # Process a batch once we have batchsize sequences
                    if len(inputs_batch) >= batchsize:
                        for output in pool.imap(byseqRunSQRNdbnseq, inputs_batch):
                            print(output, end = '')
                        inputs_batch = []
                # last batch (if anything left)
                if inputs_batch:
                    for output in pool.imap(byseqRunSQRNdbnseq, inputs_batch):
                        print(output, end = '')
                        

    else: # Running alignment-based SQUARNA

        # Get the processed sequences
        objs, fmt = ParseInput(inputseq, inputfile, inputformat,
                               fmt = fileformat, ignore = ignorewarn)

        # Get the default input lines
        defReactivities, defRestraints, defReference = ParseInput(inputseq, inputfile, inputformat,
                                                                  returndefaults = True,
                                                                  fmt = fmt,
                                                                  ignore = ignorewarn)[0]

        objs = [obj for obj in objs] # Unpack generator (in case of fasta/default format)

        # Length checks
        N = len(objs[0][1])
        assert all(len(obj[1]) == N for obj in objs),\
               'The sequences are not aligned'

        # Check reactivities for consistency and resolve them if needed
        try:
            if defReactivities:
                if len(defReactivities) != N:
                    defReactivities = list(map(float, defReactivities.split()))
                else:
                    defReactivities = [ReactDict[char] for char in defReactivities]

            assert not defReactivities or len(defReactivities) == N
        except:
            raise ValueError('Inappropriate default reactivities line:\n {}'\
                             .format(defReactivities))

        # Assert restraints and reference are of the consistent length
        # or empty line / None
        assert not defRestraints or len(defRestraints) == N,\
               'Inappropriate default restraints line:\n {}'\
               .format(defRestraints)
        assert not defReference or len(defReference) == N,\
               'Inappropriate default reference line:\n {}'\
               .format(defReference)

        # default levellimit
        if levellimit is None:
            levellimit = 3 - int(N > 500)

        # Run the alignment-based predictions
        RunSQRNdbnali(objs, defReactivities, defRestraints, defReference,
                      levellimit, freqlimit, verbose, step3,
                      paramsetnames, paramsets, threads, rankbydiff, rankby,
                      hardrest, interchainonly, toplim, outplim,
                      conslim, reactformat, poollim)
        





























        
            
