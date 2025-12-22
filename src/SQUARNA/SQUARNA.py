import os
import sys
import io
from multiprocessing import Pool

try:
    from SQRNdbnseq import RunSQRNdbnseq, ReactDict, ProcessReacts, SEPS, GAPS
    from SQRNdbnali import RunSQRNdbnali
    from SQRNrfam   import SearchRfam
except:
    from .SQRNdbnseq import RunSQRNdbnseq, ReactDict, ProcessReacts, SEPS, GAPS
    from .SQRNdbnali import RunSQRNdbnali
    from .SQRNrfam   import SearchRfam

def ParseConfig(configfile):
    """Parses the config file"""
    # Set of mandatory parameters
    params = {"algorithms",
              "bpweights",
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
            # Ignore everything after # symbol
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
                    elif key == "algorithms":
                        paramset[key] = set(val.split(','))
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
                    reactivities = ProcessReacts(list(map(float, reactivities.split())))
                else:
                    reactivities = ProcessReacts([ReactDict[char] for char in reactivities])

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
            return "stockholm", 0
        
        if line1.startswith("CLUSTAL"):
            return "clustal", 0

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
            return "fasta", (entry_lines == 1)

    return "default", (entry_lines == 1)


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
            for seqname in seqnames], len(seqnames) == 1


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

    return [('>'+name, objs[name], None, None, None) for name in names], len(names) == 1


def ParseSeq(inputseq, returndefaults):

    if returndefaults:
        return None, None, None
    return [('>inputseq',inputseq,None,None,None),]


def ParseInput(inputseq, inputname, inputformat, returndefaults = False,
               fmt = "unknown", ignore = False):
    """Parser selector"""

    if inputseq:
        return ParseSeq(inputseq, returndefaults), fmt, True
    
    if fmt == "unknown":
        fmt, single_input = GuessFormat(inputname)
        if fmt != "default":
            print("Non-default input file format is recognized: {}".format(fmt.upper()))

    if fmt == "default":
        if returndefaults:
            return next(ParseDefaultInput(inputname, inputformat, returndefaults)), fmt
        return ParseDefaultInput(inputname, inputformat,
                                 returndefaults, ignore = ignore), fmt, single_input
    elif fmt == "fasta":
        if returndefaults:
            return next(ParseFasta(inputname, returndefaults)), fmt
        return ParseFasta(inputname, returndefaults), fmt, single_input
    elif fmt == "stockholm":
        parsed, single_input = ParseStockholm(inputname, returndefaults)
        return parsed, fmt, single_input
    elif fmt == "clustal":
        parsed, single_input = ParseClustal(inputname, returndefaults)
        return parsed, fmt, single_input

def byseqRunSQRNdbnseq(args):
    """multiprocessing (single-parameter) version of RunSQRNdbnseq;
       here we parallelize over the sequences instead of parallelizing
       over alternative structures within each sequence prediction"""
    name, seq, reacts, restrs, ref, theparamsetnames,\
    theparamsets, threads, rankbydiff, rankby,\
    hardrest, interchainonly, toplim, outplim,\
    conslim, reactformat, evalonly, poollim, entropy, \
    algos, levellimit, priority, rfam = args

    # We use a printing buffer so that the output is ordered
    # instead of being mixed up due to parallelization
    with io.StringIO() as buffer:
    
        RunSQRNdbnseq(name, seq, reacts, restrs, ref, theparamsetnames,
                      theparamsets, threads, rankbydiff, rankby,
                      hardrest, interchainonly, toplim, outplim,
                      conslim, reactformat, evalonly, poollim,
                      mp = False, sink = buffer, entropy = entropy,
                      algos = algos, levellimit = levellimit,
                      priority = priority, rfam = rfam)
        return buffer.getvalue()


def Predict(inputfile = None, fileformat = "unknown", inputseq = None,
            configfile = None, inputformat = "qtrf", maxstemnum = None,
            threads = os.cpu_count(), byseq = False, algorithms = '',
            entropy = False, rankby = "r", evalonly = False, hardrest = False,
            interchainonly = False, toplim = 5, outplim = None, conslim = 1,
            poollim = 1000, reactformat = 3, alignment = False, levellimit = None,
            freqlimit = 0.35, verbose = False, step3 = "u", ignorewarn = False,
            HOME_DIR = None, write_to = None, priority = None, rfam = False,
            i = None, ff = None, c = None, config = None, s = None, seq = None,
            a = None, ali = None, algo = None, algorithm = None, rb = None,
            fl = None, freqlim = None, ll = None, levlim = None, tl = None,
            ol = None, cl = None, pl = None, pr = None, s3 = None, msn = None,
            rf = None, eo = None, hr = None, ico = None, iw = None, ignore = None,
            t = None, bs = None, v = None):
    """
        -----------------------------------------------------------------------
        Prints SQUARNA RNA secondary structure predictions from the given input

        Parameters:
            i / inputfile : string
                Path to the input file.
            ff / fileformat : unknown/fasta/default/stockholm/clustal
                Input file format. For the details of the default SQUARNA input
                format see https://github.com/febos/SQUARNA/blob/main/examples/seq_input.fas.
                "unknown"   - the format will be identified automatically.
                "default"   - default fasta-like format.
                "fasta"     - FASTA format.
                "stockholm" - STOCKHOLM format.
                "clustal"   - CLUSTAL format.
            s / seq / inputseq : string
                RNA sequence input, alternative to inputfile (higher priority).
            c / config / configfile : string
                Path to a config file or a name of a built-in config, 
                see https://github.com/febos/SQUARNA/blob/main/def.conf
                for the format details. 
                In the alignment-based mode, the default config 
                file is ali.conf. In the single-sequence mode the default
                config for sequences under 500nts is def.conf, for sequences
                between 500 and 1000nts - 500.conf, and for sequences over
                1000nts in length - 1000.conf.
                Built-in configs:
                def (def.conf) is recommended by default for RNAs under 500nts.
                alt (alt.conf) is recommended for short pseudoknotted RNAs.
                500 (500.conf) is recommended for RNAs longer 500 nts.
                1000 (1000.conf) is recommended for RNAs longer 1000 nts.
                sk (sk.conf) is recommended with SHAPE data input.
                c=nussinov (nussinov.conf) - Nussinov algorithm config.
                c=hungarian (hungarian.conf) - Hungarian algorithm config.
                c=edmonds (edmonds.conf) - Edmonds algorithm config.
                c=greedy (greedy.conf) - Greedy algorithm config.
            inputformat : string
                The order of the lines in the input file. By default, SQUARNA 
                reads the first line of the entry (among the lines after 
                the ">" line) as the seQuence (q), the second line as the 
                reacTivities (t), the third line as the Restraints (r), 
                the fourth line as the reFerence (f), and all the further lines 
                are ignored. inputformat should be a subset of qtrfx letters 
                in any order, with q being mandatory. All "x" lines will be ignored.
            msn / maxstemnum : int
                Maximum number of stems to predict in each structure. By default,
                maxstemnum is defined in a config file for each parameter set. 
                If specified explicitly it will overwrite the maxstemnum 
                values for all the parameter sets.
            t / threads : int
                Number of CPUs to use.
            bs / byseq : bool
                Parallelize the execution over the input sequences
                in the single-sequence mode. 
                By default, the execution in the single-sequence mode
                is parallelized over the structure pool within each sequence.
                Parallelizing over input sequences is recommended for 
                large input files along with fast configs.
            algo / algorithm / algorithms: string
                The algorithms to be used in single-sequence predictions.
                By default, the algorithms are derived from the config file.
                If the algorithms parameter is specified, it will overwrite the
                algorithms listed in the config file.
                The choice should be a subset of the four algorithms:
                e - Edomnds algorithm [10.6028/jres.069B.013]
                g - Greedy SQUARNA algorithm [10.1101/2023.08.28.555103]
                h - Hungarian algorithm [10.1002/nav.3800020109]
                n - Nussinov algorithm [10.1073/pnas.77.11.6309]
            rb / rankby : string
                How to rank the predicted structures. rankby should be a subset of
                letters r, s, and d in any order (r / s / rs / rd / sd / rsd).
                If both r and s are present, the structures will be ranked according
                to the total_score = structure_score * reactivity_score. If only 
                r is present, the structures will be ranked by the reactivity_score,
                and if only s is present, the structures will be ranked by the 
                structure_score. Independently, if d is present, the mutually 
                divergent structures will be put first.
            eo / evalonly : bool
                Ignored in the alignment mode.
                If specified, no predictions are made and just the reference structure
                scores are returned provided the reference is specified. 
                If non-canonical base pairs are present in the reference structure, 
                they will be considered with 0.0 weight).
            hr / hardrest : bool
                If specified, all the base pairs from the restraints line will be
                forced to be present in the predicted structures. However, it will
                not affect the structure scores, as the forced base pairs won't
                contribute to the structure score unless they were predicted without
                forcing as well.
            ico / interchainonly : bool
                Allow only inter-chain base pairs to be predicted.
            tl / toplim : int
                How many top-N structures will be subject to comparison with the reference.
            ol / outplim : int
                How many top-N structures will be printed into the stdout.
                By default, outplim = toplim.
            cl / conslim : int
                How many top-N structures will be used to derive the predicted structure consensus.
            pl / poollim : int
                Maximum number of structures allowed to populate the current structure pool (if exceeded, no bifurcation will occur anymore).
            rf / reactformat : 3/10/26
                Encoding used to output the reactivities line.
                rf=3:   0.0  <= "_" <   1/3; 
                        1/3  <= "+" <   2/3;
                        2/3  <= "#" <=  1.0;
                rf=10:  0.0  <= "0" <   0.1;
                        ....................
                        0.5  <= "5" <   0.6;
                        ....................
                        0.9  <= "9" <=  1.0;
                rf=26:  0.00 <= "a" <  0.02;
                        0.02 <= "b" <  0.06;
                        ....................
                        0.50 <= "n" <  0.54;
                        ....................
                        0.94 <= "y" <  0.98;
                        0.98 <= "z" <= 1.00.
            a / ali/ alignment : bool
                Run SQUARNA in the alignment-based mode. If specified,
                ali.conf will be used as the config file by default, 
                unless another config file is explicitly specified 
                by the user. The bpweights, minlen, and minbpscore 
                parameters for step-1 will be derived from the first 
                parameter set in the config file.
            ll / levlim / levellimit : int
                Ignored in the single-sequence mode.
                The allowed number of pseudoknot levels. All the base pairs
                of the higher levels will be removed from the structure predicted
                at step-1 and from the structure predicted at step-2. By default, 
                levellimit=3 for short alignments of no more than 500 columns, 
                and levellimit=2 for longer alignments.
            fl / freqlim / freqlimit : 0.0 <= float <= 1.0
                Ignored in the single-sequence mode.
                The percentage of sequences required to contain a base pair,
                in order for it to be added to the predicted consensus structure
                at step-2. The consensus will include all the base pairs present
                in at least "fl" share of the sequences given that the base pair
                is not in conflict (does not share a position) with a more 
                populated base pair.
            v / verbose : bool
                Run SQUARNA in the verbose mode.
                Ignored in the single-sequence mode.
            s3 / step3 : "i"/"u"/"1"/"2"
                Ignored in the single-sequence mode.
                Defines the structure that will be printed at step-3. If step3=1,
                the structure from step-1 will be printed, and the step-2 will
                be skipped completely, meaning the prediction will be super fast.
                If step3=2, the structure from step-2 will be printed. If step3=u,
                the union of base pairs of the two structures will be printed. 
                If step3=i, the intersection of base pairs of the two structures 
                will be printed.
            iw / ignore / ignorewarn : bool
                Ignore warnings.
            rfam : bool
                Enable rfam search in case of a single input sequence
            HOME_DIR : string
                Path to the folder with built-in configs.
            write_to : IO_object
                Where to write the output. By default, write_to = sys.stdout.
            pr / priority : string
                Comma-separated list of prioritized paramset names. Default: bppN,bppH1,bppH2. 
    """

    # resolve synonyms
    if i != None:
        inputfile = i
    if ff != None:
        fileformat = ff
    if config != None:
        configfile = config
    if c != None:
        configfile = c
    if seq != None:
        inputseq = seq
    if s != None:
        inputseq = s
    if ali != None:
        alignment = ali
    if a != None:
        alignment = a
    if algorithm != None:
        algorithms = algorithm
    if algo != None:
        algorithms = algo
    if rb != None:
        rankby = rb
    if freqlim != None:
        freqlimit = freqlim
    if fl != None:
        freqlimit = fl
    if levlim != None:
        levellimit = levlim
    if ll != None:
        levellimit = ll
    if tl != None:
        toplim = tl
    if ol != None:
        outplim = ol
    if cl != None:
        conslim = cl
    if pl != None:
        poollim = pl
    if pr != None:
        priority = pr
    if s3 != None:
        step3 = s3
    if msn != None:
        maxstemnum = msn
    if rf != None:
        reactformat = rf
    if eo != None:
        evalonly = eo
    if hr != None:
        hardrest = hr
    if ico != None:
        interchainonly = ico
    if ignore != None:
        ignorewarn = ignore
    if iw != None:
        ignorewarn = iw
    if t != None:
        threads = t
    if bs != None:
        byseq = bs
    if v != None:
        verbose = v
    ##################

    if HOME_DIR is None:
        HOME_DIR = os.path.dirname(os.path.abspath(__file__))

    if write_to is None:
        write_to = sys.stdout

    if inputfile != None and not os.path.exists(inputfile) and\
       os.path.exists(os.path.join(HOME_DIR, inputfile)):
        inputfile = os.path.join(HOME_DIR, inputfile)
    
    # Verifying arguments
    assert os.path.exists(str(inputfile)) or inputseq, "Input file does not exist."

    assert fileformat in {'unknown','fasta','default','stockholm','clustal'},\
           "Wrong fileformat, choose one of these: default,fasta,stockholm,clustal"

    if configfile is None:
        configfileset  = False
        configfile     = os.path.join(HOME_DIR, "def.conf")
        configfile500  = os.path.join(HOME_DIR, "500.conf")
        configfile1000 = os.path.join(HOME_DIR, "1000.conf")
        if priority is None:
            priority = set('bppN,bppH1,bppH2'.split(','))
        else:
            priority = set([x for x in priority.split(',') if x])
    else:
        configfileset = True
        if not os.path.exists(configfile):
            if os.path.exists(os.path.join(HOME_DIR, configfile+".conf")):
                configfile = os.path.join(HOME_DIR, configfile+".conf")
            elif os.path.exists(os.path.join(HOME_DIR, configfile)):
                configfile = os.path.join(HOME_DIR, configfile)
        assert os.path.exists(configfile), "Config file does not exist."
        if priority is None:
            priority = set()
        else:
            priority = set([x for x in priority.split(',') if x])

    assert ''.join(sorted(inputformat.replace('x',''))) in {"q","fq","qr","qt", "qrt",
                                                                 "fqr", "fqt", "fqrt"}, \
                   'Inappropriate inputformat value (subset of "fqrtx" with "q" being mandatory): {}'\
                   .format(inputformat)

    if maxstemnum is None:
        maxstemnum = 10**6
        maxstemnumset = False
    else:
        maxstemnumset = True
        try:
            maxstemnum = int(float(maxstemnum))
            assert maxstemnum >= 0
        except:
            raise ValueError("Inappropriate maxstemnum value (non-negative integer): {}"\
                             .format(maxstemnum))
    try:
        threads = int(float(threads))
        threads = max(1, threads)
        threads = min(threads, os.cpu_count())
    except:
        raise ValueError("Inappropriate threads value (integer): {}"\
                         .format(threads))

    try:
        algos = set(algorithms.upper())
        assert algos <= {'E','G','H','N'}
    except:
        raise ValueError('Inappropriate algorithm value (should be subset of "eghn"): {}'\
                         .format(arg.split('=', 1)[1]))

    assert rankby in {"r", "s", "rs", "dr", "ds", "drs"}, \
           'Inappropriate rankby value (r/s/rs/dr/ds/drs): {}'\
           .format(rankby)

    if outplim is None:
        outplim = toplim
        outplimset = False
    else:
        outplimset = True
        try:
            outplim = int(float(outplim))
            assert outplim > 0
            outplimset = True
        except:
            raise ValueError("Inappropriate outplim value (positive integer): {}"\
                             .format(outplim))

    try:
        toplim = int(float(toplim))
        assert toplim > 0
        if not outplimset:
            outplim = toplim
    except:
        raise ValueError("Inappropriate toplim value (positive integer): {}"\
                         .format(toplim))

    try:
        conslim = int(float(conslim))
        assert conslim > 0
    except:
        raise ValueError("Inappropriate conslim value (positive integer): {}"\
                         .format(conslim))
    try:
        poollim = int(float(poollim))
        assert poollim > 0
    except:
        raise ValueError("Inappropriate poollim value (positive integer): {}"\
                         .format(poollim))

    assert int(float(reactformat)) in {3, 10, 26},\
           "Inappropriate reactformat value (3/10/26): {}"\
           .format(reactformat)
    reactformat = int(float(reactformat))

    if not (levellimit is None):
        try:
            levellimit = int(float(levellimit))
        except:
            raise ValueError("Inappropriate levellimit value (integer): {}"\
                             .format(levellimit))

    try:
        freqlimit = float(freqlimit)
        assert 0 <= freqlimit <= 1
    except:
        raise ValueError("Inappropriate freqlimit value (float between 0.0 and 1.0): {}"\
                         .format(freqlimit))

    try:
        step3 = step3.lower()
        assert step3 in {'u', 'i', '1', '2'}
    except:
        raise ValueError("Inappropriate freqlimit value (float between 0.0 and 1.0): {}"\
                         .format(step3))

    # Process rankby
    if "d" in rankby:
        rankbydiff = True # Output diverse structures first
    else:
        rankbydiff = False
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

            inputs, fmt, single_input = ParseInput(inputseq, inputfile, inputformat,
                                                   fmt = fileformat, ignore = ignorewarn)

            if rfam:
                if not single_input:
                    print("WARNING: Found more than one sequence, rfam search disabled.",
                          file=sys.stderr)
                    rfam = False
                else:
                    inputs = list(inputs)
                    inputs[0] = list(inputs[0])
                    inputs[0][3], rfam = SearchRfam(inputs[0][1], HOME_DIR)
            
            for name, seq, reacts, restrs, ref in inputs:
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
                              conslim, reactformat, evalonly, poollim, entropy = entropy,
                              algos = algos, levellimit = levellimit, sink = write_to,
                              priority = priority, rfam = rfam)
        # Parallelizing over input sequences
        else:
            batchsize = threads*10
            with Pool(threads) as pool:
                inputs_batch = []

                inputs, fmt, single_input = ParseInput(inputseq, inputfile, inputformat,
                                                       fmt = fileformat, ignore = ignorewarn)
                if rfam:
                    if not single_input:
                        print("WARNING: Found more than one sequence, rfam search disabled.",
                              file=sys.stderr)
                        rfam = False
                    else:
                        inputs = list(inputs)
                        inputs[0] = list(inputs[0])
                        inputs[0][3], rfam = SearchRfam(inputs[0][1], HOME_DIR)
                
                for name, seq, reacts, restrs, ref in inputs:
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
                                         conslim, reactformat, evalonly, poollim,
                                         entropy, algos, levellimit, priority, rfam))
                    # Process a batch once we have batchsize sequences
                    if len(inputs_batch) >= batchsize:
                        for output in pool.imap(byseqRunSQRNdbnseq, inputs_batch):
                            print(output, end = '', file = write_to)
                        inputs_batch = []
                # last batch (if anything left)
                if inputs_batch:
                    for output in pool.imap(byseqRunSQRNdbnseq, inputs_batch):
                        print(output, end = '', file = write_to)
                        

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
                    defReactivities = ProcessReacts(list(map(float, defReactivities.split())))
                else:
                    defReactivities = ProcessReacts([ReactDict[char] for char in defReactivities])

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
                      conslim, reactformat, poollim, entropy = entropy,
                      algos = algos, sink = write_to)


def Main():

    def PrintUsage():
        print()
        print("Usage:")
        print()
        print('SQUARNA i=inputfile [OPTIONS]')
        print()
        print('SQUARNA s=ACGUGUCAC [OPTIONS]')
        print()
        print("For further details read the help message:")
        print()
        print('SQUARNA --help')
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
    configfile     = None

    inputformat = "qtrf"               # Input line order, q=seQuence,t=reacTivities,r=Restraints,f=reFerence

    maxstemnum    = None               # maximum number of stems for each structure

    threads       = os.cpu_count()     # Number of cpus to use
    byseq         = False              # Parallelize by input sequences, not by structure pool

    rankby         = "r"               # Rank by, r / s / rs / dr / ds / drs, r=reactscore,s=structscore,d=rankbydiff
    evalonly       = False             # Just evaluate the reference and do not predict anything
    hardrest       = False             # Force bp-restraints into predicted structures 
    interchainonly = False             # Forbid intra-chain base pairs

    toplim        = 5                  # Top-N to print
    outplim       = None               # Top-N structs used for metrics calculations if reference
    conslim       = 1                  # Top-N structs used for consensus
    poollim       = 100                # Maximum number of structures allowed to populate the current
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
    entropy     = False                # Calculate stem matrix entropy
    algorithms  = ""                   # Single-sequence prediction algorithms

    priority    = None                 #Comma-separated list of prioritized paramset names

    rfam        = False                # Rfam template search for structural restraints

    # Allow standard parameter input
    formatted_args = []
    cnt = 0
    while cnt < len(args):
        if args[cnt].lower() in {"-algo", "--algo", "-algorithm", "--algorithm",
                                 "-algos", "--algos", "-algorithms", "--algorithms",
                                 "-i", "--i", "-input", "--input",
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
                                 "-pr", "--pr", "-priority", "--priority",
                                 "-s3", "--s3", "-step3", "--step3",
                                 "-msn", "--msn", "-maxstemnum", "--maxstemnum",
                                 "-rf", "--rf", "-reactformat", "--reactformat",
                                 "-s", "--s", "-seq", "--seq","-sequence", "--sequence",
                                 "-t", "--t", "-threads", "--threads",}:
            formatted_args.append(args[cnt].lstrip('-')+'='+args[cnt+1])
            cnt += 1
        elif args[cnt].lower() in {"-a", "--a", "-ali", "--ali", "-alignment", "--alignment",
                                   "-bs", "--bs", "-byseq", "--byseq",
                                   "-ent", "--ent", "-entropy", "--entropy",
                                   "-eo", "--eo", "-evalonly", "--evalonly",
                                   "-hr", "--hr", "-hardrest", "--hardrest",
                                   "-iw", "--iw", "-ignore", "--ignore",
                                   "-ico", "--ico", "-interchainonly", "--interchainonly",
                                   "-rfam", "--rfam", "-v", "--v", "-verbose", "--verbose",}:
            formatted_args.append(args[cnt].lstrip('-'))
        else:
            formatted_args.append(args[cnt])
        cnt += 1

    args = formatted_args
    # Parsing arguments
    for arg in args:
        # algorithms
        if arg.lower().startswith("algo=") or\
             arg.lower().startswith("algos=") or\
             arg.lower().startswith("algorithm=") or\
             arg.lower().startswith("algorithms="):
            if not arg.split('=', 1)[1]:
                continue
            algorithms = arg.split('=', 1)[1]
        # inputseq
        if arg.lower().startswith("s=") or\
           arg.lower().startswith("seq=") or\
           arg.lower().startswith("sequence="):
            inputseq = arg.split('=', 1)[1]
        # inputfile
        elif arg.lower().startswith("i=") or\
           arg.lower().startswith("input="):
            inputfile = arg.split('=', 1)[1]
        # fileformat
        elif arg.lower().startswith("ff=") or\
           arg.lower().startswith("fileformat="):
            fileformat = arg.split('=', 1)[1].lower()
        # configfile
        elif arg.lower().startswith("c=") or\
           arg.lower().startswith("config="):
            configfile = arg.split('=', 1)[1]
        # inputformat
        elif arg.lower().startswith("if=") or\
             arg.lower().startswith("inputformat="):
            inputformat = arg.split('=', 1)[1].lower()
        # maxstemnum
        elif arg.lower().startswith("msn=") or\
             arg.lower().startswith("maxstemnum="):
            maxstemnum = arg.split('=', 1)[1]
        # threads
        elif arg.lower().startswith("t=") or\
             arg.lower().startswith("threads="):
            threads = arg.split('=', 1)[1]
        # byseq
        elif arg.lower() in {"bs", "byseq"}:
            byseq = True
        # rankby
        elif arg.lower().startswith("rb=") or\
             arg.lower().startswith("rankby="):
            rankby = ''.join(sorted(arg.split('=', 1)[1].lower()))
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
            toplim = arg.split('=', 1)[1]
        # outplim
        elif arg.lower().startswith("ol=") or\
             arg.lower().startswith("outplim="):
            outplim = arg.split('=', 1)[1]
        # conslim
        elif arg.lower().startswith("cl=") or\
             arg.lower().startswith("conslim="):
            conslim = arg.split('=', 1)[1]
        # poollim
        elif arg.lower().startswith("pl=") or\
             arg.lower().startswith("poollim="):
            poollim = arg.split('=', 1)[1]
        # priority
        elif arg.lower().startswith("pr=") or\
             arg.lower().startswith("priority="):
            priority = arg.split('=', 1)[1]
        # reactformat
        elif arg.lower().startswith("rf=") or\
             arg.lower().startswith("reactformat="):
            reactformat = arg.split('=', 1)[1]
        # alignment
        elif arg.lower() in {"a", "ali", "alignment"}:
            alignment = True
        # levellimit
        elif arg.lower().startswith("ll=") or\
             arg.lower().startswith("levlim=") or\
             arg.lower().startswith("levellim=") or\
             arg.lower().startswith("levlimit=") or\
             arg.lower().startswith("levellimit="):
            levellimit = arg.split('=', 1)[1]
        # freqlimit
        elif arg.lower().startswith("fl=") or\
             arg.lower().startswith("freqlim=") or\
             arg.lower().startswith("freqlimit=") or\
             arg.lower().startswith("frequencylim=") or\
             arg.lower().startswith("frequencylimit="):
            freqlimit = arg.split('=', 1)[1]
        # verbose
        elif arg.lower() in {"v", "verbose"}:
            verbose = True
        # ignore
        elif arg.lower() in {"iw", "ignore"}:
            ignorewarn = True
        # entropy
        elif arg.lower() in {"ent", "entropy"}:
            entropy = True
        # rfam
        elif arg.lower() == 'rfam':
            rfam = True
        # step3
        elif arg.lower().startswith("s3=") or\
             arg.lower().startswith("step3="):
            step3 = arg.split('=', 1)[1]
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

    Predict(inputfile, fileformat, inputseq,
            configfile, inputformat, maxstemnum,
            threads, byseq, algorithms, entropy,
            rankby, evalonly, hardrest,
            interchainonly, toplim, outplim, conslim,
            poollim, reactformat, alignment, levellimit,
            freqlimit, verbose, step3, ignorewarn, HOME_DIR,
            None, priority, rfam)
        


if __name__ == "__main__":

    Main()


























        
            
