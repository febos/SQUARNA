import os
import sys


from SQRNdbnseq import SQRNdbnseq, ReactDict


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
                maxstemnum = int(arg.split('=', 1)[1])
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
            reactformat = int(reactformat)

    # Process rankby
    if "d" in rankby:
        rankbydiff = True
    if "r" in rankby and "s" in rankby:
        rankby = (0, 2, 1)
    elif "r" in rankby:
        rankby = (2, 0, 1)
    elif "s"  in rankby:
        rankby = (1, 2, 0)

    # Parse the config

    # Parse the input

    # Format output (check presence for each line type)






    
